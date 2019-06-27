import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interpolate
from scipy.special import gamma
import rbdl
import pandas as pd
import matplotlib.pyplot as plt

def func(x, t, ftau, lua_model):
    c = ftau(t)
    nDof = lua_model.dof_count

    #differetnial eq: dx/dt = [qd, qdd]    
    q_int = x[:nDof]
    qd_int = x[nDof:]
    qdd =  np.zeros((2))
    dxdt = np.zeros(x.shape)

    rbdl.ForwardDynamics(lua_model,q_int,qd_int,c,qdd)

    dxdt[:nDof] = qd_int[:]
    dxdt[nDof:] = qdd[:]

    return dxdt 


if __name__ == "__main__":
    #load the lua model
    lua_model = rbdl.loadModel("/home/patrick/Documents/cartPenudlum_template/model.lua")

    #load a data file
    data = pd.read_csv("/home/patrick/Documents/cartPenudlum_template/build/RES/meshup_cart_pendulum_count_0000.csv")
    data = data.to_numpy()

    
    time = data[:,0]

    #assign the values
    q0 = data[:,1]
    q1 = data[:,2]
    qd0 = data[:,3]
    qd1 = data[:,4]
    tau0 = data[:,5]

    q = [q0,q1]
    q = np.asarray(q, dtype="double")

    qd = [qd0, qd1]
    qd = np.asarray(qd,dtype="double")

    #define state vector
    xState = np.concatenate((q,qd), axis=0)
    xState = np.asarray(xState)


    #get starting values, here iteration over all values
    tau1 = np.zeros(tau0.shape)
    tau = [tau0,tau1]
    tau = np.asarray(tau, dtype="double")

    #prepare data save array
    output = np.zeros((1,5))

    #interpolate tau
    ftau = interpolate.interp1d(time, tau, kind="linear", bounds_error=False, fill_value=(tau0[0],tau0[0]))

    for i in range(time.size-1):
        #get initial state
        xState_0 = xState[:,i]

        #integrator
        tnew = np.linspace(time[i],time[i+1], num=100, dtype="float")
               
        #integrate
        args=(ftau,lua_model)
        x = integrate.odeint(func, xState_0, tnew, args)
        
        #create control 
        control = ftau(tnew)[0]

        #save data in new array, need to append at right axis, and append control, to make one output array
        int_out = [x[:,0],x[:,1],x[:,2],x[:,3],control]
        int_out = np.asarray(int_out)
        int_out = int_out.transpose()

        output = np.append(output, int_out, axis=0)

    output = output[1:,:]
    t1 = np.linspace(time[0],time[-1], num=output.shape[0])    

    #create labeling
    labels = ["pos x","pendulum angle", "vel x", "angle vel", "control"]

    #plot the result, need to adjust time of q_int, as it is 10*999
    fig, axes = plt.subplots(1,5)
    for i in range(xState.shape[0]):
        ax = axes[i]
        ax.plot(time, xState[i,:], linestyle="--", label="original "+labels[i])
        ax.plot(t1, output[:,i],  label="integrated "+labels[i])
        if i == 0:
            # Shrink current axis by 20%
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width, box.height*0.7])
            
            # Put a legend to the right of the current axis
            ax.legend(loc='center', bbox_to_anchor=(0.5, 1.1))
        else:           
            # Shrink current axis by 20%
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width , box.height*0.7])

            # Put a legend to the right of the current axis
            ax.legend(loc='center', bbox_to_anchor=(0.5, 1.1))
    ax = axes[4]
    ax.plot(time,tau[0], label="original tau")
    ax.plot(t1,output[:,4], label="used tau")
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height*0.7])

    # Put a legend to the right of the current axis
    ax.legend(loc='center', bbox_to_anchor=(0.5, 1.1))   
    plt.savefig("/home/patrick/Documents/cart_pred/nn/visualisations/integration_origial.pdf", format="pdf")
    plt.show()


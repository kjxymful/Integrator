import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d
from scipy.special import gamma
import rbdl
import pandas as pd
import matplotlib.pyplot as plt

def func(x, t,ftau,lua_model):
    #integrator
    #differetnial eq: dx/dt = [qd, qdd]    
    q_int = np.zeros((2))
    qd_int = np.zeros((2))
    qdd = np.zeros((2))
    dxdt = np.zeros(xState_0.shape)
    
    c=ftau(t)
    for i in range(nDof):
        q_int[i] = x[i]
        qd_int[i] = x[i+nDof]
    rbdl.ForwardDynamics(lua_model,q_int,qd_int,c,qdd)
    k=0
    for j in range(nDof):
        dxdt[k] = qd_int[j]
        k=k+1
    for j in range(nDof):
        dxdt[k] = qdd[j]
        k=k+1
    return dxdt

if __name__ == "__main__":
    #load the lua model
    lua_model = rbdl.loadModel("/home/patrick/Documents/cartPenudlum_template/model.lua")
    nDof = lua_model.dof_count

    #load a data file
    data = pd.read_csv("/home/patrick/Documents/cart_pred/interpolated_data/interpolated_data_0.csv")
    data = data.to_numpy()

    time = np.linspace(0.185, 3.7, 1000)

    #assign the values
    q0 = data[:,0]
    q1 = data[:,1]
    qd0 = data[:,2]
    qd1 = data[:,3]
    tau0 = data[:,4]

    tau_time = np.linspace(0.185,3.7,tau0.shape[0])

    q = [q0,q1]
    q = np.asarray(q, dtype="double")

    qd = [qd0, qd1]
    qd = np.asarray(qd,dtype="double")

    #define state vector
    xState = np.concatenate((q,qd), axis=0)
    xState = np.asarray(xState)

    #prepare data save array
    output = np.zeros((1,5))

    #create tau
    tau1 = np.zeros(tau0.shape)
    tau = [tau0,tau1]
    tau = np.asarray(tau, dtype="double")

    #prepare data save array
    output = np.zeros((1,5))

    #interpolate tau
    ftau = interp1d(tau_time, tau, kind="linear", bounds_error=False, fill_value=(tau0[0],tau0[-1]))




    for i in range(q.shape[1]):
        #get starting values
        xState_0 = xState[:,i]
                
        args = (ftau,lua_model)
        t = np.arange(time[i],time[i+1], (time[i+1]-time[i])/10)
   
        
        #integrate
        q_int = integrate.odeint(func, xState_0, t, args)

        #create control
        control = ftau(t)[0]


        #save data in new array, need to append at right axis, and append control, to make one output array
        int_out = [q_int[:,0],q_int[:,1],q_int[:,2],q_int[:,3],control]
        int_out = np.asarray(int_out)
        int_out = int_out.transpose()

        output = np.append(output, int_out, axis=0)


    output = output[1:,:]
    t = time[:-1]
    t1 = np.linspace(0.185,3.7,10*999)

  

    #plot the result, need to adjust time of q_int, as it is 10*999
    fig, axes = plt.subplots(1,4)
    for i in range(xState.shape[0]):
        ax=axes[i]
        ax.plot(t, xState[i,:], label="original"+str(i))
        ax.plot(t1, output[:,i], linestyle="--", label="integrated"+str(i))
        ax.legend()

    plt.savefig("/home/patrick/Documents/cart_pred/nn/visualisations/integration1.pdf", format="pdf")
    plt.show()


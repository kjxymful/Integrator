/*====================================================================
 * cam-plane pendulum example
 * Copyright (c) 2015 Matthew Millard 
 * <matthew.millard@iwr.uni-heidelberg.de>
 *
 *///=================================================================


#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h> 
#include <rbdl/rbdl.h>
#include <rbdl/addons/luamodel/luamodel.h>
#include <rbdl/addons/geometry/geometry.h>
#include "csvtools.h"

#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/stepper/generation/make_controlled.hpp>

//using namespace std;
using namespace boost::numeric::odeint;

using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

//====================================================================
// Boost stuff
//====================================================================

typedef std::vector< double > state_type;
typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;

std::string outLoc;

//Multibody Variables
Model* model;
VectorNd q, qd, qdd, tau, x;  

class RightHandSide {

  public:
    RightHandSide() { };

    void operator() (const state_type &x,
                     state_type &dxdt, 
                     const double t){

        //q
        int j = 0;
        for(unsigned int i=0; i<model->q_size; i++){                
            q[i] = double(x[j]);
            j++;
        }

        //qd
        for(unsigned int i=0; i<model->qdot_size; i++){
            qd[i] = double(x[j]);
            j++;
        }

        //tau = 0
        for(unsigned int i=0; i<model->qdot_size; i++){                
            tau[i] = 0;
        }

    
        ForwardDynamics(*model,q,qd,tau,qdd);

        //populate dxdt
        //If you are using a quaternion joint, you must map wx,wy,wz to the 
        //derivatives of the 4 quaternion components using the omegaToQDot 
        //function
        j = 0;
        for(unsigned int i = 0; i < model->q_size; i++){
            dxdt[j] = double(qd[i]);
            j++;
        }
        for(unsigned int i = 0; i < model->qdot_size; i++){
            dxdt[j] = double(qdd[i]);
            j++;
        }
    }  
};


/* Problem Constants */
// Main.
int main(int argc, char** argv) {

  std::string fileName;

  // TODO lade pendulum.lua
  if (!Addons::LuaModelReadFromFile(fileName.c_str(),model)){
    std::cerr << "Error loading LuaModel: " << fileName << std::endl;
    abort();
  }

  // TODO lade meshup_cart_pendulum_count_0000.csv


  // TODO for loop über alle shooting notes
  // for ()
  // {

    q.resize(model->dof_count);
    qd.resize(model->dof_count);
    tau.resize(model->dof_count); // tau[0] = u

    // TODO Anfangswerte setzen, dh q, qd, tau (Werte aus meshup_cart_pendulum_count_0000.csv)
    // q = ...
    // qd = ...
    // tau = ...

    x.resize(model->dof_count*2);
    for(unsigned int i=0; i<q.rows();++i){
      x[i] =q[i];
      x[i+q.rows()] = qd[i];
    }

    RightHandSide rhs();

    state_type xState(x.size());
    state_type dxState(x.size());
    for(unsigned int i=0; i<x.size(); ++i){
      xState[i]   = x[i];
    }

    // TODO Anfangswerte setze, dh t0 und t1 (Werte aus meshup_cart_pendulum_count_0000.csv)
    double t;
    double t0; 
    double t1;
    unsigned int npts      = 100;

    // t0 = ...
    // t1 = ...

    double absTolVal = 1e-8;
    double relTolVal = 1e-8;

    double dt = (t1-t0)/(npts-1);
    unsigned int k=0;

    std::vector<std::vector< double > > matrixData;
    std::vector< double > rowData(model->dof_count+1);

    double a_x = 1.0 , a_dxdt = 1.0;
    controlled_stepper_type
    controlled_stepper(
        default_error_checker< double ,
                              range_algebra ,
                              default_operations >
        ( absTolVal , relTolVal , a_x , a_dxdt ) );

    double tp = 0;
    rowData[0] = 0;
    for(unsigned int z=0; z < model->dof_count; z++){
        rowData[z] = xState[z];
    }
    matrixData.push_back(rowData);

    for(unsigned int i=0; i<= npts; ++i){
      t = t0 + dt*i;

      integrate_adaptive(
          controlled_stepper ,
          rhs , xState , tp , t , (t-tp)/10 );
      tp = t;


      for(unsigned int j=0; j<x.rows();++j){
        x[j] = xState[j];
      }
      k=0;
      for(unsigned int j=0; j<model->q_size;++j){
        q[j] = xState[k];
        ++k;
      }
      for(unsigned int j=0; j<model->qdot_size;++j){
        qd[j] = xState[k];
        ++k;
      }


      rowData[0] = t;
      for(unsigned int z=0; z < model->dof_count; z++){
          rowData[z+1] = xState[z];
      }
    }
  // }

    std::string emptyHeader("");
    std::string fileNameOut(outLoc + ".csv");
    printMatrixToFile(matrixData,emptyHeader,fileNameOut);

   return 0; 
}

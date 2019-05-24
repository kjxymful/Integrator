/*====================================================================
 * cam-plane pendulum example
 * Copyright (c) 2015 Matthew Millard 
 * <matthew.millard@iwr.uni-heidelberg.de>
 *
 *///=================================================================


#include <iostream>
#include <fstream>
#include <rbdl/rbdl.h>
#include <Eigen/Core>

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


//Multibody Variables
Model* model;
VectorNd q, qd, qdd, tau, x;  

//Read a csv file
template<typename M>
M ReadCsv(const std::string& path) {
    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<double> values;
    uint rows = 0;
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }
    return Eigen::Map<const Eigen::Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, Eigen::RowMajor>>(values.data(), rows, values.size()/rows);
}

// Save data from Eigen matrix to .csv file.
const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
template<typename M> 
void WriteCsv(const std::string& path, M matrix){
    std::ofstream file(path.c_str());
    file << matrix.format(CSVFormat);
};

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

  //load lua model
  std::string filename;
  filename= "/home/patrick/Documents/cartPenudlum_template/model.lua";
  Addons::LuaModelReadFromFile(filename.c_str(), model);

  if (!Addons::LuaModelReadFromFile(filename.c_str(),model)){
    std::cerr << "Error loading LuaModel: " << filename << std::endl;
    abort();
  }

//load data in eigen matrix
  std::string path = "/home/patrick/Documents/cartPenudlum_template/build/RES/meshup_cart_pendulum_count_0000.csv";
  Eigen::MatrixXd data = ReadCsv<Eigen::MatrixXd>(path);

  // TODO for loop über alle shooting notes
  int r = data.rows();
  int c = data.cols();

  // std::vector<std::vector< double > > matrixData;
  Eigen::MatrixXd matrixData(0, c);
  // std::vector< double > rowData(model->dof_count+1);
  Eigen::MatrixXd rowData(1, model->dof_count+1);

  for (int i = 0; i<x; i++){

    Eigen::MatrixXd m(r,c-1);
    Eigen::MatrixXd qm(r,c-1);
    Eigen::MatrixXd qdm1(r,c-1);
    Eigen::MatrixXd qdm(r,c-1);
    Eigen::MatrixXd taum(r,c-1);

    m=data.rightCols(c-1);
    qm=m.leftCols(2);
    qdm1=m.rightCols(3);
    qdm=qdm1.leftCols(1);
    taum=m.rightCols(1);

    q.resize(model->dof_count);
    qd.resize(model->dof_count);
    tau.resize(model->dof_count); // tau[0] = u

      // TODO Anfangswerte setzen, dh q, qd, tau (Werte aus meshup_cart_pendulum_count_0000.csv)
    q = qm.row(i).matrix();
    qd = qdm.row(i).matrix();
    tau = taum.row(i).matrix();

    x.resize(model->dof_count*2);
    for(unsigned int i=0; i<q.rows();++i){
        x[i] =q[i];
        x[i+q.rows()] = qd[i];


    RightHandSide rhs;

    state_type xState(x.size());
    state_type dxState(x.size());
    for(unsigned int i=0; i<x.size(); ++i){
      xState[i]   = x[i];
      }

      // TODO Anfangswerte setze, dh t0 und t1 (Werte aus meshup_cart_pendulum_count_0000.csv)
      double t; //t?
      double t0=data(0,i); 
      double t1=data(0,i+1);
      unsigned int npts      = 100;


      double absTolVal = 1e-8;
      double relTolVal = 1e-8;

      double dt = (t1-t0)/(npts-1);
      unsigned int k=0;

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
      // matrixData.push_back(rowData);
      matrixData.conservativeResize(matrixData.rows()+1, Eigen::NoChange_t);
      matrixData.bottomRows(1) = rowData;

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
        matrixData.conservativeResize(matrixData.rows()+1, Eigen::NoChange_t);
        matrixData.bottomRows(1) = rowData;
      }
    }
  }
  std::string file_path = "/home/patrick/Documents/integrator/Data/integrated_data_0000.csv";
  WriteCsv(file_path, matrixData);

  return 0; 
}

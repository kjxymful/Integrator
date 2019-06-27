/*====================================================================
 * cam-plane pendulum example
 * Copyright (c) 2015 Matthew Millard 
 * <matthew.millard@iwr.uni-heidelberg.de>
 *
 *///=================================================================


#include <iostream>
#include <fstream>
#include <rbdl/rbdl.h>
#include <rbdl/addons/luamodel/luamodel.h>
#include <Eigen/Core>
#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/stepper/generation/make_controlled.hpp>

#include<sstream>
#include<string>
#include<iomanip>
#include <sys/stat.h>

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
VectorNd q, qd, qdd, tau;  

//==========================
//Read a csv file
//==========================
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
          values.push_back(std::stold(cell));
        }
        ++rows;
    }
    return Eigen::Map<const Eigen::Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, Eigen::RowMajor>>(values.data(), rows, values.size()/rows);
}
//======================
//write csv
//======================
const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
template<typename M> 
void WriteCsv(const std::string& path, M matrix){
    std::ofstream file(path.c_str());
    file << matrix.format(CSVFormat);
};

//======================
//counter
//======================
bool fileExists(const std::string& file) {
    struct stat buf;
    return (stat(file.c_str(), &buf) == 0);
}

int counter()
{
    int returnedCount = 0;
    int possibleMax = 5000; //some number you can expect.

    for (int starter = 0; starter < possibleMax; starter++){
        std::stringstream fileName_s("");
        fileName_s << "/home/patrick/Documents/integrator/Data/integrated_data_"
         << std::setfill('0') << std::setw(4) << std::to_string(starter) << ".csv";
         std::string fileName= fileName_s.str();
      
        bool status = fileExists(fileName);

        returnedCount = 1;//starter;

        if (!status)
            break;
        }

    return returnedCount;
}

//===================================
//naming the file(should be in main)
//===================================
std::string naming_file()
{
  std::string name;
  std::stringstream data_file_name_o("");
  data_file_name_o << "integrated_data_" << std::setfill('0') << std::setw(4) << std::to_string(counter()) << ".csv"; 
  
  name=data_file_name_o.str();
  return name;
}

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
  filename= "/home/patrick/Documents/integrator/model.lua";

  model = new Model();

  if (!Addons::LuaModelReadFromFile(filename.c_str(), model)){
    std::cerr << "Error loading LuaModel: " << filename << std::endl;
    abort();
  }

//load data in eigen matrix
  std::string path = "/home/patrick/Documents/cartPenudlum_template/build/RES/meshup_cart_pendulum_count_0000.csv";
  Eigen::MatrixXd data = ReadCsv<Eigen::MatrixXd>(path);
  
  // int l;  
  // int m;
  // for (l=0; l<data.rows(); l++){
  //   for (m=0; m<data.cols(); m++){
  //       if (1e-10 > fabs(data(l,m))) {
  //         data(l,m)=0;
  //       }
  //   }
  // }

  // TODO for loop über alle shooting notes
  int r = data.rows();
  int c = data.cols();

  Eigen::MatrixXd matrixData(1, c);
  Eigen::MatrixXd rowData(1, c);


  for (int j = 0; j<r-1; j++){
 
    q.resize(model->dof_count);
    qd.resize(model->dof_count);
    tau.resize(model->dof_count); // tau[0] = u
    qdd.resize(model->dof_count);

     // TODO Anfangswerte setzen, dh q, qd, tau (Werte aus meshup_cart_pendulum_count_0000.csv)
    q = data.middleCols(1, 2).row(j).transpose();
    qd = data.middleCols(3, 2).row(j).transpose();
    tau[0] = data(j, c-1);
    tau[1] = 0.;

    RightHandSide rhs;

    state_type xState(model->dof_count*2);
    for(unsigned int i=0; i<q.rows(); ++i){
        xState[i] =q[i];
        xState[i+q.rows()] = qd[i];
    }

    // TODO Anfangswerte setze, dh t0 und t1 (Werte aus meshup_cart_pendulum_count_0000.csv)
    double t0=data(j,0); 
    double t1=data(j+1,0);
    double t = t0;
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

    rowData(0,0) = t;
    for(unsigned int z=0; z < xState.size(); z++){
        rowData(0,z+1) = xState[z];
    }
    // maybe save data here


      for(unsigned int n=0; n< npts; ++n){
        t = t0 + dt*n;

        integrate_adaptive(
            controlled_stepper ,
            rhs , xState , dt*n, dt*(n+1), dt );

        rowData(0,0) = t+dt;
        for(unsigned int z=0; z < xState.size(); z++){
            rowData(0,z+1) = xState[z];
        }
      matrixData.conservativeResize(matrixData.rows()+1, matrixData.cols());
      matrixData.bottomRows(1) = rowData;
      }
  }
  std::string file_path = "/home/patrick/Documents/integrator/Data/" + naming_file();
  WriteCsv(file_path, matrixData);
  return 0; 
}

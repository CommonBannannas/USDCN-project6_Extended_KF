#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  // check inputs, vector size not zero and equal size of vectors
  if (estimations.size() != ground_truth.size() || estimations.size() == 0){
    std::cout << "Invalid estimation or ground_truth data" << std::endl;
    return rmse;
  }

  // squared residuals
  for (unsigned int i = 0; i < estimations.size(); ++i){

    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3, 4);
  // read state params
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float c0 = (px*px) + (py*py);
  float c1 = sqrt(c0);
  float c2 = (c0*c1);

  // check div by zero
  if (fabs(c0) < 0.0001){
    std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    return Hj;
  }
  //compute Jacobian
  Hj << (px/c1), (py/c1), 0, 0,
  -(py/c0), (px/c0), 0, 0,
  py*(vx*py-vy*px)/c2, px*(px*vy - py*vx)/c2, px/c1, py/c1;
  return Hj;
}

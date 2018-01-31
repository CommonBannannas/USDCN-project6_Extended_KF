#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

static VectorXd cartesian_polar(const VectorXd &x_state);
static VectorXd norm_polar(const VectorXd &z, const VectorXd &z_pred);

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = PHt * Si;

  // new estimation
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  VectorXd z_pred = cartesian_polar(x_);
  VectorXd y = norm_polar(z, z_pred);

  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

static VectorXd cartesian_polar(const VectorXd &x_state)
{

	const float eps	= 0.000001;
	float px, py, vx, vy;
	px = x_state[0];
	py = x_state[1];
	vx = x_state[2];
	vy = x_state[3];

	float rho, phi, rho_dot;
	rho = sqrt(px*px + py*py);
	phi = atan2(py, px);  // returns values between -pi and pi

	if(rho < eps)
		rho = eps;

	rho_dot = (px * vx + py * vy) / rho;

	VectorXd z_pred = VectorXd(3);
	z_pred << rho, phi, rho_dot;

	return z_pred;
}

static VectorXd norm_polar(const VectorXd &z, const VectorXd &z_pred)
{
	  VectorXd y = z - z_pred;

	  // normalize the angle
	  while(y(1) > M_PI){
	    y(1) -= M_PI*2;
	  }

	  while(y(1) < -M_PI){
	    y(1) += M_PI*2;
	  }
	  return y;
}

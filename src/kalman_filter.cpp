#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {};

KalmanFilter::~KalmanFilter() {};

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
};

void KalmanFilter::Predict() {
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
};

void KalmanFilter::Update(const VectorXd &z) {
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd K = P_ * Ht * S.inverse();
    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
};

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */
    Tools tools;
    VectorXd z_pred = tools.Cartesian2Polar(x_);
    if (z_pred[1]>M_PI){
        z_pred[1]-=(2*M_PI);
    }
    if (z_pred[1]<-M_PI){
        z_pred[1]+=(2*M_PI);
    }
    VectorXd y = z - z_pred;

    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd K = P_ * Ht * S.inverse();

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
};

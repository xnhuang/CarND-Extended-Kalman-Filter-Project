#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

    H_laser_<< 1,0,0,0,
                0,1,0,0;
    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
        0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

    //set the acceleration noise components
    noise_ax = 6;
    noise_ay = 6;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
    if (!is_initialized_) {
        /**
          * Initialize the state ekf_.x_ with the first measurement.
          * Create the covariance matrix.
          * Remember: you'll need to convert radar from polar to cartesian coordinates.
        */
        // first measurement
        cout << "EKF: " << endl;
        VectorXd x_initial_ = VectorXd(4);
        MatrixXd H_initial_;
        MatrixXd R_initial_;
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            x_initial_ = tools.Polar2Cartesian(measurement_pack.raw_measurements_);
            x_initial_[2]=0;
            x_initial_[3]=0;
            H_initial_ = tools.CalculateJacobian(x_initial_);
            R_initial_ = R_radar_;
        }
        else{
            x_initial_<< measurement_pack.raw_measurements_(0),measurement_pack.raw_measurements_(1),0,0;;
            H_initial_ = H_laser_;
            R_initial_ = R_laser_;
        }
        //the initial transition matrix F_
        MatrixXd F_initial_ = MatrixXd(4, 4);
        F_initial_ << 1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1;
        //state covariance matrix P
        MatrixXd P_initial_ = MatrixXd(4, 4);
        P_initial_ << 1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1000, 0,
                    0, 0, 0, 1000;
        MatrixXd Q_initial_ = MatrixXd(4, 4);
        Q_initial_ << 1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1;
        ekf_.Init(x_initial_,P_initial_,F_initial_,H_initial_,R_initial_,Q_initial_);
        cout << "Kalman Filter Initialization " << endl;

        previous_timestamp_ = measurement_pack.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        cout << "x_initial_ = " << ekf_.x_ << endl;
        cout << "P_initial_ = " << ekf_.P_ << endl;
        return;
    }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
    //compute the time elapsed between the current and previous measurements

    double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds

    previous_timestamp_ = measurement_pack.timestamp_;
    ekf_.F_(0,2)=dt;
    ekf_.F_(1,3)=dt;
    MatrixXd G = MatrixXd(4,2);
    G << pow(dt, 2)/2, 0,
        0, pow(dt, 2)/2,
        dt, 0,
        0, dt;

    MatrixXd a_v = MatrixXd(2,2);
    a_v << noise_ax, 0,
            0, noise_ay;

    ekf_.Q_ = G*a_v*(G.transpose());
    ekf_.Predict();

    /*****************************************************************************
    *  Update
    ****************************************************************************/

    /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
    */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

    // Radar updates
        //measurement matrix
        ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
        //measurement covariance
        ekf_.R_ = R_radar_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
    // Laser updates
        //measurement matrix
        ekf_.H_ = H_laser_;
        //measurement covariance
        ekf_.R_ = R_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }
    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}

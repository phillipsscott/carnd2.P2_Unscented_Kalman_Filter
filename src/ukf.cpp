#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <list>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0.);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5);
  P_.fill(0.);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.125;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  
  // Complete the initialization. See ukf.h for other member properties.

  is_initialized_ = false;
  //time_us_ = 0;
  
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.);
  weights_ = VectorXd(2 * n_aug_ + 1);
  double w = lambda_ / (lambda_ + n_aug_);
  weights_(0) = w;
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
	  weights_(i) = 0.5 / (lambda_ + n_aug_);
  }
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	// check sensor type
	if (!use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		cout << "Disregard reading from radar" << endl;
		return;
	}
	else if (!use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
		cout << "Disregard reading from laser" << endl;
		return;
	}

	// initialize

	if (!is_initialized_) {
		cout << "Initialize UKF\n";
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			double rh = meas_package.raw_measurements_(0);
			double ph = meas_package.raw_measurements_(1);

			double px = rh * cos(ph);
			double py = rh * sin(ph);
			x_ << px, py, 0., 0., 0.;
		}
		else {
			double px = meas_package.raw_measurements_(0);
			double py = meas_package.raw_measurements_(1);
			x_ << px, py, 0., 0., 0.;
		}

		is_initialized_ = true;

		time_us_ = meas_package.timestamp_;

		cout << "Initial measurement complete\n";
		return;
	}

	// predict

	cout << "Prediction Step\n";
	double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
	time_us_ = meas_package.timestamp_;

	// yield Sigma pts => output state mean/covar
	Prediction(dt);

	// update

	cout << "Update Step\n";

	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		cout << "RADAR\n";
		UpdateRadar(meas_package);
	}
	else {
		cout << "Lidar\n";
		UpdateLidar(meas_package);
	}
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  
	double delta_t2 = delta_t * delta_t;

	VectorXd x_aug = VectorXd(n_aug_);
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	x_aug.fill(0.);
	x_aug.head(n_x_) = x_;

	P_aug.fill(0.);
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug(5,5) = std_a_ * std_a_;
	P_aug(6, 6) = std_yawdd_ * std_yawdd_;

	MatrixXd L = P_aug.llt().matrixL();

	Xsig_aug.fill(0.);
	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i < n_aug_; i++) {
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
	}

	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		double px = Xsig_aug(0, 1);
		double py = Xsig_aug(1, i);
		double v = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yawd = Xsig_aug(4, i); 
		double nu_a = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);

		double cos_yaw = cos(yaw);
		double sin_yaw = sin(yaw);

		double px_pred, py_pred;

		if (fabs(yawd) > 0.0001) {
			px_pred = px + v / yawd * (sin(yaw + yawd * delta_t) - sin_yaw);
			py_pred = py + v / yawd * (cos_yaw - cos(yaw + yawd * delta_t));
		}
		else {
			px_pred = px + v * delta_t*cos_yaw;
			py_pred = py + v * delta_t*sin_yaw;
		}

		double v_pred = v;
		double yaw_pred = yaw + yawd * delta_t;
		double yawd_pred = yawd;

		px_pred = px_pred + 0.5*nu_a*delta_t2*cos_yaw;
		py_pred = py_pred + 0.5*nu_a*delta_t2*sin_yaw;
		v_pred = v_pred + nu_a * delta_t;
		yaw_pred = yaw_pred + 0.5*nu_yawdd*delta_t2;
		yawd_pred = yawd_pred + nu_yawdd * delta_t;

		Xsig_pred_(0, i) = px_pred;
		Xsig_pred_(1, i) = py_pred;
		Xsig_pred_(2, i) = v_pred;
		Xsig_pred_(3, i) = yaw_pred;
		Xsig_pred_(4, i) = yawd_pred;
	}

	x_.fill(0.);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		x_ = x_ + weights_(i) * Xsig_pred_.col(i);
	}

	P_.fill(0.);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

		P_ = P_ + weights_(i) * x_diff * x_diff * x_diff.transpose();
	}
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
	const int n_z = 2;
	VectorXd z = VectorXd(n_z);
	z(0) = meas_package.raw_measurements_(0);
	z(1) = meas_package.raw_measurements_(1);

	// measure

	MatrixXd Zsigma = MatrixXd(n_z, 2 * n_aug_ + 1);

	// sigma points => measurement space

	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		const double px = Xsig_pred_(0, i);
		const double py = Xsig_pred_(1, i);

		Zsigma(0, i) = px;
		Zsigma(1, i) = py;
	}

	VectorXd z_p = VectorXd(n_z);
	z_p.fill(0.);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_p = z_p + weights_(i) * Zsigma.col(i);
	}

	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		VectorXd diff = Zsigma.col(i) - z_p;

		S = S + weights_(i) * diff * diff.transpose();
	}

	MatrixXd R = MatrixXd(n_z, n_z);
	R << std_laspx_ * std_laspx_, 0,
		0, std_laspy_*std_laspy_;

	S = S + R;

	// update

	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		VectorXd z_diff = Zsigma.col(i) - z_p;
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	MatrixXd S_inv = S.inverse();
	MatrixXd K = Tc * S_inv;
	VectorXd z_diff = z - z_p;
	x_ = x_ + K * z_diff;
	P_ = P_ - K * S * K.transpose();	

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
	const int n_z = 3;
	VectorXd z = VectorXd(n_z);
	z(0) = meas_package.raw_measurements_(0);
	z(1) = meas_package.raw_measurements_(1);
	z(2) = meas_package.raw_measurements_(2);

	// measure
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		const double px = Xsig_pred_(0, i);
		const double py = Xsig_pred_(1, i);
		const double v = Xsig_pred_(2, i);
		const double yaw = Xsig_pred_(3, i);
		const double cos_v = cos(yaw)*v;
		const double sin_v = sin(yaw)*v;
		const double sqrt_p = sqrt(px*px + py * py);

		Zsig(0, i) = sqrt_p;
		Zsig(1, i) = atan2(py, px);
		Zsig(2, i) = px * cos_v + py * sin_v / sqrt_p;
	}

	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		VectorXd diff = Zsig.col(i) - z_pred;

		while (diff(1) > M_PI) diff(1) -= 2.*M_PI;
		while (diff(1) < -M_PI) diff(1) += 2.*M_PI;
		S = S + weights_(i) * diff * diff.transpose();
	}

	MatrixXd R = MatrixXd(n_z, n_z);
	R << std_radr_ * std_radr_, 0, 0,
		0, std_radphi_* std_radphi_, 0,
		0, 0, std_radrd_*std_radr_;
	S = S + R;

	// update

	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		VectorXd z_diff = Zsig.col(i) - z_pred;

		while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		while (x_diff(1) > M_PI) x_diff(1) -= 2.*M_PI;
		while (x_diff(1) < -M_PI) x_diff(1) += 2.*M_PI;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	MatrixXd S_inv = S.inverse();
	MatrixXd K = Tc * S_inv;

	VectorXd z_diff = z - z_pred;

	while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
	while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

	x_ = x_ + K * z_diff;
	P_ = P_ - K * S * K.transpose();

	//NIS_radar_ = z_diff.transpose() * S.inv * z_diff;
}

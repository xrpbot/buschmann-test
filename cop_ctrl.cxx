#include <iostream>
#include <cmath>
#include "Hermite.h"
#include <Eigen/Dense>

// Solve the equation:
// m z y_ddot - m y g = -T_x(t)
// where z is constant and T_x(t) is piecewise linear.

// The center-of-pressure p_x is defined by:
// p_y = T_x/F_z = T_x/(m*g)

// *** System parameters ***
const double sys_m = 1.0;
const double sys_g = 9.81;
const double sys_z = 1.0;

const unsigned int N_INTERVALS = 6;

// Time interval
const double h = 0.5;
const double t_end = h * N_INTERVALS;

// Reference torques T_x
const double T_x_ref[N_INTERVALS+1] = { 0.0, 0.1, 1.0, 1.1, 2.0, 2.1, 3.0};

double calc_T_x(double t)
{
    int i = floor(t / h);
    if(i < 0) i = 0;
    if(i > N_INTERVALS-1) i = N_INTERVALS-1;
    
    const double t_i = i*h;
    const double t_i1 = (i+1)*h;
    
    return ((t-t_i)*T_x_ref[i+1] + (t_i1-t)*T_x_ref[i])/h;
}

// Approximate solution
Hermite bvp_solver()
{
    const size_t n_pts = 31 + 2;
    const double dx = 0.1;
    
    Eigen::MatrixXd A(n_pts, n_pts); A.setZero();
    Eigen::MatrixXd R(n_pts, n_pts); R.setZero();
    
    for(size_t i=1; i<(n_pts-1); i++) {
        A(i, i-1) = 1./dx;
        A(i, i) = 4./dx;
        A(i, i+1) = 1./dx;
        
        R(i, i+1) = 3./(dx*dx);
        R(i, i-1) = -3./(dx*dx);
    }
    A(0, 0) = 2./dx;
    A(0, 1) = 1./dx;
    R(0, 0) = -3./(dx*dx);
    R(0, 1) = 3./(dx*dx);
    
    A(n_pts-1, n_pts-2) = 1./dx;
    A(n_pts-1, n_pts-1) = 2./dx;
    R(n_pts-1, n_pts-2) = -3./(dx*dx);
    R(n_pts-1, n_pts-1) = 3./(dx*dx);
    
    // Matrix relating k to the control points
    Eigen::MatrixXd K = A.colPivHouseholderQr().solve(R);
    
    // Matrix relating the second derivative of y to the control points
    Eigen::MatrixXd B(n_pts, n_pts);
    for(size_t i=0; i<(n_pts-1); i++) {
        for(size_t j=0; j<n_pts; j++) {
            B(i, j) = -4./dx * K(i, j) - 2./dx * K(i+1, j);
        }
        
        B(i, i) += -6./(dx*dx);
        B(i, i+1) += 6./(dx*dx);
    }
    
    // Enforce the ODE at the inner control points
    Eigen::MatrixXd D(n_pts, n_pts); D.setZero();
    Eigen::VectorXd Y(n_pts); Y.setZero();
    
    for(size_t i=1; i<(n_pts-1); i++) {
        for(size_t j=0; j<n_pts; j++) {
            D(i,j) = sys_m*sys_z*B(i, j);
        }
        D(i, i) -= sys_m*sys_g;
        Y(i) = -calc_T_x(i*dx - dx);
    }
    
    // Enforce the boundary points (y = y_cop)
    const double y_0 = 0.;
    const double y_end = calc_T_x(t_end)/(sys_m*sys_g);
    D(0, 1) = 1.;
    D(n_pts-1, n_pts-2) = 1.;
    
    Y(0) = y_0;
    Y(n_pts-1) = y_end;
    
    Eigen::VectorXd y(n_pts);
    y = D.colPivHouseholderQr().solve(Y);
    
    Hermite h;
    h.setData(n_pts, -dx, dx, y.data(), (Eigen::Matrix<double, n_pts, 1>(K*y)).data());
    return h;
}

// Analytical solution for piecewise-linear T_x
// The analytical solution is given by:
// y_i = a_i*exp(sqrt(alpha)*t) + b_i*exp(-sqrt(alpha)*t) + T_x(t)
// (typo in thesis?)
Eigen::VectorXd pwl_solver()
{
    const double alpha = sys_g/sys_z;
    const double sqrt_alpha = sqrt(alpha);
    
    // FIXME: using a dense matrix here is very inefficient.
    Eigen::MatrixXd Q(2*N_INTERVALS, 2*N_INTERVALS);
    Q.setZero();
    Eigen::VectorXd r(2*N_INTERVALS);
    r.setZero();
    
    // Left boundary condition
    // y_0 = a_0 + b_0 + T_x_0/(m*g)
    const double y_0 = 0.;
    Q(0, 0) = 1.;
    Q(0, 1) = 1.;
    
    r(0) = y_0 - T_x_ref[0]/(sys_m*sys_g);
    
    // Internal smoothness conditions
    for(unsigned int i=0; i<(N_INTERVALS-1); i++) {
        const double t_i1 = h*(i+1);
        
        // Smooth value
        // a_i*exp(sqrt(alpha)*t_{i+1}) + b_i*exp(-sqrt(alpha)*t_{i+1})
        //     - a_{i+1}*exp(sqrt(alpha)*t_{i+1}) - b_{i+1}*exp(-sqrt(alpha)*t_{i+1}) = 0
        Q(2*i+1, 2*i  ) = exp(sqrt_alpha*t_i1);
        Q(2*i+1, 2*i+1) = exp(-sqrt_alpha*t_i1);
        Q(2*i+1, 2*i+2) = -exp(sqrt_alpha*t_i1);
        Q(2*i+1, 2*i+3) = -exp(-sqrt_alpha*t_i1);
        
        // Smooth first derivative
        // a_i*exp(sqrt(alpha)*t_{i+1}) - b_i*exp(-sqrt(alpha)*t_{i+1}) + (dT_i/dt)/(m*g)
        //     - a_{i+1}*exp(sqrt(alpha)*t_{i+1}) + b_{i+1}*exp(-sqrt(alpha)*t_{i+1}) - (dT_{i+1}/dt)/(m*g) = 0
        Q(2*i+2, 2*i  ) = exp(sqrt_alpha*t_i1);
        Q(2*i+2, 2*i+1) = -exp(-sqrt_alpha*t_i1);
        Q(2*i+2, 2*i+2) = -exp(sqrt_alpha*t_i1);
        Q(2*i+2, 2*i+3) = exp(-sqrt_alpha*t_i1);
        
        r(2*i+2) = -(T_x_ref[i+1]-T_x_ref[i])/(h*sys_m*sys_g*sqrt_alpha)
            + (T_x_ref[i+2]-T_x_ref[i+1])/(h*sys_m*sys_g*sqrt_alpha);
    }
    
    // Right boundary condition
    // y_N = a_{N-1} exp(sqrt(alpha)*t_N) + b_{N-1}*exp(-sqrt(alpha)*t_N) + T_N/(m*g)
    const double y_end = T_x_ref[N_INTERVALS]/(sys_m*sys_g);
    Q(2*N_INTERVALS-1, 2*N_INTERVALS-2) = exp(sqrt_alpha*t_end);
    Q(2*N_INTERVALS-1, 2*N_INTERVALS-1) = exp(-sqrt_alpha*t_end);
    
    r(2*N_INTERVALS-1) = y_end - T_x_ref[N_INTERVALS]/(sys_m*sys_g);
    
    // Solve equations
    Eigen::VectorXd params(2*N_INTERVALS);
    params = Q.colPivHouseholderQr().solve(r);
    return params;
}

double eval_y_cog_ana(double t, const Eigen::VectorXd& params)
{
    int i = floor(t / h);
    if(i < 0) i = 0;
    if(i > N_INTERVALS-1) i = N_INTERVALS-1;
    
    const double alpha = sys_g/sys_z;
    const double sqrt_alpha = sqrt(alpha);
    
    const double a_i = params[2*i];
    const double b_i = params[2*i+1];
    const double T_x = calc_T_x(t);
    
    return a_i*exp(sqrt_alpha*t) + b_i*exp(-sqrt_alpha*t) + T_x/(sys_m*sys_g);
}

int main()
{
    // Analytical solution
    Eigen::VectorXd params = pwl_solver();
    
    // Approximate solution
    Hermite y_cog_approx = bvp_solver();
    
    std::cout << "#:1:y_cop" << std::endl;
    std::cout << "#:2:y_cog_ana" << std::endl;
    std::cout << "#:3:y_cog_approx" << std::endl;
    std::cout << "#:4:y_cop_approx" << std::endl;
    
    for(unsigned int step=0; step<=300; step++) {
        const double t = step/100.;
        const double T_x = calc_T_x(t);
        
        const double y_cop = T_x/(sys_m*sys_g);
        const double y_cop_approx = -(sys_m*sys_z*y_cog_approx.eval_d2(t) - sys_m*sys_g*y_cog_approx.eval(t))/(sys_m * sys_g);
        
        std::cout << t << " ";
        std::cout << y_cop << " ";
        std::cout << eval_y_cog_ana(t, params) << " ";
        std::cout << y_cog_approx.eval(t) << " ";
        std::cout << y_cop_approx << " ";
        
        std::cout << std::endl;
    }
    
    return 0;
}

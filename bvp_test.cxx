#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "Hermite.h"

Hermite bvp_solver();

int main()
{
    // BVP parameters
    const double omega = 3.0;
    const double t0 = 0.;
    const double t1 = 1.5;
    const double y0 = 0.7;
    const double y1 = -0.4;
    
    // Exact solution
    const double c0 = cos(omega*t0);
    const double s0 = sin(omega*t0);
    const double c1 = cos(omega*t1);
    const double s1 = sin(omega*t1);
    
    const double a = (c1*y0 - c0*y1)/(s0*c1 - s1*c0);
    const double b = (s1*y0 - s0*y1)/(s1*c0 - s0*c1);
    
    // Approximate solution
    Hermite h = bvp_solver();
    
    for(unsigned int t_idx=0; t_idx<=150; t_idx++) {
        const double t = t_idx * 0.01;
        
        std::cout << t << " ";
        std::cout << a*sin(omega*t) + b*cos(omega*t) << " ";
        std::cout << -a*omega*omega*sin(omega*t) - b*omega*omega*cos(omega*t) << " ";
        std::cout << h.eval(t) << " ";
        std::cout << h.eval_d2(t) << " ";
        std::cout << std::endl;
    }
    
    return 0;
}

Hermite bvp_solver()
{
    const size_t n_pts = 16 + 2;
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
    
    const double omega = 3.0;
    for(size_t i=1; i<(n_pts-1); i++) {
        for(size_t j=0; j<n_pts; j++) {
            D(i,j) = B(i, j);
        }
        D(i, i) += omega*omega;
    }
    
    // Enforce the boundary points
    const double y0 = 0.7;
    const double y1 = -0.4;
    D(0, 1) = 1.;
    D(n_pts-1, n_pts-2) = 1.;
    
    Eigen::VectorXd Y(n_pts); Y.setZero();
    Y(0) = y0;
    Y(n_pts-1) = y1;
    
    Eigen::VectorXd y(n_pts);
    y = D.colPivHouseholderQr().solve(Y);
    
    Hermite h;
    h.setData(n_pts, -dx, dx, y.data(), (Eigen::Matrix<double, n_pts, 1>(K*y)).data());
    return h;
}

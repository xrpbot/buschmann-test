#include <iostream>
#include <cmath>
#include <Eigen/Dense>

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

int main()
{
    const double m = 1.0;
    const double g = 9.81;
    const double z = 1.0;
    
    const double alpha = g/z;
    const double sqrt_alpha = sqrt(alpha);
    
    // FIXME: using a dense matrix here is very inefficient.
    Eigen::MatrixXd Q(2*N_INTERVALS, 2*N_INTERVALS);
    Q.setZero();
    Eigen::VectorXd r(2*N_INTERVALS);
    r.setZero();
    
    // Left boundary condition
    // y_0 = a_0 + b_0 - T_x_0/(m*g)
    const double y_0 = 0.;
    Q(0, 0) = 1.;
    Q(0, 1) = 1.;
    
    r(0) = y_0 + T_x_ref[0]/(m*g);
    
    // Internal smoothness conditions
    for(unsigned int i=0; i<N_INTERVALS-1; i++) {
        const double t_i1 = h*(i+1);
        
        // a_i*exp(sqrt(alpha)*t_{i+1}) + b_i*exp(-sqrt(alpha)*t_{i+1})
        //     - a_{i+1}*exp(sqrt(alpha)*t_{i+1}) - b_{i+1}*exp(-sqrt(alpha)*t_{i+1}) = 0
        Q(2*i+1, 2*i  ) = exp(sqrt_alpha*t_i1);
        Q(2*i+1, 2*i+1) = exp(-sqrt_alpha*t_i1);
        Q(2*i+1, 2*i+2) = -exp(sqrt_alpha*t_i1);
        Q(2*i+1, 2*i+3) = -exp(-sqrt_alpha*t_i1);
        
        // a_i*exp(sqrt(alpha)*t_{i+1}) - b_i*exp(-sqrt(alpha)*t_{i+1})
        //     - a_{i+1}*exp(sqrt(alpha)*t_{i+1}) + b_{i+1}*exp(-sqrt(alpha)*t_{i+1}) = 0
        Q(2*i+2, 2*i  ) = exp(sqrt_alpha*t_i1);
        Q(2*i+2, 2*i+1) = -exp(-sqrt_alpha*t_i1);
        Q(2*i+2, 2*i+2) = -exp(sqrt_alpha*t_i1);
        Q(2*i+2, 2*i+3) = exp(-sqrt_alpha*t_i1);
        
        r(2*i+2) = (T_x_ref[i+1]-T_x_ref[i])/(h*m*g*sqrt_alpha)
            - (T_x_ref[i+2]-T_x_ref[i+1])/(h*m*g*sqrt_alpha);
    }
    
    // Right boundary condition
    // y_N = a_{N-1} exp(sqrt(alpha)*t_N) + b_{N-1}*exp(-sqrt(alpha)*t_N) - T_N/(m*g)
    const double y_end = -T_x_ref[N_INTERVALS]/(m*g);
    Q(2*N_INTERVALS-1, 2*N_INTERVALS-2) = exp(sqrt_alpha*t_end);
    Q(2*N_INTERVALS-1, 2*N_INTERVALS-1) = exp(-sqrt_alpha*t_end);
    
    r(2*N_INTERVALS-1) = y_end + T_x_ref[N_INTERVALS]/(m*g);
    
    // Solve equations
    Eigen::VectorXd params(2*N_INTERVALS);
    
    params = Q.colPivHouseholderQr().solve(r);
    
    std::cout << "#:1:y_cop" << std::endl;
    std::cout << "#:2:y_cog" << std::endl;
    
    for(unsigned int step=0; step<=300; step++) {
        const double t = step/100.;
        
        int i = floor(t / h);
        if(i < 0) i = 0;
        if(i > N_INTERVALS-1) i = N_INTERVALS-1;
        
        const double a_i = params[2*i];
        const double b_i = params[2*i+1];
        
        const double T_x = calc_T_x(t);
        const double y = a_i*exp(sqrt_alpha*t) + b_i*exp(-sqrt_alpha*t) - T_x/(m*g);
        
        std::cout << t << " " << T_x/(m*g) << " " << -y;
        
        std::cout << std::endl;
    }
    
    return 0;
}

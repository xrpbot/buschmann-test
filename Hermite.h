#ifndef CW_HERMITE_H
#define CW_HERMITE_H

#include <vector>

class Hermite
{
  public:
    Hermite();
    void setData(unsigned int n, double x0, double dx, const double* y, const double *k);
    
    double eval(double x);
    double eval_d(double x);
    double eval_d2(double x);
    
    unsigned int index(double x);
    
  private:
    unsigned int fn;   // number of *points*
    double fx0, fdx;
    
    std::vector<double> fy;   // value at the interpolation points
    std::vector<double> fk;   // (scaled) derivative at the interpolation points
};

#endif

#include "Hermite.h"
#include <cmath>

Hermite::Hermite()
    : fn(0) {}

void Hermite::setData(unsigned int n, double x0, double dx, const double* y, const double *k)
{
    fn = n;
    fx0 = x0;
    fdx = dx;
    
    if(fn == 1) fn = 0;
    
    fy.resize(fn);
    fk.resize(fn);
    
    for(unsigned int i=0; i<fn; i++) {
        fy[i] = y[i];
        fk[i] = fdx*k[i];
    }
}

unsigned int Hermite::index(double x)
{
    if(fn == 0) return 0;
    
    int i = floor((x-fx0)/fdx);
    if(i < 0) i = 0;
    if(i > fn-2) i = fn-2;
    
    return i;
}

double Hermite::eval(double x)
{
    if(fn == 0) return 0.;
    
    int i = index(x);
    const double t = (x-fx0)/fdx - i;
    
    const double y1 = fy[i];
    const double y2 = fy[i+1];
    
    const double a = fk[i] - (y2 - y1);
    const double b = (y2 - y1) - fk[i+1];
    
    return (1.-t)*y1 + t*y2 + t*(1.-t)*(a*(1.-t) + b*t);
}

double Hermite::eval_d(double x)
{
    if(fn == 0) return 0.;
    
    int i = index(x);
    const double t = (x-fx0)/fdx - i;
    
    const double y1 = fy[i];
    const double y2 = fy[i+1];
    
    const double a = fk[i] - (y2 - y1);
    const double b = (y2 - y1) - fk[i+1];
    
    return (-y1 + y2 + a - 4.*a*t + 2.*b*t + 3.*a*t*t - 3.*b*t*t)/fdx;
}

double Hermite::eval_d2(double x)
{
    if(fn == 0) return 0.;
    
    int i = index(x);
    const double t = (x-fx0)/fdx - i;
    
    const double y1 = fy[i];
    const double y2 = fy[i+1];
    
    const double a = fk[i] - (y2 - y1);
    const double b = (y2 - y1) - fk[i+1];
    
    return (-4*a + 2.*b + 6.*a*t - 6.*b*t)/(fdx*fdx);
}


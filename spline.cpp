/*
Abraham Flores
Notre Dame Physics REU 2016
6/9/2016
Language C++


SPLINE INTERPOLATION:
*see header for details
*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <utility>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_laguerre.h>
#include "spline.h"
/*
Enter function required for evaluation here:
*/
namespace spline{
    double CubicIntegrate(double *x,double *y, int n){

        gsl_interp_accel *acc = gsl_interp_accel_alloc (); //Returns pointer to accelerator object, tracks state of lookups
        const gsl_interp_type *t = gsl_interp_cspline;//Cubic Spline with periodic boundary conditions, result: piecewise cubic on each interval
        gsl_spline *spline = gsl_spline_alloc (t,n);//Returns a pointer to a newly allocated interpolation object of type t for n size data-points

        gsl_spline_init(spline, x, y, n);

        double integral_value;
        integral_value = gsl_spline_eval_integ(spline,x[0],x[n-1],acc);//Evaluate Integral

        gsl_spline_free (spline);
        gsl_interp_accel_free (acc);

        return integral_value;
    }

        double AkimaIntegrate(double *x,double *y, int n){

        gsl_interp_accel *acc = gsl_interp_accel_alloc (); //Returns pointer to accelerator object, tracks state of lookups
        const gsl_interp_type *t = gsl_interp_akima;//Cubic Spline with periodic boundary conditions, result: piecewise cubic on each interval
        gsl_spline *spline = gsl_spline_alloc (t,n);//Returns a pointer to a newly allocated interpolation object of type t for n size data-points

        gsl_spline_init(spline, x, y, n);

        double integral_value;
        integral_value = gsl_spline_eval_integ(spline,x[0],x[n-1],acc);//Evaluate Integral

        gsl_spline_free (spline);
        gsl_interp_accel_free (acc);

        return integral_value;
    }

        double SteffenIntegrate(double *x,double *y, int n){

        gsl_interp_accel *acc = gsl_interp_accel_alloc (); //Returns pointer to accelerator object, tracks state of lookups
        const gsl_interp_type *t = gsl_interp_steffen;//Cubic Spline with periodic boundary conditions, result: piecewise cubic on each interval
        gsl_spline *spline = gsl_spline_alloc (t,n);//Returns a pointer to a newly allocated interpolation object of type t for n size data-points

        gsl_spline_init(spline, x, y, n);

        double integral_value;
        integral_value = gsl_spline_eval_integ(spline,x[0],x[n-1],acc);//Evaluate Integral

        gsl_spline_free (spline);
        gsl_interp_accel_free (acc);

        return integral_value;
    }
}

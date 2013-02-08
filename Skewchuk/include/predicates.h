#ifndef __predicates_h__
#define __predicates_h__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double orient2d( double*, double*, double* );
double orient2dfast( double*, double*, double* );

double incircle( double*, double*, double*, double* );
double incirclefast( double*, double*, double*, double* );

#endif //__predicates_h__
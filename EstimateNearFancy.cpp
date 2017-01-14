#include"EstimateNear.h"

#include<math.h>



double xm=-1;

double ym=-1;

double sr2=1.41421356;



int EstimateNear(double *x, double*y, double xo, double yo, double w1, double w2)

{

	if (1)

	{

		double a = w2/sr2-w1/w2*w1/sr2;

		xm = (-a*w1/w2+ sqrt( (a*w1/w2)*(a*w1/w2)-4*(w1/w2*w1/w2+w2/w1*w2/w1)*(a*a-w2*w2) ))

			/ 2 / (w1/w2*w1/w2+w2/w1*w2/w1);

		ym = w2/sr2+w1/w2*(xm-w1/sr2);

	}



	double xi,yi;

	yi = (yo-xo*w1/w2+w2*(w1/w2*w1/w2))/(1+w1/w2*w1/w2);

	xi = (w2-yi)*w1/w2;



	if (yi<0)

	{

		*x = w1;

		*y = w2;

		return(0);

	}

	if (yi < w2/sr2)

	{

		*y = ( ym*w1/(w1-xm)-ym/(w1-xm)*(xo-yo*w2/w1) )/ (1+ym/(w1-xm)*w2/w1);

		*x = xo + (*y-yo)*w2/w1;

		return(0);

	}

	if (xi<xm)

	{

		*x = xm;

		*y = ym;

		return(0);

	}

	// Must be the last case:

	*x = xo;

	*y = w2-xo/xm*(ym-w2);

	return(0);

}
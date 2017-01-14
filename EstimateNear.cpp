#include"EstimateNear.h"

#include<math.h>



int EstimateNear(double *x, double*y, double xo, double yo, double w1, double w2)

{

	*x = w1;

	*y = 0;

	if (yo> .00001) *y=.00001;

	return(0);

}
#include<math.h> // For sqrt() and abs()



/*  Function FindNear()

 *

 *  This function takes x,y as input for the initial guess of the closest point on

 *    the ellipse to the point of interest, xo,yo.

 *  Returns x,y which will then hold the true closest point to xo,yo on the ellipse.

 *

 *   Note briefly that the solver works by finding the intersection of a line

 *  from the center of curvature for the old x,y to the point of interest, xo,yo

 *  with the ellipse (or a tangent to it).  This is fast, and it can be proved not

 *   to EVER go unstable.  For ANY ellipse.  The proof is at the bottom of this file.

 *

 *   Note that as xo,yo gets far away from the ellipse compared to say the magnitude

 *      of w2, the algorithm slows down, as discussed at the foot of the stability proof.

 *

 *  This function depends on the requirement that w1>=1pixel and w2>=1pixel.

 */

int scount=0;



int FindNear(double *x, double *y, //Pass x,y by reference

			  double xo, double yo,//xo,yo is the point whose distance from ellipse is needed.

 			  double w1, double w2) // w1, w2 are the major and minor axis ordinates, rspv.

{

/*  x,y starting values.

 */

	double xold, yold;			



/* Variables for computing the center of curvature

 *	for a point on the ellipse.	

 *		 The curvature is calculated by finding 2 points near to each other on

 *		the ellipse, then intersecting the lines perpendicular to the ellipse 

 *		that pass through them.

 */

	double epsilon;				// Determines how close the two points used are.

	double xoffset, yoffset;	// (xoffset, yoffset) and (x,y) are the two points used.

	double slope1, slope2;		// Used for calculating the intersection of the

								//   perpendicular lines.

	double xc, yc;				// The resulting center of curvature.



/* Variables for the logic of iterative solver.

 */

	int iterate_forceflag;		

	int slope_vertical_flag;	// Set if the line from xo,yo to xc,yc is vertical; used to

								//   avoid divide by zero but still deal with the subcase.

	int exactx_flag;			// Set if need to find exact intersection with ellipse of

								//   the line from xo,yo to xc,yc.



/* Variables used for calculating new x,y

 */

	// Variables for moving tangentially along the ellipse

	double slope_cut;			// Slope of line xo,yo through xc,yc.

	double slope_tangent;		// Slope of line tangent to x,y.

	double a,b;					// For solving quadratic equation

	double hold1;				// For solving quadratic equation

	double x1,x2;				// For solving quadratic equation

	double y1,y2;				// For solving quadratic equation



	// Variable for moving radially from the ellipse

	double ratio;				// Ratio of actual to required magnitude of (x,y);

								//  it is required that x,y be on the ellipse,

								//   ie., that it not wander radially.



	xold = *x; yold = *y;



	iterate_forceflag=1;		// Force entry into while loop



	

	if ( (w1/w2)>=600.0 )			// The ellipse is too eccentric, thus hard to calculate,

	{							//  so treat it as a line segment from (0,0) to (w1,0).

		if (xo < w1)

		{

			*x = xo;

			*y = 0.0;

		}

		else

		{

			*x = w1;

			*y = 0.0;

		}

	}

	else

	if (yo < 0.000001)			// Since xo,yo is on x axis, so is x,y.

	{

		*x=w1;

		*y=0;

	}

	else						// Iteratively solve for nearest x,y on ellipse to xo,yo.

	while ( iterate_forceflag

	||	((xold-(*x))*(xold-(*x))+(yold-(*y))*(yold-(*y))>.000001)//more than <>th of a pixel of movement last time.

	// Note: For high eccentricity ellipses with big border widths, you basically

	// get  Width as the radius of a circle and tiny tiny pixel movements have huge effects;

	//  I experimented to come up with .000001 and have now set the maximum w1/w2=600,

	// which is enforced above by treating such ellipses as a line segment from (0,0) to (w1,0).

	||  ((*x)/w1*(*x)/w1+(*y)/w2*(*y)/w2> 1.00001) || ((*x)/w1*(*x)/w1+(*y)/w2*(*y)/w2< 0.99999)

		)

	{

		iterate_forceflag=0;

		xold=*x; yold=*y;



/* Move tangent to ellipse

 *

 */

		//Obtain center of curvature for xold,yold to cut the ellipse for new x,y.

		//A circle centered at xc,yc would share the same curvature at xold,yold as the

		//ellipse.

		//

		//The line xo,yo to xc,yc will cut the ellipse near to the ideal x,y sought.

		epsilon = 0.000001;

		if (xold<0.001)

		{

			xc = xold;			// put the center in place to wedge x,y toward (0,w2).

			yc = w2 - w1*w1/w2;

		}

		else

		{

			if (xold<epsilon) epsilon = -epsilon;



			xoffset = xold - w1/w2*yold*epsilon;

			yoffset = yold + w2/w1*xold*epsilon;

			if (xoffset<0.001)

			{

				xc = xoffset;  	// put the center in place to wedge x,y toward (0,w2).

				yc = w2-w1*w1/w2;

			}

			else

			{

				slope1=(w1/w2*yold/(w2/w1*xold));  // for a line off of (x,y)

				slope2=yoffset/xoffset*w1/w2*w1/w2;  // for a line through (xoffset, yoffset)

				xc = (yoffset-yold+slope1*xold-slope2*xoffset)/(slope1-slope2);

				yc = yold+slope1*(xc-xold);

			}

		}



		// Determine if the line from xo,yo to xc,yc is vertical;

		//   set slope_cut as slope.

		if (abs(xo-xc)<0.000001)

		{

			slope_vertical_flag=1;

		}

		else

		{

			slope_vertical_flag=0;

			slope_cut=(yo-yc)/(xo-xc);

		}



		// Determine whether linear or exact intersection with the ellipse

		//	will be required, then move x,y to the intersection of the (ellipse)

		//  with the line segment from xc,yc to xo,yo.

		exactx_flag=0;

		if (( (*y)<0.001)&&(xo<w1))   // if tangent to ellipse is vertical line

			exactx_flag=1;		  //  and doesn't intersect xc,yc to xo,yo.

		if (exactx_flag==0)

			if (-w1/w2*yold*(-(yo-yold))+w2/w1*xold*(xo-xold)<0) // if tangent to 

											// ellipse is not vertical

				exactx_flag=1;				// and doesn't intersect xc,yc to xo,yo.

									// Note the test used shows that the

		    //vector (-w1*yold/w2, w2*xold/w1) is clockwise of the vector

		    //       (xo-xold, yo-yold).  Thus the tangent passes "above" xo,yo.

		if (exactx_flag)

		{	// Since tangent to ellipse is too coarse, find new x,y by:

			//  intersect ellipse exactly with line through xo,yo and xc,yc.

			if (!slope_vertical_flag)

			{

				a=w2*w2+slope_cut*slope_cut*w1*w1;

				b=2*(yc-xc*slope_cut)*slope_cut*w1*w1;

				hold1=sqrt(b*b-4*a*

					(

						(yc-xc*slope_cut)*(yc-xc*slope_cut)*w1*w1

						- w1*w1*w2*w2

					)

						);

				x1= (-b+hold1)/2/a;

				x2= (-b-hold1)/2/a;

				y1= yc+(x1-xc)*slope_cut;

				if (y1>0)

				{

					*x=x1;

					*y=y1;

				}

				else

				{

					y2 = yc+(x2-xc)*slope_cut;

					*y=y2;

					*x=x2;

				}

			}

			else // xo,yo through xc,yc runs vertical, so just grab ellipse at xo.

			{

				*x=(xo+xc)/2;  // Just take the average; it's too small to matter

				*y=w2*sqrt(1-xo/w1*xo/w1);

			}

		}

		else

		{	// Since close enough to ideal x,y find the next iteration's x,y by:

			//  intersect tangent to ellipse from x,y with line through xo,yo and xc,yc.

			slope_tangent=-xold/yold*w2/w1*w2/w1;

			if (slope_vertical_flag)

			{

				//vertical line xo,yo - xc,yc case:

				*x = (xo+xc)/2;  // Just take the average; it's too small to matter.

				*y = yold+((*x)-xold)*slope_tangent;

				//maybe x belongs better as something else weighted average maybe but

				// error seems small 1e-3 so why bother.

			}

			else

			{

				// normal slope case:

				*x = (yc-yold+xold*slope_tangent-xc*slope_cut)/(slope_tangent-slope_cut);

				*y = yc + ((*x)-xc)*slope_cut;

			}

		}

	

		

	/*

	 * Move radial to ellipse

	 */

		ratio = 1/(1+0.5*((*x)/w1*(*x)/w1+(*y)/w2*(*y)/w2-1));  // Basically this estimates

		//       (1 + (r^2 -1))^0.5 as a first order MacLauren polynomial (Taylor)

		//       and uses its reciprocal to recenter x,y on the ellipse.

		//       That is r = sqrt((x/w1)^2+(y/w2)^2) is the radius if you project

		//			the ellipse onto a unit circle.

		(*x)*=ratio;

		(*y)*=ratio;



	}  // A sufficiently accurate x,y has been found.



	/*  Note about the Taylor approximation: some tested limits for using this are:

	This is a basic program:



	'Using a first order Taylor approximation to sqrt(radius^2)

	' to divide a vector (The guess x,y) will not "ring",or go unstable,

	' for radius between 1 and 15 times the circle radius.

	CLS

	FOR i = 1.00001 TO 15 STEP .0001

	ne = i * 1 / (1 + .5 * (i ^ 2 - 1))

	min = 1 / i

	PRINT i, ne, min

	IF ne <= min THEN PRINT ne, min

	NEXT i

	*/





	return(0);

}





/*  Proof of the stability of this algorithm for finding the closest point

 *  on the ellipse to xo,yo:

 *

 *  The proof relies on a geometric construction.

 *

 *  Draw an ellipse, centered on the origin, with major axis (w1, 0), minor axis (0,w2).

 *  Consider only the portion of the ellipse x>0, y>0;

 *   here create a graph of the centers of curvature for all the points of the ellipse.

 *   Note that the center of curvature for the point (w1,0) is w2*w2/w1,

 *      and the center of curvature for the point (0,w2) is w1*w1/w2.

 *     By eyeballing it for the intermediate points you can figure out where the

 *      center of curvature is for any given point on the ellipse.

 *   Now my proof is this: My algorithm finds new points x,y by moving only

 *    counterclockwise along the ellipse, and it never overshoots the ideal x,y position.

 *    Therefore, since it always gets a good deal closer, but never overshoots, it's stable.

 *   This proof relies on this theorem: Any line segment (x,y) to its center of curvature

 *       (xc,yc) will never intersect with a line generated by finding any point (x2,y2)

 *       counterclockwise on the ellipse of (x,y) through which a line is drawn to the

 *		 center of curvature of (x2,y2).

 *		

 *		 Consider the movement from (x,y)'s center (xc,yc) of the locus of points

 *         which are centers of curvature.  The new xc,yc if the ellipse were a circle

 *         would not move.  However, it is elliptical the center will move in a 2 step:

 *          because moving counterclockwise along the ellipse opens up the center of

 *          curvature the new xc,yc will move in a second quadrant direction proportional

 *          to the vector (x2-x, y2-y); it will also move radially as the center of curvature

 *          becomes a greater distance from the point on the ellipse.  Therefore it can

 *          be seen that the old line segment (x,y)-(xc,yc) will not intersect with the

 *          line through (x2,y2) and its center of curvature.  The theorem is proved:

 *			A line through a point and its center of curvature counterclockwise to a given

 *				point will not intersect with the segment through that point and its center

 *				of curvature.

 *   The proof requires it to be shown that the new (x,y) will never overshoot

 *       the ideal (x,y).  Now draw an xo,yo and its corresponding ideal (x,y).

 *		 Now add its center of curvature, which is of course colinear with these 2 points.

 *       Draw the line through these three points.  

 *		 Note that any point clockwise of the ideal (x,y) will have a center of curvature

 *         which is clockwise of the ideal (x,y)'s center of curvature, and closer

 *         in to the ellipse than it is.  Considering the theorem it can be noted

 *         that the center of curvature of the (working) (x,y) will be clockwise

 *		   of the line through the ideal (x,y) and xo,yo.  This follows because otherwise

 *         that line would intersect with the line from (x,y) through its center of 

 *		   curvature.  Therefore a line drawn from (x,y)'s center of curvature through

 *         xo,yo will pass through the ellipse clockwise of the ideal (x,y).  Thus the

 *         intersection of (xc,yc)-(xo,yo) with the ellipse will not overshoot the

 *         ideal (x,y).  Therefore it cannot go unstable by overshoot.

 *			   Also note that in the case of the circle, there's no overshoot.  Draw it.

 *   The proof of convergence also requires it be shown that the new (x,y)

 *		 will be a good deal near

 *       to the ideal (x,y).  This can be examined by constructing it geometrically

 *       to prove it for yourself.  However, it is important to note that this

 *       algorithm slows down as xo,yo becomes far off of the ellipse.

 *       To see this, iterate through the algorithm pencil and ruler style moving

 *         xo,yo outward.  The fraction of the distance which is covered between

 *         the last (x,y) value and the ideal (x,y) value comes out to 

 *         the ratio     distance from xc,yc to x,y/ distance from xc,yc to xo,yo.

 *          So approximately 

 *              the fraction toward completion = curvature/(distance off ellipse+curvature)

 */


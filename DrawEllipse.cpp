/*

 * DrawEllipse.cpp

 * 

 * Render an elliptical frame of constant width and

 *  place the frame in a .BMP file called "frame.bmp".

 *

 * functions used:

 *    MkBList		Makes a list for guessing the nearest point on the ellipse.

 *    ListGuess     Determines for a point the nearest point in the list, 

 *                    travelling along the ellipse toward the major axis.

 *					The point must be outside the ellipse or on it.

 *    FindNear		Iterates to find the closest point on the ellipse.

 */



#include<math.h>     // For the abs function.

#include<stdio.h>	 // For the output to .bmp file

#include<WINDOWS.H>  // For the output to .bmp file



/* The geometric algorithm for iterating to find the nearest point on ellipse:

 */

#include"FindNear.h"

/* int FindNear(double *x, double *y, //Pass x,y by reference

 *				 double xo, double yo,

 *				 double w1, double w2);

 */

#include"EstimateNear.h"

/* int EstimateNear(double *x, double*y, double xo, double yo, double w1, double w2);

 */



/* Information about the elliptical frame */

#define W1 600.0					// Width of the frame

#define W2 200.0					// Height of the frame (Must be <= WIDTH of the frame)

#define WIDTH 70.0					// Thickness of the frame in pixels



/* For the image output to .bmp file: */

#define ALLIGN_DWORD(x) ((((x)+3)/4)*4)

#define IWIDTH 720

#define IHEIGHT 486



int main()

{



/* Open up the output bitmap file and enter in the header information

 */

	DWORD BytesPerLine;

	BytesPerLine= IWIDTH * 3;

	BytesPerLine= ALLIGN_DWORD(BytesPerLine);

	DWORD RasterSize = BytesPerLine * (DWORD) IHEIGHT;



	BITMAPFILEHEADER myheader;



	myheader.bfType = (UINT) (('M'<<8)|'B');

	myheader.bfSize = (DWORD) ( sizeof(BITMAPFILEHEADER)

		+ sizeof(BITMAPINFO)

		+ RasterSize );

	myheader.bfReserved1 = (WORD)0;

	myheader.bfReserved2 = (WORD)0;

	myheader.bfOffBits = (DWORD)(sizeof(BITMAPFILEHEADER)

			+sizeof(BITMAPINFO));// Follow the header directly

	// by the BITMAPINFO struct, directly followed by the image info.



	BITMAPINFOHEADER mybmih;

	mybmih.biSize=(DWORD)sizeof(BITMAPINFOHEADER);

	mybmih.biWidth=(LONG)IWIDTH;

	mybmih.biHeight=(LONG)IHEIGHT;

	mybmih.biPlanes=(WORD)1;

	mybmih.biBitCount=(WORD)24;

	mybmih.biCompression=(DWORD)BI_RGB;

	mybmih.biSizeImage=(DWORD)0;// since in BI_RGB format.

	mybmih.biXPelsPerMeter=(LONG)0;

	mybmih.biYPelsPerMeter=(LONG)0;

	mybmih.biClrUsed=(DWORD)0;

	mybmih.biClrImportant=(DWORD)0;



	BITMAPINFO mybmi;

	mybmi.bmiHeader=mybmih;



	FILE *fp;

	fp=fopen ("frame.bmp","wb");

	fwrite (&myheader, sizeof(BITMAPFILEHEADER),1,fp);



	fwrite (&mybmi, sizeof(BITMAPINFO),1,fp);

	

/* Construct the image here:

 */  

	struct pixel{

		BYTE blue;

		BYTE green;

		BYTE red;

	};

	struct pixel image;



	double x,y;				// Coordinates turned into the closest point on the ellipse.

	int holder;

	int i,j;

	BYTE val;

	double xo,yo;

	double w1, w2;

	double width;

	w1= (double) (W1);

	w2= (double) (W2);

	width = (double) WIDTH;



	for (j=0; j<IHEIGHT;j++)

	{

		for (i=0;i<IWIDTH;i++)

		{

			// bgr order.

			// Give val the grayscale 8bit output intensity.

			xo = (double)i;

			yo = (double)j;



			EstimateNear(&x, &y, xo, yo, w1, w2);

			if (xo/w1*xo/w1+yo/w2*yo/w2<1 && (w1/w2<600.0))

			{

				val = (BYTE)0;

			}

			else

			{

				if (xo/( w1*(width+w2)/w2 )*xo/( w1*(width+w2)/w2 )+yo/(width+w2)*yo/(width+w2)

			     > 1)

				{

				 val = (BYTE)0;

				}

				else

				{

					FindNear(&x, &y, //Pass x,y by reference

						xo, yo,

						w1, w2);				//Iterate to solve for the

															// nearest point on the ellipse

					if ( (xo-x)*(xo-x) + (yo-y)*(yo-y)

					  < width*width

					)

					{

						val = (BYTE)255;						//If distance from the point

					}				// to the nearest point on the ellipse (x,y) is less

									// than the width of the frame, then it's in the frame

									//   so plot it.

					else

						val = (BYTE)0;

				}

			}

			image.blue=(BYTE)val; // blue

			image.green=(BYTE)val; // green

			image.red=(BYTE)val; // red

			fwrite (&image, sizeof(pixel), 1,fp);

			// output it:

		}

		

		// now put in filler so it's a multiple of 4 bytes per row.

		holder = (4-(3*IWIDTH)%4)%4;

		while (holder>0)

		{

			putc((BYTE)0, fp);

			holder--;

		}

	}

	fclose (fp);



	return(0);

}


/*
	Source Code for Part III Essay - "Flast" Flow in Curved Pipes (ha ha)
	Author: Marshall Ward
	Date: 21 January 2002

	Uses a compact finite difference method (FDM), at the moment, to
solve the boundary layer equations for high Reynolds number fluid flow in
a curved pipe. The equations were expanded in x (the angle) and the
coefficents of the velocities are determined by solving differential
equations in y (the distance from the wall, in the boundary layer).

	Later, the differential equation solver should be a high order
compact FDM (HOC FDM), but for now I'm jsut experimenting with the idea of
numerically solving ODEs.

Program notes:
1) Note use of "<="s on order loops, as well as ord+1 in the initialization -
	this is so that "zeroth order" actually executes something, namely the
	solving of the lowest order problem.

*/

#include "pipe.h"
#include <iomanip.h>

#define WC0 1.763	// Leading order core flow (determined self-consistently)

// Boundary conditions
#define UA 0.0		// u(0) = 0
#define UB 0.0		// u(infinity) = 0
#define VA 0.0		// v(0) = 0
#define WA 0.0		// w(0) = 0
#define WB WC0		// w(infinity) = wc0

// Range of Banded Matrix
// reordering {u,v,w} may reduce it
#define KL	5		// Lower range of banded matrix
#define KU	5		// Upper range of banded matrix

///////////////////////////////////////////////////////////////////////////////

// Constructor - creates the "pipe", described by the velocities

Pipe::Pipe(int size, int ord, double y_start, double y_end)
{
	int i;	// Counter

	ord++;	// To work in computer language - but remember that zero order
			// means that ord=1

	// Initialize variables
	M = size;
	N = M - 2;
	order = ord;
	a = y_start;
	b = y_end;
	h = (b-a)/(size-1);

	// Intialize velocities
	u  = new velocity [order];
	v  = new velocity [order];
	w  = new velocity [order];
	wc = new double [order];

	for(i = 0; i < order; i++)
	{
		u[i] = new double [M];
		v[i] = new double [M];
		w[i] = new double [M];
	}

	// Set up the arrays which the CLAPACK solver will use
	// Dimensions = (2*KL + KU + 1, 3*size)
	AB = new vector[2*KL+KU+1];
	for(i = 0; i < 2*KL+KU+1; i++)
		AB[i] = new double[3*N];
	B = new double[3*N];

	// A lot of entries in AB are never declared, so it's a good idea to zero
	//	them out:
	zeroAB();

	initializeVelocities();

}	// End constructor

///////////////////////////////////////////////////////////////////////////////

// Destructor

Pipe::~Pipe(void)
{
	int i;	// Counter

	for(i = 0; i < order; i++)
	{
		delete [] u[i];
		delete [] v[i];
		delete [] w[i];
	}

	delete [] u;
	delete [] v;
	delete [] w;
	delete [] wc;

	for(i = 0; i < 2*KL+KU+1; i++)
		delete [] AB[i];
	delete [] AB;
	delete [] B;

}	// End destructor

///////////////////////////////////////////////////////////////////////////////

// initializeVelocities - determine initial guesses for the velocities
//		This is done by a linear relationship from A to B for each component
//	For the moment, only zeroth order is initialized

void Pipe::initializeVelocities(void)
{
	int i;	// Counter

//	u[0][0] = UA;
//	v[0][0] = VA;
//	w[0][0] = WA;
	wc[0] = WC0;

	for(i = 0; i < M; i++)
	{
		u[0][i] = UA + i*(UB - UA)/(b - a)*h;
		v[0][i] = VA;
		w[0][i] = WA + i*(WB - WA)/(b - a)*h;
	}
}	// End initializeVelocities

///////////////////////////////////////////////////////////////////////////////

// zeroAB - fills AB with zeros, since I don't intialize all of the values in
//	setupMatrices.
//	NB - this is actually because I don't implement the matrix properly...

void Pipe::zeroAB(void)
{
	int i, j;

	for(i = 0; i < (2*KL+KU+1); i++)
	{
		for(j = 0; j < 3*N; j++)
		{
			AB[i][j] = 0;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

// setupMatrices - Set up AB and B to pass them to DGBSV (banded matrix solver)

// AB is properly initialized!

// ***In the future I think it may be better to create a matrix
//	so that a lot of these explicit steps can be written using loops.

void Pipe::setupMatrices(void)
{
	int i, j;		// Counters
	int start, end;		// Ranges for internal loop
	double u_i, v_i, w_i;	// Values of u and v and point[i];
	double t_u, t_w;	// Evaluated derivatives of u and w (v not needed)

	// i = 0 (indexes uvw values at i+1)
	u_i = u[0][1];
	v_i = v[0][1];
	w_i = w[0][1];

	t_u = (u[0][2] - u[0][0])/(2*h);
	t_w = (w[0][2] - w[0][0])/(2*h);

	// Fragmented start anti-diagonal
	// Loop is from [KL+KU...KL] [0..KU]
	AB[KL+KU+0-0][0] = 2+2*(h*h)*u_i;
	AB[KL+KU+0-1][1] = (h*h)*t_u;
	AB[KL+KU+0-2][2] = (h*h)*w_i;
	AB[KL+KU+0-3][3] = -1+h/2*v_i;
	AB[KL+KU+0-4][4] = 0;
	AB[KL+KU+0-5][5] = 0;

	B[0] = -u[0][2] + 2*u_i - u[0][0]
		+ (h*h)*(v_i*t_u + u_i*u_i + (w_i*w_i - WC0*WC0)/2);

	AB[KL+KU+1-0][0] = -2*h;
	AB[KL+KU+1-1][1] = 0;
	AB[KL+KU+1-2][2] = 0;
	AB[KL+KU+1-3][3] = 0;
	AB[KL+KU+1-4][4] = -1;
	AB[KL+KU+1-5][5] = 0;

	B[1] = -u_i;

	AB[KL+KU+2-0][0] = 0;
	AB[KL+KU+2-1][1] = (h*h)*t_w;
	AB[KL+KU+2-2][2] = 2;
	AB[KL+KU+2-3][3] = 0;
	AB[KL+KU+2-4][4] = 0;
	AB[KL+KU+2-5][5] = -1+h/2*v_i;

	B[2] = v_i*t_w;

	for(i = 1; i < (N-1); i++)
	{
		u_i = u[0][i+1];
		v_i = v[0][i+1];
		w_i = w[0][i+1];

		t_u = (u[0][i+2] - u[0][i])/(2*h);
		t_w = (w[0][i+2] - w[0][i])/(2*h);

		// Unfragmented middle anti-diagonals
		AB[2*KL+KU-2-0][3*i-3] = -1-h/2*v_i;
		AB[2*KL+KU-2-1][3*i-2] = 0;
		AB[2*KL+KU-2-2][3*i-1] = 0;
		AB[2*KL+KU-2-3][3*i+0] = 2+2*(h*h)*u_i;
		AB[2*KL+KU-2-4][3*i+1] = (h*h)*t_u;
		AB[2*KL+KU-2-5][3*i+2] = (h*h)*w_i;
		AB[2*KL+KU-2-6][3*i+3] = -1+h/2*v_i;
		AB[2*KL+KU-2-7][3*i+4] = 0;
		AB[2*KL+KU-2-8][3*i+5] = 0;

		B[3*i] = -u[0][i] + 2*u_i - u[0][i+2]
			+ (h*h)*(v_i*t_u + u_i*u_i + (w_i*w_i - WC0*WC0)/2);

		AB[2*KL+KU-1-0][3*i-3] = 0;
		AB[2*KL+KU-1-1][3*i-2] = 1;
		AB[2*KL+KU-1-2][3*i-1] = 0;
		AB[2*KL+KU-1-3][3*i+0] = -2*h;
		AB[2*KL+KU-1-4][3*i+1] = 0;
		AB[2*KL+KU-1-5][3*i+2] = 0;
		AB[2*KL+KU-1-6][3*i+3] = 0;
		AB[2*KL+KU-1-7][3*i+4] = -1;
		AB[2*KL+KU-1-8][3*i+5] = 0;

		B[3*i+1] = -u_i;

		AB[2*KL+KU+0-0][3*i-3] = 0;
		AB[2*KL+KU+0-1][3*i-2] = 0;
		AB[2*KL+KU+0-2][3*i-1] = -1-h/2*v_i;
		AB[2*KL+KU+0-3][3*i+0] = 0;
		AB[2*KL+KU+0-4][3*i+1] = (h*h)*t_w;
		AB[2*KL+KU+0-5][3*i+2] = 2;
		AB[2*KL+KU+0-6][3*i+3] = 0;
		AB[2*KL+KU+0-7][3*i+4] = 0;
		AB[2*KL+KU+0-8][3*i+5] = -1+h/2*v_i;

		B[3*i+2] = v_i*t_w;

		// Assign B values
	}

	// i = N-1 (w/ upper bounds)
	u_i = u[0][N];	// u[N] = u[M-2] = next to last pt.
	v_i = v[0][N];
	w_i = w[0][N];

	t_u = (u[0][N+1] - u[0][N-1])/(2*h);
	t_w = (w[0][N+1] - w[0][N-1])/(2*h);

	// Fragmented end anti-diagonal
	AB[2*KL+KU-2-0][3*N-1-KL+0] = -1-h/2*v_i;
	AB[2*KL+KU-2-1][3*N-1-KL+1] = 0;
	AB[2*KL+KU-2-2][3*N-1-KL+2] = 0;
	AB[2*KL+KU-2-3][3*N-1-KL+3] = 2+2*(h*h)*u_i;
	AB[2*KL+KU-2-4][3*N-1-KL+4] = (h*h)*t_u;
	AB[2*KL+KU-2-5][3*N-1-KL+5] = (h*h)*w_i;

	// B[3(N-1)] = ??
	B[3*(N-1)] = -u[0][N-1] + 2*u_i - u[0][N+1]
		+ (h*h)*(v_i*t_u + u_i*u_i + (w_i*w_i - WC0*WC0)/2);

	AB[2*KL+KU-1-0][3*N-1-KL+0] = 0;
	AB[2*KL+KU-1-1][3*N-1-KL+1] = 1;
	AB[2*KL+KU-1-2][3*N-1-KL+2] = 0;
	AB[2*KL+KU-1-3][3*N-1-KL+3] = -2*h;
	AB[2*KL+KU-1-4][3*N-1-KL+4] = 0;
	AB[2*KL+KU-1-5][3*N-1-KL+5] = 0;

	B[3*(N-1)+1] = -u_i;

	AB[2*KL+KU-0-0][3*N-1-KL+0] = 0;
	AB[2*KL+KU-0-1][3*N-1-KL+1] = 0;
	AB[2*KL+KU-0-2][3*N-1-KL+2] = -1-h/2*v_i;
	AB[2*KL+KU-0-3][3*N-1-KL+3] = 0;
	AB[2*KL+KU-0-4][3*N-1-KL+4] = (h*h)*t_w;
	AB[2*KL+KU-0-5][3*N-1-KL+5] = 2;

	B[3*(N-1)+2] = v_i*t_w;


}	// End setupMatrices

///////////////////////////////////////////////////////////////////////////////

// solveAB - solves the matrix eqn AB*x=B for x, puts answer in B
//	(This uses a CLAPACK routine. Or, it will if everything is working)
//
// Oh GOD, this is so annoying. I hate FORTRAN, and I hate myself for not
//	knowing it

void Pipe::solveAB(void)
{
	int n, kl, ku, nrhs, ldab, *ipiv, ldb, info;
	double ab, b;

	n = 3*N;
	kl = KL;
	ku = KU;
	nrhs = 1;
	ldb = n;

	ipiv = new int [n];

	// GRRRR.. ab and b aren't in the right format

	dgbsv_(&n, &kl, &ku, &nrhs, AB, &ldab, ipiv, B, &ldb, &info);
}

///////////////////////////////////////////////////////////////////////////////

// displayAB - this is entirely for debugging, to see if I did setupM ok.

void Pipe::displayAB(void)
{
	int i, j;

	cout << endl;
	cout << "M: " << M << endl;
	cout << "N: " << N << endl;
	cout << "h: " << h << endl;

	for(i = 0; i < M; i++)
	{
		cout << "u[" << i << "]: " << setprecision(3) << u[0][i] << " ¦ ";
		cout << "v[" << i << "]: " << setprecision(3) << v[0][i] << " ¦ ";
		cout << "w[" << i << "]: " << setprecision(3) << w[0][i];
		cout << endl;
	}

	cout << endl << "Contents of AB:" << endl;

	for(i = KL; i < 2*KL+KU+1; i++)
	{
		for(j = 0; j < 3*N; j++)
		{
			cout << AB[i][j] << " \t ";
		}
		cout << endl;
	}

	cout << "Contents of B:" << endl;

	for(i = 0; i < 3*N; i++)
	{
		cout << B[i] << endl;
	}
}

///////////////////////////////////////////////////////////////////////////////

// min - returns minimum value of two numbers

int Pipe::min(int i, int j)
{
	if(i <= j) return i;
	else return j;

}	// End min

///////////////////////////////////////////////////////////////////////////////

// max - returns maximum value of two numbers

int Pipe::max(int i, int j)
{
	if(i >= j) return i;
	else return j;

}	// End max

///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

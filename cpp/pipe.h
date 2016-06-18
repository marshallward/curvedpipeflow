#ifndef PIPE_H
#define PIPE_H

typedef double * velocity;
typedef velocity * velocity_ptr;

typedef double * vector;
typedef vector * matrix;

class Pipe
{
///////////////////////////////////////////////////////////////////////////////

	public:
		Pipe(int, int, double, double);		// Creates pipe and sets up uvw
		~Pipe();							// Destructs pipe
		void initializeVelocities();		// Sets initial guesses
		void zeroAB();						// Fills AB with zeros
		void setupMatrices();				// Prepares matrices for CLAPACK
		void solveAB();						// Solves AB*x=b (in my dreams)
		void displayAB();					// DEBUG - display AB to check
		int min(int, int);					// returns minimum
		int max(int, int);					// returns maximum

	private:
		velocity_ptr u;		// radial velocity in BL, towards center positive
		velocity_ptr v;		// tangential velocity in BL, CCW positive
		velocity_ptr w;		// cross-sectional flow in BL, into page positive
		velocity wc;		// cross-sectional core flow, into page positive

		double h;			// Step size
		double a;			// Lower boundary (y=0)
		double b;			// Upper boundary (y->infinity)

		int M;				// number of nodes
		int N;				// number of internal nodes (M-2)
		int order;			// Order or calculation + 1 (for convenience)

		matrix AB;			// For CLAPACK Ax=B solver, banded storage
		vector B;			// RHS of matrix eqn Ax = B.

};	// End Pipe

#endif

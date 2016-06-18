/* runpipe.cpp - client program for Pipe

	This program uses the pipe for analysis, and is the interface between the
	user and the pipe.
*/

#include "pipe.h"

int main(void)
{
	Pipe pipe(10, 0, 0, 10);

	pipe.setupMatrices();

	return 0;
}

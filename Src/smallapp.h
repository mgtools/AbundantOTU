#ifndef __SMALL_APP_
#define __SMALL_APP_

#define MBIG 1000000000
#define MSEED 16183398
#define MZ 0
#define FAC (1.0/MBIG)

#include "lib.h"

namespace smallapp
{
	int  	ran_number(int n, int *idum);
	double 	random1(int *idum);
	char 	char2intgen(char c);
	char 	char2int(char c);
	FILE*	ckopen(char *name, char *mode);
	int	maxlength(int *index, int size);
	long	trans_seq(char *seq, int len);
	double	PoissonDist(int ave, int mincount);
	int     getnumread(FILE *fp);

};

#endif


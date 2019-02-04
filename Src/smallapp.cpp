#include <limits.h>
#include "smallapp.h"

#define MBIG 1000000000
#define MSEED 16183398
#define MZ 0
#define FAC (1.0/MBIG)

double smallapp::random1(int *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj, mk;
	int i, ii,k;

	if(*idum < 0 || iff == 0)	{ /* initialization */
		iff = 1;
		mj = MSEED - (*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55] = mj;
		mk = 1;
		for(i = 1; i < 54; i ++)	{
			ii = (21 + i) % 55;
			ma[ii] = mk;
			mk = mj - mk;
			if(mk < MZ)	mk += MBIG;
			mj = ma[ii];
		}
		for(k = 1; k <= 4; k ++)	{
			for(i = 1; i <= 55; i ++)	{
				ma[i] -= ma[1 + (i + 30) % 55];
				if(ma[i] < MZ)	ma[i] += MBIG;
			}
		}
		inext = 0;
		inextp = 31;
		*idum = 1;
	}
	if( ++inext == 56) inext = 1;
	if( ++inextp == 56) inextp = 1;
	mj = ma[inext] - ma[inextp];
	if(mj < MZ) mj += MBIG;
	ma[inext] = mj;
	return(mj * FAC);
}

int smallapp::ran_number(int n, int *idum)
{
	double	     t;
	int          p;

	t = random1(idum);
	if(t == 1.0)	{
		t = t - 1.0e-10;
	}
	p = ( int ) (n * t);
	return( p );
}

char smallapp::char2int(char c)
{
	int	i, k;
	//char	con[4] = {'0', '1', '2', '3'};

	if(c >= 'A' && c <= 'Z')	c = c - 'A' + 'a';

	for(i = 0; i < total_nuc; i ++)		{
		if(c == na_name[i])	{
			k = i;
			break;
		}
	}

	if((i == total_nuc) || (k > 3)) {
		k = rand() % 4;
	}
	/*
	int	idum = 1; //to be checked!!
	if(i == total_nuc)	{
		printf("Not found %c\n", c);
		k = ran_number(4, &idum);
	}

	if(k > 3)	{
		k = ran_number(4, &idum);
	}
	*/
	return(k);
}

char smallapp::char2intgen(char c)
{
	int	i, k;

	if(c >= 'A' && c <= 'Z')	c = c - 'A' + 'a';
	for(i = 0; i < total_nuc; i ++)		{
		if(c == na_name[i])	{
			k = i;
			break;
		}
	}

	int	idum = 1; //to be checked!!

	if(i == total_nuc)	{
		printf("Not found %c\n", c);
		k = ran_number(4, &idum);
	}

	return(k);
}

FILE* smallapp::ckopen(char *name, char *mode)
{
	FILE *fp;

	if ((fp = fopen(name, mode)) == NULL)   {
		printf("Cannot open %s.\n", name);
		exit(-1);
	}
	return(fp);
}

int smallapp::maxlength(int *index, int size)
{
	int     i, j;

	j = 0;
	for(i = 0; i < size; i ++)      {
		if(index[i] > j)        j = index[i];
	}
	return(j);
}

long smallapp::trans_seq(char *seq, int len)
{               
	int     i, j, k;
	long     res;

	res = 0;
	//printf("len = %d\n", len);
	for(i = 0; i < len; i ++)       {
		if(res > (LONG_MAX - seq[i]) / 4) {
			printf("res exceeds LONG_MAX\n");
			printf("consider smaller k-mer, len %d, seq %c\n", len, na_name[seq[i]] - 'a' + 'A');
			exit(0);
		}
		res = res * 4 + seq[i];
		//printf("i %d seq %c %d res %ld\n", i, seq[i], seq[i], res);
	}

	return(res);
}       

double smallapp::PoissonDist(int ave, int mincount)
{
	int     i, j, k;
	double  logsum, logprob;

	logsum = 0;
	for(i = 1; i <= mincount; i ++) {
		logsum += log(i);
	}
	logprob = mincount * log(ave) - ave - logsum;
	return(logprob);
}

int smallapp::getnumread(FILE *fp) 
{
	int     n;
	char    str[MAX_SEQLENGTH + 2];

	n = 0;  
	while(fgets(str, MAX_SEQLENGTH, fp))    {
		if(str[0] == '>')       {
			n ++;
		}
	}
	return(n);
}               

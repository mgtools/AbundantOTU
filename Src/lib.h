#ifndef __LIB_H__
#define __LIB_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

/* NOW() */
#include <time.h>
#define NOW(fp) \
	{ \
	time_t now = time(NULL); \
	fputs(ctime(&now), fp); \
	}
#define total_nuc 16
#define na_name "gactnrywsmkhbvdx"

#define pow1(x) ((x) * (x))
#define rev(x) ((2 + (x)) % 4)
#define max(x,y) ((x) > (y) ? (x) : (y))
#define min(x,y) ((x) < (y) ? (x) : (y))
#define numc(i, j) ((i) < (j)? (j)*((j)+1)/2+(i): (i)*((i)+1)/2+(j))
#define MAX_SEQLENGTH 1000000
#define LARGENUMBER 10000000

typedef struct smallhash {
	long    hashval;
	int     count;
	int     pos;
	struct smallhash *left;
	struct smallhash *right;
} SMALLTABLE;

using namespace std;

#endif

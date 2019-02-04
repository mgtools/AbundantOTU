#include "hashapp.h"

SMALLTABLE* hashapp::ins_larg_hash(SMALLTABLE *hash_table, long ctl)
{
	return ins_larg_hash(hash_table, ctl, -1);
}

SMALLTABLE* hashapp::ins_larg_hash(SMALLTABLE *hash_table, long ctl, int pos0)
//YY, Aug 3, 2011, add pos0
{
	int	i, j, k, l;

	if(!hash_table)	{
		//printf("create new hash_table\n");
		hash_table = new SMALLTABLE;
		hash_table -> hashval = ctl;
		hash_table -> count = 1;
		hash_table -> left = NULL;
		hash_table -> right = NULL;
		hash_table -> pos = pos0; //YY, Aug 3, 2011
	} else	{
		//printf("to insert, ctl %d, hashval %d\n", ctl, hash_table->hashval);
		if(hash_table -> hashval == ctl)	{
			hash_table -> count ++;
		} else if(hash_table -> hashval < ctl)	{
			hash_table -> right = ins_larg_hash(hash_table -> right, ctl, pos0);
		} else	{
			hash_table -> left = ins_larg_hash(hash_table -> left, ctl, pos0);
		}
	}
	return(hash_table);
}

SMALLTABLE* hashapp::free_table(SMALLTABLE *hash)
{
	if(hash -> left)	{
		hash -> left = free_table(hash -> left);
	}
	if(hash -> right)	{
		hash -> right = free_table(hash -> right);
	}
	delete hash;
	return(NULL);
}

int hashapp::TopHash(SMALLTABLE **hash_table, int top, int coverthresh, long int n_ban, int *num_cover, long int **ctl)
{
	int     j, k, x, t;
	long int        i, ctltmp[2];
	int     maxn = 0;
	SMALLTABLE *hash_temp;

	for(t = 0; t < top; t ++) {
		j = 0;
		for(i = 0; i < n_ban; i ++)     {
			hash_temp = hash_table[i];
			k = MaxOneHash(hash_temp, ctltmp);
			//printf("t = %d, k = %d, ctltmp %ld %ld\n", t, k, ctltmp[0], ctltmp[1]);
			//getchar();
			if((k > j) && (t == 0 || k < maxn))     {
				j = k;
				ctl[t][0] = i;
				ctl[t][1] = ctltmp[1];
				//printf("ctl %ld %ld\n", ctl[t][0], ctl[t][1]);
			}
		}
		if(j < coverthresh) break;
		num_cover[t] = j;
		//maxn = int(j * 0.9);
		maxn = int(j * 0.8);
		//printf("top %d, num_cover %d\n", t, num_cover[t]);
	}
	return t;
}

int hashapp::SmartHash2(SMALLTABLE **hash_table, int top, int guesscov, int coverthresh, long int n_ban, int *num_cover, long int **ctl)
{
	int     j, k, x, t;
	long int        i, ctltmp[2];
	SMALLTABLE *hash_temp;
	int     *maxn = new int[top];

	num_cover[0] = MaxHash(hash_table, ctl[0], n_ban);
	maxn[0] = num_cover[0];
	printf("cover %d ", maxn[0]);
	for(t = 1; t < top; t ++) {
		maxn[t] = int(maxn[t - 1] * 0.5); 
		num_cover[t] = 0;
		printf("%d ", maxn[t]);
	}
	printf("\n");

	j = 0;
	for(i = 0; i < n_ban; i ++)     {
		hash_temp = hash_table[i];
		k = MaxOneHash(hash_temp, ctltmp);
		for(t = 1; t < top; t ++) {
			if(num_cover[t] < k && k <= maxn[t]) {
				num_cover[t] = k;
				ctl[t][0] = i;
				ctl[t][1] = ctltmp[1];
			}
		}
	}

	int	actual = 1;
	for(t = 1; t < top; t ++) {
		if(num_cover[t] < coverthresh) break;
		if(num_cover[t] != num_cover[actual - 1]) {
			num_cover[actual] = num_cover[t];
			ctl[actual][0] = ctl[t][0]; 
			ctl[actual][1] = ctl[t][1];
			actual ++;
		}
	}
	return actual;
}

int hashapp::SmartHashS2B(SMALLTABLE **hash_table, int top, int guesscov, int coverthresh, long int n_ban, int *num_cover, long int **ctl)
{
	int     j, k, x, t;
	long int        i, ctltmp[2];
	SMALLTABLE *hash_temp;
	int     *maxn = new int[top];

	num_cover[top - 1] = MaxHash(hash_table, ctl[top - 1], n_ban);
	int	maxcover = num_cover[top - 1];

	if(guesscov == 0 || guesscov > maxcover) {
		maxn[top - 1] = maxcover;
		for(t = top - 2; t >= 0; t --) {
			maxn[t] = int(maxn[t + 1] * 0.5); 
			num_cover[t] = 0;
		}
		for(i = 0; i < n_ban; i ++)     {
			hash_temp = hash_table[i];
			k = MaxOneHash(hash_temp, ctltmp);
			for(t = 0; t < top - 1; t ++) {
				if(num_cover[t] < k && k <= maxn[t]) {
					num_cover[t] = k;
					ctl[t][0] = i;
					ctl[t][1] = ctltmp[1];
				}
			}
		}
	}
	else {
		double	ratio = pow(guesscov * 1.0 / maxcover, 1.0 / (top - 2));
		if(guesscov * ratio < 100) {
			ratio = pow(guesscov * 1.0 / maxcover, 1.0 / (top - 1));
		}
		maxn[top - 1] = maxcover;
		for(t = top - 2; t >= 0; t --) {
			maxn[t] = int(maxn[t + 1] * ratio);
			num_cover[t] = 0;
		}
		//printf(" --maxcover %d, guesscov %d, lastcover %d, ratio %.2f\n", maxcover, guesscov, maxn[top - 1], ratio);
		for(i = 0; i < n_ban; i ++)     {
			hash_temp = hash_table[i];
			k = MaxOneHash(hash_temp, ctltmp);
			for(t = 0; t < top - 1; t ++) {
				if(num_cover[t] < k && k <= maxn[t]) {
					num_cover[t] = k;
					ctl[t][0] = i;
					ctl[t][1] = ctltmp[1];
				}
			}
		}
	}

	int	actual = 0;
	for(t = 0; t < top; t ++) {
		if(num_cover[t] < coverthresh) break;
		if((actual == 0) || (num_cover[t] != num_cover[actual - 1])) {
			num_cover[actual] = num_cover[t];
			ctl[actual][0] = ctl[t][0]; 
			ctl[actual][1] = ctl[t][1];
			actual ++;
		}
	}
	return actual;
}

int hashapp::SmartHash(SMALLTABLE **hash_table, int top, int guesscov, int coverthresh, long int n_ban, int *num_cover, long int **ctl)
{
	int     j, k, x, t;
	long int        i, ctltmp[2];
	SMALLTABLE *hash_temp;
	int     *maxn = new int[top];

	num_cover[0] = MaxHash(hash_table, ctl[0], n_ban);
	int	maxcover = num_cover[0];

	if(guesscov == 0 || guesscov > maxcover) {
		maxn[0] = maxcover;
		for(t = 1; t < top; t ++) {
			maxn[t] = int(maxn[t - 1] * 0.5); 
			num_cover[t] = 0;
		}
		for(i = 0; i < n_ban; i ++)     {
			hash_temp = hash_table[i];
			k = MaxOneHash(hash_temp, ctltmp);
			for(t = 1; t < top; t ++) {
				if(num_cover[t] < k && k <= maxn[t]) {
					num_cover[t] = k;
					ctl[t][0] = i;
					ctl[t][1] = ctltmp[1];
				}
			}
		}
	}
	else {
		double	ratio = pow(guesscov * 1.0 / maxcover, 1.0 / (top - 2));
		if(guesscov * ratio < 100) {
			ratio = pow(guesscov * 1.0 / maxcover, 1.0 / (top - 1));
		}
		maxn[0] = maxcover;
		for(t = 1; t < top; t ++) {
			maxn[t] = int(maxn[t - 1] * ratio);
			num_cover[t] = 0;
		}
		//printf(" --maxcover %d, guesscov %d, lastcover %d, ratio %.2f\n", maxcover, guesscov, maxn[top - 1], ratio);
		for(i = 0; i < n_ban; i ++)     {
			hash_temp = hash_table[i];
			k = MaxOneHash(hash_temp, ctltmp);
			for(t = 1; t < top; t ++) {
				if(num_cover[t] < k && k <= maxn[t]) {
					num_cover[t] = k;
					ctl[t][0] = i;
					ctl[t][1] = ctltmp[1];
				}
			}
			//if(k > maxn[top - 1]) break;
		}
	}

	int	actual = 1;
	for(t = 1; t < top; t ++) {
		if(num_cover[t] < coverthresh) break;
		if(num_cover[t] != num_cover[actual - 1]) {
			num_cover[actual] = num_cover[t];
			ctl[actual][0] = ctl[t][0]; 
			ctl[actual][1] = ctl[t][1];
			actual ++;
		}
	}
	return actual;
}

/*
int hashapp::SmartHash(SMALLTABLE **hash_table, int top, int guesscov, int coverthresh, long int n_ban, int *num_cover, long int **ctl)
{
int     j, k, x, t;
long int        i, ctltmp[2];
SMALLTABLE *hash_temp;
int     *maxn = new int[top];

num_cover[0] = MaxHash(hash_table, ctl[0], n_ban);
int	maxcover = num_cover[0];

if(guesscov == 0) {
maxn[0] = maxcover;
for(t = 1; t < top; t ++) {
//maxn[t] = int(maxn[t - 1] * 0.9); 
maxn[t] = int(maxn[t - 1] * 0.5); 
num_cover[t] = 0;
}
for(i = 0; i < n_ban; i ++)     {
hash_temp = hash_table[i];
k = MaxOneHash(hash_temp, ctltmp);
for(t = 1; t < top; t ++) {
if(num_cover[t] < k && k <= maxn[t]) {
num_cover[t] = k;
ctl[t][0] = i;
ctl[t][1] = ctltmp[1];
}
}
}
}
else {
//double ratio = 2.0;
double	ratio = 4.0;
double	maxguess = guesscov;
for(t = 1; t < top; t ++) maxguess *= ratio;
if(maxguess > maxcover) {
ratio = pow(maxcover * 1.0 / guesscov, 1.0 / (top - 1));
}
printf(" --gusscover %d, max-cover %d, ratio changed to %.2f\n", guesscov, num_cover[0], ratio);
maxn[0] = guesscov;
num_cover[0] = 0;
for(t = 1; t < top; t ++) {
maxn[t] = int(maxn[t - 1] * ratio);
num_cover[t] = 0;
}
if(maxguess > maxcover) maxn[top - 1] = maxcover;
for(i = 0; i < n_ban; i ++)     {
hash_temp = hash_table[i];
k = MaxOneHash(hash_temp, ctltmp);
for(t = 0; t < top; t ++) {
if(num_cover[t] < k && k <= maxn[t]) {
num_cover[t] = k;
ctl[t][0] = i;
ctl[t][1] = ctltmp[1];
}
}
//if(k > maxn[top - 1]) break;
}
}

int	actual = 1;
for(t = 1; t < top; t ++) {
if(num_cover[t] < coverthresh) break;
if(num_cover[t] != num_cover[actual - 1]) {
num_cover[actual] = num_cover[t];
ctl[actual][0] = ctl[t][0]; 
ctl[actual][1] = ctl[t][1];
actual ++;
}
}
return actual;
}
*/

int hashapp::MaxHash(SMALLTABLE **hash_table, long *ctl, int n_ban)
{
	int	j, k;
	long	i, ctltmp[2];
	SMALLTABLE *hash_temp;

	j = 0;
	for(i = 0; i < n_ban; i ++)	{
		hash_temp = hash_table[i];
		k = MaxOneHash(hash_temp, ctltmp);
		if(k > j)	{
			j = k;
			ctl[0] = i;
			ctl[1] = ctltmp[1];
		}
	}
	return(j);
}


int hashapp::MaxOneHash(SMALLTABLE *hash_temp, long *ctl)
{
	int	i, j, k; 
	long	ctltmp[2];

	if(!hash_temp)	return(0);
	ctl[1] = hash_temp -> hashval;
	j = hash_temp -> count;

	k = MaxOneHash(hash_temp -> left, ctltmp);
	if(k > j)	{
		j = k;
		ctl[1] = ctltmp[1];
	}
	k = MaxOneHash(hash_temp -> right, ctltmp);
	if(k > j)	{
		j = k;
		ctl[1] = ctltmp[1];
	}
	return(j);
}

int hashapp::getcount(SMALLTABLE *hash_table, long ctl)
{
	int	c;
	if(!hash_table) {
		return(0);
	} else  {
		if(hash_table -> hashval == ctl)        {
			return(hash_table -> count);
		} else if(hash_table -> hashval < ctl)  {
			c = getcount(hash_table -> right, ctl);
		} else  {
			c = getcount(hash_table -> left, ctl);
		}
	}
	return(c);
}

SMALLTABLE* hashapp::FindHash(SMALLTABLE *hash_table, long ctl)
{
	SMALLTABLE *hash_temp, *hash_temp1;

	hash_temp = hash_table;
	if(hash_temp)   {
		if(hash_temp -> hashval == ctl) {
			return(hash_temp);
		}
		hash_temp1 = FindHash(hash_temp -> right, ctl);
		if(!hash_temp1) {
			hash_temp = FindHash(hash_temp -> left, ctl);
		} else  {
			hash_temp = hash_temp1;
		}
	}
	return(hash_temp);
}

SMALLTABLE* hashapp::ins_seq_hash(SMALLTABLE *hash_table, long ctl, int pos)
{
	if(!hash_table) {
		hash_table = new SMALLTABLE;
		hash_table -> hashval = ctl;
		hash_table -> count = 1;
		hash_table -> pos = pos;
		hash_table -> left = NULL;
		hash_table -> right = NULL;
	} else  {
		if(hash_table -> hashval == ctl)        {
			hash_table -> count ++;
		} else if(hash_table -> hashval < ctl)  {
			hash_table -> right = ins_larg_hash(hash_table -> right, ctl, pos);
			//YY, Aug 3, 2011 add pos
		} else  {
			hash_table -> left = ins_larg_hash(hash_table -> left, ctl, pos);
			//YY, Aug 3, 2011 add pos
		}
	}
	return(hash_table);
}

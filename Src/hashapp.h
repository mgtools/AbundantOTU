#ifndef __HASHAPP_H_
#define __HASHAPP_H_

#include "lib.h"

namespace hashapp
{
	SMALLTABLE* ins_larg_hash(SMALLTABLE *hash_table, long ctl);
	SMALLTABLE* ins_larg_hash(SMALLTABLE *hash_table, long ctl, int pos0); //YY, Aug 3, 2011
	SMALLTABLE* free_table(SMALLTABLE *hash);
	int	TopHash(SMALLTABLE **hash_table, int top, int coverthresh, long int n_ban, int *num_cover, long int **ctl);
	int	SmartHash(SMALLTABLE **hash_table, int top, int guesscov, int coverthresh, long int n_ban, int *num_cover, long int **ctl);
	int	SmartHashS2B(SMALLTABLE **hash_table, int top, int guesscov, int coverthresh, long int n_ban, int *num_cover, long int **ctl);
	int	SmartHash2(SMALLTABLE **hash_table, int top, int guesscov, int coverthresh, long int n_ban, int *num_cover, long int **ctl);
	int	MaxHash(SMALLTABLE **hash_table, long *ctl, int n_ban);
	int	MaxOneHash(SMALLTABLE *hash_temp, long *ctl);
	int	getcount(SMALLTABLE *hash_table, long int ctl);
	SMALLTABLE* FindHash(SMALLTABLE *hash_table, long ctl);
	SMALLTABLE* ins_seq_hash(SMALLTABLE *hash_table, long ctl, int pos);
};

#endif

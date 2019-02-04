#ifndef __ALIGN_APP_
#define __ALIGN_APP_

#include "lib.h"
#include "smallapp.h"
#include "hashapp.h"
#include <vector>

class alignapp 
{
protected:
	int 	bandwidth;
	char 	strandtag;
	double 	alnperc;
	int 	word_len;
	int	kmer;
	long 	n_ban, n_ban1;
	int	kmer_recruit;
	long	n_ban_recruit, n_ban1_recruit;
	int 	idum;
	int 	mu1, mu2, gamma0;
	int 	mgsize;
	//double 	PercWeight;
	double 	penalty;
	int	num_con;
	int	num_con_valid;
	int	verbose;
	int	abundantonly;

public:
	alignapp(void);
	void	initenv(void);
	void	compute_n_ban_recruit(void);
	void	compute_n_ban(void);
	int 	CAlign(char *consseq, char **seq, int *len_seq, int num_seq, int *num_mismatch, int maxconlen, int mingsize, int ic0, double PW);
	int 	align2seq(SMALLTABLE **hash_table, char *Conseq, int len_Conseq, char *src_seq, int len_seq, int speedup, double alnperc, vector<int> &aln, int &alnlen);
	int 	global_align(char *seq1, char *seq2, int len1, int len2, int b, vector<int> &aln, int &alnlen);
	int 	global_align_eff(char *seq1, char *seq2, int len1, int len2, int b, vector<int> &aln, int &alnlen);

	void	set_kmer(int value) { kmer = value; compute_n_ban(); }
	void	set_bandwidth(int value) { bandwidth = value; }
	//void	set_PercWeight(float value) { PercWeight = value; }
	void	set_alnperc(float value) { alnperc = value; printf("alnperc changed %f\n", alnperc); }
	void	set_mgsize(int value) { mgsize = value; }
	void	set_strandtag(void) { strandtag = 1; }
	void    setverbose(void) { verbose = 1; }
        void    setabundantonly(void) { abundantonly = 1; printf("set abundantonly to 1\n"); }
	int     buildfullhash(SMALLTABLE **hash_table, FILE *fp);
	void    OneHashBuild(SMALLTABLE **hash_table, char *src_seq, int len_seq);
	void    OneHashCover(SMALLTABLE **hash_table, char *src_seq, int len_seq);

};

#endif

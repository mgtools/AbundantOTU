#ifndef __CAAlign_H_
#define __CAAlign_H_

#include "alignapp.h"
#include "smallapp.h"

class calign:public alignapp
{
private:
	int	top;
	int     hash_ban; 
	int     coverthresh;
	int	seq_sym;
	int     num_seq, *len_seq, *len_seq_pe, totseq;
	bool	paired;
	bool	ifswap;
	char    **src_name;
	char    **src_seq;
	char    **src_name_pe;
	char    **src_seq_pe;
	int	**src_qual;
	char    *index;
	int     *pos;
	int	*seq2clust;
	double	*seq2clustdiff;
	int     *num_pref, *num_suf;
	int     *kmerseq, kmerpos;
	char    **prefseq, **sufseq;
	char    **Conseq, **Conseq_pe;
	int	*mem_Conseq;
	int     *len_prefseq, *len_sufseq;
	double  cov;
	double	P;
	char    *tmpseq;
	SMALLTABLE **hash_table, **hash_table_pe;
	int     *len_Conseq, *len_Conseq_pe;
	int	*con_sort;
	int	*con2sort;
	int	ifrefine;
	int	max_seqlen;
	int	speedup;
	int	seqnum_limit;
	int	start_idx;
	int	unassigned_tot;
	int	unassigned_valid;
	int	*unassigned_idx;
	bool	checkalignlen;

	char	fileprefix[1000];
	FILE	*log;

public:
	~calign(void) { cleanup(); }
	void 	ca2otu(char *prefix);
	void 	ca2otu(int smallclust, char *prefix);
	void	ini_calign(int, int);
	void	copyreads(int, int *, char**, char**, int*);
	void    getreads(int tot, char **seq, char **name, int* len);
	void	getreads(char *seqfile);
	void	getreads_pe(char *seqfile, char *seqfile_pe);
	void	getreads(char *seqfile, int minlen);
	void	getreads_pe(char *seqfile, char *seqfile_pe, int minlen);
	int	readseq(char **seq, int *len_seq, char **src_name, FILE *fp, int minlen);
	void	getqual(char *qualityfile);
	void	iteralign(void);
	void	greedyrecruit(void);
	char*	oneconsensus(int &len);
	char*	pickconsensus(int &conseqlen, int &seqi); //pick longest unassigned sequence as consensus
	int	getconsensus(int num_cover, long *ctl, char *conseq); //consensus by consensus alignment
	void 	map2consensus(int &conseqlen, char *conseq, int& assignedtot, int *assignedseq, double *assignedcov);
	void 	map2consensus_pe(int &conseqlen, char *conseq, int &conseqlen_pe, char *conseq_pe, int& assignedtot, int *assignedseq, double *assignedcov);
	void	map2consensus_pe2(int &len_conseq, char *conseq, int &len_conseq_pe, char *conseq_pe, int& assignedtot, int *assignedseq, double *assignedcov);
	int	prefsufseq(int num_cover, long int* ctl, char **src_seq, int *len_seq, int num_seq, char *index,
		char **prefseq, int *len_prefseq, char **sufseq, int *len_sufseq, int *kmerseq, int *kmerpos);
	void	cleanup(void);
	void	removesmall(int small_clust_size, bool ifempty);

	void	saveresults(void);
	char*   int2seq(char *src_seq, int len_seq);
	char*   seq2int(char *src_seq, int len_seq);
	void    writeseq(char *src_seq, char *src_name, int len_seq, FILE *fp);
	void	printaln(char *seq1, int len1, char *seq2, int len2, vector<int> aln);
	void    writesingleton(void);
	void    writeconseq(void);
	void    writeconseq_pe(void);
	void    writeclust(void);
	void    writeclustsize(void);
	void	writelog(int conid, int conlen, char *conseq, int assignedtot, int *assigned);
	void	writelog_pe(int conid, int conlen, char *conseq, int conlen_pe, char *conseq_pe, int assignedtot, int *assigned);

	void	assignseq(int num_con, int tot, int *assignedseq, double *assignedcov);
	void	getunassigned(int ifreshuffle);
	int	hashcover(char **src_seq, int *len_seq, int num_seq, char *index, int *num_cover_list, long **ctl_list);
	void    OneHash(SMALLTABLE **hash_table, char *src_seq, int len_seq);
	void    OneSeqHash(SMALLTABLE **hash_table, char *src_seq, int len_seq);
	void	diffstat(int num0, int id, int *diff2id, double *diff, int &min_i, int &max_i, int &median_i);
	void	sortcons(void);

	int     align2intseq(char *seq1, int len1, char *seq2, int len2, int &alnlen);
	int     align2charseq(char *seq1, char *seq2);

	//functions for subotu
	int     getclustsize(int clust);
	int     getclustnum(void);
	int     getclustinfo(int clust, int *mem, double *diff, char *conseq);
	void    subotu(int minclustsize, char *fileprefix);
	void 	setfileprefix(char *tmp);
	void	setswap(void) { ifswap = true; }
	void	chimeracheck(int conseq_len, char *conseq, int num, int *seqlist);
};

#endif

//Program Name: calign
//Latest modification on April 22th, 2014
//Developer: Yuzhen Ye <yye@indiana.edu> and Yongan Zhao <yongzhao@indiana.edu>
//Affiliation: School of Informatics and Computing, Indiana University, Bloomington

#include <limits.h>
#include <algorithm>
#include "calign.h"

using namespace smallapp;
using namespace hashapp;

int ifprint = 0;

void calign::ca2otu(char *prefix)
{
	int	smallclust = 5;
	ca2otu(smallclust, prefix);
}

void calign::ca2otu(int smallclust, char *prefix)
{

	strcpy(fileprefix, prefix);

	int	modelist[2] = {1, 0};
	int	i, t, mode;
	time_t  jobstart = time(NULL);
	printf("\nnow get consensus by consensus alignment\n\n");
	iteralign();

        if(abundantonly == 1) {
                sortcons();
                removesmall(smallclust, true);
        }
        else {
                time_t jobfinished = time(NULL);
                double  timeused = difftime(jobfinished, jobstart);
                if(timeused < 60) printf("Time used for abundantOTU: %d sec\n", int(timeused));
                else printf("Time used for abundantOTU: %d min (%.1f hours)\n", int(timeused / 60), timeused / 3600.0);
                removesmall(3, true); //add by YY
                greedyrecruit(); //modified by YY
                sortcons();
                removesmall(1, true); //modified by YY
        }
	
	//removesmall(3, true); //add by YY
	//greedyrecruit(); //modified by YY
	//sortcons();
	//removesmall(1, true); //modified by YY
	cout<<"Total clusters: "<<num_con_valid<<endl;
}

void calign::cleanup(void)
{
	int	i;
	delete[] sufseq;
	delete[] prefseq;
	delete[] len_prefseq;
	delete[] len_sufseq;
	delete[] tmpseq;
	delete[] pos;
	delete[] num_pref;
	delete[] num_suf;
	delete[] kmerseq;
	delete[] index;
	delete[] unassigned_idx;
	int	dumy = 10; 
	//to avoid segmentation in extreme cases where all input sequences are singletons
	//fixed on June 6, 2011 by YY
	for(i = 0; i < num_seq + dumy; i ++)	{
		delete[] Conseq[i];
	}
	delete[] Conseq;
	delete[] len_Conseq;
	delete[] mem_Conseq;
	delete[] seq2clust;
	delete[] seq2clustdiff;
	for(i = 0; i < num_seq; i ++)   {
		if(seq_sym == 0) {
			delete[] src_seq[i];
			delete[] src_name[i];
		}
	}
	delete[] src_seq;
	delete[] src_name;
	delete[] len_seq;
}

void calign::ini_calign(int n, int nc)
{
	src_name = new char*[n];
	src_seq = new char*[nc];
	len_seq = new int[nc]; 
	if(paired) { //paired-end reads
		src_name_pe = new char*[n];
		src_seq_pe = new char*[nc];
		len_seq_pe = new int[nc]; 
	}
	index = new char[nc];
	unassigned_valid = 0;
	unassigned_tot = 0;
	unassigned_idx = new int[nc];
	pos = new int[nc];
	seq2clust = new int[nc];
	seq2clustdiff = new double[nc];
	num_pref = new int[nc];
	num_suf = new int[nc];
	kmerseq = new int[nc];
	prefseq = new char*[nc];
	sufseq = new char*[nc];
	len_prefseq = new int[nc];
	len_sufseq = new int[nc];
	int	i;
	for(i = 0; i < nc; i ++) {
		index[i] = pos[i] = num_pref[i] = num_suf[i] = 0;
		kmerseq[i] = len_prefseq[i] = len_sufseq[i] = 0;
		seq2clust[i] = -1;
		prefseq[i] = NULL;
		sufseq[i] = NULL;
		len_seq[i] = 0; //Apr 16, 2011
	}
	int	dumy = 10; //to avoid segmentation in extreme cases where all input sequences are singletons
	Conseq = new char*[nc + dumy];
	len_Conseq = new int[nc + dumy];
	mem_Conseq = new int[nc + dumy];
	if(paired) {
		Conseq_pe = new char*[nc + dumy];
		len_Conseq_pe = new int[nc + dumy];
	}
	for(i = 0; i < nc + dumy; i ++) {
		len_Conseq[i] = mem_Conseq[i] = 0;
		Conseq[i] = NULL;
		if(paired) {
			len_Conseq_pe[i] = 0;
			Conseq_pe[i] = NULL;
		}
	}

	tmpseq = new char[MAX_SEQLENGTH];
	coverthresh = 5;
	num_con = 0;
	start_idx = 0;
	con_sort = con2sort = NULL;
	checkalignlen = true;
}

//get reads from "old" reads
void calign::copyreads(int tot, int *idx, char **seq0, char **name0, int *seqlen0)
{
	seq_sym = 1;
	//paired = false;
	//strcpy(fileprefix, "tmp");

	ini_calign(tot, tot);

	int	i, j;
	max_seqlen = 0;
	for(i = 0; i < tot; i ++) {
		j = idx[i];
		len_seq[i] = seqlen0[j];
		src_seq[i] = seq0[j]; //assign address only! 
		src_name[i] = name0[j]; //v1.1
		if(max_seqlen < len_seq[i]) max_seqlen = len_seq[i];
	}
	num_seq = totseq = tot;
	//printf("%d seq copied\n", tot);
}

//get reads from reads (create new reads)
void calign::getreads(int tot, char **seq, char **name, int* len)
{
	ini_calign(tot, tot);
	int     i, j;
	for(i = 0; i < tot; i ++) {
		len_seq[i] = len[i];// added by yongzhao. strlen(seq[i]);
		src_seq[i] = seq[i];// added by yongzhao. seq2int(seq[i], len_seq[i]);
		src_name[i] = name[i]; // added by yongzhao. new char[strlen(name[i]) + 1];
		//strcpy(src_name[i], name[i]);// added by yongzhao.
		writeseq(src_seq[i], src_name[i], len_seq[i], stdout);
	}
	num_seq = totseq = tot;
	max_seqlen = len_seq[0];
	int     min_seqlen = MAX_SEQLENGTH;
	for(int i = 1; i < num_seq; i ++) {
		if(len_seq[i] > max_seqlen) max_seqlen = len_seq[i];
		if(min_seqlen > len_seq[i]) min_seqlen = len_seq[i];
	}
	printf("totseq %d, max_seqlen %d, min_seqlen %d\n", num_seq, max_seqlen, min_seqlen);

	src_qual = NULL;
}

//get reads from file
void calign::getreads(char *seqfile)
{
	getreads(seqfile, kmer);
}

void calign::getreads_pe(char *seqfile, char *seqfile_pe)
{
	getreads_pe(seqfile, seqfile_pe, kmer);
}

//get reads from file
void calign::getreads(char *seqfile, int minlen)
{
	seq_sym = 0;
	paired = false;

	//get read num first
	FILE	*fp0 = smallapp::ckopen(seqfile, "r");
	int	n = getnumread(fp0);
	fclose(fp0);
	ini_calign(n, n);
	//now get reads
	FILE* 	fp = smallapp::ckopen(seqfile, "r");
	printf("now read seqfile: %s..\n", seqfile);
	num_seq = readseq(src_seq, len_seq, src_name, fp, minlen);
	fclose(fp);

	totseq = num_seq;

	max_seqlen = len_seq[0];
	int	min_seqlen = MAX_SEQLENGTH;
	for(int i = 1; i < num_seq; i ++) {
		if(len_seq[i] > max_seqlen) max_seqlen = len_seq[i];
		if(min_seqlen > len_seq[i]) min_seqlen = len_seq[i];
	}
	printf("totseq %d, max_seqlen %d, min_seqlen %d\n", num_seq, max_seqlen, min_seqlen);
	fflush(stdout);

	src_qual = NULL;
}

//get reads from file
void calign::getreads_pe(char *seqfile, char *seqfile_pe, int minlen)
{
	seq_sym = 0;
	paired = true;
	ifswap = false;

	//check/get read num first
	FILE	*fp0 = smallapp::ckopen(seqfile, "r");
	int	n = getnumread(fp0);
	fclose(fp0);
	FILE	*fp0_pe = smallapp::ckopen(seqfile_pe, "r");
	int	n_pe = getnumread(fp0_pe);
	fclose(fp0_pe);
	if(n != n_pe) {
		printf("Input sequences not paired!! %d vs %d\n", n, n_pe);
		exit(0);
	}
	ini_calign(n, n);

	//now get reads
	FILE* 	fp = smallapp::ckopen(seqfile, "r");
	printf("now read seqfile: %s..\n", seqfile);
	num_seq = readseq(src_seq, len_seq, src_name, fp, minlen);
	fclose(fp);
	FILE* 	fp_pe = smallapp::ckopen(seqfile_pe, "r");
	printf("now read seqfile: %s..\n", seqfile_pe);
	num_seq = readseq(src_seq_pe, len_seq_pe, src_name_pe, fp_pe, minlen);
	fclose(fp_pe);

	totseq = num_seq;

	max_seqlen = len_seq[0];
	int	max_seqlen_pe = len_seq_pe[0];
	int	min_seqlen = MAX_SEQLENGTH;
	int	min_seqlen_pe = MAX_SEQLENGTH;
	for(int i = 1; i < num_seq; i ++) {
		if(len_seq[i] > max_seqlen) max_seqlen = len_seq[i];
		if(min_seqlen > len_seq[i]) min_seqlen = len_seq[i];
		if(len_seq_pe[i] > max_seqlen_pe) max_seqlen_pe = len_seq_pe[i];
		if(min_seqlen_pe > len_seq[i]) min_seqlen_pe = len_seq_pe[i];
	}
	printf("totseq-pair %d, max_seqlen %d, min_seqlen %d paired-ends maxlen %d minlen %d\n", num_seq, max_seqlen, min_seqlen, max_seqlen_pe, min_seqlen_pe);
	fflush(stdout);

	src_qual = NULL;
}

int calign::readseq(char **seq, int *len_seq, char **src_name, FILE *fp, int minlen)
{
	int	i, j, num_seq0, currlen;
	char	str[MAX_SEQLENGTH + 2], currname[1000];
	char 	*src_seq;
	src_seq = new char[MAX_SEQLENGTH];
	int	shortremove = 0;
	num_seq0 = 0;
	currlen = 0;
	currname[0] = 0;
	while(fgets(str, MAX_SEQLENGTH, fp))	{
		if(str[0] == '>')	{
			if(currlen >= minlen) {
				seq[num_seq0] = new char[currlen];
				for(i = 0; i < currlen; i ++)	{
					seq[num_seq0][i] = src_seq[i];
				}
				src_name[num_seq0] = new char[strlen(currname) + 1];
				strcpy(src_name[num_seq0], currname);
				len_seq[num_seq0] = currlen;
				num_seq0 ++;
			}
			else if(currlen > 0) {
				shortremove ++;
			}
			currlen = 0;
			sscanf(&str[1], "%s", currname);
		}
		else {
			if(currlen + strlen(str) > MAX_SEQLENGTH) {
				printf("please increase MAX_SEQLENGTH (%d)\n", MAX_SEQLENGTH);
				exit(0);
			} //Apr 21, 2011
			for(i = 0; i < strlen(str); i ++)	{
				if(str[i] >= 'A' && str[i] <= 'Z' || str[i] >= 'a' && str[i] <= 'z')	{
					src_seq[currlen ++] = smallapp::char2int(str[i]);
				}
			}
		}
	}
	if(currlen >= minlen)	{
		seq[num_seq0] = new char[currlen];
		for(i = 0; i < currlen; i ++)	{
			seq[num_seq0][i] = src_seq[i];
		}
		src_name[num_seq0] = new char[strlen(currname) + 1];
		strcpy(src_name[num_seq0], currname);
		len_seq[num_seq0] = currlen;
		num_seq0 ++;
	}
	else if(currlen > 0) {
		shortremove ++;
	}

	if(shortremove) printf("!!!Note: %d short sequences (< %d bp) were removed from the calculation; remain %d\n", shortremove, minlen, num_seq0);
	delete[] src_seq;
	return(num_seq0);
}

void calign::getqual(char *qualfile)
{
	FILE* 	fp = smallapp::ckopen(qualfile, "r");
	src_qual = new int*[num_seq];
	int	num = -1;
	int	posadd = 0;
	char	str[MAX_SEQLENGTH + 2];
	char	name[1000];
	int	i, j, qual;
	while(fgets(str, MAX_SEQLENGTH, fp))	{
		if(str[0] == '>') {
			num ++;
			sscanf(str + 1, "%s", name);
			if(strcmp(name, src_name[num])) {
				printf("quality file error seq %s vs %s\n", name, src_name[num]);
				exit(0);
			}
			src_qual[num] = new int[len_seq[num]];
			posadd = 0;
		}
		else {
			for(i = 0; i < strlen(str); i ++) {
				j = i - 2;
				if(j < 0) j = 0;
				if(str[i] == ' ' || i == strlen(str) - 1) {
					sscanf(str + j, "%d", &qual);
					if(posadd < len_seq[num]) {
						src_qual[num][posadd] = qual;
					}
					posadd ++;	
				}
			}
		}
	}
	fclose(fp);
	num ++;
	if(num != num_seq) {
		printf("quality file is not consistent with sequence file\n");
		exit(0);
	}
}

char* calign::oneconsensus(int &conseqlen)
{
	top = 1;
	int     *num_cover_list = new int[top];
	long    **ctl_list = new long*[top];
	int	t;
	for(t = 0; t < top; t ++) ctl_list[t] = new long[2];

	getunassigned(0);
	int num_cover_tot = hashcover(src_seq, len_seq, num_seq, index, num_cover_list, ctl_list);
	//printf("num_cover_tot %d, highest-one %ld thresh %d\n", num_cover_tot, num_cover_list[0], coverthresh);
	int num_cover_real = prefsufseq(num_cover_list[0], ctl_list[0], src_seq, len_seq, num_seq, index,
		prefseq, len_prefseq, sufseq, len_sufseq, kmerseq, &kmerpos);
	//introduce num_cover_real
	char	*conseq = new char[2 * max_seqlen + 1];

	//conseqlen = getconsensus(num_cover_list[0], ctl_list[0], conseq);
	conseqlen = getconsensus(num_cover_real, ctl_list[0], conseq);

	for(t = 0; t < top; t ++) delete[] ctl_list[t];
	delete[] ctl_list;
	delete[] num_cover_list;

	if(checkalignlen && conseqlen > 1.5 * max_seqlen) {
		conseqlen = max_seqlen;
	}  //YY April 25, 2014 to be evaluated

	return conseq;
}

void calign::iteralign(void)
{
	top = 4;

	int	maxmax = 500000; //or 2000000
	if(num_seq < maxmax) {
		speedup = 0;
		seqnum_limit = num_seq;
	}
	else {
		speedup = 1;
		seqnum_limit = maxmax; //v1.1
	}

	int	l1, l2, i, j, k, n;
	int     t, num_cover_tot;
	int     *num_cover_list = new int[top];
	long    **ctl_list = new long*[top];
	for(t = 0; t < top; t ++) ctl_list[t] = new long[2];
	getunassigned(speedup); //v1.1
	num_cover_tot = hashcover(src_seq, len_seq, num_seq, index, num_cover_list, ctl_list);
	printf("num_seq %d, num_cover_tot %d, highest-one %ld thresh %d\n", num_seq, num_cover_tot, num_cover_list[0], coverthresh);

	int	convergence = 0;
	int	num_cover_real = 0;
	char	temp[1000];
	sprintf(temp, "%s.log", fileprefix);
	log = smallapp::ckopen(temp, "w");
	while(num_cover_tot > 0 && totseq > 0)	{
		//printf("mark 1\n");
		//fflush(stdout);
		int	maxassigned = 0;
		int	assigned, maxconseqlen, maxconseqlen_pe, conseqlen, conseqlen_pe;
		int	*assignedseq, *maxassignedseq;
		double	*assignedcov, *maxassignedcov;
		char	*conseq, *maxconseq, *conseq_pe, *maxconseq_pe;
		for(t = 0; t < num_cover_tot; t ++) {//iterate over several k-mers
			if(verbose) printf(" Consens %d, try k-mer %d, cover %ld ctl %ld\n", num_con + 1, t, num_cover_list[t], ctl_list[t]);
			fflush(stdout);
			num_cover_real = prefsufseq(num_cover_list[t], ctl_list[t], src_seq, len_seq, num_seq, index,
				prefseq, len_prefseq, sufseq, len_sufseq, kmerseq, &kmerpos);
			//printf("    after prefsufseq\n"); 
			//fflush(stdout);
			conseq = new char[2 * max_seqlen + 1];
			if(paired) conseq_pe = new char[2 * max_seqlen + 1];
			assignedseq = new int[totseq];
			assignedcov = new double[totseq];
			//printf("    before getconsensus\n"); 
			//fflush(stdout);
			conseqlen = getconsensus(num_cover_real, ctl_list[t], conseq); //note num_cover_real, April 21, 2011
			//printf("    after getconsensus\n"); 
			if(paired) map2consensus_pe(conseqlen, conseq, conseqlen_pe, conseq_pe, assigned, assignedseq, assignedcov);
			else map2consensus(conseqlen, conseq, assigned, assignedseq, assignedcov);
			if(verbose) printf("  --conseq-length=%d, seq-included=%d\n", conseqlen, assigned);
			fflush(stdout);
			if(t == 0 || assigned > maxassigned) {
				maxconseqlen = conseqlen;
				if(paired) maxconseqlen_pe = conseqlen_pe;
				maxassigned = assigned;
				if(t > 0) {
					delete[] maxassignedseq;
					delete[] maxassignedcov;
					delete[] maxconseq;
					if(paired) delete[] maxconseq_pe;
				}
				maxassignedseq = assignedseq;
				maxassignedcov = assignedcov;
				maxconseq = conseq;
				if(paired) maxconseq_pe = conseq_pe; 
			}
			else {
				delete[] assignedseq;
				delete[] assignedcov;
				delete[] conseq;
				if(paired) delete[] conseq_pe;
			}
			if ((assigned > 1000) || (assigned > 100 && assigned > num_cover_list[t] * 0.3)) break; 
			//the k-mer leads good consensus, otherwise try other k-mer
		}
		totseq -= maxassigned;

		printf("#Consensus %d len %d, seq-included %d (tot-assigned %d %.2f%% remain %d)\n", 
			num_con + 1, conseqlen, maxassigned, num_seq - totseq, (num_seq - totseq)*100.0/num_seq, totseq);
		if(num_con + 1 == 559) ifprint = 1;
		else ifprint = 0;
		if(paired) writelog_pe(num_con, maxconseqlen, maxconseq, maxconseqlen_pe, maxconseq_pe, maxassigned, maxassignedseq);
		else writelog(num_con, maxconseqlen, maxconseq, maxassigned, maxassignedseq);
		fflush(stdout);

		assignseq(num_con, maxassigned, maxassignedseq, maxassignedcov);
		fflush(stdout);

		if(maxassigned < 5) convergence += 1;
		else convergence = 0;

		Conseq[num_con] = maxconseq;
		len_Conseq[num_con] = maxconseqlen;
		if(paired) {
			Conseq_pe[num_con] = maxconseq_pe;
			len_Conseq_pe[num_con] = maxconseqlen_pe;
		}
		mem_Conseq[num_con] = maxassigned;
		//if(len_Conseq[num_con] > 500) {
		//	chimeracheck(len_Conseq[num_con], Conseq[num_con], maxassigned, maxassignedseq);
		//}
		num_con ++; //modified by YY
		fflush(stdout);
		getunassigned(speedup); //v1.1
		delete[] maxassignedseq;
		delete[] maxassignedcov;
		if(convergence >= 2) break;
		fflush(stdout);
		num_cover_tot = hashcover(src_seq, len_seq, num_seq, index, num_cover_list, ctl_list);
		fflush(stdout);
	}
	for(t = 0; t < top; t ++) delete[] ctl_list[t];
	delete[] ctl_list;
	delete[] num_cover_list;
	fclose(log);
}

void calign::chimeracheck(int conseq_len, char *conseq, int num, int *seqlist) 
{
	hash_table = new SMALLTABLE*[n_ban_recruit];
	int	i, j, a;
	for(i = 0; i < n_ban_recruit; i ++) hash_table[i] = NULL;
	OneSeqHash(hash_table, conseq, conseq_len);
	vector<int> aln;
	int k = 0;
	int *cov = new int[conseq_len];
	int *mid = new int[conseq_len];
	for(i = 0; i < conseq_len; i ++) cov[i] = mid[i] = 0;
	for(int i0 = 0; i0 < num; i0 ++) { //v1.1
		i = seqlist[i0];
		int	alnspeedup = 1; //v1.1
		int	alnlen = 0; //v1.1
		int mis = align2seq(hash_table, conseq, conseq_len, src_seq[i], len_seq[i], alnspeedup, alnperc, aln, alnlen); //v1.1
		printaln(conseq, conseq_len, src_seq[i], len_seq[i], aln);
		int	p1 = 0;
		int	start, end;
		start = end = -1;
		for (j = 0; j < aln.size(); j ++) {
			a = aln[j];
			if(a != 3) p1 ++; //not gap in conseq
			if(a != 2) cov[p1 - 1] ++; //not gap in seq
			if(start == -1 && a != 2) start = p1;
			if(a != 2) end = p1;
		}
		mid[(end - start) / 2 + start] ++;
	}
	for(i = 0; i < conseq_len; i ++) printf("%d %d %d\n", i, cov[i], mid[i]); 
	for(i = 0; i < n_ban_recruit; i ++)	{
		if(hash_table[i])	{
			hash_table[i] = hashapp::free_table(hash_table[i]);
		}
	}
	delete[] hash_table;
	delete[] mid;
	delete[] cov;
}

//remaining sequences
void calign::greedyrecruit(void)
{
	int	*assignedseq;
	double	*assignedcov;
	int	assigned, conseqlen;

	int	convergence = 0;
	int	seqi;
	int	smallclust = 5;
	int	step = 0;
	getunassigned(0);
	while((unassigned_tot > 0) && (step <= 10))	{ //use step to terminate,YY Dec 31, 2015
		assignedseq = new int[unassigned_tot];
		assignedcov = new double[unassigned_tot];

		Conseq[num_con] = pickconsensus(conseqlen, seqi);
		printf("pickconsensus: seq %d %s, len %d\n", seqi + 1, src_name[seqi], conseqlen);
		if(!paired) map2consensus(conseqlen, Conseq[num_con], assigned, assignedseq, assignedcov);
		else {
			Conseq_pe[num_con] = new char[len_seq_pe[seqi] + 1];
			for(int j = 0; j < len_seq_pe[seqi]; j ++) {
				Conseq_pe[num_con][j] = src_seq_pe[seqi][j];
			}
			len_Conseq_pe[num_con] = len_seq_pe[seqi];
			Conseq_pe[num_con][len_Conseq_pe[num_con]] = 0;
			map2consensus_pe2(conseqlen, Conseq[num_con], len_Conseq_pe[num_con], Conseq_pe[num_con], assigned, assignedseq, assignedcov);
		}
		len_Conseq[num_con] = conseqlen;

		if(assigned == 0) {
			printf("Warning: map2consensus error, NO sequences recruited to seq %d\n", seqi);
			//exit(0);
			assigned = 1; //just assign itself, YY, Jan9, 2013
			assignedseq[0] = seqi;
			assignedcov[0] = 0;
		}

		totseq -= assigned;
		printf("#greedyrecruit %d len %d seq-included %d (tot-assigned %d %.2f%% remain %d)\n", num_con + 1, len_Conseq[num_con], assigned, num_seq-totseq, (num_seq-totseq)*100.0/num_seq, totseq);
		fflush(stdout);

		assignseq(num_con, assigned, assignedseq, assignedcov);

		delete[] assignedseq;
		delete[] assignedcov;

		mem_Conseq[num_con] = assigned;
		num_con ++;

		getunassigned(0);
	
		if(assigned >= smallclust) step = 0;
		step += 1;
	}
}

//select the longest sequence that remain unassigned as consensus
char* calign::pickconsensus(int &conseqlen, int &seq)
{
	int	thisnumseq = num_seq;
	thisnumseq = num_seq;

	int	maxlen = 0;
	int	maxi = 0;
	int	thislen;
	//for(int i = 0; i < thisnumseq; i ++)	{
	//	if(index[i] == 1)	continue;	
	for(int i0 = 0; i0 < unassigned_tot; i0 ++) {
		int i = unassigned_idx[i0];
		if(paired) thislen = len_seq[i] + len_seq_pe[i]; //use combined length for paired-end reads
		else thislen = len_seq[i];
		if(thislen > maxlen) {
			maxlen = thislen;
			maxi = i;
		}
	}
	seq = maxi;
	conseqlen = len_seq[maxi];
	char* conseq = new char[conseqlen + 1];
	int	j;
	for(j = 0; j < conseqlen; j ++) {
		conseq[j] = src_seq[maxi][j];
	}
	conseq[conseqlen] = 0;
	return conseq;
}

//get consensus sequence by consensus alignment
int calign::getconsensus(int num_cover, long *ctl, char *conseq)
{
	int	i, j, n;
	int	l1 = smallapp::maxlength(len_prefseq, num_cover);
	int	l2 = smallapp::maxlength(len_sufseq, num_cover);
	//derive consense first
	//conseq = new char[l1 + l2 + kmer + 2 * bandwidth + 1];
	conseq[0] = 0;
	for(j = 0; j < num_cover; j ++)	{
		num_pref[j] = 0;
	}
	//backward extension
	n = CAlign(tmpseq, prefseq, len_prefseq, num_cover, num_pref, l1 + bandwidth, mgsize, 0, alnperc);
	if(verbose) printf("  --backward n = %d\n", n);
	for(j = 0; j < num_cover; j ++)	{
		num_suf[j] = num_pref[j];
	}
	for(i = 0; i < n; i ++)	{
		conseq[i] = tmpseq[n - 1 - i];
	} //assign prefix to Conseq
	for(i = 0; i < kmer; i ++)	{
		conseq[n + i] = src_seq[kmerseq[0]][kmerpos + i];
	} //assign kmer sequence to Conseq

	int	len_conseq = kmer + n;

	//forward extension
	n = CAlign(tmpseq, sufseq, len_sufseq, num_cover, num_suf, l2 + bandwidth, mgsize, n, alnperc * 2);
	if(verbose) printf("  --forward n = %d\n", n);
	for(i = 0; i < n; i ++)	{
		conseq[i + len_conseq] = tmpseq[i];
	} //assign sufix to Conseq
	len_conseq += n;
	conseq[len_conseq] = 0;

	for(i = 0; i < num_cover; i ++)	{
		//printf("i %d, num_cover %d prelen %d sufseq %d\n", i, num_cover, len_prefseq[i], len_sufseq[i]);
		delete[] prefseq[i];
		delete[] sufseq[i];
	}

	if(checkalignlen && len_conseq > 1.5 * max_seqlen) {
		len_conseq = max_seqlen;
	}  //YY April 25, 2014 to be evaluated

	return len_conseq;
}

int calign::align2charseq(char *seq1, char *seq2)
{
	char    *seq1_int = seq2int(seq1, strlen(seq1));
	char    *seq2_int = seq2int(seq2, strlen(seq2));
	int	alnlen = 0;
	int     mis = align2intseq(seq1_int, strlen(seq1), seq2_int, strlen(seq2), alnlen);
	delete[] seq1_int;
	delete[] seq2_int;
	return mis;
}

//add alnlen by YY
int calign::align2intseq(char *seq1, int len1, char *seq2, int len2, int &alnlen)
{
	hash_table = new SMALLTABLE*[n_ban_recruit];
	int     i;
	for(i = 0; i < n_ban_recruit; i ++) hash_table[i] = NULL;
	OneSeqHash(hash_table, seq1, len1);
	vector<int> aln;
	alnlen = 0;
	int mis = align2seq(hash_table, seq1, len1, seq2, len2, 0, 0, aln, alnlen);
	for(i = 0; i < n_ban_recruit; i ++)     {
		if(hash_table[i])       {
			hash_table[i] = hashapp::free_table(hash_table[i]);
		}
	}
	return mis;
}

//-----------------------------------------------------------------------------------------------
void calign::map2consensus(int &len_conseq, char *conseq, int& assignedtot, int *assignedseq, double *assignedcov)
{
	hash_table = new SMALLTABLE*[n_ban_recruit];
	int	i, tryn, t;
	for(i = 0; i < n_ban_recruit; i ++) hash_table[i] = NULL;

	//then define which sequences can be included for the current consensus seq
	int	k = 0;
	int miscut = int(alnperc * len_conseq + 0.5); //round up
	vector<int> aln;

	tryn = 2;

	for(t = 0; t < tryn; t ++) {
		OneSeqHash(hash_table, conseq, len_conseq);
		k = 0;
		for(int i0 = 0; i0 < unassigned_tot; i0 ++) { //v1.1
			i = unassigned_idx[i0];
			int	alnspeedup = 1; //v1.1
			int	alnlen = 0; //v1.1
			int mis = align2seq(hash_table, conseq, len_conseq, src_seq[i], len_seq[i], alnspeedup, alnperc, aln, alnlen); //v1.1
			int	miscut1 = int(alnlen * alnperc + 0.5); //v1.1
			if(mis <= miscut1) {
				//printaln(conseq, len_conseq, src_seq[i], len_seq[i], aln);
				cov = double(mis) / alnlen; //v1.1
				assignedseq[k] = i;
				assignedcov[k ++] = cov;
			}
		}
		for(i = 0; i < n_ban_recruit; i ++)	{
			if(hash_table[i])	{
				hash_table[i] = hashapp::free_table(hash_table[i]);
			}
		}
		if(k < 10 || k >= 10000) break; //v1.1
		if(t == 0 && tryn > 1) {
			calign *ca = new calign;
			ca -> copyreads(k, assignedseq, src_seq, src_name, len_seq);
			int	len_newcons;
			char	*newcons = ca -> oneconsensus(len_newcons);
			delete ca;
			int	changed = 1;
			if(len_newcons == len_conseq) {
				for(i = 0; i < len_conseq; i ++) {
					if(newcons[i] != conseq[i]) break;
				}
				if(i == len_conseq) changed = 0;
			}
			if(!changed) {
				break; //no change to consensus
			}
			else {
				len_conseq = len_newcons;
				for(i = 0; i < len_conseq; i ++) conseq[i] = newcons[i];
				delete[] newcons;
			}
		}
	}

	assignedtot = k;

	delete[] hash_table;
}

//-----------------------------------------------------------------------------------------------
//recruitment using both ends; but only the consensens sequence of reads 1 are given for this function -- different from map2consensus_pe2() 
void calign::map2consensus_pe(int &len_conseq, char *conseq, int &len_conseq_pe, char *conseq_pe, int& assignedtot, int *assignedseq, double *assignedcov)
{
	hash_table = new SMALLTABLE*[n_ban_recruit];
	int	i, tryn, t;
	for(i = 0; i < n_ban_recruit; i ++) hash_table[i] = NULL;

	//then define which sequences can be included for the current consensus seq
	int	k = 0;
	vector<int> aln_1;
	vector<int> aln_2;
	int	alnspeedup = 1;
	int	alnlen_1, alnlen_2, alnlen, mis, miscut1;

	//read 1
	OneSeqHash(hash_table, conseq, len_conseq);
	k = 0;
	int	recheck = 0;
	int	*rechecklist = new int[unassigned_tot];
	int	*alenlist = new int[unassigned_tot];
	int 	*mislist = new int[unassigned_tot];	
	int	*goodlist = new int[unassigned_tot];
	int	goodnum = 0;
	for(int i0 = 0; i0 < unassigned_tot; i0 ++) { //v1.1
		i = unassigned_idx[i0];
		alnlen = 0; //v1.1
		mis = align2seq(hash_table, conseq, len_conseq, src_seq[i], len_seq[i], alnspeedup, alnperc, aln_1, alnlen); //v1.1
		miscut1 = int(alnlen * alnperc + 0.5); //v1.1
		if(mis <= miscut1) {
			//printaln(conseq, len_conseq, src_seq[i], len_seq[i], aln_1);
			cov = double(mis) / alnlen; //v1.1
			assignedseq[k] = i;
			assignedcov[k ++] = cov;
		}
		//if(mis < 100) {
		if(mis <= 2 * miscut1) { //YY, May 27, 2015
			rechecklist[recheck] = i;
			alenlist[recheck] = alnlen;
			mislist[recheck] = mis;
			recheck ++;
		}
		if(mis < 2) {
			goodlist[goodnum ++] = i;
		}
	}
	int	assigned_read1 = k;
	for(i = 0; i < n_ban_recruit; i ++)	{
		if(hash_table[i])	{
			hash_table[i] = hashapp::free_table(hash_table[i]);
		}
	}
	delete[] hash_table;

	printf("  Recruited to read1-consensus %d, good ones %d\n", k, goodnum);
	fflush(stdout);
	if(assigned_read1 <= 0) {
		assignedtot = 0;
		delete[] rechecklist;
		delete[] alenlist;
		delete[] mislist;
		delete[] goodlist;
		return;
	}

	//read 2 -- consensus
	calign *ca = new calign;
	ca -> copyreads(assigned_read1, assignedseq, src_seq_pe, src_name_pe, len_seq_pe); //consensus for the paired ends using all pairs (not only good ones)
	char *cons_pe = ca -> oneconsensus(len_conseq_pe);
	for(int i0 = 0; i0 < len_conseq_pe; i0++) conseq_pe[i0] = cons_pe[i0];
	conseq_pe[len_conseq_pe] = 0;
	delete[] cons_pe;
	delete ca;

	//redo the recruitment using both ends
	hash_table_pe = new SMALLTABLE*[n_ban_recruit];
	for(i = 0; i < n_ban_recruit; i ++) hash_table_pe[i] = NULL;
	OneSeqHash(hash_table_pe, conseq_pe, len_conseq_pe);
	k = 0;
	for(int i0 = 0; i0 < recheck; i0 ++) { //v1.1
		i = rechecklist[i0];
		int mis_1 = mislist[i0];
		int alnlen_1 = alenlist[i0];
		int mis_2 = align2seq(hash_table_pe, conseq_pe, len_conseq_pe, src_seq_pe[i], len_seq_pe[i], alnspeedup, alnperc, aln_2, alnlen_2); //v1.1
		int mis_combined = mis_1 + mis_2;
		int alnlen_combined = alnlen_1 + alnlen_2;
		int miscut_combined = int(alnlen_combined * alnperc + 0.5);
		if(mis_combined <= miscut_combined) {
			//printaln(conseq, len_conseq, src_seq[i], len_seq[i], aln);
			cov = double(mis_combined) / alnlen_combined; //v1.1
			assignedseq[k] = i;
			assignedcov[k ++] = cov;
		}
	}
	assignedtot = k;
	printf("Reads recruited based on one read %d, both reads %d\n", assigned_read1, k);
	fflush(stdout);

	for(i = 0; i < n_ban_recruit; i ++)	{
		if(hash_table_pe[i])	{
			hash_table_pe[i] = hashapp::free_table(hash_table_pe[i]);
		}	
	}
	delete[] hash_table_pe;
	delete[] rechecklist;
	delete[] alenlist;
	delete[] mislist;
	delete[] goodlist;
	if(verbose) printf("end of map2consensus_pe\n");
	fflush(stdout);
}

//-----------------------------------------------------------------------------------------------
//recruitment using both ends; both consensens sequences are given for this function -- different from map2consensus_pe() 
void calign::map2consensus_pe2(int &len_conseq, char *conseq, int &len_conseq_pe, char *conseq_pe, int& assignedtot, int *assignedseq, double *assignedcov)
{
	hash_table = new SMALLTABLE*[n_ban_recruit];
	hash_table_pe = new SMALLTABLE*[n_ban_recruit];
	for(int i = 0; i < n_ban_recruit; i ++) hash_table[i] = hash_table_pe[i] = NULL;
	int	k = 0;
	vector<int> aln_1;
	vector<int> aln_2;
	int	alnspeedup = 1;
	int	alnlen_1, alnlen_2, alnlen, mis, miscut1;
	OneSeqHash(hash_table, conseq, len_conseq);
	OneSeqHash(hash_table_pe, conseq_pe, len_conseq_pe);
	for(int i0 = 0; i0 < unassigned_tot; i0 ++) { //v1.1
		int i = unassigned_idx[i0];
		//printf("i %d %d\n (out of %d)", i0, i, unassigned_tot);
		alnlen_1 = alnlen_2 = 0;
		int mis_1 = align2seq(hash_table, conseq, len_conseq, src_seq[i], len_seq[i], alnspeedup, alnperc, aln_1, alnlen_1); //v1.1
		int mis_2 = align2seq(hash_table_pe, conseq_pe, len_conseq_pe, src_seq_pe[i], len_seq_pe[i], alnspeedup, alnperc, aln_2, alnlen_2); //v1.1
		int mis_combined = mis_1 + mis_2;
		int alnlen_combined = alnlen_1 + alnlen_2;
		int miscut_combined = int(alnlen_combined * alnperc + 0.5);
		if(mis_combined <= miscut_combined) {
			//printaln(conseq, len_conseq, src_seq[i], len_seq[i], aln);
			cov = double(mis_combined) / alnlen_combined; //v1.1
			assignedseq[k] = i;
			assignedcov[k ++] = cov;
		}
	}
	assignedtot = k;
	printf("Reads recruited based on paired-end reads %d\n", k);

	for(int i = 0; i < n_ban_recruit; i ++)	{
		if(hash_table[i])	{
			hash_table[i] = hashapp::free_table(hash_table[i]);
		}
		if(hash_table_pe[i])	{
			hash_table_pe[i] = hashapp::free_table(hash_table_pe[i]);
		}	
	}
	delete[] hash_table;
	delete[] hash_table_pe;
}

void calign::printaln(char *seq1, int len1, char *seq2, int len2, vector<int> aln)
{
	char	*alnseq1 = new char[len1 + len2 + 1];
	char	*alnseq2 = new char[len1 + len2 + 1];
	char	*alnstr = new char[len1 + len2 + 1];
	alnseq1[0] = alnseq2[0] = alnstr[0]= 0;
	int	p1, p2, i, a;
	p1 = p2 = 0;
	//writeseq(seq1, "seq1", len1, stdout);
	//writeseq(seq2, "seq2", len2, stdout);
	for (i = 0; i < aln.size(); i ++) {
		a = aln[i];
		if(a == 0 || a == 1) {
			sprintf(alnseq1, "%s%c", alnseq1, na_name[seq1[p1 ++]] - 'a' + 'A');
			sprintf(alnseq2, "%s%c", alnseq2, na_name[seq2[p2 ++]] - 'a' + 'A');
			if(a == 0) strcat(alnstr, " ");
			else strcat(alnstr, "x");
		}
		else if(a == 2) {
			sprintf(alnseq1, "%s%c", alnseq1, na_name[seq1[p1 ++]] - 'a' + 'A');
			strcat(alnseq2, "-");
			strcat(alnstr, "x");
		}
		else {
			strcat(alnseq1, "-");
			sprintf(alnseq2, "%s%c", alnseq2, na_name[seq2[p2 ++]] - 'a' + 'A');
			strcat(alnstr, "x");
		}
	}
	printf("%-8s%s\n", "aln1", alnseq1);
	printf("%-8s%s\n", " ", alnstr);
	printf("%-8s%s\n", "aln2", alnseq2);
	delete[] alnseq1;
	delete[] alnseq2;
	delete[] alnstr;
}

//-----------------------------------------------------------------------------------------------
void calign::assignseq(int num_con, int tot, int *assignedseq, double *assignedcov)
{
	int     i0, i;
	for(i0 = 0; i0 < tot; i0 ++) {
		i = assignedseq[i0];
		index[i] = 1;
		seq2clust[i] = num_con;
		seq2clustdiff[i] = assignedcov[i0];
	}
}


//get the most frequent k-mer by using hash
int calign::hashcover(char **src_seq, int *len_seq, int num_seq, char *index, int *num_cover_list, long **ctl_list)
{
	//printf("in hashcover mark 1, unassigned_valid %d\n", unassigned_valid);
	fflush(stdout);
	SMALLTABLE **hash_table;
	long    i, ctl1, ctl2;
	hash_table = new SMALLTABLE*[n_ban];
	for(i = 0; i < n_ban; i ++) hash_table[i] = NULL;

	int	n = 0;
	for(int i0 = 0; i0 < unassigned_valid; i0 ++) { //v1.1
		i = unassigned_idx[i0];
		OneHash(hash_table, src_seq[i], len_seq[i]);
		n ++;
	}
	//printf("in hashcover mark 2, n = %d\n", n);
	fflush(stdout);

	int	howmany;
	if(top == 1) {
		num_cover_list[0] = MaxHash(hash_table, ctl_list[0], n_ban);
		howmany = 1;
		//printf("in hashcover mark 3-1\n");
		//fflush(stdout);
	}
	else if(num_con == 0 || mem_Conseq[num_con - 1] < 10) {
		howmany = SmartHash(hash_table, top, 0, coverthresh, n_ban, num_cover_list, ctl_list);
		//printf("in hashcover mark 3-2, num_con %d, mem_Conseq %d\n", num_con, mem_Conseq[num_con - 1]);
		//fflush(stdout);
	}
	else {
		howmany = SmartHash(hash_table, top, mem_Conseq[num_con - 1], coverthresh, n_ban, num_cover_list, ctl_list);
		//printf("in hashcover mark 3-3, num_con %d, mem_Conseq %d\n", num_con, mem_Conseq[num_con - 1]);
		//fflush(stdout);
	}

	for(i = 0; i < n_ban; i ++)	{
		if(hash_table[i])	{
			hash_table[i] = hashapp::free_table(hash_table[i]);
		}
	}
	//printf("in hashcover mark 4\n");
	//fflush(stdout);
	delete[] hash_table;

	//printf("in hashcover mark 5\n");
	//fflush(stdout);

	return howmany;
}

void calign::getunassigned(int ifspeedup)
{ //v1.1
	int     unassigned0 = 0;
	int     i;
	for(i = 0; i < num_seq; i ++) {
		if(index[i] == 0) unassigned0 += 1;
	}
	int     unassigned = 0;
	if(ifspeedup && (unassigned0 > seqnum_limit)) {
		int     sseq = start_idx;
		if(verbose) printf("Reshuffled, start at sseq %d\n", sseq);
		for(i = sseq; i < num_seq; i ++) {
			if(index[i] == 0) unassigned_idx[unassigned ++] = i;
		}
		for(i = 0; i < sseq; i ++) {
			if(index[i] == 0) unassigned_idx[unassigned ++] = i;
		}
		start_idx += seqnum_limit;
		if(start_idx > num_seq) start_idx -= num_seq;
	}
	else {
		for(i = 0; i < num_seq; i ++) {
			if(index[i] == 0) unassigned_idx[unassigned ++] = i;
		}
	}
	unassigned_tot = unassigned;
	if(ifspeedup && (unassigned > seqnum_limit)) unassigned_valid = seqnum_limit;
	else unassigned_valid = unassigned;
}

int calign::prefsufseq(int num_cover, long int* ctl, char **src_seq, int *len_seq, int num_seq, char *index,
						char **prefseq, int *len_prefseq, char **sufseq, int *len_sufseq, int *kmerseq, int *kmerpos)
{
	int	i, j, n, n0, m, k;
	long    hash_ban1, hash_ban2, ctl1, ctl2;
	int	seqwithmul = 0;
	//int	pos[100];
	n = n0 = 0;
	for(int i0 = 0; i0 < unassigned_valid; i0 ++) { //v1.1
		i = unassigned_idx[i0];
		hash_ban1 = smallapp::trans_seq(&src_seq[i][0], word_len);
		ctl1 = smallapp::trans_seq(&src_seq[i][word_len], kmer - word_len);
		int	findseed = 0; //April 21, 2011 to deal with multiple seeds in one seq--pick the first one
		if(hash_ban1 == ctl[0] && ctl1 == ctl[1])	{
			kmerseq[n] = i;
			if(n == 0)	{
				*kmerpos = 0;
			}
			prefseq[n] = new char[1];
			len_prefseq[n] = 0;
			len_sufseq[n] = len_seq[i] - kmer;
			sufseq[n] = new char[len_sufseq[n]];
			for(m = 0, k = kmer; k < len_seq[i]; k ++, m ++)	{
				sufseq[n][m] = src_seq[i][k];
			}
			n ++;
			n0 ++;
			//pos[findseed] = 0;
			findseed ++; //April 21, 2011 pick the first one
		}
		for(j = 1; j <= len_seq[i] - kmer; j ++)   {
			hash_ban1 = (hash_ban1 * 4 + src_seq[i][j + word_len - 1]) % n_ban;
			ctl1 = (ctl1 * 4 + src_seq[i][j + kmer - 1]) % n_ban1;
			if(hash_ban1 == ctl[0] && ctl1 == ctl[1])	{
				if(n == 0)	{
					*kmerpos = j;
				}
				if(findseed < 1) { //April 21, 2011 pick the first one
					kmerseq[n] = i; //Aug 8, 2011, by YY, move into loop!
					prefseq[n] = new char[j];
					len_prefseq[n] = j;
					for(m = j - 1, k = 0; k < j; k ++, m --)	{
						prefseq[n][k] = src_seq[i][m];
					}
					len_sufseq[n] = len_seq[i] - kmer - j;
					sufseq[n] = new char[len_sufseq[n]];
					for(m = 0, k = j + kmer; k < len_seq[i]; k ++, m ++)	{
						sufseq[n][m] = src_seq[i][k];
					}
					n ++;
				}
				//pos[findseed] = j;
				findseed ++;
				n0 ++; //April 21, 2011 pick the first one
			}
		}
		if(findseed > 1) {
			seqwithmul += 1;
			//printf("seq %d seed-locations %d\n", i, findseed);
			//for(j = 0; j < findseed; j ++) {
			//	printf("pos %d %d\n", j, pos[j]);
			//}
		}
	}
	//printf("n = %d, num_cover %d\n", n, num_cover);
	if(n0 != num_cover)	{
		printf("Coverage not equal: num_cover %d n %d\n", num_cover, n);
		exit(0);
	}
	if(n0 > n) {
		printf("!!!Warning: %d sequence(s) have multiple seeds\n", seqwithmul);
	}
	return n;
}

void calign::OneHash(SMALLTABLE **hash_table, char *src_seq, int len_seq)
{
	int     j;
	long    hash_ban1, hash_ban2;
	long    ctl1, ctl2;

	hash_ban1 = smallapp::trans_seq(&src_seq[0], word_len);
	ctl1 = smallapp::trans_seq(&src_seq[word_len], kmer - word_len);
	if(hash_ban1 >= LONG_MAX || ctl1 >= LONG_MAX) {
		printf("hash_ban1 or ctl1 exceed LONG_MAX\n");
		exit(0);
	}
	hash_table[hash_ban1] = hashapp::ins_larg_hash(hash_table[hash_ban1], ctl1);
	for(j = 1; j <= len_seq - kmer; j ++)   {
		hash_ban1 = (hash_ban1 * 4 + src_seq[j + word_len - 1]) % n_ban;
		ctl1 = (ctl1 * 4 + src_seq[j + kmer - 1]) % n_ban1;
		hash_table[hash_ban1] = hashapp::ins_larg_hash(hash_table[hash_ban1], ctl1);
	}
}

void calign::OneSeqHash(SMALLTABLE **hash_table, char *src_seq, int len_seq)
{
	int     i, j;
	long    hash_ban1, hash_ban2;
	long     ctl1, ctl2;

	hash_ban1 = smallapp::trans_seq(&src_seq[0], word_len);
	ctl1 = smallapp::trans_seq(&src_seq[word_len], kmer_recruit - word_len);
	hash_table[hash_ban1] = hashapp::ins_seq_hash(hash_table[hash_ban1], ctl1, 0);
	for(j = 1; j <= len_seq - kmer_recruit; j ++)   {
		hash_ban1 = (hash_ban1 * 4 + src_seq[j + word_len - 1]) % n_ban_recruit;
		ctl1 = (ctl1 * 4 + src_seq[j + kmer_recruit - 1]) % n_ban1_recruit;
		hash_table[hash_ban1] = hashapp::ins_seq_hash(hash_table[hash_ban1], ctl1, j);
	}
}

using namespace std;
void calign::diffstat(int num0, int id, int *diff2id, double *diff, int &min_i, int &max_i, int &median_i)
{
	vector<pair<double, int> > myvector;
	int	i, num, index;
	num = index = 0;
	for(i = 0; i < num0; i ++) {
		if(diff2id[i] != id) continue;	
		myvector.push_back(std::pair<double, int>(diff[i], i));
		num ++;
	}
	if(num == 0) {
		printf("warning: empty cluster\n");
		exit(0);
	}
	sort(myvector.begin(), myvector.end()); 
	//sort by first element, i.e., the num of sequences included in each consensus

	min_i = myvector[0].second;
	max_i = myvector[num - 1].second;
	median_i = myvector[num / 2].second;

	myvector.clear();
}

void calign::sortcons(void)
{
	vector<pair<int, int> > myvector;
	int	index = 0;
	int	i;
	for(i = 0; i < num_con; i ++) {
		myvector.push_back(std::pair<int, int>(mem_Conseq[i], index++));
	}

	sort(myvector.begin(), myvector.end()); 
	//sort by first element, i.e., the num of sequences included in each consensus
	if(con_sort != NULL) {
		delete[] con_sort;
		delete[] con2sort;
	}
	con_sort = new int[num_con];
	con2sort = new int[num_con];
	index = num_con - 1;
	for(vector<pair<int, int> >::iterator it = myvector.begin(); it != myvector.end(); it ++) {
		con2sort[(*it).second] = index;
		con_sort[index --] = (*it).second;
	} //rbegin: reverse begin
	myvector.clear();
}

//"release" small clusters 
void calign::removesmall(int small_clust_size, bool ifempty)
{
	num_con_valid = 0;

	int     i, j;
	for(i = 0; i < num_con; i ++) {
		if(mem_Conseq[i] < small_clust_size) {
			for(j = 0; j < num_seq; j ++) {
				if(seq2clust[j] == i) index[j] = 0;
			}
			totseq += mem_Conseq[i];
			if(ifempty) mem_Conseq[i] = 0; //to be checked!!
		}
		else {
			num_con_valid ++;
		}
	}
	//printf("cluster %d, valid %d\n", num_con, num_con_valid);
}

void calign::saveresults(void)
{
	printf("start writing reports\n");

	writeclust();
	writeclustsize();

	if(paired) writeconseq_pe();
	else writeconseq();

	if((abundantonly == 1) && (totseq > 0)) {
		writesingleton();
	}
}

char* calign::int2seq(char *src_seq, int len)
{
	char    *seq = new char[len + 1];
	for(int i = 0; i < len; i ++)
		seq[i] = na_name[src_seq[i]] - 'a' + 'A';
	seq[len] = 0;
	return seq;
}

char* calign::seq2int(char *src_seq, int len)
{
	char    *seq = new char[len + 1];
	for(int i = 0; i < len; i ++)
		seq[i] = smallapp::char2int(src_seq[i]);
	seq[len] = 0;
	return seq;
}

void calign::writeseq(char *src_seq, char *src_name, int len, FILE *fp)
{
	int	i, j, k;

	fprintf(fp, ">%s\n", src_name);
	for(i = 0; i < len; i ++)	{
		fprintf(fp, "%c", na_name[src_seq[i]] - 'a' + 'A');
		if(i % 60 == 59)	{
			fprintf(fp, "\n");
		}
	}
	if(i % 60 != 0)	{
		fprintf(fp, "\n");
	}
}

void calign::writesingleton(void)
{
	FILE	*fp;
	char	temp[1000];
	sprintf(temp, "%s.single", fileprefix);
	fp = smallapp::ckopen(temp, "w");

	int	thisnumseq = num_seq;

	int	i, j, k;
	int	add = 0;
	for(i = 0; i < thisnumseq; i ++)	{
		if(index[i] == 1)	{
			continue;
		}
		add ++;
		writeseq(src_seq[i], src_name[i], len_seq[i], fp);
	}
	fclose(fp);

	printf("%d singleton (unclustered) sequences saved to %s\n", add, temp);
}

void calign::writeconseq(void)
{
	int	i0, i, j, k; 
	char	temp[1000], str[1000];
	sprintf(temp, "%s.cons", fileprefix);
	FILE	*fp = smallapp::ckopen(temp, "w");

	for(i0 = 0; i0 < num_con_valid; i0 ++)	{
		i = con_sort[i0];
		sprintf(str, "Consensus%d len=%d tot-seq=%d", i0 + 1, len_Conseq[i], mem_Conseq[i]);
		writeseq(Conseq[i], str, len_Conseq[i], fp);
	}
	fclose(fp);

	printf("%d consensus sequences saved to %s\n", num_con_valid, temp);
}

void calign::writeconseq_pe(void)
{
	int	i0, i, j, k; 
	char	temp[1000], temp_pe[1000], str[1000];
	if(ifswap) {
		sprintf(temp, "%s.1.cons", fileprefix);
		sprintf(temp_pe, "%s.2.cons", fileprefix);
	}
	else {
		sprintf(temp, "%s.2.cons", fileprefix);
		sprintf(temp_pe, "%s.1.cons", fileprefix);
	}
	printf("write consensus ..\n");
	FILE	*fp = smallapp::ckopen(temp, "w");
	FILE	*fp_pe = smallapp::ckopen(temp_pe, "w");
	for(i0 = 0; i0 < num_con_valid; i0 ++)	{
		i = con_sort[i0];
		sprintf(str, "Consensus%d len=%d tot-seq=%d", i0 + 1, len_Conseq[i], mem_Conseq[i]);
		writeseq(Conseq[i], str, len_Conseq[i], fp);
		sprintf(str, "Consensus%d len=%d tot-seq=%d", i0 + 1, len_Conseq_pe[i], mem_Conseq[i]);
		writeseq(Conseq_pe[i], str, len_Conseq_pe[i], fp_pe);
	}
	fclose(fp);
	fclose(fp_pe);

	printf("%d consensus sequences saved to %s\n", num_con_valid, temp);
}

void calign::writeclust(void)
{
	char	temp[1000];
	sprintf(temp, "%s.clust", fileprefix);
	FILE	*fp = smallapp::ckopen(temp, "w");

	int	i0, i, j, min_i, max_i, median_i;
	for(i0 = 0; i0 < num_con_valid; i0 ++) {
		i = con_sort[i0];
		diffstat(num_seq, i, seq2clust, seq2clustdiff, min_i, max_i, median_i);
		fprintf(fp, "#Consensus %d include seq=%d (%s:%.1f%% %s:%.1f%% %s:%.1f%%)\n", 
			i0 + 1, mem_Conseq[i], src_name[min_i], seq2clustdiff[min_i]*100, src_name[max_i], 
			seq2clustdiff[max_i]*100, src_name[median_i], seq2clustdiff[median_i]*100);
		for(j = 0; j < num_seq; j ++) {
			if(seq2clust[j] != i)	continue;
			if(paired) {
				if(ifswap) fprintf(fp, "  seq %d %s len:%d:%d diff:%.1f%%\n", j + 1, src_name[j], len_seq_pe[j], len_seq[j], seq2clustdiff[j] * 100);
				else fprintf(fp, "  seq %d %s len:%d:%d diff:%.1f%%\n", j + 1, src_name[j], len_seq[j], len_seq_pe[j], seq2clustdiff[j] * 100);
			}
			else fprintf(fp, "  seq %d %s len:%d diff:%.1f%%\n", j + 1, src_name[j], len_seq[j], seq2clustdiff[j] * 100);
		}
	}

	fclose(fp);
}

void calign::writelog(int conid, int conlen, char *conseq, int assignedtot, int *assigned)
{
	char str[1000];	
	sprintf(str, "Consensus%d len=%d tot-seq=%d", conid + 1, conlen, assignedtot);
	writeseq(conseq, str, conlen, log);
	int	i0, i;
	for(i0 = 0; i0 < assignedtot; i0 ++) {
		i = assigned[i0];
		fprintf(log, " - %s\n", src_name[i]);	
	}
	fflush(log);
}

void calign::writelog_pe(int conid, int conlen, char *conseq, int conlen_pe, char *conseq_pe, int assignedtot, int *assigned)
{
	char str[1000];	
	sprintf(str, "Consensus%d-R1 len=%d tot-seq=%d", conid + 1, conlen, assignedtot);
	writeseq(conseq, str, conlen, log);
	sprintf(str, "Consensus%d-R2 len=%d tot-seq=%d", conid + 1, conlen_pe, assignedtot);
	writeseq(conseq_pe, str, conlen_pe, log);
	int	i0, i;
	for(i0 = 0; i0 < assignedtot; i0 ++) {
		i = assigned[i0];
		//fprintf(log, " i0 %d i %d\n", i0, i);	
		fprintf(log, " - %s\n", src_name[i]);	
	}
	fflush(log);
	if(verbose) printf("end of writelog_pe\n");
}

void calign::writeclustsize(void)
{
	char	temp[1000], tmp[100];
	sprintf(temp, "%s.clustsize", fileprefix);
	FILE	*fp = smallapp::ckopen(temp, "w");

	int	i0, i, j;
	if(paired) fprintf(fp, "#id description length length-pe include-seq accumulated accumulated%\n");
	else fprintf(fp, "#id description length include-seq accumulated accumulated%\n");
	int	add = 0;
	for(i0 = 0; i0 < num_con_valid; i0 ++) {
		i = con_sort[i0];
		add += mem_Conseq[i];
		sprintf(tmp, "Consensus%d", i0 + 1);
		if(paired) {
			if(ifswap) fprintf(fp, "%-4d %-14s %-4d %-4d %-6d %-6d %8.2f%%\n", i0 + 1, tmp, len_Conseq_pe[i], len_Conseq[i], mem_Conseq[i], add, add * 100.0 / num_seq);
			else fprintf(fp, "%-4d %-14s %-4d %-4d %-6d %-6d %8.2f%%\n", i0 + 1, tmp, len_Conseq[i], len_Conseq_pe[i], mem_Conseq[i], add, add * 100.0 / num_seq);
		}
		else fprintf(fp, "%-4d %-14s %-4d %-6d %-6d %8.2f%%\n", i0 + 1, tmp, len_Conseq[i], mem_Conseq[i], add, add * 100.0 / num_seq);
	}

	fclose(fp);
	printf("Cluster size saved to file %s\n", temp);
}

//subotu added by YY
int calign::getclustnum(void)
{
	return num_con_valid;
}

int calign::getclustsize(int clust)
{
	return mem_Conseq[con_sort[clust]];
}

int calign::getclustinfo(int clust0, int *mem, double *dif, char *con)
{
	int	clust = con_sort[clust0]; //note input index: unsorted 
	int	add = 0;
	int	i;
	for(i = 0; i < num_seq; i ++) {
		if(seq2clust[i] == clust) {
			dif[add] = seq2clustdiff[i];
			mem[add] = i;
			add ++;
		}
	}
	for(i = 0; i < len_Conseq[clust]; i ++) {
		con[i] = Conseq[clust][i];
	}
	con[i] = 0;
	return len_Conseq[clust];
}

void calign::subotu(int minclustsize, char *fileprefix)
{
	char	temp[1000], temp2[1000], str[100];
	sprintf(temp, "%s.subotu-clust", fileprefix);
	FILE	*fp = smallapp::ckopen(temp, "w");
	sprintf(temp2, "%s.subotu-cons", fileprefix);
	FILE	*fp2 = smallapp::ckopen(temp2, "w");

	int	i0, i, j, k, min_i, max_i, median_i, index_ori;
	char	*conseq;
	for(i0 = 0; i0 < num_con_valid; i0 ++) {
		i = con_sort[i0];
		if(mem_Conseq[i] < minclustsize) break;
		diffstat(num_seq, i, seq2clust, seq2clustdiff, min_i, max_i, median_i);
		fprintf(fp, "#Consensus %d include seq=%d (%s:%.1f%% %s:%.1f%% %s:%.1f%%)\n", 
			i0 + 1, mem_Conseq[i], src_name[min_i], seq2clustdiff[min_i]*100, src_name[max_i], 
			seq2clustdiff[max_i]*100, src_name[median_i], seq2clustdiff[median_i]*100);
		sprintf(str, "Consensus%d len=%d tot-seq=%d", i0 + 1, len_Conseq[i], mem_Conseq[i]);
		writeseq(Conseq[i], str, len_Conseq[i], fp2);
		int	*assignedseq = new int[mem_Conseq[i]];
		int	*assigned2sub = new int[mem_Conseq[i]];
		printf("\n\nNow infer subotu %d with seq %d\n", i0, mem_Conseq[i]);
		int	add = 0;
		for(j = 0; j < num_seq; j ++) {
			if(seq2clust[j] == i)	{
				assigned2sub[add] = -1;
				assignedseq[add ++] = j;
			}
		}
                calign *casub = new calign;
		casub -> copyreads(add, assignedseq, src_seq, src_name, len_seq);
		//casub -> ca2otu();
		casub -> set_alnperc(0.01);
		casub -> iteralign();
		casub -> sortcons();
		casub -> removesmall(0, true);
		int	clustnum = casub -> getclustnum();
		printf("clust %d %d seq %d, subotu %d\n", i0, i, add, clustnum);
		int	assigned2subtot = 0;
		int	maxsize = 0;
		int	lastclustsize = 0;
		for(j = 0; j < clustnum; j ++) {
			int	clustsize = casub -> getclustsize(j);
			if(j == 0) { 
				maxsize = clustsize;
			}
			//else if((clustsize < maxsize * 0.05) || clustsize <= 5) {
			//else if((clustsize < lastclustsize * 0.1) || clustsize <= 5) {
			else if(clustsize < lastclustsize * 0.01 || clustsize < 10) {
				break;
			}
			printf("subotu %d size %d\n", j, clustsize);
			int 	*clustmem = new int[clustsize];
			double	*seqdif = new double[clustsize];
			char	*conseq = new char[len_Conseq[i] * 2];
			int	conlen = casub -> getclustinfo(j, clustmem, seqdif, conseq);
			for(k = 0; k < clustsize; k ++) {
				assigned2sub[clustmem[k]] = j;
			}
			assigned2subtot += clustsize;
			fprintf(fp, "#subotu %d.%d include seq=%d\n", i0 + 1, j + 1, clustsize);
			for(k = 0; k < clustsize; k ++) {
				index_ori = assignedseq[clustmem[k]];
				fprintf(fp, "  seq %d %s len:%d diff:%.1f%% diff-subotu:%.1f%%\n", index_ori + 1, 
					src_name[index_ori], len_seq[index_ori], seq2clustdiff[index_ori] * 100, seqdif[k] * 100);
			}
			sprintf(str, "subotu%d.%d len=%d tot-seq=%d", i0 + 1, j + 1, conlen, clustsize);
			writeseq(conseq, str, conlen, fp2);
			delete[] clustmem;
			delete[] seqdif;
			delete[] conseq;
			printf("  sub-otu %d include %d\n", j, clustsize);
			lastclustsize = clustsize;
		}
		if(assigned2subtot < mem_Conseq[i]) {
			fprintf(fp, "#subotu %d.- include seq=%d\n", i0 + 1, mem_Conseq[i] - assigned2subtot); 
			for(k = 0; k < add; k ++) {
				if(assigned2sub[k] != -1) continue;
				index_ori = assignedseq[k];
				fprintf(fp, "  seq %d %s len:%d diff:%.1f%%\n", index_ori + 1, 
					src_name[index_ori], len_seq[index_ori], seq2clustdiff[index_ori] * 100);
			}
		}
		fflush(stdout);
		delete casub;
		delete[] assignedseq;
		delete[] assigned2sub;
	}

	fclose(fp);
	printf("subotu saved to %s and %s\n", temp, temp2);
}

void calign::setfileprefix(char *tmp) {
	strcpy(fileprefix, tmp);
}

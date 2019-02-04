//Program Name: alignapp
//Latest modification on Oct 29th, 2010
//Developer: Yuzhen Ye <yye@indiana.edu>
//Affiliation: School of Informatics and Computing, Indiana University, Bloomington

#include "alignapp.h"
#include <limits.h>

using namespace hashapp;

alignapp::alignapp(void)
{
	initenv();
}

void alignapp::initenv(void)
{
	kmer = 40;
	bandwidth = 5;
	strandtag = 0;
	alnperc = 0.03;
	penalty = 0.5;
	mgsize = 0;
	word_len = 12;

	abundantonly = 0;
	//Score for match, mismatch and indel are all 1 -- edit distance is used 
	mu1 = 1;
	mu2 = 1;
	gamma0 = 3;

	idum = 0; //Bug reported by Raad, Feb 8, 2013 
	smallapp::random1(&idum);

	num_con = 0;
	compute_n_ban_recruit();
	compute_n_ban();

	verbose = 0;
}

void alignapp::compute_n_ban_recruit(void)
{
	int	i;
	kmer_recruit = 20;
	n_ban_recruit = 1;
	for(i = 0; i < word_len; i ++)  {
		if(n_ban_recruit > LONG_MAX / 4) {
			printf("n_ban larger than LONG_MAX; stop\n");
			exit(0);
		}
		n_ban_recruit *= 4;
	}
	n_ban1_recruit = 1;
	for(i = 0; i < kmer_recruit - word_len; i ++)   {
		if(n_ban1_recruit > LONG_MAX / 4) {
			printf("n_ban1_recruit larger than LONG_MAX; stop\n");
			exit(0);
		}
		n_ban1_recruit *= 4;
	}
}

void alignapp::compute_n_ban(void)
{
	int	i;
	n_ban = 1;
	for(i = 0; i < word_len; i ++)  {
		if(n_ban > LONG_MAX / 4) {
			printf("n_ban larger than LONG_MAX; stop\n");
			exit(0);
		}
		n_ban *= 4;
	}
	n_ban1 = 1;
	for(i = 0; i < kmer - word_len; i ++)   {
		if(n_ban1 > LONG_MAX / 4) {
			printf("n_ban1 larger than LONG_MAX; stop\n");
			exit(0);
		}
		n_ban1 *= 4;
	}
}

int alignapp::CAlign(char *consseq, char **seq, int *len_seq, int num_seq, int *num_mismatch, int maxconlen, int mingsize, int ic0, double PW)
{
	int	i, j, k, l, m, n, m1, q;
	int	len, ic, num;
	int	score, ovscore, conscore, **maxscore, totscore, lastscore, currscore;
	int	**scoreline, **scoretmp, **mutcount, *mscore, **muttmp, *num_mm, *num_mm0;

	mscore = new int[num_seq];
	num_mm = new int[num_seq];
	num_mm0 = new int[num_seq];
	scoretmp = new int*[num_seq];
	muttmp = new int*[num_seq];
	scoreline = new int*[num_seq];
	maxscore = new int*[num_seq];
	mutcount = new int*[num_seq];
	for(k = 0; k < num_seq; k ++)	{
		scoretmp[k] = new int[8 * bandwidth + 4];
		muttmp[k] = new int[8 * bandwidth + 4];
		for(j = 0; j < 8 * bandwidth + 4; j ++) scoretmp[k][j] = muttmp[k][j] = 0;
		scoreline[k] = new int[2 * bandwidth + 1];
		mutcount[k] = new int[2 * bandwidth + 1];
		for(j = 0; j < 2 * bandwidth + 1; j ++) scoreline[k][j] = mutcount[k][j] = 0;
		maxscore[k] = new int[4];
		for(j = 0; j < 4; j ++) maxscore[k][j] = 0;
		num_mm0[k] = num_mismatch[k];
		num_mm[k] = 0; //YY, Feb 8, 2013
	}

	//Initialization for scoreline
	for(k = 0; k < num_seq; k ++)	{
		for(m = 0; m <= 2 * bandwidth; m ++)	{
			j = m - bandwidth;
			if(j >= 0)	{
				scoreline[k][m] = -gamma0 * j;
				mutcount[k][m] = j;
			} else	{
				scoreline[k][m] = gamma0 * j;
				mutcount[k][m] = -j;
			}
		}
	}


	//printf("PW %f\n", PW);
	double	maxscoreper = -LARGENUMBER;
	int	maxscorenum = 0;
	len = 0;
	ovscore = 0;
	num = num_seq;
	for(i = 0; i < maxconlen; i ++)	{
		ic = (int) (((double) i + kmer + ic0) * PW);
		conscore = -LARGENUMBER;
		for(j = 0; j < 4; j ++)	{
			totscore = 0;
			for(k = 0; k < num_seq; k ++)	{
				score = 0;
				lastscore = -LARGENUMBER;
				maxscore[k][j] = -LARGENUMBER;
				for(m1 = max(-i, -bandwidth), m = m1 + bandwidth; m1 <= min(len_seq[k] - 1 - i, bandwidth); m ++, m1 ++)	{
					currscore = scoreline[k][m];
					// Replacement/matching at position i
					if(j == seq[k][i + m1])	{
						score = currscore + mu1;
						muttmp[k][4 * m + j] = mutcount[k][m];
					} else	{
						score = currscore - mu2;
						muttmp[k][4 * m + j] = mutcount[k][m] + 1;
					}
					if(score < lastscore - gamma0)	{
						//Deletion in the reference sequence k at position i + m1
						score = lastscore - gamma0;
						muttmp[k][4 * m + j] = mutcount[k][m - 1] + 1;
					} else if(m1 < bandwidth && i + m1 < len_seq[k] - 1 && score < scoreline[k][m + 1] - gamma0)	{
						//Deletion in the consensus sequence at position i
						score = scoreline[k][m + 1] - gamma0;
						muttmp[k][4 * m + j] = mutcount[k][m + 1] + 1;
					}
					scoretmp[k][4 * m + j] = score;
					if(score > maxscore[k][j])	{
						maxscore[k][j] = score;
					}
					lastscore = score;
				}

				//The extension of the consensus sequence goes beyond sequence k, then use last score stored for sequence k
				if(maxscore[k][j] == -LARGENUMBER)	{
					maxscore[k][j] = 0;
				} else {
					mscore[k] = maxscore[k][j];
				}

				if((num_mm[k] + num_mm0[k] <= ic) && (i - len_seq[k] <= bandwidth))	{
					totscore += maxscore[k][j];
				}
			}
			if(totscore > conscore)	{
				conscore = totscore;
				consseq[i] = j;
			}
		}

		num = 0;
		double	scorenew = 0;
		for(k = 0; k < num_seq; k ++)   {
			score = -LARGENUMBER;
			for(m1 = max(-i, -bandwidth), m = m1 + bandwidth; m1 <= min(len_seq[k] - 1 - i, bandwidth); m ++, m1 ++)	{
				scoreline[k][m] = scoretmp[k][4 * m + consseq[i]];
				mutcount[k][m] = muttmp[k][4 * m + consseq[i]];
				if(scoreline[k][m] > score)	{
					score = scoreline[k][m];
					num_mm[k] = muttmp[k][4 * m + consseq[i]];
				}
			}
			if(num_mm[k] + num_mm0[k] <= ic && i - len_seq[k] <= bandwidth)	{
				num ++;
				scorenew += score;
			}
		}

		//if(conscore <= 0 && num == 0) break; //to-be-tested YY April 24, 2014
		if(conscore <= 0 && (!(num >= 5 || num > 0.5 * num_seq))) break; //to-be-tested YY April 24, 2014

		if((scorenew > 0)&& (num > mingsize) && ((scorenew / num) > maxscoreper) && (num > maxscorenum / 10)) { //to-be-tested, like local alignment
			ovscore = conscore;
			//len = i;
			len = i + 1;  //v1.1
			for(k = 0; k < num_seq; k ++)   {
				num_mismatch[k] = num_mm[k];
			}
			maxscoreper = scorenew / num;
			maxscorenum = num;
		}
	}

	for(k = 0; k < num_seq; k ++)	{
		delete[] scoretmp[k];
		delete[] muttmp[k];
		delete[] scoreline[k];
		delete[] maxscore[k];
		delete[] mutcount[k];
	}
	delete[]  muttmp;
	delete[]  scoretmp;
	delete[]  scoreline;
	delete[] mscore;
	delete[] num_mm;
	delete[] num_mm0;
	delete[]  maxscore;
	delete[]  mutcount;

	return len;
}

int alignapp::align2seq(SMALLTABLE **hash_table, char *Conseq, int len_Conseq, char *src_seq, int len_seq, int speedup, double alnperc, vector<int> &aln, int &alnlen)
{
	int	i, j, k, l, m, n, max_i, max_j;
	SMALLTABLE *hash_temp;
	int     hash_ban1, hash_ban2;
	long     ctl1, ctl2;
	int	*count;
	double	id;
	alnlen = 0;
	count = new int[len_seq + len_Conseq];
	for(i = 0; i < len_seq + len_Conseq; i ++) count[i] = 0;

	hash_ban1 = smallapp::trans_seq(&src_seq[0], word_len);
	ctl1 = smallapp::trans_seq(&src_seq[word_len], kmer_recruit - word_len);
	hash_temp = FindHash(hash_table[hash_ban1], ctl1);
	if(hash_temp)	{
		m = hash_temp -> pos;
		count[m + len_seq] ++;
	}
	for(j = 1; j <= len_seq - kmer_recruit; j ++)   {
		hash_ban1 = (hash_ban1 * 4 + src_seq[j + word_len - 1]) % n_ban_recruit;
		ctl1 = (ctl1 * 4 + src_seq[j + kmer_recruit - 1]) % n_ban1_recruit;
		hash_temp = FindHash(hash_table[hash_ban1], ctl1);
		if(hash_temp)	{
			//if(hash_temp -> pos > len_Conseq) {
			//	printf("Error: hash_temp pos %d vs %d\n", hash_temp -> pos, len_Conseq);
			//	exit(0);
			//}
			m = hash_temp -> pos - j;
			count[m + len_seq] ++;
		}
	}
	max_j = 0;
	int	tot_count = 0; //v1.1
	for(i = 0; i < len_seq + len_Conseq; i ++)	{
		if(max_j < count[i])	{
			max_i = i;
			max_j = count[i];
		}
		tot_count += count[i];
	}
	delete[] count;

	m = max_i - len_seq;

	int	tag = 0;
	if(speedup) {//v1.1
		int     slen = len_Conseq < len_seq?len_Conseq:len_seq;
		int     miscut = int(slen * alnperc + 0.5);
		int     min_diag0 = int((slen * 1.0 / (miscut + 1) - kmer_recruit));
		int     min_tot = min_diag0 * (miscut + 1);
		int     min_diag = (miscut + 1) > 4?min_diag0 * 2:min_diag0;
		if(max_j < min_diag || tot_count < min_tot) tag = 1;
	}
	if(max_j == 0 || tag) {
		k = 100;
	}
	else if(m >= 0)	{
		k = global_align_eff(Conseq + m, src_seq, len_Conseq - m, len_seq, bandwidth, aln, alnlen);
		for(i = 0; i < m; i ++) aln.insert(aln.begin(), 2); //deletion (of the second seq) at the beginning
	} else	{ //m < 0
		k = global_align_eff(Conseq, src_seq - m, len_Conseq, len_seq + m, bandwidth, aln, alnlen);
		for(i = 0; i < -m; i ++) aln.insert(aln.begin(), 3); //insertion (of the second seq) at the beginning
	}

	return(k);
}

//Compute edit distance between sequence 1 and 2, with no end gap penalt
int alignapp::global_align(char *seq1, char *seq2, int len1, int len2, int b, vector<int> &aln, int &alnlen)
{
	int	i, j, k, l, m, n, b2;
	int	*score, curr_score, lastscore, bestscore;
	int	bestscore_i, bestscore_k;
	int	**srecord = new int*[len1];
	for(i = 0; i < len1; i ++) {
		srecord[i] = new int[len2];
		for(j = 0; j < len2; j ++) srecord[i][j] = LARGENUMBER;
	}
	b2 = b * 2;
	score = new int[b2 + 1];
	for(i = 0; i <= b2; i ++)	{
		score[i] = 0;
	}
	bestscore = LARGENUMBER;
	for(i = 0; i < len1; i ++)	{
		if(i <= b)	{
			lastscore = 0;
		} else	{
			lastscore = LARGENUMBER;
		}
		for(k = max(0, i - b), j = k - i + b; k <= min(i + b, len2 - 1); j ++, k ++)	{
			curr_score = lastscore + 1;
			//printf("seq1[%d] %c seq2[%d] %c\n", i, seq1[i], k, seq2[k]);
			if(seq1[i] == seq2[k])	{
				if(score[j] < curr_score)	{
					curr_score = score[j];
				}
			} else	{
				if(score[j] + 1 < curr_score)	{
					curr_score = score[j] + 1;
				}
			}
			if(j < b2)	{
				if(score[j + 1] + 1 < curr_score)	{
					curr_score = score[j + 1] + 1;
				}
			}
			score[j] = curr_score;
			lastscore = curr_score;
			srecord[i][k] = curr_score;
			if(curr_score < bestscore && (i == len1 - 1))	{
				bestscore = curr_score;
				bestscore_i = i;
				bestscore_k = k;
			}
		}
		if(curr_score < bestscore && (k == len2))	{
			bestscore = curr_score;
			bestscore_i = i;
			bestscore_k = k - 1;
		}
	}
	delete[] score;

	//get alignment
	bool	ifcheck = true;
	if(ifcheck && (bestscore != LARGENUMBER)) {
		aln.clear();
		int curr_i = bestscore_i;
		int curr_k = bestscore_k;
		int	checkscore = 0;
		int	mod;
		alnlen = 0; //v1.1
		while(curr_i >= 0 && curr_k >= 0) {
			mod = 0;
			if((curr_i > 0) && (srecord[curr_i][curr_k] == srecord[curr_i - 1][curr_k] + 1)) mod = 2; //deletion
			else if((curr_k > 0) && (srecord[curr_i][curr_k] == srecord[curr_i][curr_k - 1] + 1)) mod = 3; //insertion
			else { 
				if(seq1[curr_i] == seq2[curr_k]) mod = 0; //match
				else mod = 1; //mismatch
			}
			if(mod != 0) 	checkscore += 1;
			if(mod != 2)	curr_k -= 1; 
			if(mod != 3) 	curr_i -= 1;
			aln.insert(aln.begin(), mod);
			alnlen += 1;
		}
		while(curr_k >= 0) {
			curr_k -= 1;
			aln.insert(aln.begin(), 3);
			alnlen += 1;
		}
		while(curr_i >= 0) {
			curr_i -= 1;
			aln.insert(aln.begin(), 2);
			alnlen += 1;
		}

		if(checkscore != bestscore) {
			printf("len1 %d len2 %d, bestscore (%d %d) score %d\n", len1, len2, bestscore_i, bestscore_k, bestscore);
			printf("!!Warning: bestscore do not match %d vs %d\n", checkscore, bestscore);
			//exit(0);
		}
	}

	for(i = 0; i < len1; i ++) {
		delete[] srecord[i];
	}
	delete[] srecord;

	return(bestscore);
}


//Compute edit distance between sequence 1 and 2, with no end gap penalt, memory efficient
int alignapp::global_align_eff(char *seq1, char *seq2, int len1, int len2, int b, vector<int> &aln, int &alnlen)
{
	int	i, j, k, l, m, n, b2;
	int	*score, curr_score, lastscore, bestscore;
	int	bestscore_i, bestscore_k;
	int	**srecord = new int*[len1];
	b2 = b * 2;
	for(i = 0; i < len1; i ++) {
		//srecord[i] = new int[len2];
		//for(j = 0; j < len2; j ++) srecord[i][j] = 1000;
		srecord[i] = new int[b2 + 1];
		for(j = 0; j < b2 + 1; j ++) {
			srecord[i][j] = LARGENUMBER;
		}
	}
	score = new int[b2 + 1];
	for(i = 0; i <= b2; i ++)	{
		score[i] = 0;
	}
	bestscore = LARGENUMBER;
	for(i = 0; i < len1; i ++)	{
		if(i <= b)	{
			lastscore = 0;
		} else	{
			lastscore = LARGENUMBER;
		}
		for(k = max(0, i - b), j = k - i + b; k <= min(i + b, len2 - 1); j ++, k ++)	{
			curr_score = lastscore + 1;
			//printf("seq1[%d] %c seq2[%d] %c\n", i, seq1[i], k, seq2[k]);
			if(seq1[i] == seq2[k])	{
				if(score[j] < curr_score)	{
					curr_score = score[j];
				}
			} else	{
				if(score[j] + 1 < curr_score)	{
					curr_score = score[j] + 1;
				}
			}
			if(j < b2)	{
				if(score[j + 1] + 1 < curr_score)	{
					curr_score = score[j + 1] + 1;
				}
			}
			score[j] = curr_score;
			lastscore = curr_score;
			//srecord[i][k] = curr_score;
			srecord[i][k - i + b] = curr_score;
			if(curr_score < bestscore && (i == len1 - 1))	{
				bestscore = curr_score;
				bestscore_i = i;
				bestscore_k = k;
			}
		}
		if(curr_score < bestscore && (k == len2))	{
			bestscore = curr_score;
			bestscore_i = i;
			bestscore_k = k - 1;
		}
	}
	delete[] score;

	//get alignment
	bool ifcheck = true;
	if(ifcheck && (bestscore != LARGENUMBER)) {
		aln.clear();
		int curr_i = bestscore_i;
		int curr_k = bestscore_k;
		int	checkscore = 0;
		int	mod;
		alnlen = 0; //v1.1
		while(curr_i >= 0 && curr_k >= 0) {
			//if((curr_i > 0) && (srecord[curr_i][curr_k] == srecord[curr_i - 1][curr_k] + 1)) mod = 2; //deletion
			if((curr_i > 0) && (curr_k <= curr_i - 1 + b) && (curr_k >= curr_i - 1 - b) && (srecord[curr_i][curr_k - curr_i + b] == srecord[curr_i - 1][curr_k - (curr_i - 1) + b] + 1)) mod = 2; //deletion
			//else if((curr_k > 0) && (srecord[curr_i][curr_k] == srecord[curr_i][curr_k - 1] + 1)) mod = 3; //insertion
			else if((curr_k > 0) && (curr_k - 1 <= curr_i + b) && (curr_k - 1 >= curr_i - b) && (srecord[curr_i][curr_k - curr_i + b] == srecord[curr_i][curr_k - 1 - curr_i + b] + 1)) mod = 3; //insertion
			//note check bandwidth here, fixed June 6, 2011, Ye
			else { 
				if(seq1[curr_i] == seq2[curr_k]) mod = 0; //match
				else mod = 1; //mismatch
			}
			if(mod != 0) 	checkscore += 1;
			if(mod != 2)	curr_k -= 1; 
			if(mod != 3) 	curr_i -= 1;
			aln.insert(aln.begin(), mod);
			alnlen += 1;
		}
		while(curr_k >= 0) {
			curr_k -= 1;
			aln.insert(aln.begin(), 3);
			alnlen += 1;
		}
		while(curr_i >= 0) {
			curr_i -= 1;
			aln.insert(aln.begin(), 2);
			alnlen += 1;
		}

		if(checkscore != bestscore) {
			printf("len1 %d len2 %d, bestscore (%d %d) score %d\n", len1, len2, bestscore_i, bestscore_k, bestscore);
			printf("!!Warning: bestscore do not match %d vs %d\n", checkscore, bestscore);
			//exit(0);
		}
	}

	for(i = 0; i < len1; i ++) {
		delete[] srecord[i];
	}
	delete[] srecord;

	return(bestscore);
}

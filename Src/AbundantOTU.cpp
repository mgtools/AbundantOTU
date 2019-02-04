//Program Name: AbundantOTU
//Latest modification on April 22, 2014
//Yuzhen Ye <yye@indiana.edu> and Yongan Zhao <yongzhao@indiana.edu>
//Affiliation: School of Informatics and Computing, Indiana University, Bloomington
#include "calign.h"
#include <time.h>

void printusage(char *errorinfo);

main(int argc, char **argv)
{
	char	inpfile[100], inpfile_pe[100], outfile[100];
	inpfile[0] = inpfile_pe[0] = outfile[0] = 0;
	int	num;
	double	numf;
	int	anchor = 1;
	int 	subotu = 0;
	int 	mlen = 0;

	calign	*aligneng = new calign;

	for(int i = 0; i < argc; i ++) {
		if((!strcmp(argv[i], "-i")) && argc > i + 1) {
			strcpy(inpfile, argv[i + 1]);
		}
		else if((!strcmp(argv[i], "-j")) && argc > i + 1) {
			strcpy(inpfile_pe, argv[i + 1]);
		}
		else if((!strcmp(argv[i], "-o")) && argc > i + 1) {
			strcpy(outfile, argv[i + 1]);
		}
		else if((!strcmp(argv[i], "-a")) && argc > i + 1) {
			if(!strcmp(argv[i + 1], "2")) anchor = 2;
		}
		else if((!strcmp(argv[i], "-k")) && argc > i + 1) {
			sscanf(argv[i + 1],"%d", &num);
			aligneng -> set_kmer(num);
		}
		else if((!strcmp(argv[i], "-b")) && argc > i + 1) {
			sscanf(argv[++ i],"%d", &num);
			aligneng -> set_bandwidth(num);
		}
		else if(!strcmp(argv[i], "-s")) {
			aligneng -> set_strandtag();
		}
		else if((!strcmp(argv[i], "-d")) && argc > i + 1) {
			sscanf(argv[i + 1],"%lf", &numf);
			aligneng -> set_alnperc(numf);
		}
		else if((!strcmp(argv[i], "-g")) && argc > i + 1) {
			sscanf(argv[i + 1],"%d", &num);
			aligneng -> set_mgsize(num);
		}
		else if(!strcmp(argv[i], "-v")) {
			aligneng -> setverbose();
		}
		else if(!strcmp(argv[i], "-mlen") && argc > i + 1) {
			sscanf(argv[i + 1], "%d", &mlen); 
		}
		else if (!strcmp(argv[i], "-sb")) {
			subotu = 1;	
		}
                else if (!strcmp(argv[i], "-abundantonly")) {
                        aligneng -> setabundantonly();
                }
	}

	if(inpfile[0] == 0) printusage("Error: input sequence is not given\n");
	if(outfile[0] == 0) {
		strcpy(outfile, inpfile);
	} 

	time_t jobstart = time(NULL);

	if(inpfile_pe[0] != 0) {
		if(anchor == 1) { 
			if(mlen != 0) aligneng -> getreads_pe(inpfile, inpfile_pe, mlen);
			else aligneng -> getreads_pe(inpfile, inpfile_pe);
		}
		else {
			if(mlen != 0) aligneng -> getreads_pe(inpfile_pe, inpfile, mlen);
			else aligneng -> getreads_pe(inpfile_pe, inpfile);
			aligneng -> setswap();
		}
	}
	else {
		if(mlen != 0) aligneng -> getreads(inpfile, mlen);
		else aligneng -> getreads(inpfile);
	}
/*
	if(mlen != 0) {
		if(inpfile_pe[0] != 0) aligneng -> getreads_pe(inpfile, inpfile_pe, mlen);
		else aligneng -> getreads(inpfile, mlen);
	}
	else {
		if(inpfile_pe[0] != 0) aligneng -> getreads_pe(inpfile, inpfile_pe);
		else aligneng -> getreads(inpfile);
	}
*/
	aligneng -> ca2otu(outfile);

	aligneng -> saveresults();

        //if(subotu) aligneng -> subotu(50, outfile);
        if(subotu) aligneng -> subotu(10, outfile);

	delete aligneng;

	time_t jobfinished = time(NULL);
	double  timeused = difftime(jobfinished, jobstart);
	if(timeused < 60) {
		printf("Time used: %d sec\n", int(timeused));
	}
	else {
		printf("Time used: %d min (%.1f hours)\n", int(timeused / 60), timeused / 3600.0);
	}
}

void printusage(char *error)
{
	printf("%s\n", error);
	printf("AbundantOTU+ -i SourceSeqFile(fasta) -o OutputSeqFile(base-name-only)\n");
	printf("   -i SourceSeqFile (FASTA format): The input sequence file name\n");
	printf("   -j SourceSeqFile (FASTA format): The input sequence file name (paired end!)\n");
	printf("   -o OutputSeqfile (FASTA format): The output file name of consensus sequences\n\n");
	printf("  Other optinal parameters:\n");
        printf("   -abundantonly (default: off)\n");
	printf("   -mlen shortest-reads-to-be-included (default: the same as k-tuple)\n");
	printf("   -a 2 (default off): to use r2 as anchors; only work when -j is given\n");
	printf("   -d dissimilarity-cutoff: default 0.03 (3%% difference)\n");
	printf("   -k k-tuple: length of seed; default 40\n");
	printf("   -b bandwidth: band size for the banded alignment\n");
	printf("   -sb : do subotu inference; default off\n");

	exit(0);
}

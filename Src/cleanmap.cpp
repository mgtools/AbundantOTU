//Program Name: cleanseq
//Latest modification on Dec 1st, 2014
//Yuzhen Ye <yye@indiana.edu>
//Affiliation: School of Informatics and Computing, Indiana University, Bloomington
#include "cleanseq.h"
#include <time.h>

void printusage(char *errorinfo);

main(int argc, char **argv)
{
	std::vector <string> inpfile(2), newfasta(2), newfastq(2);
	std::vector <int> clip(2, 0);
	string mapfile, barcodefile;
	bool ifraw = false;

	CleanSeq	*eng = new CleanSeq;
	int	seq = 0;

	for(int i = 0; i < argc; i ++) {
		if((!strcmp(argv[i], "-i")) && argc > i + 3) {
			inpfile[0] = argv[i + 1];
			newfastq[0] = argv[i + 2];
			newfasta[0] = argv[i + 3];
			seq ++;
		}
		else if((!strcmp(argv[i], "-j")) && argc > i + 3) {
			inpfile[1] = argv[i + 1];
			newfastq[1] = argv[i + 2];
			newfasta[1] = argv[i + 3];
			seq ++;
		}
		else if((!strcmp(argv[i], "-b")) && argc > i + 1) {
			barcodefile = argv[i + 1];
		}
		else if((!strcmp(argv[i], "-m")) && argc > i + 1) {
			mapfile = argv[i + 1];
		}
		else if((!strcmp(argv[i], "-c1")) && argc > i + 1) {
			int tmp = 0;
			sscanf(argv[i + 1], "%d", &tmp);
			clip[0] = tmp;
		}
		else if((!strcmp(argv[i], "-c2")) && argc > i + 1) {
			int tmp = 0;
			sscanf(argv[i + 1], "%d", &tmp);
			clip[1] = tmp;
		}
		else if(!strcmp(argv[i], "-raw")) ifraw = true;
	}

	if(seq < 1) printusage("Error: input sequence is not given (-i)\n");
	if(barcodefile.empty()) printusage("Error: barcode is not given (-b)\n");
	if(mapfile.empty()) printusage("Error: mapfile is not given (-m)\n");

	time_t jobstart = time(NULL);


	eng -> GetBarcode(barcodefile);
	if(seq == 2) 
		if(ifraw)
			eng -> ProcessPairedRaw(inpfile, newfastq, newfasta, mapfile, clip);
		else
			eng -> ProcessPaired(inpfile, newfastq, newfasta, mapfile, clip);
		
	else 
		eng -> ProcessAssembled(inpfile[0], newfastq[0], newfasta[0], mapfile);

	delete eng;

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
	cout<<error<<endl;
	cout<<"CleanMap -i input-fastq-file output-fastq-file output-fasta-file -b barcode.txt -m seq-sample-map-file"<<endl;
	cout<<" Optional:"<<endl; 
	cout<<" -j input-fastq-file output-fastq-file output-fasta-file"<<endl;
	exit(0);
}

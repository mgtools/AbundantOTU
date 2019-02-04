#ifndef _CLEANSEQ_
#define _CLEANSEQ_

#define MAXS 1000
#define BARL 12
#define PRIL 30 

#include "lib.h"

class CleanSeq {
	private:
		char	**sample;
		char 	**barcode_f, **barcode_r;
		char 	**primer_f, **primer_r;
		char 	**space_f, **space_r;
		int	barcodeMis;
		int	numN;	
		int	ifAssembled;
		int	totalSample;
		int	scoreCut;
		int	checkEndLen;
	public:
		CleanSeq(void);
		~CleanSeq(void);
		char* GetWord(const char *str, char *sep, int& beg);
		void RevComp(char *seq);
		std::string RevComp(string seq);
		void GetBarcode(string barcodefile);
		bool NPass(string seq);
		bool QualityPass(string seq);
		bool FindSampleAssembled(string seq, int &sid, int &beg, int &end);
		bool FindSamplePaired(std::vector <string> header, std::vector <string> seq, int &sid, int &beg, int &end);
		void ProcessAssembled(string seqfile, string newfastq, string newfasta, string sampleidfile);
		void ProcessPaired(std::vector <string> seqfile, std::vector <string> newfastq, std::vector <string> newfasta, string sampleidfile, std::vector <int> clip);
		void ProcessPairedRaw(std::vector <string> seqfile, std::vector <string> newfastq, std::vector <string> newfasta, string sampleidfile, std::vector <int> clip);
};

#endif

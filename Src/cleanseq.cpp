#include "cleanseq.h"
using namespace std;

int ifprint = 0;

CleanSeq::CleanSeq(void)
{
	barcodeMis = 1;
	totalSample = 0;
	numN = 1;
	scoreCut = 20;
	checkEndLen = 20;
}

CleanSeq::~CleanSeq(void)
{
	for(int i = 0; i < totalSample; i ++) {
		delete[] sample[i];
		delete[] barcode_f[i];
		delete[] barcode_r[i];
		delete[] primer_f[i];
		delete[] primer_r[i];
		delete[] space_f[i];
		delete[] space_r[i];
	}	
	delete[] sample;
	delete[] barcode_f;
	delete[] barcode_r;
	delete[] primer_f;
	delete[] primer_r;
	delete[] space_f;
	delete[] space_r;
}

char* CleanSeq::GetWord(const char *str, char *sep, int &beg)
{
	int i;
	for(i = beg; i < strlen(str); i ++) {
		if(!strncmp(str + i, sep, strlen(sep))) break; 
	} 
	char *tmp = new char[i - beg + 1];
	if(i > beg)  {
		strncpy(tmp, str + beg, i - beg); 
	}
	tmp[i - beg] = '\0';
	if(beg > strlen(str)) beg = strlen(str) - 1;
	else beg = i + strlen(sep);
	return tmp;
}

void CleanSeq::RevComp(char *seq)
{
	int slen = strlen(seq);
	for(int i = 0; i < slen; i ++) {
		if(seq[i] == 'A') seq[i] = 'T';
		else if(seq[i] == 'T') seq[i] = 'A';
		else if(seq[i] == 'G') seq[i] = 'C';
		else if(seq[i] == 'C') seq[i] = 'G';
		else seq[i] = 'N';
	}
	//reverse
	for(int i = 0; i < slen / 2; i ++) {
		char tmp = seq[i];
		seq[i] = seq[slen - i - 1];
		seq[slen - i - 1] = tmp;
	}
}

std::string CleanSeq::RevComp(string seq)
{
	if(ifprint) cout<<"Ori-Seq "<<seq<<endl;
	int slen = seq.length();
	for(int i = 0; i < slen; i ++) {
		if(seq[i] == 'A') seq[i] = 'T';
		else if(seq[i] == 'T') seq[i] = 'A';
		else if(seq[i] == 'G') seq[i] = 'C';
		else if(seq[i] == 'C') seq[i] = 'G';
		else seq[i] = 'N';
	}
	//reverse
	for(int i = 0; i < slen / 2; i ++) {
		char tmp = seq[i];
		seq[i] = seq[slen - i - 1];
		seq[slen - i - 1] = tmp;
	}
	if(ifprint) cout<<"Rev-Seq "<<seq<<endl;
	return seq;
}

void CleanSeq::GetBarcode(string barcodefile)
{
	sample = new char*[MAXS];
	barcode_f = new char*[MAXS];
	barcode_r = new char*[MAXS];
	primer_f = new char*[MAXS];
	primer_r = new char*[MAXS];
	space_f = new char*[MAXS];
	space_r = new char*[MAXS];
	std::ifstream inf(barcodefile.c_str());
	cout<<"barcodefile "<<barcodefile<<endl;
	if(!inf) {
		cout<<"cannot open file "<<barcodefile<<endl;
		exit(0);
	}
	std::string s;
	std::getline(inf, s);
	int	t = 0;
	while(!inf.eof()) {
		std::getline(inf, s);
		if(s.empty()) continue;
		int b = 0;
		sample[t] = GetWord(s.c_str(), "\t", b);
		barcode_f[t] = GetWord(s.c_str(), "\t", b); 
		barcode_r[t] = GetWord(s.c_str(), "\t", b);
		RevComp(barcode_r[t]);
		primer_f[t] = GetWord(s.c_str(), "\t", b); 
		primer_r[t] = GetWord(s.c_str(), "\t", b); 
		RevComp(primer_r[t]);
		space_f[t] = GetWord(s.c_str(), "\t", b);
		space_r[t] = GetWord(s.c_str(), "\t", b);
		RevComp(space_r[t]);
		printf("sample %d, barcode R1#%s# R2#%s# primer #%s# #%s# space #%s# #%s#\n", t, barcode_f[t], barcode_r[t], primer_f[t], primer_r[t], space_f[t], space_r[t]);
		t ++;
	}
	totalSample = t;
	inf.close();
	cout<<"Total Sample "<<totalSample<<endl;
}
/*
bool CleanSeq::NPass(const char *seq)
{
	int	n = 0;
	for(int i = 0; i < strlen(seq); i ++) {
		if(seq[i] == 'N') {
			n++;
			if(n > numN) return false;
		}
	}
	return true;
}
*/
bool CleanSeq::NPass(string seq)
{
	int	n = 0;
	for(int i = 0; i < seq.length(); i ++) {
		if(seq[i] == 'N') {
			n++;
			if(n > numN) return false;
		}
	}
	return true;
}

//PHRAP 33
bool CleanSeq::QualityPass(string str)
{
	int totScore = 0;
	for(int i = 0; i < str.length(); i ++) {
		int score = str[i] - '!';
		totScore += score;
	}
	int aveScore = totScore / str.length();
	return (aveScore >= scoreCut);
}

bool CleanSeq::FindSampleAssembled(string seq, int &sid, int &beg, int &end)
{
	int	s, i, ipos, jpos, mis_f0, mis_r0, mis_f, mis_r, p;
	int	slen = seq.length();
	for(s = 0; s < totalSample; s ++) {
		ipos = jpos = -1;
		mis_f = mis_r = 0;
		for(i = 0; i < checkEndLen; i ++) {
			mis_f0 = 0;
			for(p = 0; p < BARL; p ++) {
				if(seq[i + p] != barcode_f[s][p]) {
					mis_f0 += 1;
					if(mis_f0 > barcodeMis) break; 
				}
			}
			if(p == BARL) {
				ipos = i;
				mis_f = mis_f0;
				break;
			}
		}
		for(i = 0; i < checkEndLen; i ++) {
			mis_r0 = 0;
			for(p = 0; p < BARL; p ++) {
				if(seq[slen - i - p - 1] != barcode_r[s][BARL - p - 1]) {
					mis_r0 += 1;
					if(mis_r0 > barcodeMis) break;
				}
			}
			if(p == BARL) {
				jpos = slen - i - 1;
				mis_r = mis_r0;
				break;
			}
		}
		//printf("sample %d, %s %s, ipos %d, jpos %d\n", s, barcode_f[s], barcode_r[s], ipos, jpos);
		if((ipos != -1) && (jpos != -1) && (mis_f + mis_r <= barcodeMis)) {
			sid = s;
			beg = ipos;
			end = jpos;
			return true;
		}
	}
	return false;
}

bool CleanSeq::FindSamplePaired(std::vector <string> header, std::vector <string> seq, int &sid, int &beg, int &end)
{
	int	s, i, ipos, jpos, mis_f0, mis_r0, mis_f, mis_r, p;
	//cout<<"header0 "<<header[0]<<endl;
	string str("M00532:71:000000000-AB1LY:1:1101:19234:1465");
	std::size_t found = header[0].find(str);
	if(ifprint) cout<<"header "<<header[0]<<" target "<<str<<" found "<<found<<endl;
	if(found!=std::string::npos) {
		cout<<"r1 "<<seq[0]<<endl;
		cout<<"r2 "<<seq[1]<<endl;
		getchar();
	}
	for(s = 0; s < totalSample; s ++) {
		ipos = jpos = -1;
		mis_f = mis_r = 0;
		for(i = 0; i < checkEndLen; i ++) {
			mis_f0 = 0;
			for(p = 0; p < BARL; p ++) {
				if(seq[0][i + p] != barcode_f[s][p]) {
					mis_f0 += 1;
					if(mis_f0 > barcodeMis) break; 
				}
			}
			if(p == BARL) {
				ipos = i;
				mis_f = mis_f0;
				break;
			}
		}
		for(i = 0; i < checkEndLen; i ++) {
			mis_r0 = 0;
			for(p = 0; p < BARL; p ++) {
				if(seq[1][seq[1].length() - i - p - 1] != barcode_r[s][BARL - p - 1]) {
					mis_r0 += 1;
					if(mis_r0 > barcodeMis) break;
				}
			}
			if(p == BARL) {
				jpos = seq[1].length() - i - 1;
				mis_r = mis_r0;
				break;
			}
		}
		//if(ifprint && ipos != -1) {
		//	printf("sample %d, %s %s, ipos %d, jpos %d\n", s, barcode_f[s], barcode_r[s], ipos, jpos);
		//	getchar();
		//}
		if((ipos != -1) && (jpos != -1) && (mis_f + mis_r <= barcodeMis)) {
			//printf("sample %d, ipos %d, jpos %d\n", s, ipos, jpos);
			sid = s;
			beg = ipos;
			end = jpos;
			return true;
		}
	}
	return false;
}



void CleanSeq::ProcessAssembled(string seqfile, string newfastq, string newfasta, string sampleidfile)
{
	std::ifstream inf(seqfile.c_str());
	if(!inf) {
		cout<<"open file error "<<seqfile<<endl;
		exit(0);
	}
	std::ofstream out(newfastq.c_str());
	if(!out) {
		cout<<"open file error "<<newfastq<<endl;
		exit(0);
	}
	std::ofstream outfa(newfasta.c_str());
	if(!outfa) {
		cout<<"open file error "<<newfasta<<endl;
		exit(0);
	}
	std::ofstream map(sampleidfile.c_str());
	if(!map) {
		cout<<"open file error "<<sampleidfile<<endl;
		exit(0);
	}
	int totalRead = 0;
	int validRead = 0;

	std::string header, seqseq, header_q, quality;
	int	sid, beg, end;
	int	badN = 0;
	int	badA = 0;
	int	notAssigned = 0;
	while(!inf.eof()) {
		std::getline(inf, header);
		std::getline(inf, seqseq);
		std::getline(inf, header_q);
		std::getline(inf, quality);
		totalRead ++;
		if (totalRead % 100000 == 0) {
			cout<<"checked "<<totalRead<<" valid "<<validRead<<endl;
		}
		if(!NPass(seqseq)) {
			badN ++;
			continue;
		}
		if(!QualityPass(quality)) {
			badA ++;
			continue; 
		}
		if(!FindSampleAssembled(seqseq, sid, beg, end)) {
			notAssigned ++;
			continue;
		}
		out<<header<<endl<<seqseq<<endl<<header_q<<endl<<quality<<endl;
		header.erase(0, 1);
		outfa<<">"<<header<<endl;
		int b = beg + strlen(barcode_f[sid]) + strlen(space_f[sid]) + strlen(primer_f[sid]);
		int e = end - strlen(barcode_r[sid]) - strlen(space_r[sid]) - strlen(primer_r[sid]);
		int a = 0;
		//printf("sample %d barcode %s %s, spacer %s %s, primer %s %s\n", sid, barcode_f[sid], barcode_r[sid], space_f[sid], space_r[sid], primer_f[sid], primer_r[sid]);
		//printf("seq len %d, b %d, e %d\n", seqseq.length(), b, e);
		//getchar();
		for(int p = b; p < e; p ++) {
			outfa<<seqseq[p];
			a ++;
			if((a % 70 == 0) && (p != e - 1)) outfa<<endl;
		}
		outfa<<endl;
		map<<validRead<<" "<<sid<<endl;
		validRead ++;
	}
	inf.close();
	out.close();
	outfa.close();
	map.close();
	cout<<"filtered because of Ns "<<badN<<endl;
	cout<<"filtered because of low average score "<<badA<<endl;
	cout<<"filtered because adaptor seqs not found "<<notAssigned<<endl;
	cout<<"checked "<<totalRead<<" valid "<<validRead<<endl;
}


void CleanSeq::ProcessPaired(std::vector <string> seqfile, std::vector <string> newfastq, std::vector <string> newfasta, string sampleidfile, std::vector <int> clip)
{
	std::vector <std::ifstream *> ifs;
	std::vector <std::ofstream *> out;
	std::vector <std::ofstream *> outfa;
	int i;
	for(i = 0; i < 2; i ++) {
		std::ifstream *f = new std::ifstream(seqfile[i].c_str(), std::ios::in);
		cout<<"file "<<i<<" "<<seqfile[i]<<endl;
		if(!f) {
			cout<<"open file error "<<seqfile[i]<<endl;
			exit(0);
		}
		ifs.push_back(f);

		std::ofstream *out0 = new std::ofstream(newfastq[i].c_str(), std::ios::out);
		if(!out0) {
			cout<<"open file error "<<newfastq[i]<<endl;
			exit(0);
		}
		out.push_back(out0);

		std::ofstream *outfa0 = new std::ofstream(newfasta[i].c_str(), std::ios::out);
		if(!outfa0) {
			cout<<"open file error "<<newfasta[i]<<endl;
			exit(0);
		}
		outfa.push_back(outfa0);
	}
	std::ofstream map(sampleidfile.c_str());
	if(!map) {
		cout<<"open file error "<<sampleidfile<<endl;
		exit(0);
	}
	int totalRead = 0;
	int validRead = 0;

	std::vector <std::string> header(2), seqseq(2), header_q(2), quality(2);
	int	beg, end, b, e;
	int	sid;
	int	badN = 0;
	int	badA = 0;
	int	notAssigned = 0;
	while(!ifs[0]->eof()) {
		int valid = 0;
		totalRead ++;
		if (totalRead % 100000 == 0) {
			cout<<"checked "<<totalRead<<" valid "<<validRead<<endl;
		}
		for(i = 0; i < 2; i ++) {
			std::getline(*(ifs[i]), header[i]);
			std::getline(*(ifs[i]), seqseq[i]);
			std::getline(*(ifs[i]), header_q[i]);
			std::getline(*(ifs[i]), quality[i]);
			if(!NPass(seqseq[i])) {
				badN ++;
				continue;
			}
			if(!QualityPass(quality[i])) {
				badA ++;
				continue; 
			}
			valid ++;
		}
		if(valid < 2) continue;
		seqseq[1] = RevComp(seqseq[1]); //reverse complement r2
		if(!FindSamplePaired(header, seqseq, sid, beg, end)) {
			notAssigned ++;
			continue;
		}
		for(i = 0; i < 2; i ++) {
			*out[i]<<header[i]<<endl<<seqseq[i]<<endl<<header_q[i]<<endl<<quality[i]<<endl;
			header[i].erase(0, 1);
			*outfa[i]<<">"<<header[i]<<endl;
			if(i == 0) { 
				b = beg + strlen(barcode_f[sid]) + strlen(space_f[sid]) + strlen(primer_f[sid]);
				e = seqseq[i].length() - clip[0];
			}
			else { //rev-complement
				b = clip[1];
				e = end - strlen(barcode_r[sid]) - strlen(space_r[sid]) - strlen(primer_r[sid]);
			}
			int a = 0;
			for(int p = b; p < e; p ++) {
				*outfa[i]<<seqseq[i][p];
				a ++;
				if((a % 70 == 0) && (p != e - 1)) *outfa[i]<<endl;
			}
			*outfa[i]<<endl;
		}
		validRead ++;
		map<<validRead<<" "<<sid<<endl;
	}
	for(i = 0; i < 2; i ++) {
		(*ifs[i]).close();
		(*out[i]).close();
		(*outfa[i]).close();
	}
	map.close();
	cout<<"filtered because of Ns "<<badN<<endl;
	cout<<"filtered because of low average score "<<badA<<endl;
	cout<<"filtered because adaptor seqs not found "<<notAssigned<<endl;
	cout<<"checked "<<totalRead<<" valid "<<validRead<<endl;
}

void CleanSeq::ProcessPairedRaw(std::vector <string> seqfile, std::vector <string> newfastq, std::vector <string> newfasta, string sampleidfile, std::vector <int> clip)
{
	std::vector <std::ifstream *> ifs;
	std::vector <std::ofstream *> out;
	std::vector <std::ofstream *> outfa;
	int i;
	for(i = 0; i < 2; i ++) {
		std::ifstream *f = new std::ifstream(seqfile[i].c_str(), std::ios::in);
		cout<<"file "<<i<<" "<<seqfile[i]<<endl;
		if(!f) {
			cout<<"open file error "<<seqfile[i]<<endl;
			exit(0);
		}
		ifs.push_back(f);

		std::ofstream *out0 = new std::ofstream(newfastq[i].c_str(), std::ios::out);
		if(!out0) {
			cout<<"open file error "<<newfastq[i]<<endl;
			exit(0);
		}
		out.push_back(out0);

		std::ofstream *outfa0 = new std::ofstream(newfasta[i].c_str(), std::ios::out);
		if(!outfa0) {
			cout<<"open file error "<<newfasta[i]<<endl;
			exit(0);
		}
		outfa.push_back(outfa0);
	}
	std::ofstream map(sampleidfile.c_str());
	if(!map) {
		cout<<"open file error "<<sampleidfile<<endl;
		exit(0);
	}
	int totalRead = 0;
	int validRead = 0;

	std::vector <std::string> header(2), seqseq(2), header_q(2), quality(2);
	int	beg, end, b, e;
	int	sid;
	int	badN = 0;
	int	badA = 0;
	int	notAssigned = 0;
	while(!ifs[0]->eof()) {
		int valid = 0;
		totalRead ++;
		if (totalRead % 100000 == 0) {
			cout<<"checked "<<totalRead<<" valid "<<validRead<<endl;
		}
		for(i = 0; i < 2; i ++) {
			std::getline(*(ifs[i]), header[i]);
			std::getline(*(ifs[i]), seqseq[i]);
			std::getline(*(ifs[i]), header_q[i]);
			std::getline(*(ifs[i]), quality[i]);
			if(!NPass(seqseq[i])) {
				badN ++;
				continue;
			}
			if(!QualityPass(quality[i])) {
				badA ++;
				continue; 
			}
			valid ++;
		}
		//if(valid < 2) continue; //no quality control
		seqseq[1] = RevComp(seqseq[1]); //reverse complement r2
		if(!FindSamplePaired(header, seqseq, sid, beg, end)) {
			notAssigned ++;
			continue;
		}
		for(i = 0; i < 2; i ++) {
			*out[i]<<header[i]<<endl<<seqseq[i]<<endl<<header_q[i]<<endl<<quality[i]<<endl;
			header[i].erase(0, 1);
			*outfa[i]<<">"<<header[i]<<endl;
			if(i == 0) { 
				b = beg + strlen(barcode_f[sid]) + strlen(space_f[sid]) + strlen(primer_f[sid]);
				e = seqseq[i].length() - clip[0];
			}
			else { //rev-complement
				b = clip[1];
				e = end - strlen(barcode_r[sid]) - strlen(space_r[sid]) - strlen(primer_r[sid]);
			}
			int a = 0;
			//for(int p = b; p < e; p ++) {
			for(int p = 0; p < end; p ++) { //no trimming
				*outfa[i]<<seqseq[i][p];
				a ++;
				//if((a % 70 == 0) && (p != e - 1)) *outfa[i]<<endl;
				if((a % 70 == 0) && (p != end - 1)) *outfa[i]<<endl;
			}
			*outfa[i]<<endl;
		}
		validRead ++;
		map<<validRead<<" "<<sid<<endl;
	}
	for(i = 0; i < 2; i ++) {
		(*ifs[i]).close();
		(*out[i]).close();
		(*outfa[i]).close();
	}
	map.close();
	//cout<<"filtered because of Ns "<<badN<<endl;
	//cout<<"filtered because of low average score "<<badA<<endl;
	cout<<"filtered because adaptor seqs not found "<<notAssigned<<endl;
	cout<<"checked "<<totalRead<<" valid "<<validRead<<endl;
}

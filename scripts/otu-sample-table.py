#!/usr/bin/env python
#this script accepts two types of files with reads-sample information

import os
import sys

infofile, clustfile, minseqsup, taxfile, outfile, normalize = "", "", 1, "", "", "yes"
for idx in range(len(sys.argv)):
	if sys.argv[idx] == '-i' and len(sys.argv) > idx + 1:
		infofile = sys.argv[idx + 1]	
	elif sys.argv[idx] == '-c' and len(sys.argv) > idx + 1:
		clustfile = sys.argv[idx + 1]	
	elif sys.argv[idx] == '-m' and len(sys.argv) > idx + 1:
		minseqsup = int(sys.argv[idx + 1])	
	elif sys.argv[idx] == '-t' and len(sys.argv) > idx + 1:
		taxfile = sys.argv[idx + 1]	
	elif sys.argv[idx] == '-o' and len(sys.argv) > idx + 1:
		outfile = sys.argv[idx + 1]	
	elif sys.argv[idx] == '-n' and len(sys.argv) > idx + 1:
		normalize = sys.argv[idx + 1]
if not (infofile and clustfile and outfile):
	print sys.argv[0] + " -i info-file -c clust-file -o output-file <-m minimum-seq-support> <-t taxonomic-assignment> <-n no/yes>"
	sys.exit()

print "process", infofile
ifile = open(infofile, "r")
sampleid, sampleidx, sampleotu = [], [], []
id2sample = [0,]  #real data starts from 1
beg, end = 0, 0
for aline in ifile:
	subs = aline[:-1].split()
	if len(subs) != 3:
		sys.exit(infofile + " incorrect -- which can be prepared using merge.sh script")
	sample, beg, end = subs[0], int(subs[1]), int(subs[2])
	if subs[1] not in sampleid:
		sampleid.append(sample)
		sampleidx.append(beg)
sampleidx.append(end)
ifile.close()
	
print "process", clustfile
ifile = open(clustfile, "r")
con = 0
otuid0 = []
otuid = []
otureads = []
for aline in ifile:
	if aline[0] == '#':
		con += 1
		sidx = 0
		subs = aline.split()
		size = int(subs[3][4:]) 
		if size < minseqsup:
			break 
		otuid0.append("Consensus" + subs[1])
		otuid.append("otu" + subs[1])
		otureads.append([0] * len(sampleid))
	else:
		subs = aline.strip().split()
		seqid = int(subs[1])
		sidx = -1
		for idx in range(sidx, len(sampleid)):
			if (seqid >= sampleidx[idx]) and (seqid < sampleidx[idx + 1]):
				sidx = idx
				break
		otureads[len(otuid)-1][sidx] += 1
ifile.close()

if normalize == "yes":
	for sidx in range(len(sampleid)):
		tot = 0
		for oidx in range(len(otuid)):
			tot += otureads[oidx][sidx]
		for oidx in range(len(otuid)):
			otureads[oidx][sidx] = "%e" % (float(otureads[oidx][sidx]) / tot)

if taxfile:
	otu2tax = {}
	inf = open(taxfile, "r")
	for aline in inf:
		subs = aline.strip().split("\t")
		if len(subs) < 2:
			sys.exit(taxfile + " incorrect -- one line per otu (otu-id and taxonomic-lable seperated by a tab) is expected")
		otu2tax[subs[0]] = subs[1]
	inf.close()
	for oidx in range(len(otuid)):
		oldid = otuid[oidx]
		if otuid0[oidx] in otu2tax:
			otuid[oidx] = otu2tax[otuid0[oidx]] + "|" + oldid

ofile = open(outfile, "w")
print >>ofile, "\t" + "\t".join(sampleid)

for oidx in range(len(otuid)):
	tmp = [otuid[oidx],]
	for sidx in range(len(sampleid) - 1):
		tmp.append(str(otureads[oidx][sidx]))
	print >>ofile, "\t".join(tmp)
ofile.close()

print "otu-sample table saved to", outfile

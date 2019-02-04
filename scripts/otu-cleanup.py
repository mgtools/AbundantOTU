#!/usr/bin/env python

import sys

if len(sys.argv) < 5:
	print "usage:", sys.argv[0], "AbundantOTU-output-basename new-file-basename minsize-to-keep other-keep-info-file(s)"
	sys.exit()

basenameold, basenamenew = sys.argv[1], sys.argv[2]
minsize = int(sys.argv[3])
infofiles = sys.argv[4:]

valid = []

#first check size
try:
	inf = open(basenameold + ".clustsize", "r")
except IOError:
	print "open", basenameold + ".clustsize", "error"	
	sys.exit()
for aline in inf:
	if aline[0] == '#':
		if "length-pe" in aline:
			paired = True
		else:
			paired = False
		continue
	subs = aline[:-1].split()
	if paired:
		size = int(subs[4])
		totalreads = int(subs[5])
	else:
		size = int(subs[3])
		totalreads = int(subs[4])
	if size < minsize: 	
		valid.append(0)
	else:
		valid.append(1)
total = len(valid)

#then check other flags
for n in range(len(infofiles)):
	inf = open(infofiles[n], "r")
	m = 0
	for aline in inf:
		subs = aline.split()
		flag = subs[1]
		if flag == "discard":
			valid[m] = 0
		m += 1
	inf.close()
	if m != total:
		print "error in file", infofiles[n]
		sys.exit()

totalvalid = sum(valid) 
print "total consensus", len(valid), "valid", totalvalid

#output cleaned files		
##clustsize file
inf = open(basenameold + ".clustsize", "r")
aline = inf.readline()
out = open(basenamenew + ".clustsize", "w")
print >>out, aline[:-1]
add = 0
accu = 0
addvalid = 0
for aline in inf:
	if not valid[add]:
		add += 1
		continue
	addvalid += 1
	subs = aline[:-1].split()
	if paired:
		accu += int(subs[4])
		accuper = accu * 100.0 / totalreads
		print >>out, "%-4d %-14s %-4s %-4s %-6s %-6d %8.2f%%" % (addvalid, subs[1], subs[2], subs[3], subs[4], accu, accuper)
	else:
		accu += int(subs[3])
		accuper = accu * 100.0 / totalreads
		print >>out, "%-4d %-14s %-4s %-6s %-6d %8.2f%%" % (addvalid, subs[1], subs[2], subs[3], accu, accuper)
	add += 1
inf.close()
out.close()

##clust file
inf = open(basenameold + ".clust", "r")
out = open(basenamenew + ".clust", "w")
add = -1
for aline in inf:
	if aline[0] == '#':
		add += 1		
		if valid[add]:
			print >>out, aline[:-1]
	elif valid[add]:
		print >>out, aline[:-1]
inf.close()
out.close()

#consensus file
oldfiles, newfiles = [], []
if paired:
	oldfiles = [basenameold + ".1.cons", basenameold + ".2.cons"]	
	newfiles = [basenamenew + ".1.cons", basenamenew + ".2.cons"]	
else:
	oldfiles = [basenameold + ".cons",]
	newfiles = [basenamenew + ".cons",]
for n in range(len(oldfiles)):
	inf = open(oldfiles[n], "r")
	out = open(newfiles[n], "w")
	add = -1
	for aline in inf:
		if aline[0] == ">":
			add += 1
			if valid[add]:
				print >>out, aline[:-1]
		elif valid[add]:
			print >>out, aline[:-1]
	inf.close()
	out.close()

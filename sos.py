#!/usr/bin/env python

import os
import re

items = [
["mul_m4_m4m4", "mult_m4_m4m4", [0, 2, 1]],
["mul_m3_m3m4", "mult_m3_m3m4", [0, 2, 1]],
]

# create regular expressions
subs = []

sdetect = "("

for item in items:
	old, new, shuffle = item

	sfrom = r'([^A-Za-z_:])%s\s*\(' % old

	if len(shuffle):
		for s in shuffle:
			sfrom += r'(.*)'

			if s != shuffle[-1]:
				sfrom += r','
		
		sfrom += r'\)'
	
	sto = r'\1%s(' % new

	if len(shuffle):
		for s in shuffle:
			sto += r'\%d' % (s+2)

			if s != shuffle[-1]:
				sto += r','

		sto += r')'
	
	subs += [[sfrom, sto]]

	if sdetect != "(":
		sdetect += "|"
	sdetect += old

sdetect += ")"
rdetect = re.compile(sdetect)

for i in range(0, len(subs)):
	print(subs[i])
	subs[i][0] = re.compile(subs[i][0])
	#subs[i][1] = re.compile(subs[i][1])

# apply to files
files = []

for dirpath, dirnames, filenames in os.walk("source"):
	for filename in filenames:
		if filename.endswith(".c") or filename.endswith(".cpp") or filename.endswith(".h"):
			files += [dirpath + os.sep + filename]

# default cache is too small, and there doesn't seem
# to be a way to compile patterns for sub()
re._MAXCACHE = len(items)*2 + 50

for file in files:
	f = open(file, "r")
	changed =False

	lines = f.readlines()
	newlines = []
	for line in lines:
		newline = line
		if rdetect.search(line):
			for sub in subs:
				newline = sub[0].sub(sub[1], newline)
				if line != newline:
					changed = True
		newlines += [newline]
	
	f.close()

	if changed:
		print("modified " + file)
		f = open(file, "w")
		f.writelines(newlines)
		f.close()

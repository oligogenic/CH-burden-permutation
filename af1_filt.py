import sys

f = open(sys.argv[1])

h = f.readline()

print(h.strip())


for l in f:
	sl = l.split('\t')
	if sl[33] == '' or float(sl[33]) <= 0.01:
		print(l.strip())



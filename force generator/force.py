import numpy as np
import xlwings as xw

forcebook = xw.Book('reactions.xlsx')
forcesheet = forcebook.sheets('Joint Reactions')

range1 = forcesheet.range('A4')
range2 = forcesheet.range('E4')

lables = {'DEAD': 13, "LIVE": 3, "W30": 71, "W90": 72, "W120": 73, \
		"qv": 91, "qh0": 81, "qh90": 82, "qh30": 83, "qh120": 84, \
		"qh60": 85, "qh150": 86, "wv": 74, "SNOW": 5}

def rotate(a, x0, y0):
	a = a*np.pi/180
	i = np.sin(a)
	j = np.cos(a)
	return [x0*j-y0*i, x0*i+y0*j]

members = set()
for i in range1.expand().value:
	members.add(i[0])

colorall = set()
for i in range(4, len(range1.expand().value)+4):
	colori = forcesheet.range('A'+str(i)).color
	if colori is not None:
		colorall.add(colori)

coltypes = {}
for i in colorall:
	colorimember = set()
	for j in range(4, len(range1.expand().value)+4):
		memberi = forcesheet.range('A'+str(j))
		if memberi.color == i:
			colorimember.add(memberi.value)
	coltypes[i] = colorimember

for i in colorall:
	members = members - coltypes[i]
	
coltypes['white'] = members

NMALL = []
for i, j in zip(range1.expand().value, range2.expand().value):
	NMi = set()
	N, M33, M22 = j[2], -(j[0]*1.5+j[4]), -j[3]+j[1]*1.5
	NMALL.append([i[0], i[1], N, M33, M22])
	
def dicNM(st, lit=NMALL):
	NMdic = {}
	for i in lit:
		if i[0] in st:
			NMdic[i[0]] = []
	for i in NMdic.keys():
		for j in lit:
			if j[0] == i:
				NMdic[i].append(j[1:])
	return NMdic

def represent(dic):
	NALL = []
	M3ALL = []
	M2ALL = []
	for i, j in dic.items():
		for k in j:
			NALL.append(k[1])
			M3ALL.append(k[2])
			M2ALL.append(k[3])
	
	rep = set()
		
	for i, j in dic.items():
		for k in j:
			if k[1] == max(NALL) or k[1] == min(NALL) or \
				k[2] == max(M3ALL) or k[2] == min(M3ALL) or \
				k[3] == max(M3ALL) or k[3] == min(M3ALL):
				rep.add(i)
				
	return rep

colreps = {}
for i in coltypes.keys():	
	colreps[i] = represent(dicNM(coltypes[i]))

outfilescodes = []	
for i, j in colreps.items():
	for k in j:
		outfilescodes.append([i, k])
		
print(outfilescodes)
	

def repsout(outfilescode, lst=NMALL):
	with open(str(outfilescode[0])+'_'+str(outfilescode[1])+'_forces.txt', 'w') as f:
		for k in lst:
			if k[0] == outfilescode[1]:
				f.write('{:},{:},{:.1f},{:.1f},{:.1f},{:.1f},'.\
				format(k[1], lables[k[1]], k[2], k[3], k[4], 0))
				f.write('\n')
		
for i in outfilescodes:					
	repsout(i)








		






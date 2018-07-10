import numpy as np
import os

def trans_comma_str(string):
	'''to transform comma strings into lists of numbers'''
	item = []
	lst = []
	for i in string[:-1]:
		if i != ",":
			item.append(i)
		else:
			lst.append(''.join(item))
			item = []
			continue

	return lst


def rotate(a, x0, y0):
	a = a*np.pi/180
	i = np.sin(a)
	j = np.cos(a)
	return [x0*j-y0*i, x0*i+y0*j]


class Readfile():
	'''create the database'''
	def __init__(s, f):
		s.data = {}
		s.forces = {}
		with open(f) as s.f:
			for i in s.f.readlines():
				s.data[trans_comma_str(i)[1]] = [float(i) for i in trans_comma_str(i)[-4:]]
			s.f.close()
		for i, j in s.data.items():
			try:
				x = j[2]/j[0]
			except ZeroDivisionError:
				x = 9999
			try:
				y = j[1]/j[0]
			except ZeroDivisionError:
				y = 9999
			a = rotate(j[3], x, y)
			s.forces[i] = [j[0], j[0]*a[1], j[0]*a[0]]


def force_inputfile(note):
	while True:
		try:
			forcefile = input(note)
			if forcefile:
				return Readfile(forcefile)
			else:
				print("use the default sample file")
				return Readfile("f_default_in.txt")
		except Exception:
			print("can't find the file or the default, try again!")


def f_input(note):
	while True:
		try:
			a = float(input(note))
		except Exception:
			print("input floats or ints!")
		else:
			return a


def n_input(note):
	while True:
		try:
			a = int(input(note))
		except Exception:
			print("input ints!")
		else:
			return a


def s_input(note):
	return str(input(note))


def forceinput():
	return list(map(f_input, ["N:", "M3:", "M2:"]))


def combos(pairs):
	len_p = len(pairs)
	len_each = []

	combos_num = 1
	for i in pairs:
		combos_num *= len(i)
		len_each.append(len(i))

	combos = np.zeros((combos_num, len_p)).T

	repeat_num = 1
	repeat = [0]
	for i in len_each[:-1]:
		repeat_num *= i
		repeat.append(repeat_num)

	group_num = 1
	for i in range(len_p):
		n = int(combos_num/(len_each[i]*group_num))
		group_num *= len_each[i]

		for j in range(len_each[i]):
			combos[i, (j)*n:(j+1)*n] = pairs[i][j]

	combolistraw = combos.tolist()

	combos = []
	combos.append(combolistraw[0])
	for i,j in zip(combolistraw[1:], repeat[1:]):
		xi = [i[a:a+int(len(i)/j)] for a in range(0, len(i), int(len(i)/j))]

		for b in xi[1:]:
			b[:] = xi[0][:]

		xi = [j for i in xi for j in i]
		combos.append(xi)

	return(np.array(combos).T)


class Dead():
	def __init__(s, tag):
		s.factors = [1, 1.2]
		s.factors2 = [1.35]
		s.t = tag

	def livemain(s):
		return s.factors

	def deadmain(s):
		return s.factors2

	def quake(s, gamma_RE):
		return [i*gamma_RE for i in s.factors]

	def tag(s):
		return s.t


class Live():
	factors = [0, 1.4]
	factors2 = [1, 1.2]
	def __init__(s, tag, combof, quakef):
		s.c = combof
		s.q = quakef
		s.t = tag

	def major(s):
		return s.factors

	def combo(s):
		return [i*s.c for i in s.factors]

	def quakew(s, gamma_RE):
		return [i*s.q*gamma_RE for i in s.factors]

	def quake(s, gamma_RE):
		return [i*s.q*gamma_RE for i in s.factors2]

	def tag(s):
		return s.t


class Quake():
	factors = [1.3, 0.5]
	def __init__(s, tag, gamma_RE):
		s.g = gamma_RE
		s.t = tag

	def major(s):
		return [s.factors[0]*s.g]

	def combo(s):
		return [s.factors[1]*s.g]

	def tag(s):
		return s.t


def comboedF(forces, forcepatterns):
	c1 = combos(list(forcepatterns.values()))

	f1 = []
	for i, j in forcepatterns.items():
		if i in forces.keys():
			f1.append(forces[i])

	f1A = np.array(f1).T

	results = []
	for i in c1:
		results.append(((i*f1A).sum(axis=1)).tolist())

	fout1.write("\n")
	print("load pattern tags:".title(), file=fout1)
	print(list(forcepatterns.keys()), file=fout1)
	fout1.write("combo coefficients accordingly:".title())
	fout1.write("\n")
	print(repr(c1), file=fout1)

	return results, c1.shape[0]


def forceformattedout(forcepairs):
	n = forcepairs[0]
	mx = forcepairs[1]
	my = forcepairs[2]
	try:
		x = my/n
	except ZeroDivisionError:
		x = 9999
	try:
		y = mx/n
	except ZeroDivisionError:
		y = 9999

	if (x > 0) and (y >= 0):
		try:
			thitar = np.arctan(x/y)*180/np.pi
		except ZeroDivisionError:
			thitar = 90
	elif (x >= 0) and (y < 0):
		try:
			thitar = np.arctan(y/x)*180/np.pi+90
		except ZeroDivisionError:
			thitar = 180
	elif (x < 0) and (y <= 0):
		try:
			thitar = np.arctan(x/y)*180/np.pi+180
		except ZeroDivisionError:
			thitar = 270
	elif  (x <= 0) and (y > 0):
		thitar =360-np.arctan(x/y)*180/np.pi

	length = (x**2+y**2)**0.5
	m3 = n*length
	print('{:.1f},{:.1f},{:.1f},{:.1f},{:.1f},'.format(thitar, n, mx, my, m3), file=fout2)


def multiple_combo_generator(forces, forcepatterns):

	tags_count = [0, 0, 0]
	tags_record = [[], [], []]
	for i in forces.keys():
		for j in tags_need_loop:
			if i[0] is j:
				tags_count[tags_need_loop.index(j)] += 1
				tags_record[tags_need_loop.index(j)].append(i)

	value_removed = {}
	for i,k in list(forcepatterns.items()):
		for j in tags_need_loop:
			if i is j:
				value_removed[i] = k
				del forcepatterns[i]

	tags_required = []
	for i, k in value_removed.items():
		for j in tags_need_loop:
			if i is j:
				tags_required.append(tags_record[tags_need_loop.index(i)])

	possible_combos = combos(tags_required)

	for i in possible_combos:
		new_fp = {}
		value_add = {}
		for k, l in value_removed.items():
			for j in i:
				if str(int(j))[0] is k:
					value_add[str(int(j))] = l
		new_fp.update(forcepatterns)
		new_fp.update(value_add)
		yield new_fp


def outpath(note):
	while True:
		try:
			outpath = input(note)
			if outpath:
				return outpath
			else:
				print("use the default path")
				return "out"
		except Exception:
			print("errors, try again!")


dc = Dead('11')
lc = Live('2', 1, 0)
d = Dead('12')
sd = Dead('13')
f = Live('3', 0.7, 0.5)
l = Live('4', 0.7, 0.5)
w = Live('7', 0.6, 0.2)
s = Live('5', 0.7, 0.5)
r = Live('6', 0.7, 0.5)
g_RE = 0.75
qh = Quake('8', g_RE)
qv = Quake('9', g_RE)

tags_need_loop = ['7', '8', '9']

cba = {}
cba["construction_cases"] = {dc.tag():dc.livemain(), lc.tag():lc.major()}
cba["deadmajor_wind_cases"] = {d.tag():d.deadmain(), d.tag():sd.deadmain(), f.tag():f.combo(), l.tag():l.combo(), w.tag():w.combo()}
cba["deadmajor_snow_cases"] = {d.tag():d.deadmain(), d.tag():sd.deadmain(), f.tag():f.combo(), l.tag():l.combo(), s.tag():s.combo()}
cba["wind_cases"] = {d.tag():d.livemain(), sd.tag():sd.livemain(), f.tag():f.combo(), l.tag():l.combo(), w.tag():w.major()}
cba["snow_cases"] = {d.tag():d.livemain(), sd.tag():sd.livemain(), f.tag():f.combo(), l.tag():l.combo(), s.tag():s.major()}
cba["roof_cases"] = {d.tag():d.livemain(), sd.tag():sd.livemain(), f.tag():f.combo(), l.tag():l.combo(), r.tag():r.major()}
cba["floor_wind_cases"] = {d.tag():d.livemain(), sd.tag():sd.livemain(), f.tag():f.major(), l.tag():l.combo(), w.tag():w.combo()}
cba["floor_snow_cases"] = {d.tag():d.livemain(), sd.tag():sd.livemain(), f.tag():f.major(), l.tag():l.combo(), s.tag():s.combo()}
cba["floor_roof_cases"] = {d.tag():d.livemain(), sd.tag():sd.livemain(), f.tag():f.major(), l.tag():l.combo(), r.tag():r.combo()}
cba["live_wind_cases"] = {d.tag():d.livemain(), sd.tag():sd.livemain(), f.tag():f.combo(), l.tag():l.major(), w.tag():w.combo()}
cba["live_snow_cases"] = {d.tag():d.livemain(), sd.tag():sd.livemain(), f.tag():f.combo(), l.tag():l.major(), s.tag():s.combo()}
cba["live_roof_cases"] = {d.tag():d.livemain(), sd.tag():sd.livemain(), f.tag():f.combo(), l.tag():l.major(), r.tag():r.combo()}
cba["quakeh_wind_cases"] = {d.tag():d.quake(g_RE), sd.tag():sd.quake(g_RE), f.tag():f.quake(g_RE), qh.tag():qh.major(), qv.tag():qv.combo(), w.tag():w.quakew(g_RE)}
cba["quakeh_snow_cases"] = {d.tag():d.quake(g_RE), sd.tag():sd.quake(g_RE), f.tag():f.quake(g_RE), qh.tag():qh.major(), qv.tag():qv.combo(), s.tag():s.quakew(g_RE)}
cba["quakev_wind_cases"] = {d.tag():d.quake(g_RE), sd.tag():sd.quake(g_RE), f.tag():f.quake(g_RE), qh.tag():qh.combo(), qv.tag():qh.major(), w.tag():w.quakew(g_RE)}
cba["quakev_snow_cases"] = {d.tag():d.quake(g_RE), sd.tag():sd.quake(g_RE), f.tag():f.quake(g_RE), qh.tag():qh.combo(), qv.tag():qh.major(), s.tag():s.quakew(g_RE)}
loadcasetypes = len(cba.items())

f1 = force_inputfile('the force file name:').forces
loadpatterns = len(f1)


foldname = outpath("foldername for files:")
outfolder = os.getcwd()+"\\"+foldname
if not os.path.exists(outfolder):
	os.makedirs(outfolder)

fout1 = open(outfolder+"\\combos.txt", "w")
fout1.write("Note:")
fout1.write("\n")
fout1.write("Dead loads for construction: {}".format(dc.tag()))
fout1.write("\n")
fout1.write("Dead loads for usage: {}".format(d.tag()))
fout1.write("\n")
fout1.write("Super dead loads for usage i.e. the deadload of the building: {}".format(sd.tag()))
fout1.write("\n")
fout1.write("Live loads for construction: {}".format(lc.tag()))
fout1.write("\n")
fout1.write("Live loads for usage i.e. the outdoor live load: {}".format(l.tag()))
fout1.write("\n")
fout1.write("Live loads on floors for usage: {}".format(f.tag()))
fout1.write("\n")
fout1.write("Snow loads on the roof: {}".format(s.tag()))
fout1.write("\n")
fout1.write("Roof live loads for usage: {}".format(r.tag()))
fout1.write("\n")
fout1.write("Wind loads from multiple directions: {}".format(w.tag()))
fout1.write("\n")
fout1.write("Horizontal quake loads from multiple directions: {}".format(qh.tag()))
fout1.write("\n")
fout1.write("Vertical quake loads for up and down: {}".format(qv.tag()))
fout1.write("\n")
fout1.write("\n")
fout1.write("Loadcase names:".title())
fout1.write("\n")

fout2 = open("combined_forces.txt", "w")


for i, j in cba.items():
	print(i, file=fout1)
	for a, b in j.items():
		print(str(a).title(), [int(round(m*100))/100 for m in b], file=fout1)
	fout1.write("\n")


t1 = tags_need_loop[0]
t2 = tags_need_loop[1]
t3 = tags_need_loop[2]

loadcases = 0
loadcombos = 0
for cbn, cb in cba.items():
	print(cb.keys())
	if (t1 in cb.keys()) or (t2 in cb.keys()) or (t3 in cb.keys()):
		for i in multiple_combo_generator(f1, cb):
			ccc1 = comboedF(f1, i)
			for j in ccc1[0]:
				forceformattedout(j)
			loadcases += 1
			loadcombos += ccc1[1]
	else:
		ccc2 = comboedF(f1, cb)
		for j in ccc2[0]:
			forceformattedout(j)
		loadcases += 1
		loadcombos += ccc2[1]


fout1.write("\n")
fout1.write("\n")
fout1.write("The total number of  load pattern is {}.".format(loadpatterns))
fout1.write("\n")
fout1.write("The total number of  loadcase types is {}.".format(loadcasetypes))
fout1.write("\n")
fout1.write("The total number of  loadcases is {}.".format(loadcases))
fout1.write("\n")
fout1.write("\n")
fout1.write("The total number of  combinations is {}.".format(loadcombos))
fout1.write("\n")
fout1.close()
fout2.close()








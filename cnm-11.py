import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import os


def rotate(a, x0, y0):
	a = a*np.pi/180
	i = np.sin(a)
	j = np.cos(a)
	return [x0*j-y0*i, x0*i+y0*j]


def trans_comma_str(string):
	'''to transform comma strings into lists of numbers'''
	item = []
	lst = []
	for i in string[:-1]:
		if i != ",":
			item.append(i)
		else:
			lst.append(float(''.join(item)))
			item = []
			continue
	return lst


class Readfile():
	'''create the database'''
	def __init__(s, f):
		with open(f) as s.f:
			s.conc = trans_comma_str(s.f.readline())
			s.rebar = trans_comma_str(s.f.readline())
			s.profile = trans_comma_str(s.f.readline())
			s.conc_sec = []
			s.rebar_sec = []
			s.profile_sec = []
			for i in s.f.readlines():
				if trans_comma_str(i)[0] == 1:
					s.conc_sec.append(trans_comma_str(i)[1:])
				elif trans_comma_str(i)[0] == 2:
					s.rebar_sec.append(trans_comma_str(i)[1:])
				elif trans_comma_str(i)[0] == 3:
					s.profile_sec.append(trans_comma_str(i)[1:])
				else:
					s.f.close()


def Readforcefile(note):
	'''create the forces database'''

	while True:
		try:
			forcefile = input(note)
			f = open(forcefile)
		except Exception:
			print("use the sample file.")
			f =  open("combined_forces.txt")
			break
		else:
			break

	forces = []
	for i in f.readlines():
		forces.append(trans_comma_str(i))
	return forces



def sec_post(sec, thitar, x0, y0):
	'''post-process the section in order to 1, move the global
	coordinate system to the a point, rotate the section
	(counterclockwisely)'''
	new_sec = []
	for i in sec:
		k = [0,0,0]
		k[0] = i[0]
		k[-2:] = rotate(thitar, i[-2]-x0, i[-1]-y0)
		new_sec.append(k)
	return new_sec


class Concrete():
	def __init__(s, c):
		if 15 <= c <= 80:
			s.c = c
		else:
			s.c = False

	def n(s):
		if s.c:
			return min(2, 2-(s.c-50)/60)

	def fc(s):
		if s.c:
			alphac1 = max(0.76, 0.76+(0.82-0.76)*(s.c-50)/30)
			alphac2 = min(1, 1-(1-0.87)*(s.c-40)/40)
			return alphac1*alphac2*0.88*s.c/1.4

	def strain_0(s):
		if s.c:
			return max(0.002+0.5*(s.c-50)*1e-5, 0.002)

	def strain_cu(s):
		if s.c:
			return min(0.0033-(s.c-50)*1e-5, 0.0033)


def stressc(strain, fc, n, strain_0, strain_cu):
	if 0 < strain <= strain_0:
		return fc*(1-(1-strain/strain_0)**n)
	elif strain_0 < strain <= strain_cu:
		return fc
	elif strain <= 0:
		return 0


def stressst(strain, E, fy, strainmax):
	if strain <= -fy/E:
		return -fy
	elif -fy/E < strain <= fy/E:
		return strain*E
	else:
		return fy


def strain_pattern(x, strain_cu, strainmax):
	'''to determine the stain pattern of a section. x is the relative
	height of concrete under compression'''
	x0 = strain_cu/(strain_cu+abs(strainmax))
	if  x <= x0:
		b = abs(strainmax)
		a = x*b/(1-x)
	else:
		a = strain_cu
		b = (1-x)*a/x
	return [a, b]


def stain_value(x, a, b):
	'''the stain value at any point under a certain strain pattern'''
	return a-(a+b)*x


def h_sec(sec):
	j = []
	for i in sec:
		j.append(i[-1])
	y1 = sec[j.index(max(j[:]))][-1]
	y2 = -sec[j.index(min(j[:]))][-1]
	return [y1, y1+y2]


def h0_sec(sec_c, sec_s):
	j = []
	k = []
	for i in sec_c:
		j.append(i[-1])
	for i in sec_s:
		k.append(i[-1])
	y1 = sec_c[j.index(max(j[:]))][-1]
	y2 = -sec_s[k.index(min(k[:]))][-1]
	return y1+y2


def stc(sec_cr, y1, h0, sp):
	st = []
	for i in sec_cr:
		svi = stain_value((y1-i[-1])/h0, *sp)
		sti = stressc(svi, fc, n, strain_0, strain_cu)
		st.append(sti)
	return st


def sts(sec_sr, y1, h0, sp, E, fy, strainmax):
	st = []
	for i in sec_sr:
		svi = stain_value((y1-i[-1])/h0, *sp)
		sti = stressst(svi, E, fy, strainmax)
		st.append(sti)
	return st


def nm(sec, st):
	n = []
	m3 = []
	m2 = []
	for i, j in zip(sec, st):
		ni = j*i[0]
		n.append(ni)
		m3.append(ni*i[-1])
		m2.append(ni*i[-2])
	return [sum(n)/1e3, sum(m3)/1e6, sum(m2)/1e6]


def angleinput(note):
	while True:
		try:
			ang_in = input(note)
			if ang_in :
				if 0 <= float(ang_in) <= 359:
					return float(ang_in)
				else:
					return float(ang_in)%360
			else:
				print("use the default value 0")
				return 0
		except Exception:
			print("input a number or nothing!")


def axes_num_input():
	while True:
		try:
			ani = input("input the number of main axis:")
			if ani:
				return int(float(ani))
			else:
				print("use the default value 0 for 360 view")
				return 0
		except Exception:
			print("input a number or nothing!")


def sec_input(note):
	while True:
		try:
			secfilename = input(note)
			if secfilename:
				return Readfile(secfilename)
			else:
				print("use the default sample file")
				return Readfile("default.txt")
		except Exception:
			print("can't find the file or the default, try again!")


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

'''
def range_input():
	while True:
		try:
			rng = int(input("Input searching range now:"))
			if rng <= 0:
				continue
			elif rng >90:
				print("the seaching range shall be less than 90!")
				continue
		except Exception:
			print("numbers only!")
		else:
			return rng
'''

def possible_range(axes_num, FA):
	note = "you must assign a searching range now."
	if axes_num >= 1:
		interval = 360/(axes_num*2)
		r_range = int(interval/2)
		collection = []
		for i in range(axes_num*2):
			collection.append((2*i+1)*interval/2)
		for i, j in zip(collection[:-1], collection[1:]):
			if i <= FA <j:
				ro = (i+j)/2
				break
		else:
			ro = 0
	else:
		ro = 0
		r_range = 90
	return [ro, r_range]


def curve3(AN, plt3D=False):

	FA = AN-AA

	conc_sec_AN = sec_post(s1.conc_sec, AN, *cendroid)
	rebar_sec_AN = sec_post(s1.rebar_sec, AN, *cendroid)
	profile_sec_AN = sec_post(s1.profile_sec, AN, *cendroid)

	y1 = h_sec(conc_sec_AN)[0]
	h = h_sec(conc_sec_AN)[1]
	h0 = h0_sec(conc_sec_AN, rebar_sec_AN)

	major = possible_range(axes_num, FA)[0]
	r_range = possible_range(axes_num, FA)[1]
	NR = []
	MXR = []
	MYR = []
	thitark = []

	N3R = []
	M3RX = []
	M3RY = []

	for k in p_zone:

		r = y1 - h0*k
		thitar = 1
		signal = True

		while signal:

			for i in range(-thitar*10, thitar*10):

				rotation = i/10
				rr = r*np.cos(rotation*np.pi/180)
				conc_sec_rotate = sec_post(conc_sec_AN, rotation, 0, 0)
				rebar_sec_rotate = sec_post(rebar_sec_AN, rotation, 0, 0)
				profile_sec_rotate = sec_post(profile_sec_AN, rotation, 0, 0)
				y1r = h_sec(conc_sec_rotate)[0]
				hr = h_sec(conc_sec_rotate)[1]
				h0r = h0_sec(conc_sec_rotate, rebar_sec_rotate)
				kr = (y1r-rr)/h0r
				spr = strain_pattern(kr, strain_cu, rebar_strainmax)

				stccr = stc(conc_sec_rotate, y1r, h0r, spr)
				nm3r = nm(conc_sec_AN, stccr)

				stcpr = stc(profile_sec_rotate, y1r, h0r, spr)
				nm4r = nm(profile_sec_AN, stcpr)

				stsr0 = sts(rebar_sec_rotate, y1r, h0r, spr, *s1.rebar)
				stsrd = dict(enumerate(stsr0))
				rebar_sec_p = []
				rebar_sec_t = []
				stsrp = []
				stsrt = []
				for a1, a2 in stsrd.items():
					if a2 >= 0:
						rebar_sec_p.append(rebar_sec_AN[a1])
						stsrp.append(stsr0[a1])
					else:
						rebar_sec_t.append(rebar_sec_AN[a1])
						stsrt.append(stsr0[a1])
				nm1rp = nm(rebar_sec_p, stsrp)
				nm1rt = nm(rebar_sec_t, stsrt)

				stspr0 = sts(profile_sec_rotate, y1r, h0r, spr, *s1.profile)
				stsprd = dict(enumerate(stspr0))
				profile_sec_p = []
				profile_sec_t = []
				stsprp = []
				stsprt = []
				for b1, b2 in stsprd.items():
					if b2 >= 0:
						profile_sec_p.append(profile_sec_AN[b1])
						stsprp.append(stspr0[b1])
					else:
						profile_sec_t.append(profile_sec_AN[b1])
						stsprt.append(stspr0[b1])
				nm2rp = nm(profile_sec_p, stsprp)
				nm2rt = nm(profile_sec_t, stsprt)

				M2rp = nm1rp[2]+nm2rp[2]+nm3r[2]-nm4r[2]
				N2rp = nm1rp[0]+nm2rp[0]+nm3r[0]-nm4r[0]
				try:
					xrp = M2rp/N2rp
				except Exception:
					xrp = 0

				M2rt = nm1rt[2]+nm2rt[2]
				N2rt = nm1rt[0]+nm2rt[0]
				try:
					xrt = M2rt/N2rt
				except Exception:
					xrt = 0

				#print("{}/{}={}(xrp), {}/{}={}(xrt)".format(M2rp, N2rp, xrp, M2rt, N2rt, xrt))

				if abs(xrp-xrt) <= 2:
					nm3 = nm(conc_sec, stccr)
					nm4 = nm(profile_sec, stcpr)
					nm1 = nm(rebar_sec, stsr0)
					nm2 = nm(profile_sec, stspr0)

					nm3_AN = nm(conc_sec_AN, stccr)
					nm4_AN = nm(profile_sec_AN, stcpr)
					nm1_AN = nm(rebar_sec_AN, stsr0)
					nm2_AN = nm(profile_sec_AN, stspr0)

					signal = False

			else:
				thitar = thitar + 1
				if FA-major >= 0:
					if thitar > r_range-(major-FA):
						print("can't find a proper range")
						signal = False
				if FA-major < 0:
					if thitar < (major-FA)-r_range:
						print("can't find a proper range")
						signal = False


		N = nm1[0]+nm2[0]+nm3[0]-nm4[0]
		MX = nm1[1]+nm2[1]+nm3[1]-nm4[1]
		MY = nm1[2]+nm2[2]+nm3[2]-nm4[2]

		N3 = nm1_AN[0]+nm2_AN[0]+nm3_AN[0]-nm4_AN[0]
		M3 = nm1_AN[1]+nm2_AN[1]+nm3_AN[1]-nm4_AN[1]

		NR.append(N)
		MXR.append(MX)
		MYR.append(MY)
		thitark.append(rotation)

		if plt3D:
			M3RX.append(M3*np.sin(AN*np.pi/180))
			M3RY.append(M3*np.cos(AN*np.pi/180))
			N3R.append(N3)
			#ax.plot(M3RX, M3RY, zs=N3R, zdir='z', c='r')


		else:
			print("{0:5.1f}{1:5.1f}{2:5.1f}".format(rotation, k, kr))
			print('{0:20.2f}{1:20.2f}{2:20.2f}'.format(N, MX, MY))

	if plt3D:
		return [M3RX, M3RY, N3R]
	else:
		return [NR, MXR, MYR, thitark, AN]


def outputfile(curve_AN):
	foldname = outpath("foldername for files:")
	outfolder = os.getcwd()+"\\"+foldname
	if not os.path.exists(outfolder):
		os.makedirs(outfolder)

	with open(outfolder+"\\"+str(AA)+"_"+str(curve_AN[-1])+".txt", "w") as fout:
		fout.write("main axis rotation(clockwisely):")
		fout.write(str(AA))
		fout.write("\n")
		fout.write("force angle input:")
		fout.write(str(curve_AN[-1]))
		fout.write("\n")
		fout.write('{0:>20}{1:>20}{2:>20}'.format("N", "Mxr", "Myr"))
		fout.write("\n")
		for i, j, k in zip(*curve_AN[0:3]):
			fout.write('{0:20.4f}{1:20.4f}{2:20.4f}'.format(i, j, k))
			fout.write("\n")
		fout.close()


def plot_2D(curve_AN):

	fig1, ax = plt.subplots()
	ax.plot(curve_AN[1], curve_AN[0], 'k')
	ax.plot(curve_AN[2], curve_AN[0], 'r-.')
	ax.grid(True, linestyle='-.')
	ax.tick_params(labelcolor='b', labelsize='medium', width=3)
	ax.set(title='N-M',xlabel='M', ylabel='N')

	fig2, ax = plt.subplots()
	ax.plot(p_zone, curve_AN[3])
	xlab = 'relative height under compression'
	ylab = 'the angle between local 33 axes and the neutral axes'
	ax.set(xlabel=xlab, ylabel=ylab)

	fig3, ax = plt.subplots()
	ax.plot(xcoordsc, ycoordsc, 'x')
	ax.plot(xcoordsr, ycoordsr, 'o')
	ax.plot(xcoordsp, ycoordsp, '^')

	plt.show()


def plot_3D():
	try:
		range_need = int(360/(2*axes_num))
	except ZeroDivisionError:
		range_need = 360
	for j in range(0, range_need+10 , 10):
		print(j)
		data = curve3(j, plt3D=True)
	plt.show()


def plot_3D2():

	fig = plt.figure()
	ax = Axes3D(fig)

	forces = Readforcefile("input the section inner forces file:")

	XF, YF, ZF = [], [], []
	for i in forces:
		XF.append(i[-1]*np.sin(i[0]*np.pi/180))
		YF.append(i[-1]*np.cos(i[0]*np.pi/180))
		ZF.append(i[1])
	ax.scatter(XF, YF, ZF, zdir="z",s=20, c="g", marker="^")

	try:
		range_need = int(360/(2*axes_num))
	except ZeroDivisionError:
		range_need = 360



	X, Y, Z = [], [], []
	for j in range(0, range_need+10 , 10):
		print(j)
		curvej = curve3(j, plt3D=True)
		X.append(curvej[0])
		Y.append(curvej[1])
		Z.append(curvej[2])

	X0 = [j for i in X for j in i]
	Y0 = [j for i in Y for j in i]
	Z0 = [j for i in Z for j in i]
	ax.scatter(X0, Y0, Z0, zdir="z",s=20, c="r", marker="o")


	pairs = []
	for i, j, k in zip(X, Y, Z):
		for a, b, c in zip(i, j, k):
			pairs.append([a, b, c])

	pairs = sorted(pairs, key=lambda x:x[2])

	MXs = []
	for i in pairs:
		MXs.append(i[0])

	thepoint = MXs.index(max(MXs))

	X11, Y11, Z11 = [], [], []
	X12, Y12, Z12 = [], [], []

	for i,j in enumerate(pairs):
		if i<= thepoint:
			X11.append(j[0])
			Y11.append(j[1])
			Z11.append(j[2])
		else:
			X12.append(j[0])
			Y12.append(j[1])
			Z12.append(j[2])

#	ax.plot_trisurf(X11, Y11, Z11)
#	ax.plot_trisurf(X12, Y12, Z12, color="r")

	plt.show()


def whattodo():
	note1 = "3D, 2D, Text, Quit:"
	note2 = "input the N angle:"
	while True:
		key = input(note1)
		if key.lower()[:1] == "3":
			plot_3D2()
			break
		elif key.lower()[:1] == "2":
			AN = angleinput(note2)
			curve_AN = curve3(AN)
			plot_2D(curve_AN)
			break
		elif key.lower()[:1] == "t":
			AN = angleinput(note2)
			curve_AN = curve3(AN)
			outputfile(curve_AN)
		elif key.lower()[:1] == "q":
			break
		else:
			print("do 2D")
			AN = angleinput(note2)
			curve_AN = curve3(AN)
			plot_2D(curve_AN)


#data pre-process
s1 = sec_input("the section file name:")

c = Concrete(*s1.conc)
strain_cu = c.strain_cu()
fc = c.fc()
n = c.n()
strain_0 = c.strain_0()
strain_cu = c.strain_cu()
rebar_strainmax = s1.rebar[2]

#to obtain the cendroid of the section axial resistance
nnn = []
mm2 = []
mm3 = []
cendroid = []

for i in s1.conc_sec:
	ni = fc*i[0]
	m2i = ni*i[1]
	m3i = ni*i[2]
	nnn.append(ni)
	mm2.append(m2i)
	mm3.append(m3i)

for i in s1.rebar_sec:
	nri = s1.rebar[1]*i[0]
	m2ri = nri*i[1]
	m3ri = nri*i[2]
	nnn.append(nri)
	mm2.append(m2ri)
	mm3.append(m3ri)

for i in s1.profile_sec:
	npi = (s1.profile[1]-fc)*i[0]
	m2pi = npi*i[1]
	m3pi = npi*i[2]
	nnn.append(npi)
	mm2.append(m2pi)
	mm3.append(m3pi)

cendroid[:] = [sum(mm2)/sum(nnn), sum(mm3)/sum(nnn)]

#the relative compression zone
p_zone = [i/10 for i in range(0, 21)]

#to input the number of main axes
axes_num = axes_num_input()

AA = angleinput("input the first main axis angle from the gloabl y clockwisely:")

#data for 2D plot and the text output
conc_sec = sec_post(s1.conc_sec, 0, *cendroid)
rebar_sec = sec_post(s1.rebar_sec, 0, *cendroid)
profile_sec = sec_post(s1.profile_sec, 0, *cendroid)

xcoordsc = []
ycoordsc = []
for i in conc_sec:
	xcoordsc.append(i[-2]/10)
	ycoordsc.append(i[-1]/10)

xcoordsr = []
ycoordsr = []
for i in rebar_sec:
	xcoordsr.append(i[-2]/10)
	ycoordsr.append(i[-1]/10)

xcoordsp = []
ycoordsp = []
for i in profile_sec:
	xcoordsp.append(i[-2]/10)
	ycoordsp.append(i[-1]/10)

#start she shit!!
whattodo()

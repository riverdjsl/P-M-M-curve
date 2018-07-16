import numpy as np

def number_input(note):
	while True:
		try:
			x = float(input(note))
		except Exception:
			print("try again!")
		else:
			return x

b = number_input("input the outer width:")
h = number_input("input the outer height:")
ds = number_input("input the outer diameter:")
tw = number_input("input the inner thickness:")
h_r_n = number_input("input the number of rebars vertically: ")
h_r_a = number_input("input the area of rebars vertically: ")	
b_r_n = number_input("input the number of rebars horizontally: ")
b_r_a = number_input("input the area of rebars horizontally: ")	
c_r_a = number_input("input the area of corner rebars: ")	
a_s = number_input("input the cover thickness of rebars: ")	
element_num = number_input("input the element num in both directions: ")
eles_num = number_input("input the element num of pipecore: ")

element_area = b*h/element_num**2
ele_h = h/element_num
ele_b = b/element_num
r_sp_h = (h-2*a_s)/(h_r_n+1)
r_sp_b = (b-2*a_s)/(b_r_n+1)
xs0 = (b-ds)/2+ds/2
ys0 = (h-ds)/2+ds/2
rr = ds/2-tw/2
els_s = 360/eles_num

f = open("rec_sec_w_pipe.txt", "w")
f.write("30,")
f.write("\n")
f.write("200000, 360, -0.01,")
f.write("\n")
f.write("200600, 295, -0.01,")
f.write("\n")

for i in [k*ele_b+ele_b/2 for k in range(int(element_num))]:
	for j in [k*ele_h+ele_h/2 for k in range(int(element_num))]:
		f.write('1,{:.1f},{:.1f},{:.1f},'.format(element_area, i, j))
		f.write("\n")

for i in [a_s, h-a_s]:
	for j in [a_s, b-a_s]:
		f.write('2,{:.1f},{:.1f},{:.1f},'.format(c_r_a, j, i))
		f.write("\n")

for i in [(k+1)*r_sp_b+a_s for k in range(int(b_r_n))]:
	for j in [a_s, h-a_s]:
		f.write('2,{:.1f},{:.1f},{:.1f},'.format(b_r_a, i, j))
		f.write("\n")

for i in [(k+1)*r_sp_h+a_s for k in range(int(h_r_n))]:
	for j in [a_s, b-a_s]:
		f.write('2,{:.1f},{:.1f},{:.1f},'.format(h_r_a, j, i))
		f.write("\n")

for i in [k for k in range(int(eles_num))]:
	x = xs0+rr*np.cos(els_s/2+i*els_s)
	y = ys0+rr*np.sin(els_s/2+i*els_s)
	f.write('3,{:.1f},{:.1f},{:.1f},'.format(els_s*np.pi/180*rr*tw, x, y))
	f.write("\n")
	


f.close()


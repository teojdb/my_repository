#!/usr/bin/python

from sys import argv

dm = open(argv[1])

logo = """+------------------------------------+ 
| CACULATE TPA TENSORS AND DELTA_TLM | 
|		teojdb		     |
+------------------------------------+ """

print logo
print "Program usage: ./gfsm_multi.py QR.out"

dm_lines = dm.readlines()

##### EXTRACT DATA FROM QR.OUT #####

s = {}
e={}
w={}
w['w0']=0
for i in range(len(dm_lines)):
	if dm_lines[i].find("Excitation energies are calculated for symmetry no.")>=0:
		nr_ex=int(dm_lines[i].split()[0])
		states=range(1,nr_ex+1) # 1,2,3,4,...,nr_ex
		states_0=range(nr_ex+1)
	if dm_lines[i].startswith("     Final DFT energy:"):
		e_ground = float(dm_lines[i].split()[3])
	if dm_lines[i].find("Dipole moment components")>=0:
			s['x00'] = float(dm_lines[i+5].split()[1])
			s['y00'] = float(dm_lines[i+6].split()[1])
			s['z00'] = float(dm_lines[i+7].split()[1])
for i in range(len(dm_lines)):
	for st in states:
		if dm_lines[i].startswith(" @ Excited state no:    %d in symmetry"%st):
			
			e['ex%d'%st] = float(dm_lines[i+8].split()[4])
			w['w%d'%st] = e['ex%d'%st] - e_ground
			s['x0%d'%st] = float(dm_lines[i+11].split()[9])
			s['y0%d'%st] = float(dm_lines[i+14].split()[9])
			s['z0%d'%st] = float(dm_lines[i+17].split()[9])
print nr_ex
vec = []
for a in range(len(dm_lines)):
	if dm_lines[a].startswith("@Transition moment <B | A | C> in a.u. for"):
		vec.append(dm_lines[a:a+7])

for i in range(len(vec)):
	liczba = float(vec[i][5].split()[8])
	s[vec[i][1][43].lower()+vec[i][2][50]+vec[i][3][50]]= liczba

c = ['x','y','z']
for st in states:
	for i in range(len(c)):
		
		s["%s%d%d"%(c[i],st,st)] = -s["%s%d%d"%(c[i],st,st)]+s["%s00"%c[i]]


for i in range(nr_ex+1):
	for j in range(nr_ex+1):
		if i>j:		
			s["x%d%d"%(i,j)] = s["x%d%d"%(j,i)]
			s["y%d%d"%(i,j)] = s["y%d%d"%(j,i)]
			s["z%d%d"%(i,j)] = s["z%d%d"%(j,i)]

##### OPEN OUTPUT FILES #####

output = open('%s_gfsmmulti_GEN.out'%argv[1].split(".")[0], 'w')
output2 = open('%s_gfsmmulti_3SMdetails_GEN.out'%argv[1].split(".")[0], 'w')
output3 = open('%s_gfsmmulti_2SMdetails_GEN.out'%argv[1].split(".")[0], 'w')
output4 = open('%s_gfsmmulti_FSM_GEN.out'%argv[1].split(".")[0], 'w')

##### WRITE XDIPLEN, YDIPLEN, ZDIPLEN #####

print>>output2, "*** XDIPLEN ***"
print>>output2, " ", " ".join([x.rjust(10) for x in map(str,range(nr_ex+1))])
for i in range(nr_ex+1):
  print>>output2, i,
  for j in range(nr_ex+1):
    print>>output2, str(round(float(s["x%d%d"%(i,j)]),3)).rjust(10),
  print>>output2
print>>output2

print>>output2, "*** YDIPLEN ***"
print>>output2, " ", " ".join([x.rjust(10) for x in map(str,range(nr_ex+1))])
for i in range(nr_ex+1):
  print>>output2, i,
  for j in range(nr_ex+1):
    print>>output2, str(round(float(s["y%d%d"%(i,j)]),3)).rjust(10),
  print>>output2
print>>output2

print>>output2, "*** ZDIPLEN ***"
print>>output2, " ", " ".join([x.rjust(10) for x in map(str,range(nr_ex+1))])
for i in range(nr_ex+1):
  print>>output2, i,
  for j in range(nr_ex+1):
    print>>output2, str(round(float(s["z%d%d"%(i,j)]),3)).rjust(10),
  print>>output2
print>>output2


##### SOS MODEL #####

no_states = states
dE ={}
print>>output, 95*'-'

def calculate_tensors(wf,n):
	for i in range(nr_ex+1):
		dE['dE%d'%i] = w['w%d'%i] - wf/2
	sxx={}
	syy={}
	szz={}
	sxy={}
	sxz={}
	syz={}

	sxx['sxx_1'] = 2*( (s['x00']*s['x0%d'%n]/dE['dE0']) + (s['x01']*s['x1%d'%n]/dE['dE1']))
	syy['syy_1'] = 2*( (s['y00']*s['y0%d'%n]/dE['dE0']) + (s['y01']*s['y1%d'%n]/dE['dE1']))
	szz['szz_1'] = 2*( (s['z00']*s['z0%d'%n]/dE['dE0']) + (s['z01']*s['z1%d'%n]/dE['dE1']))
	sxy['sxy_1'] = ((s['x00']*s['y0%d'%n] + s['y00']*s['x0%d'%n])/dE['dE0']) + ((s['x01']*s['y1%d'%n]+s['y01']*s['x1%d'%n])/dE['dE1'])
	syz['syz_1'] = ((s['y00']*s['z0%d'%n] + s['z00']*s['y0%d'%n])/dE['dE0']) + ((s['y01']*s['z1%d'%n]+s['z01']*s['y1%d'%n])/dE['dE1'])
	sxz['sxz_1'] = ((s['x00']*s['z0%d'%n] + s['z00']*s['x0%d'%n])/dE['dE0']) + ((s['x01']*s['z1%d'%n]+s['z01']*s['x1%d'%n])/dE['dE1'])

#	sxx['sxx_1'] = 2*((s['x01']*s['x1%d'%n]/dE['dE1']))
#	syy['syy_1'] = 2*((s['y01']*s['y1%d'%n]/dE['dE1']))
#	szz['szz_1'] = 2*((s['z01']*s['z1%d'%n]/dE['dE1']))
#	sxy['sxy_1'] = ((s['x01']*s['y1%d'%n]+s['y01']*s['x1%d'%n])/dE['dE1'])
#	syz['syz_1'] = ((s['y01']*s['z1%d'%n]+s['z01']*s['y1%d'%n])/dE['dE1'])
#	sxz['sxz_1'] = ((s['x01']*s['z1%d'%n]+s['z01']*s['x1%d'%n])/dE['dE1'])


	a = 2
	while a <= len(states):
		sxx['sxx_%d'%a] = sxx['sxx_%d'%(a-1)] + 2*s['x0%d'%a]*s['x%d%d'%(a,n)]/dE['dE%d'%a]
		syy['syy_%d'%a] = syy['syy_%d'%(a-1)] + 2*s['y0%d'%a]*s['y%d%d'%(a,n)]/dE['dE%d'%a]
		szz['szz_%d'%a] = szz['szz_%d'%(a-1)] + 2*s['z0%d'%a]*s['z%d%d'%(a,n)]/dE['dE%d'%a]
		sxy['sxy_%d'%a] = sxy['sxy_%d'%(a-1)] + ((s['x0%d'%a]*s['y%d%d'%(a,n)] + s['y0%d'%a]*s['x%d%d'%(a,n)])/dE['dE%d'%a])
		syz['syz_%d'%a] = syz['syz_%d'%(a-1)] + ((s['y0%d'%a]*s['z%d%d'%(a,n)] + s['z0%d'%a]*s['y%d%d'%(a,n)])/dE['dE%d'%a])
		sxz['sxz_%d'%a] = sxz['sxz_%d'%(a-1)] + ((s['x0%d'%a]*s['z%d%d'%(a,n)] + s['z0%d'%a]*s['x%d%d'%(a,n)])/dE['dE%d'%a])
		a+=1

        deltas = {}
        for i in states:
		deltas["delta_%d"%i] = ((6*(sxx['sxx_%d'%i])**2 + syy['syy_%d'%i]**2 + szz['szz_%d'%i]**2) + 8*(sxy['sxy_%d'%i]**2+sxz['sxz_%d'%i]**2 + syz['syz_%d'%i]**2) + 4*(sxx['sxx_%d'%i]*syy['syy_%d'%i] + syy['syy_%d'%i]*szz['szz_%d'%i] + sxx['sxx_%d'%i]*szz['szz_%d'%i]))/30


	lt={}
	for i in states:	
		lt['list_tensors_%d'%i] = (sxx['sxx_%d'%i],syy['syy_%d'%i], szz['szz_%d'%i], sxy['sxy_%d'%i], syz['syz_%d'%i], sxz['sxz_%d'%i], deltas["delta_%d"%i])
	
	return lt



for i in states:
	print 'transition 0-->%d'%i
	print>>output, "Transition 0-->%d"%i
	print>>output, 95*'-'
	print>>output, "%s & %10s & %10s & %10s & %10s & %10s & %10s & %10s \\\ "%('i','Sxx', 'Syy', 'Szz', 'Sxy', 'Syz', 'Sxz', 'Delta')
	for j in states:
		print>>output, "%d"%j, "& % 10.5f & % 10.5f & % 10.5f & % 10.5f & % 10.5f & % 10.5f & % 10.1f \\\ "%calculate_tensors(w['w%d'%i],i)['list_tensors_%d'%j]
	print>>output, 95*'-'

############ 3 STATE MODEL ##############

def l(vec):
	l = (vec[0]**2+vec[1]**2+vec[2]**2)**0.5
	return l
def cos(vec1,vec2):
	l1 = (vec1[0]**2+vec1[1]**2+vec1[2]**2)**0.5	
	l2 = (vec2[0]**2+vec2[1]**2+vec2[2]**2)**0.5
	cos = (vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2])/(l1*l2)
	return cos
def diff(vec1,vec2):
	vec_diff = [vec1[0]-vec2[0],vec1[1]-vec2[1],vec1[2]-vec2[2]]
	return vec_diff

a={}
# Slownik a zawiera wektory

for i in states_0:
	for j in states_0:
		a["mu%d%d"%(i,j)] = [s['x%d%d'%(i,j)],s['y%d%d'%(i,j)],s['z%d%d'%(i,j)]]

print>>output
print>>output, '3SM with different intermediates'
print>>output
print>>output, '%s & %15s & %15s & %15s & %15s \\\ '%('i', 'delta_ii', 'delta_tt', 'delta_it', 'delta_tlm')

transitions = range(1,nr_ex+1)
intermediates = range(1,nr_ex+1)


d={}
for f in transitions:

	print>>output, 80*'-'
	print>>output, "Transition 0-->%d"%f
	print>>output, 80*'-'

	print>>output2, 115*'-'
	print>>output2, "Transition 0-->%d"%f
        print>>output2, 115*'-'
        print>>output2, " i & %7s & %7s & %7s & %7s & %7s & %12s & %7s & %7s & %7s & %12s & %12s & %12s & %12s & %12s & %12s & %12s\\\ "%("|mu_00|", "|mu_0i|", "|mu_ii|", "|mu_if|", "wi", "|cos^if_0i|", "|mu_0f|", "|mu_ff|", "wf", "|cos^ff_0f|", "|cos^if_0i|","|cos^ff_0f|","|cos^0f_0i|","|cos^ff_if|","|cos^ff_0i|","|cos^if_0f|")
	print>>output2, 115*'-'

	for i in intermediates:
		if i!= f:
			vec_diff = diff(a['mu%d%d'%(f,f)],a['mu%d%d'%(0,0)])
# DELTA^ii and DELTA^jj :
			d['d^ii_i%d_f%d'%(i,f)] = 8*((l(a['mu%d%d'%(0,i)])*l(a['mu%d%d'%(i,f)])/(w['w%d'%i]-0.5*w['w%d'%f]))**2)*(2*cos(a['mu%d%d'%(i,f)],a['mu%d%d'%(0,i)])**2+1)/30

			d['d^ff_i%d_f%d'%(i,f)] = 8*((l(a['mu%d%d'%(0,f)])*l(vec_diff)/(w['w%d'%f]-0.5*w['w%d'%f]))**2)*(2*cos(vec_diff,a['mu%d%d'%(0,f)])**2+1)/30

# DELTA^ij:

			d['d^if_i%d_f%d'%(i,f)] = 8*(l(a['mu%d%d'%(0,i)])*l(a['mu%d%d'%(0,f)])*l(a['mu%d%d'%(i,f)])*l(vec_diff)/((w['w%d'%i]-0.5*w['w%d'%f])*(w['w%d'%f]-0.5*w['w%d'%f])))*\
(cos(a['mu%d%d'%(0,i)],a['mu%d%d'%(i,f)])*cos(a['mu%d%d'%(0,f)],vec_diff) + \
 cos(a['mu%d%d'%(0,i)],a['mu%d%d'%(0,f)])*cos(a['mu%d%d'%(i,f)],vec_diff)+ \
 cos(a['mu%d%d'%(0,i)],vec_diff)*cos(a['mu%d%d'%(0,f)],a['mu%d%d'%(i,f)]))/30
#			print 8*(l(a['mu0%d'%f])*l(a['mu0%d'%i])*l(a['mu%d%d'%(i,f)])*l(a['mu%d%d'%(i,i)])/((w['w%d'%i]-0.5*w['w%d'%f])*(w['w%d'%f]-0.5*w['w%d'%f])))

# DELTA_3SM:
			d['d_3SM_i%d_f%d'%(i,f)] =  d['d^ii_i%d_f%d'%(i,f)]+ d['d^ff_i%d_f%d'%(i,f)] + 2*d['d^if_i%d_f%d'%(i,f)]

# WRITING DATA IN FILE1:
			print>>output, "%d & % 15.6f & % 15.6f & % 15.6f & % 15.6f \\\ "%(i, d['d^ii_i%d_f%d'%(i,f)], d['d^ff_i%d_f%d'%(i,f)], d['d^if_i%d_f%d'%(i,f)], d['d_3SM_i%d_f%d'%(i,f)])


# WRITING DATA IN FILE2:
                        print>>output2, "%2d & %7.3f & %7.3f & %7.3f & %7.3f & %7.3f & %12.3f & %7.3f & %7.3f & %7.3f & %12.3f & %12.3f & %12.3f & %12.3f & %12.3f & %12.3f & %12.3f \\\ "%(i, l(a['mu00']) , l(a['mu0%d'%i]), l(a['mu%d%d'%(i,i)]), l(a['mu%d%d'%(i,f)]), w['w%d'%i], cos(a['mu%d%d'%(i,f)],a['mu0%d'%i]), l(a['mu0%d'%f]), l(a['mu%d%d'%(f,f)]),  w['w%d'%f], cos(a['mu%d%d'%(f,f)],a['mu0%d'%f]),
cos(a['mu%d%d'%(0,i)],a['mu%d%d'%(i,f)]),cos(a['mu%d%d'%(0,f)],a['mu%d%d'%(f,f)]),cos(a['mu%d%d'%(0,i)],a['mu%d%d'%(0,f)]), cos(a['mu%d%d'%(i,f)],a['mu%d%d'%(f,f)]), cos(a['mu%d%d'%(0,i)],a['mu%d%d'%(f,f)]), cos(a['mu%d%d'%(0,f)],a['mu%d%d'%(i,f)]))


print>>output, 80*'-'
print>>output2, 115*'-'


############## GFSM #################


def gfsm_equation(i,j,f):
	nom = 8*l(a["mu%d%d"%(0,i)])*l(a["mu%d%d"%(0,j)])*l(a["mu%d%d"%(i,f)])*l(a["mu%d%d"%(j,f)])
	denom = (w['w%d'%i]-0.5*w['w%d'%f])*(w['w%d'%j]-0.5*w['w%d'%f])
	X = cos(a['mu%d%d'%(0,i)],a['mu%d%d'%(i,f)])*cos(a['mu%d%d'%(0,j)],a['mu%d%d'%(j,f)]) + \
	cos(a['mu%d%d'%(0,i)],a['mu%d%d'%(0,j)])*cos(a['mu%d%d'%(i,f)],a['mu%d%d'%(j,f)])+ \
 	cos(a['mu%d%d'%(0,j)],a['mu%d%d'%(i,f)])*cos(a['mu%d%d'%(0,i)],a['mu%d%d'%(j,f)])
	gfsm = nom*X/(denom*30)
	return gfsm


def gfsm_summation(list_states):
	import numpy
#	list_states = [0] + inter_states + [f]
	f = list_states[-1]
	dim = len(list_states) 
	mat = numpy.zeros((dim,dim))
	for a in range(dim):
		i = list_states[a]
		for b in range(dim):
			j = list_states[b]
			mat[a][b] += gfsm_equation(i,j,f)

	delta_FSM = sum(sum(mat))
	return delta_FSM


def unique_combinations(n,k):
	# n - set (list)
	# k - number of combinations
	from itertools import combinations
	sets = []
	for comb in combinations(n, k):
		sets.append(comb)
	return sets

def domain(f):
	GFSM_sets = []
	import copy 
	inter = copy.deepcopy(transitions)
	inter.remove(f)
	for k in range(len(inter)+1):
		comb = unique_combinations(inter,k)
		for c in comb: 
			FSM = [0] + list(c) + [f]
			GFSM_sets.append(FSM)
	return GFSM_sets



###### WRITE GFSM DATA IN FILE4 ########

print>>output4, "%30s   & %30s \\\ "%("i,j \in", "delta_GFSM")
print>>output4, 70*'-'

for f in transitions:
	f_set = domain(f)
	print>>output4, ("S_%d <--- S_0"%f).center(70)
        print>>output4, 70*'-'
	for s in f_set:
		print>>output4, "%30s   & %30.3f \\\ "%(str(s), gfsm_summation(s))

	print>>output4, 70*'-'


###### OPTIONAL CHECK TEST ##########

def check_test():
	print>>output4
	print>>output4, "CHECK TEST"
	print>>output4
	test_set = [[0,2,2,1],[0,3,3,1],[0,3,3,2],[0,2,3,3,1]]
	for i in test_set:
		print>>output4, "%30s   & %30.3f \\\ "%(str(i), gfsm_summation(i))


########## 2SM DETAILS ##############

print>>output3, "%7s & %7s & %7s & %12s & %7s %12s %12s \\\ "%("|mu_0f|", "|mu_ff|", "|mu_00|", "|del_mu_ff|", "wf", "|cos^0f_ff|", "delta_2SM")


for f in transitions:
        print>>output3, 70*'-'
        print>>output3, "Transition 0-->%d"%f
        print>>output3, 70*'-'
	vec_diff = diff(a['mu%d%d'%(f,f)],a['mu%d%d'%(0,0)])
	delta_2SM = 8*((2*l(a['mu%d%d'%(0,f)])*l(vec_diff)/w['w%d'%f])**2)*(2*cos(vec_diff,a['mu0%d'%f])**2 + 1)/30
 	print>>output3, "%7.3f & %7.3f & %7.3f & %12.3f & %7.3f %12.3f & %12.3f \\\ "%(l(a['mu%d%d'%(0,f)]), l(a['mu%d%d'%(f,f)]),  l(a['mu%d%d'%(0,0)]), l(vec_diff), w['w%d'%f], cos(vec_diff,a['mu0%d'%f]), delta_2SM)

######### OUTPUTS #################


print "Output written in file: %s"%('%s_gfsmmulti_GEN.out'%argv[1].split(".")[0])
print "Details in output: %s"%('%s_gfsmmulti_3SMdetails_GEN.out'%argv[1].split(".")[0])
print "2SM details in output: %s"%('%s_gfsmmulti_2SMdetails_GEN.out'%argv[1].split(".")[0])
print "FSM in output: %s"%('%s_gfsmmulti_FSM_GEN.out'%argv[1].split(".")[0])




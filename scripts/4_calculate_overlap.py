#!/usr/bin/python

helptext = """***********************
|    CALCULATE BETA 2  |
|- an overlap parameter|
|   Joanna Bednarska   |
|        08.2017       |
************************
Usage: ./overlap_beta.py [option] prefix_INHOMO_norca.inp prefix2_exper_scaled.dat
Files "ref00_sim00.dat", "BROADENING" in working directory 
OPTIONS:
-b       PRODUCE NORCA INP FILES WITH DIFFERENT BROADENING
-fit     FIT THE SET OF SPECTRA TO REFERENCE WITH DEFAULT START AND STOP
-fit2    FIT THE SET OF SPECTRA TO REFERENCE WITH DEFINED START AND STOP
-o       CALCULATE SINGLE OVERLAP (BETA PARAMETER) BETWEEN TWO SPECTRA WITH DEFAULT START AND STOP
-o2 	 CALCULATE OVERLAP (BETA PARAMETER) BETWEEN TWO SPECTRA WITH DEFINED START AND STOP
"""

from sys import argv

def read(ref_file, sim_file):
  ref_lines = open(ref_file).readlines()
  sim_lines = open(sim_file).readlines()
  ref_data = {}
  sim_data = {}
  for l in ref_lines:
    data = map(float, l.strip().split())
    ref_data[data[0]] = data[3]
  for l in sim_lines:
    data = map(float, l.strip().split())
    sim_data[data[0]] = data[3]
  return ref_data, sim_data

def ref00_sim00():
  inp = open("ref00_sim00.dat").readlines()
  ref00,sim00,start,stop = map(float, inp[1].split())
  return ref00, sim00,start,stop

def shift_00(ref00, sim00):
  shift = ref00 - sim00
  print "Shift: nm", shift
  simul_data_shift = {}
  for k in sim_data.keys():
    simul_data_shift[k + shift] = sim_data[k]
  sim_ready = {}
  sort = simul_data_shift.keys()
  sort.sort()
  for s in sort:
    sim_ready[s] = simul_data_shift[s]
  start_ref = min(ref_data.keys())
  stop_ref = max(ref_data.keys())
  start_sim = min(sim_ready.keys())
  stop_sim = max(sim_ready.keys())
  start = max(start_sim,start_ref)
  stop = min(stop_sim,stop_ref)
  return sim_ready, start, stop


def extrapol((x1, y1),(x2,y2)):
  a = (y1 - y2)/(x1-x2)                                                         
  b = y2 - (y1 - y2)*x2/(x1 - x2)
  return a, b


def func_sim(sim_ready, x):
  wave = sim_ready.keys()
  wave.sort() 
  for i in range(len(wave)-1):
    if wave[i] <= x and x < wave[i+1]:
      a,b = extrapol((wave[i],sim_ready[wave[i]]),(wave[i+1],sim_ready[wave[i+1]])) 
      y = a*x +b
  return y


def func_ref(ref_data, x):
  wave = ref_data.keys()
  wave.sort()
  for i in range(len(wave)-1):
    if wave[i] <= x and x < wave[i+1]:
      a,b = extrapol((wave[i],ref_data[wave[i]]),(wave[i+1],ref_data[wave[i+1]])) 
      y = a*x +b
  return y


def integral(dx):
  num_rec = int((stop - start)/dx)
  print num_rec
  area_diff = 0
  area_ref = 0
  for n in range(num_rec):
#    print "cal. step", start+n*dx
    if start+n*dx >= stop: break
    a = abs(func_ref(ref_data, start+n*dx) - func_sim(sim_ready, start+n*dx))*dx
    area_diff += a
    f = abs(func_ref(ref_data, start+n*dx))*dx
    area_ref += f
  beta = (area_diff**2/area_ref**2)**0.5
  return beta

def orca_inp(list_inhomo, norca_inp):
  inp = open(norca_inp).readlines()
  for h in list_inhomo:
    h = str(h)
    norca_out0 = norca_inp.replace("norca.inp", "0_norca.inp")
    norca_out = norca_out0.replace("INHOMO",h)
    out = open(norca_out, "w")
    for l in range(len(inp)):
      if "1.0 0.0 0.0" in inp[l]:
        print "ok"
        print>>out, inp[l].strip().replace("BROADENING",h)
      else:
        print>>out, inp[l].strip()


def broadening():
  broad = map(float, open("BROADENING").readline().split())
  return broad

def convert_cm(cm):
	nm = 10**7/cm
	eV = 1.23981*10**(-4)*cm
	return nm,eV


###### MAIN MODULE ########

keywords = argv
print keywords


if len(keywords) == 1:
	print helptext
	exit()
else:
	option = argv[1]

# PRODUCE NORCA INP FILES WITH DIFFERENT BROADENING
if option == "-b":
	  broadenings = broadening()
	  norca_inp = argv[2]
	  orca_inp(broadenings,norca_inp)

 
# FIT THE SET OF SPECTRA TO REFERENCE WITH DEFAULT START AND STOP
elif option == "-fit":
	broadenings = broadening()
	ref = argv[3]
	r00 = ref00_sim00()[0]
	s00 = ref00_sim00()[1]
	prefix_sim = argv[2].split("INHOMO")[0]
	for b in broadenings:
		sim = "%s%s_0_norca_abs_scaled.dat"%(prefix_sim,b)
		ref_data, sim_data = read(ref,sim)
		sim_ready, start, stop = shift_00(r00, s00)
		out = open("%s_shift.dat"%sim[:-4], "w")
		out2 = open("BETA.OUT", "a")
		print>>out2, "%10s  %10.5f"%(b, integral(50))	
		keys = sim_ready.keys()
		keys.sort()
		for k in keys:
 		  nm, eV = convert_cm(k)
		  print>>out, "%15.5f   %15.5f  %15.5f  %15.5f"%(k, eV, nm, sim_ready[k])


# FIT THE SET OF SPECTRA TO REFERENCE WITH DEFINED START AND STOP
elif option == "-fit2":
	broadenings = broadening()
	ref = argv[3]
	r00 = ref00_sim00()[0]
	s00 = ref00_sim00()[1]
	start = ref00_sim00()[2]
	stop  = ref00_sim00()[3]
	prefix_sim = argv[2].split("INHOMO")[0]
	for b in broadenings:
		sim = "%s%s_0_norca_abs_scaled.dat"%(prefix_sim,b)
		ref_data, sim_data = read(ref,sim)
		sim_ready = shift_00(r00, s00)[0]
		out = open("%s_shift.dat"%sim[:-4], "w")
		out2 = open("BETA.OUT", "a")
		print>>out2, "%10s  %10.5f"%(b, integral(50))	
		keys = sim_ready.keys()
		keys.sort()
		for k in keys:
 		  nm, eV = convert_cm(k)
		  print>>out, "%15.5f   %15.5f  %15.5f  %15.5f"%(k, eV, nm, sim_ready[k])



# CALCULATE OVERLAP (BETA PARAMETER) BETWEEN TWO SPECTRA WITH DEFAULT START AND STOP
elif option == "-o":
	ref = argv[2]
        sim = argv[3]
        r00 = ref00_sim00()[0]
        s00 = ref00_sim00()[1]
        ref_data, sim_data = read(ref,sim)
        sim_ready, start, stop = shift_00(r00, s00)
        print start, stop
        print "%15s   beta = %10.5f"%(sim, integral(50))
        out = open("%s_shift.dat"%sim[:-4], "w")
        keys = sim_ready.keys()
        keys.sort()
        for k in keys:
          nm, eV = convert_cm(k)
          print>>out, "%15.5f   %15.5f  %15.5f  %15.5f"%(k, eV, nm, sim_ready[k])


# CALCULATE OVERLAP (BETA PARAMETER) BETWEEN TWO SPECTRA WITH DEFINED START AND STOP 
elif option == "-o2":
	ref = argv[2]
        sim = argv[3]
        r00 = ref00_sim00()[0]
        s00 = ref00_sim00()[1]
        start = ref00_sim00()[2]
        stop  = ref00_sim00()[3]
        ref_data, sim_data = read(ref,sim)
        sim_ready = shift_00(r00, s00)[0]
        print "%15s   beta = %10.5f"%(sim, integral(50))
        out = open("%s_shift.dat"%sim[:-4], "w")
        keys = sim_ready.keys()
        keys.sort()
        for k in keys:
          nm, eV = convert_cm(k)
          print>>out, "%15.5f   %15.5f  %15.5f  %15.5f"%(k, eV, nm, sim_ready[k])

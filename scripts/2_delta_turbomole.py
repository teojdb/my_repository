#!/usr/bin/python2.7

logo = """
 *******************************
 | DISPLACEMENT FROM TURBO CC2 |
 |      Joanna Bednarska       |
 |          12.2016            |
 ******************************* """

helptext = ''' Program usage ./delta_turbo_orca.py -a --> produce simple output 
./delta_turbo_orca.py -c --> check data
./delta_turbo_orca.py -v prefix_mol homo inhomo --> produce orca inp '''

print logo
print helptext

import numpy
from sys import argv

def extract_freq(vib):
  lines_freq = open(vib).readlines() #OTWORZ PLIK vibspectrum
  freq_all_cm = []
  for l in lines_freq:
    if not l.startswith("$") and not l.startswith("#"):
      freq_all_cm.append(float(l[26:35].strip()))
  natoms = (len(freq_all_cm))/3
  return freq_all_cm, natoms

def extract_nmodes(vib_nmodes): 
  lines_modes = open(vib_nmodes).readlines() #OTWORZ PLIK vib_nmodes_all
  nmodes_list = []
  for l in lines_modes:
    if not l.startswith("$") and not l.startswith("#"):
      nmodes_list.extend(map(float, l[7:].split()))
  return nmodes_list

def extract_redmass(contr):
  lines_redmass = open(contr).readlines() #OTWORZ PLIK control 
  red_masses = []
  for i in range(len(lines_redmass)):
    if lines_redmass[i].startswith("$vibrational reduced masses"):
      start = i+1
    if lines_redmass[i].startswith("$nvibro"):
      stop = i
  a = start
  while a < stop:
    red_masses.extend(map(float,lines_redmass[a].split()))
    a+=1
  return red_masses

def extract_grad(ricc2):
  lines_grad = open(ricc2).readlines() #OTWORZ PLIK ricc2.out
  grad = []
  for i in range(len(lines_grad)):
    if lines_grad[i].find("cartesian gradient of the energy (hartree/bohr)")!=-1:
      start = i+3
    if lines_grad[i].find("resulting FORCE  (fx,fy,fz) =")!=-1:   
      stop = i-2
    if lines_grad[i].startswith(" | sym | multi | state |          CC2 excitation energies"):
      d = lines_grad[i+4].split()
      eV = float(d[7])
      au = float(d[9])
      cm = float(d[11])
  j = start
  while j <= stop:
    if not lines_grad[j].startswith("  ATOM") and not lines_grad[j].isspace():
      grad.append([x.replace("D", "E") for x in lines_grad[j].split()[1:]])
    j+=1
  natoms = extract_freq("vibspectrum")[1]
  if natoms%5 == 0:
    Grad = []
    c = 0
    while c < len(grad):
      b = 0
      while b <= 4:
        a = 0
        while a <= 2:
          Grad.append(grad[a+c][b])
          a+=1
        b+=1
      c+=3
    GRAD = map(float, Grad)
  else:
    GRAD_head = []
    c = 0
#    grad[0][0], grad[1][0], grad[2][0], grad[0][1], grad[1][1], grad[2][1] .. grad[2][4]
#    grad[3][0],
    while c < len(grad)-3:
      b = 0
      while b <= 4:
        a = 0
        while a <= 2:
          GRAD_head.append(grad[a+c][b])
          a+=1
        b+=1
      c+=3
    last = natoms%5
    GRAD_tail = []                     
    c = 0
    b = 0
    while b <= last-1:
      a = len(grad)-3
      while a < len(grad):
        GRAD_tail.append(grad[a][b])
        a+=1
      b+=1
    GRAD = map(float, GRAD_head + GRAD_tail)
  return GRAD, eV, au, cm
  

def calc_displ():
  output = open("output", "w")
  freq_all_cm = extract_freq("vibspectrum")[0]
  redmasses = extract_redmass("control")
  natoms = extract_freq("vibspectrum")[1]
  grad = extract_grad("ricc2.out")[0]
  nmodes_list = extract_nmodes("vib_normal_modes")

  ### FREQ ###
  freq_all_au = [] #WSZYSTKIE FREQ W A.U.
  for i in freq_all_cm:
  	freq_all_au.append(4.556327599*10**(-6)*i)
  
  freq_bezzer_cm = [] #FREQ BEZ ROTACJI I TRANSLACJI W CM-1
  freq_bezzer_au = [] #FREQ BEZ ROTACJI I TRANSLACJI W A.U.
  #USUN ZERA Z FREQ
  for i in freq_all_cm:
  	if i != 0:	
  		freq_bezzer_cm.append(i) 
  
  for i in freq_all_au:	
  	if i != 0:
  		freq_bezzer_au.append(i)
  ###
  
  nm = len(freq_bezzer_cm) #LICZBA MODOW NORMALNYCH
  frequencies = freq_bezzer_au #FREQ NM
  
  # WYDRUKUJ W OUTPUCIE FREQ W CM-1 I W AU
  print>>output, "%12s & %14s \\\ "%('Frequencies [cm-1]', 'Frequencies [a.u.]')
  i = 0
  while i < len(frequencies):
  	print>>output, '%18.2f & %18.4e \\\ '%(freq_bezzer_cm[i], freq_bezzer_au[i])
  	i+=1
  print>>output
  print>>output, "Numer of normal modes = %d"%nm
  print>>output
  ###
  
  trans_rot = len(freq_all_cm) - len(freq_bezzer_cm) #LICZBA TRANSLACJE+ROTACJE
  
   ### GRAD ###
  
  grad_text = map(str,grad) # LISTA GRAD STRINGOW
  
  # WYDRUKUJ GRADIENTY W OUTPUCIE
  a = 0
  b = 3
  grad_output = []
  while a < len(grad):
  	c = []
  	for i in range(a, b):
  		c.append(grad_text[i])
  	grad_output.append(c)
  	a+=3
  	b+=3	
  print>>output, "Excited energy gradient [Hartree/Bohr]"
  print>>output, "   % 15s % 15s % 15s"%('X', 'Y', 'Z')
  i = 0
  while i < len(grad_output):
  	print>>output, '%3d % 15s % 15s % 15s'%(i+1, grad_output[i][0], grad_output[i][1], grad_output[i][2])
  	i+=1
  #####
  
  # grad_mat TO WEKTOR GRADIENTOW:
  # grad_mat(1x3N) = [[ Fx1 Fy1 Fz1 Fx2 Fy2 Fz3 .... FxN FyN FzN]]
  
  grad_mat0 = numpy.zeros((1,3*natoms))
  grad_mat = grad_mat0 + grad
  
  # W PLIKU NMODES WSPOLRZEDNE WEKTOROW WLASNYCH ATOMOW SA DRUKOWANE KOLUMNAMI !!!
  
  # ODCZYTAJ PLIK NMODES I ZAPISZ JE W POSTACI MACIERZY nmodes_mat(3N,3N)
  # nmodes(3N, 3N)
  # nmodes = [[x1m1 x1m2 x1m3 ... x1m(3N-6) ] 
  #	  [y1m1 y1m2 z1m3] ....          ]
  #	  [z1m1 z1m2 z1m3  ....          ]]
  #          ^
  #          |
  #	   WEKTORY WLASNE DLA DANEGO MODU
 
#  for i in lines_nmodes:
#  	nmodes_list.extend(map(float, i.split()))


  ### NMODES ###
  nmode = []
  a = 0
  b = 3*natoms
  while a < len(nmodes_list):
  	c = []
  	for i in range(a,b):
  		c.append(nmodes_list[i])
  	nmode.append(c) 
  	a += 3*natoms
  	b += 3*natoms
  nmodes = nmode
  m = numpy.zeros((3*natoms, 3*natoms)) #MACIERZ ZEROWA 3Nx3N
  mat = m + nmodes #mat 
  nmodes_mat = numpy.array(nmodes)
  ###
  
  #SPRAWDZENIE MNOZENIA MACIERZY
  #d = 0
  #e = 0
  #list = []
  #while d < len(nmodes[0]):
  #	list.append(grad[d]*nmodes[0][d])
  #	d +=1
  
  #ZROB MACIERZ WAZONA - KAZDA KOLUMNA PODZIELONA PRZEZ PIERWIASTEK Z REDMASS (A.U.) DANEGO MODU
  from math import sqrt
  i = 0
  while i < len(redmasses):
  	nmodes_mat[:,i] /= sqrt(redmasses[i])
  	i += 1                         
  ###
  
  #POLICZ MASY ZREDUKOWANE - 1 / SUMA KWADRATOW ELEMENTOW KOLUMNY
  red_mat = nmodes_mat[:,:]
  redmass_calc = list(1/sum(red_mat**2))
  ###
  
  #WYDRUKUJ W OUTPUCIE MASY ZREDUKOWANE ODCZYTANE I OBLICZONE
  print>>output
  print>>output, "REDUCED MASSES"
  print>>output, "% 8s & % 12s \\\ "%('Read', 'Calculated') 
  i = 0
  while i < len(redmasses): 
  	print>>output, '%8.5f & %12.5f \\\ '%(redmasses[i], redmass_calc[i])
  	i+=1
  print>>output
  ###
  
  #POMNOZ MACIERZE GRADIENTU(1,3N) I NMODES(3N,3N) I WYDRUKUJ WYNIK W OUTPUCIE
  grad_norm = numpy.dot(grad_mat, nmodes_mat)
  grad_norm_text = list(grad_norm)
  map(str, grad_norm_text)
  print>>output
  print>>output, 'GRADIENT OF EXCITED STATE ENERGY dE/dQ [Hartree/sqrt(amu)] nQ = 3N - 6 '
  i = 0
  gradyy = []
  while i < len(grad_norm[0][trans_rot:]):
        gradyy.append(float(grad_norm_text[0][trans_rot:][i]))
  	print>>output, 'dE/dQ(%d)= % 10s'%(i+1, grad_norm_text[0][trans_rot:][i])
  	i+=1
  print>>output
  ###
  
  # OBLICZ DISPLACEMENTY, FREQUENCIES = FREQ_BEZZER_AU
  # DIS_AU[I] = GRAD[I]/((FREQ_AU**2)*SQRT(REDMASS*1.8228*10**3))
  # DIS_DIM[I] = GRAD/SQRT(FREQ_AU**3*1.8228*10**3)
  
  print>>output, 'DISPLACEMENTS'
  displacements =[]
  a = 0
  while a < len(frequencies):
  	dis = grad_norm[0][a+trans_rot]/(frequencies[a]**2)
  	displacements.append(dis)
  	a+=1
  i = 0
  dis_au =[]
  dis_dim = []
  while i < len(frequencies):
  	dis_au.append(grad_norm[0][i+trans_rot]/(frequencies[i]**2*sqrt(redmasses[i+trans_rot])*1.82288848*10**(3)))
  	dis_dim.append(-grad_norm[0][i+trans_rot]/(sqrt(1.82288848*10**(3)*frequencies[i]**3)))
  	i+=1
  ###
  
  #WYDRUKUJ DISPLACEMENTY
  i = 0
  print>>output, 'Frequencies[cm-1] & Displacements[dimless] & Displacements[a.u] \\\ '
  while i < len(dis_dim):
  	print>>output, '%-17.2f & % 22.5e & % 19.5e \\\ '%(freq_bezzer_cm[i], dis_dim[i], dis_au[i])
  	i+=1
  ###
  reorg = []
  for i in range(len(freq_bezzer_cm)):
  	r = 0.5*(dis_dim[i]**2)*freq_bezzer_cm[i]
  	reorg.append(r)
  
  ereorg = sum(reorg)
  print>>output, "Ereorg = ", ereorg
  outt = open("CC2_HR.OUT", "w")
  for i in range(len(dis_dim)):
    print>>outt, "& %10.3f & %10.3e & %10.3e & %10.3f &"%(freq_bezzer_cm[i],gradyy[i],freq_bezzer_au[i],dis_dim[i])

  return freq_bezzer_cm, dis_dim, ereorg

def orca_inp(prefix, start, stop, points, homo, inhomo):
  atoms = extract_freq("vibspectrum")[1]
  frequencies_cm, dis_dim, en_reorg = calc_displ()
  vertical_energy_cm = extract_grad("ricc2.out")[3] 
  output = open("%s_%s_%s_orca.inp"%(prefix,homo,inhomo),'w')
  print>>output, "%sim"
  print>>output, "Model IMDHO"
  print>>output, "AbsRange %d, %d"%(start,stop)
  print>>output, "NAbsPoints %d"%points
  print>>output, "AbsScaleMode Rel"
  print>>output, "EnInput=EV"
  print>>output, "end"
  print>>output
  print>>output, "$el_states"
  print>>output, "1" 
  print>>output, "1 %.1f %f %f 1.0 0.0 0.0"%(vertical_energy_cm, float(homo), float(inhomo))
  print>>output
  print>>output, "$ss" 
  print>>output, "1" 
  print>>output, "1 0" 
  print>>output
  print>>output, "$vib_freq_gs"
  print>>output, "%d"%(3*atoms -6)
  a = 0
  while a < len(frequencies_cm):
  	print>>output, '% 5d  % f'%(a+1, frequencies_cm[a])
  	a+=1
  print>>output, "$sdnc"
  print>>output, "%d 1"%(3*atoms -6)
  print>>output, "	1"
  a = 0
  while a < len(dis_dim):
  	print>>output, '% 5d  % f'%(a+1, dis_dim[a])
  	a+=1
  print>>output
  output.close()
  print "OUTPUT WRITTEN IN FILE: ", output


####################### RUN PROGRAM #########################

if len(argv) == 2 and argv[1] == "-a":
  calc_displ()
elif len(argv) == 2 and argv[1] == "-c":
  print "Frequencies [cm-1]: "
  print extract_freq("vibspectrum")[0]
  print "No of atoms: "
  print extract_freq("vibspectrum")[1]
  print "Normal modes: "
  print extract_nmodes("vib_normal_modes")[:10]
  print extract_nmodes("vib_normal_modes")[-10:]
  print "Reduced masses: "
  print extract_redmass("control")
  print "Ex. st. gradient: "
  print extract_grad("ricc2.out")
elif len(argv)>2 and argv[1] == "-v":
  prefix = argv[2]
  homo = argv[3]
  inhomo = argv[4]
  start = 15000 
  stop = 40000
  points = 20000
  orca_inp(prefix, start, stop, points, homo, inhomo)


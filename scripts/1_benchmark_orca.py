#!/usr/bin/python

from math import sqrt
from os import listdir
from sys import argv

print "****************************************"
print "* BENCHMARK OF REORGANISATION ENERGIES *"
print "*          Joanna Bednarska            *"
print "****************************************"
print "Do not use slash at the end of directory name of geometries!"
print "****************************************"
print "OPTIONS: -h produce HESSIAN inputs"
print "         -g produce GRADIENT inputs"
print " 	-a make ANALYSIS"

if len(argv) > 4:
	helptext = """Usage: ./benchmark.py [opcja] directory_geometries
OPTIONS: -h HESSIAN
         -g GRADIENT
	 -a ANALYSIS
*** file 'basis.inp' should contain basis functions ***
*** file 'functionals.inp' should contain functionals/other methods ***
DIRECTORY TREE:
..\\benchmark:
..\\..\\basis.inp
..\\..\\functionals.inp
..\\..\\xyz
..\\..\\..\\czasteczka.xyz
..\\..\\optgeo
..\\..\\..\\ czasteczka_metoda_baza.xyz
..\\..\\hessian_inp
..\\..\\..\\hessian_czasteczka_metoda_baza.inp
..\\..\\hessian_log
..\\..\\..\\hessian_czasteczka_metoda_baza.log
..\\..\\hessian_fchk
..\\..\\..\\hessian_czasteczka_metoda_baza.fchk
..\\..\\gradient_inp
..\\..\\..\\gradient_czasteczka_metoda_baza.inp
..\\..\\gradient_log
..\\..\\..\\gradient_czasteczka_metoda_baza.log
..\\..\\gradient_fchk
..\\..\\..\\gradient_czasteczka_metoda_baza.fchk
..\\..\\vib
..\\..\\..\\orca
..\\..\\..\\..\\vib_czasteczka_metoda_baza_orca.inp
..\\..\\..\\norca
..\\..\\..\\..\\vib_czasteczka_metoda_baza_norca.inp
..\\..\\ereorg.tex ---> OUTPUT
"""
	print helptext
	exit(1)


def read_files():
	directory_xyz = argv[2]
	basis = open('basis.inp').readline().strip().split()
	functionals = open('functionals.inp').readline().strip().split()
	geometries = listdir(directory_xyz)
	molecules = [g[:-4] for g in geometries]
	directory_base =  '/'.join(argv[2].split('/')[:-1])
	return [basis,functionals,geometries,directory_xyz,directory_base, molecules]



def create_hessian_inputs(opt):
  basis, functionals, geometries, directory_xyz, directory_base = read_files()[:-1]
  for geom in geometries:
    coord = open('%s/%s'%(directory_xyz,geom),'r').readlines()
    for f in functionals:
      for b in basis:
	inp = open('%s/hessian_inp/hessian_%s_%s_%s.inp'%(directory_base, geom[:-4], f, b),'w')
	if opt[0] == "chk":
	  print>>inp, '%%chk=hessian_%s_%s_%s'%(geom[:-4],f,b)
	print>>inp, '#p %s/%s '%(f,b) + opt[1]
	print>>inp
        print>>inp, "%s %s"%(geom, opt[1])
        print>>inp
	print>>inp, "0 1"
	for line in coord:
	  print>>inp, line.strip()
        print>>inp
	inp.close()


def geo_symbol(geo_opt, output):
  for i in geo_opt:
    if i.split()[1] == '1':
      print>>output, "H   %  5f  % 5f  % 5f"%(float(i.split()[3]), float(i.split()[4]), float(i.split()[5]))	
    if i.split()[1] == '6':
      print>>output, "C   %  5f  % 5f  % 5f"%(float(i.split()[3]), float(i.split()[4]), float(i.split()[5]))
    if i.split()[1] == '7':
      print>>output, "N   %  5f  % 5f  % 5f"%(float(i.split()[3]), float(i.split()[4]), float(i.split()[5]))
    if i.split()[1] == '8':
      print>>output, "O   %  5f  % 5f  % 5f"%(float(i.split()[3]), float(i.split()[4]), float(i.split()[5]))
    if i.split()[1] == '9':
      print>>output, "F   %  5f  % 5f  % 5f"%(float(i.split()[3]), float(i.split()[4]), float(i.split()[5]))
    if i.split()[1] == '5':
      print>>output, "B   %  5f  % 5f  % 5f"%(float(i.split()[3]), float(i.split()[4]), float(i.split()[5]))
    if i.split()[1] == '16':
      print>>output, "S   %  5f  % 5f  % 5f"%(float(i.split()[3]), float(i.split()[4]), float(i.split()[5]))
    if i.split()[1] == '30':
      print>>output, "Zn  %  5f  % 5f  % 5f"%(float(i.split()[3]), float(i.split()[4]), float(i.split()[5]))
    if i.split()[1] == '17':
      print>>output, "Cl  %  5f  % 5f  % 5f"%(float(i.split()[3]), float(i.split()[4]), float(i.split()[5]))
    if i.split()[1] == '35':
      print>>output, "Br  %  5f  % 5f  % 5f"%(float(i.split()[3]), float(i.split()[4]), float(i.split()[5]))



def extract_optgeo():
  directory_base = read_files()[4]
  directory_hessian_log = "%s/hessian_inp"%directory_base
  print directory_hessian_log
  files_log = listdir(directory_hessian_log)
  opt_files = []
  opt_chk = []
  for file in files_log: 
    if file.endswith("-CHCl3.log"):
      print "ok"
      read_log = open("%s/%s"%(directory_hessian_log, file))
      log_lines = read_log.readlines()
      stat = 0
      freq = 0
      stand = 0
      for i in log_lines:                         
        if i.find('NAtoms')>=0:
          natoms = int(i.split()[1])
          break
      for i in log_lines:
        if i.find('Stationary')>=0:
          stat = log_lines.index(i)
      for i in log_lines:
        if i.find('NImag=0')!=-1:
          freq = 1
      if stat == 0 or freq == 0:
        print 'Not optimized', "stat:", stat, "freq:", freq, file
      else:
        print 'Optimized correctly', file
        opt_files.append(file)
        opt_chk.append(file.replace("log","chk"))
        for a in log_lines[stat:]:
          if a.find('Standard orientation')>=0:                               
            stand = log_lines[stat:].index(a)
            break
        if stand == 0:
          print "Nosymm" 
          for a in log_lines[stat:]:
            if a.find('Input orientation')>=0:                               
              stand = log_lines[stat:].index(a)
              break
        geo_optimized = [x.strip() for x in log_lines[stat+stand+5:stat+stand+5+natoms]]
        out = open('%s/optgeo/%s_optgeo.xyz'%(directory_base, file[8:-4]), 'w')
        geo_symbol(geo_optimized,out)
        read_log.close()
  print "Correctly optimized molecules (%d): "%len(opt_files)
  print " ".join(opt_files)
  print "Correct chk files: "
  print " ".join(opt_chk)


# Depending on the file name different variables: g, f, b!
def create_gradient_inputs(opt):
  directory_base = read_files()[4]
  files_optgeo = listdir("%s/optgeo"%directory_base)
  for geom in files_optgeo:
    coord = open('%s/optgeo/%s'%(directory_base,geom),'r').readlines()
    g = geom.split("_")[0]
    f = geom.split("_")[1]
    b = geom.split("_")[2]
#    g = "_".join(geom.split("_")[:2])
#    f = geom.split("_")[2]
#    b = geom.split("_")[3]
    inp = open('%s/gradient_inp/gradient_%s_%s_%s.inp'%(directory_base, g, f, b),'w')
    print g,f,b
    if opt[0] == "chk":
      print>>inp, '%%chk=gradient_%s_%s_%s'%(g,f,b)
    print>>inp, '#p %s/%s '%(f,b) + opt[1]
    print>>inp
    print>>inp, "%s %s"%(geom, opt[1])
    print>>inp
    print>>inp, "0 1"
    for line in coord:
      print>>inp, line.strip()
    print>>inp
    inp.close()




def ereorg(hess_fchk, grad_fchk):
  import numpy
  ### READ FILES ###
  directory_base = read_files()[4]
  hessian_fchk = open("%s/%s"%(directory_base,hess_fchk)) 
  gradient_fchk = open("%s/%s"%(directory_base,grad_fchk)) 
  hessian_lines = hessian_fchk.readlines() 
  gradient_lines = gradient_fchk.readlines() 
  
  ### READ NATOMS AND FREQUENCIES AND REDUCED MASSES ###
  freq_reduced = []
  i =0 
  while i < len(hessian_lines):
    if hessian_lines[i].startswith("Number of atoms"):
      natoms = int(hessian_lines[i].split()[4])
    if hessian_lines[i].startswith("Number of Normal Modes"):
      nm = int(hessian_lines[i].split()[5])
      linijki = 2*nm/5 +1 # lines freq + reduced mases 5 columns
    if hessian_lines[i].startswith("Vib-E2"):
      for a in range(linijki):
        f = map(float,hessian_lines[i+a+1].split()) #podzielona linijka
        freq_reduced.extend(f)
    i+=1
  frequencies_cm = freq_reduced[:nm] # frequencies [cm-1]
  reduced = freq_reduced[nm:2*nm] # reduced masses [au]   
  frequencies_au = [4.556327599*10**(-6)*i for i in frequencies_cm]
  
  ### READ NORMAL MODES ### 
  vib_modes = []
  a = 0
  while a < len(hessian_lines):
    if hessian_lines[a].startswith('Vib-Modes'):
      nm_natom_3 = int(hessian_lines[a].split()[3])
      if nm_natom_3%5 != 0:
        vib_lines = nm_natom_3/5 +1 
      else:
        vib_lines = nm_natom_3/5
      for i in range(vib_lines):
        linijka = hessian_lines[a+i+1].split()
	for b in linijka:
	  vib_modes.append(float(b))
    a+=1
  
  ### CREATE LIST 3natoms x 3natoms OF NORMAL MODES
  trzyN = 3*natoms
  vib_list_3N3N = []
  a = 0
  b = 0
  while a < len(vib_modes):
    mac = []
    for c in range(a, b+trzyN):
        mac.append(float(vib_modes[c]))
    vib_list_3N3N.append(mac)
    a+=trzyN
    b+=trzyN

  ### CONVERT LIST TO MATRIX ###
  # mat(3N-6, 3N) --> 3N-6 rows i 3N columns
  # mat_tr(3N, 3N-6) 
  # mat_tr = [[x1m1 x1m2 x1m3  ... x1m(3N-6) ] 
  #           [y1m1 y1m2 z1m3  ....          ]
  #           [z1m1 z1m2 z1m3  ....          ]]
  #            ^
  #            |
  #            EIGEN VECTORS FOR FIRST MODE
  # mat[0,:] --> show 1 row
  # mat[:,0] --> show 1 column
  m = numpy.zeros((len(vib_list_3N3N),len(vib_list_3N3N[0])))
  mat = m + vib_list_3N3N
  mat_tr = numpy.transpose(mat)

  ### DIVIDE VIB MODES COLUMNS BY SQRT(RED_MASS[A.U.]) ###
  i = 0
  while i < len(reduced):
    mat_tr[:,i] /= sqrt(reduced[i])
    i += 1

  ### READ GRADIENTS ###
  forces = []
  i = 0
  while i < len(gradient_lines):
    if gradient_lines[i].startswith("Cartesian Gradient"):
      nr = int(gradient_lines[i].split()[4])
      if nr%5 !=0:
        linijki = nr/5 +1
      elif nr%5 ==0:
        linijki = nr/5
      for a in range(linijki):
         force = gradient_lines[i+a+1].split()
         for b in force:
           forces.append(float(b))
    i+=1
    if gradient_lines[i].startswith("Dipole Moment") or gradient_lines[i].startswith("ETran NETS"): break
  force_matrix0 = numpy.zeros((1,trzyN))
  force_matrix = force_matrix0 + forces
  
  ### MULTIPLY GRADIENT(1,3N) x VIB(3N,3N) ###
  ### PRODUCT IS DE/DQ GRADIENT OF THE EXCITED STATE ENERGY W.R.T. NORMAL MODE [HARTREE/SQRT(AMU)]
  grad_norm = numpy.dot(force_matrix,mat_tr)

  ### CALCULATE DISPLACEMENT ###
  # DIS_AU[I] = GRAD[I]/((FREQ_AU**2)*SQRT(REDMASS*1.8228*10**3))
  # DIS_DIM[I] = GRAD/SQRT(FREQ_AU**3*1.8228*10**3)
  displacements = []
  displ = []
  a = 0
  while a < len(frequencies_au):
    dis = grad_norm[0][a]/(-frequencies_au[a]**2)
    displacements.append(dis)
    displ.append(dis*reduced[a])
    a+=1
  i = 0
  dis_au =[]
  dis_dim = []
  while i < len(displacements):
#    dis_au.append(-displacements[i]/(sqrt(reduced[i])*1.82288848*10**(3)))
    dis_dim.append(-grad_norm[0][i]/((sqrt(1.82288848*10**(3)))*(sqrt(frequencies_au[i]**3))))
    dis_au.append(-grad_norm[0][i]/((sqrt(1.82288848*10**(3)))*(sqrt(frequencies_au[i]**3)))/sqrt(frequencies_au[i]))   
    i+=1
#  print '%10s %15s %15s'%("freq[cm-1]", "dis[au]", "dis[]")
#  for i in range(len(dis_au)):
#    print '%10.2f  %15.5E  %15.5E'%(frequencies_cm[i], dis_au[i], dis_dim[i])
  adiabatic = []
  for i in range(len(frequencies_au)):
    adiabatic.append(frequencies_cm[i]*dis_dim[i]**2)
  ereorg = 0.5*sum(adiabatic)
  
  i = 0
  while i < len(gradient_lines):
    if gradient_lines[i].startswith("SCF Energy"):
      ground = float(gradient_lines[i].split()[3])
#      print "ground"
    if gradient_lines[i].startswith("Total Energy"):
      excited = float(gradient_lines[i].split()[3])
#      print "excited"
    if gradient_lines[i].startswith("ETran state values"):
      t = gradient_lines[i+1].split()
#      print "trans"
      trans_moment = map(float, [t[1], t[2], t[3]])
    i+=1

  vertical_energy_au = excited-ground
  vertical_energy_cm = 219474.63*vertical_energy_au
  adiabatic_energy_cm = vertical_energy_cm-ereorg
#  print "G/(R^0.5*W^2)  G/W^2     DIS_AU     DIS_DIM  "
#  for i in range(len(displacements)):
#    print "%13.3f  %10.3f %10.3f %10.3f"%(displacements[i], displ[i], dis_au[i],dis_dim[i])

  return ereorg, frequencies_cm, dis_dim, natoms, vertical_energy_cm, adiabatic_energy_cm, trans_moment



def orca_inp(start, stop, points, norm, adiab, homo, inhomo, location):
  directory_base = read_files()[4]
  basis = read_files()[0]
  functionals = read_files()[1]
  molecules = read_files()[5]
  if location == "specific": 
    hess_vib = listdir("%s/vib/hess_vib"%directory_base) 
    grad_vib = listdir("%s/vib/grad_vib"%directory_base)
  elif location == "standard":
    hess_vib = listdir("%s/hessian_fchk"%directory_base)
    grad_vib = listdir("%s/gradient_fchk"%directory_base)
  for b in basis:
    for mol in molecules:
      print "Analysing molecule:", mol
      results = []
      for f in functionals: 
        if "hessian_%s_%s_%s.fchk"%(mol,f,b) in hess_vib and 'gradient_%s_%s_%s.fchk'%(mol,f,b) in grad_vib:
	  print "Analysing functional:", f
	  if norm == 0:
	    output = open("%s/vib/orca/%s_%s_%s_%s_%s_orca.inp"%(directory_base,mol,f,b,homo,inhomo),'w')
	  elif norm == 1:
	    output = open("%s/vib/norca/%s_%s_%s_%s_%s_norca.inp"%(directory_base,mol,f,b,homo,inhomo),'w') 
          if location == "specific":
            en_reorg_module = ereorg('vib/hess_vib/hessian_%s_%s_%s.fchk'%(mol,f,b),'vib/grad_vib/gradient_%s_%s_%s.fchk'%(mol,f,b))
          elif location == "standard":
            en_reorg_module = ereorg('hessian_fchk/hessian_%s_%s_%s.fchk'%(mol,f,b),'gradient_fchk/gradient_%s_%s_%s.fchk'%(mol,f,b))
	  en_reorg = en_reorg_module[0]
	  frequencies_cm = en_reorg_module[1]
	  dis_dim = en_reorg_module[2] 
          atoms = en_reorg_module[3]
          vertical_energy_cm = en_reorg_module[4]
	  adiabatic_energy_cm = en_reorg_module[5]
          trans_moment = en_reorg_module[6]

          print>>output, "%sim"
          print>>output, "Model IMDHO"
          print>>output, "AbsRange %d, %d"%(start,stop)
          print>>output, "NAbsPoints %d"%points
          print>>output, "AbsScaleMode Rel"
          if adiab == 0:
       	    print>>output, "EnInput=EV"
          print>>output, "end"
          print>>output
          print>>output, "$el_states"
          print>>output, "1" 
	  if adiab == 0 and norm == 1:
            print>>output, "1 %.1f %f %f 1.0 0.0 0.0"%(vertical_energy_cm, homo, inhomo)
          elif adiab == 0 and norm == 0:
            print>>output, "1 %.1f %f %f %f %f %f"%(vertical_energy_cm, homo, inhomo, trans_moment[0], trans_moment[1], trans_moment[2])
          elif adiab == 1 and norm == 1:
            print>>output, "1 %.1f %f %f 1.0 0.0 0.0"%(adiabatic_energy_cm, homo, inhomo)
          elif adiab == 1 and norm == 0:
            print>>output, "1 %.1f %f %f %f %f %f"%(adiabatic_energy_cm, homo, inhomo, trans_moment[0], trans_moment[1], trans_moment[2])
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

#        else:
#          print "CAN'T FIND hessian_%s_%s_%s.fchk OR gradient_%s_%s_%s.fchk FILE!"%(mol,f,b,mol,f,b)


def table(hess_fchk, grad_fchk):
#  basis = ["cc-pVTZ"]
#  functionals = ["CAM-B3LYP"]
#  molecules = ["dans","pna" ]
  basis = read_files()[0]
  functionals = read_files()[1]
  molecules = read_files()[5] 
  if len(hess_fchk) != len(grad_fchk):
    print "FILES DON'T MATCH"
  tabs = []
  for b in basis:
    tab = open("ereorg_%s.tex"%b, "w")
    tabs.append("ereorg_%s.tex"%b)
    print>>tab, "\\begin{table}"
    print>>tab, "\\begin{tabular}{ | c | " + len(functionals)*"c" + "| } \\hline\\hline"
    print>>tab, "&" + "  &  ".join(functionals) + "\\\\"
    for mol in molecules:
      print "Analysing molecule:", mol
      results = []
      for f in functionals: 
        if "hessian_%s_%s_%s.fchk"%(mol,f,b) in hess_fchk and 'gradient_%s_%s_%s.fchk'%(mol,f,b) in grad_fchk:
          result = round(ereorg('hessian_fchk/hessian_%s_%s_%s.fchk'%(mol,f,b),'gradient_fchk/gradient_%s_%s_%s.fchk'%(mol,f,b))[0],2)
	else:
          result = "-"
        results.append(str(result).center(10)) 
      print>>tab, "\_".join(mol.split("_")), '&', '&'.join(results), '\\\\'
    print>>tab, "\\hline\\hline"
    print>>tab, "\\end{tabular}"
    print>>tab, "\\caption{\label{tab: %s} Benchmark of reorganisation energy of investigated molecules. }"%b
    print>>tab, "\\end{table}"
    tab.close()
  return tabs
  print "OUTPUT FILES: ", tabs

def ereorg2(hess_fchk, grad_fchk):
  import numpy
  ### READ FILES ###
  hessian_lines = open(hess_fchk).readlines() 
  gradient_lines = open(grad_fchk).readlines() 
  
  ### READ NATOMS AND FREQUENCIES AND REDUCED MASSES ###
  freq_reduced = []
  i =0 
  while i < len(hessian_lines):
    if hessian_lines[i].startswith("Number of atoms"):
      natoms = int(hessian_lines[i].split()[4])
    if hessian_lines[i].startswith("Number of Normal Modes"):
      nm = int(hessian_lines[i].split()[5])
      linijki = 2*nm/5 +1 # lines freq + reduced mases 5 columns
    if hessian_lines[i].startswith("Vib-E2"):
      for a in range(linijki):
        f = map(float,hessian_lines[i+a+1].split()) #podzielona linijka
        freq_reduced.extend(f)
    i+=1
  frequencies_cm = freq_reduced[:nm] # frequencies [cm-1]
  reduced = freq_reduced[nm:2*nm] # reduced masses [au]   
  frequencies_au = [4.556327599*10**(-6)*i for i in frequencies_cm]
  
  ### READ NORMAL MODES ### 
  vib_modes = []
  a = 0
  while a < len(hessian_lines):
    if hessian_lines[a].startswith('Vib-Modes'):
      nm_natom_3 = int(hessian_lines[a].split()[3])
      if nm_natom_3%5 != 0:
        vib_lines = nm_natom_3/5 +1 
      else:
        vib_lines = nm_natom_3/5
      for i in range(vib_lines):
        linijka = hessian_lines[a+i+1].split()
	for b in linijka:
	  vib_modes.append(float(b))
    a+=1
  
  ### CREATE LIST 3natoms x 3natoms OF NORMAL MODES
  trzyN = 3*natoms
  vib_list_3N3N = []
  a = 0
  b = 0
  while a < len(vib_modes):
    mac = []
    for c in range(a, b+trzyN):
        mac.append(float(vib_modes[c]))
    vib_list_3N3N.append(mac)
    a+=trzyN
    b+=trzyN

  ### CONVERT LIST TO MATRIX ###
  # mat(3N-6, 3N) --> 3N-6 rows i 3N columns
  # mat_tr(3N, 3N-6) 
  # mat_tr = [[x1m1 x1m2 x1m3  ... x1m(3N-6) ] 
  #           [y1m1 y1m2 z1m3  ....          ]
  #           [z1m1 z1m2 z1m3  ....          ]]
  #            ^
  #            |
  #            EIGEN VECTORS FOR FIRST MODE
  # mat[0,:] --> show 1 row
  # mat[:,0] --> show 1 column
  m = numpy.zeros((len(vib_list_3N3N),len(vib_list_3N3N[0])))
  mat = m + vib_list_3N3N
  mat_tr = numpy.transpose(mat)

  ### DIVIDE VIB MODES COLUMNS BY SQRT(RED_MASS[A.U.]) ###
  i = 0
  while i < len(reduced):
    mat_tr[:,i] /= sqrt(reduced[i])
    i += 1

  ### READ GRADIENTS ###
  forces = []
  i = 0
  while i < len(gradient_lines):
    if gradient_lines[i].startswith("Cartesian Gradient"):
      nr = int(gradient_lines[i].split()[4])
      if nr%5 !=0:
        linijki = nr/5 +1
      elif nr%5 ==0:
        linijki = nr/5
      for a in range(linijki):
         force = gradient_lines[i+a+1].split()
         for b in force:
           forces.append(float(b))
    i+=1
    if gradient_lines[i].startswith("Dipole Moment") or gradient_lines[i].startswith("ETran NETS"): break
  force_matrix0 = numpy.zeros((1,trzyN))
  force_matrix = force_matrix0 + forces
  ### MULTIPLY GRADIENT(1,3N) x VIB(3N,3N) ###
  ### PRODUCT IS DE/DQ GRADIENT OF THE EXCITED STATE ENERGY W.R.T. NORMAL MODE [HARTREE/SQRT(AMU)]
  grad_norm = numpy.dot(force_matrix,mat_tr)

  ### CALCULATE DISPLACEMENT ###
  # DIS_AU[I] = GRAD[I]/((FREQ_AU**2)*SQRT(REDMASS*1.8228*10**3))
  # DIS_DIM[I] = GRAD/SQRT(FREQ_AU**3*1.8228*10**3)
  displacements = []
  displ = []
  a = 0
  while a < len(frequencies_au):
    dis = grad_norm[0][a]/(-frequencies_au[a]**2)
    displacements.append(dis)
    displ.append(dis*reduced[a])
    a+=1
  i = 0
  dis_au =[]
  dis_dim = []
  while i < len(displacements):
#    dis_au.append(-displacements[i]/(sqrt(reduced[i])*1.82288848*10**(3)))
    dis_dim.append(-grad_norm[0][i]/((sqrt(1.82288848*10**(3)))*(sqrt(frequencies_au[i]**3))))
    dis_au.append(-grad_norm[0][i]/((sqrt(1.82288848*10**(3)))*(sqrt(frequencies_au[i]**3)))/sqrt(frequencies_au[i]))   
    i+=1
#  print '%10s %15s %15s'%("freq[cm-1]", "dis[au]", "dis[]")
#  for i in range(len(dis_au)):
#    print '%10.2f  %15.5E  %15.5E'%(frequencies_cm[i], dis_au[i], dis_dim[i])
  adiabatic = []
  for i in range(len(frequencies_au)):
    adiabatic.append(frequencies_cm[i]*dis_dim[i]**2)
  ereorg = 0.5*sum(adiabatic)
  
  i = 0
  while i < len(gradient_lines):
    if gradient_lines[i].startswith("SCF Energy"):
      ground = float(gradient_lines[i].split()[3])
#      print "ground"
    if gradient_lines[i].startswith("Total Energy"):
      excited = float(gradient_lines[i].split()[3])
#      print "excited"
    if gradient_lines[i].startswith("ETran state values"):
      t = gradient_lines[i+1].split()
#      print "trans"
      trans_moment = map(float, [t[1], t[2], t[3]])
    i+=1

  vertical_energy_au = excited-ground
  vertical_energy_cm = 219474.63*vertical_energy_au
  adiabatic_energy_cm = vertical_energy_cm-ereorg
#  print "G/(R^0.5*W^2)  G/W^2     DIS_AU     DIS_DIM  "
#  for i in range(len(displacements)):
#    print "%13.3f  %10.3f %10.3f %10.3f"%(displacements[i], displ[i], dis_au[i],dis_dim[i])

#  return ereorg, frequencies_cm, dis_dim, natoms, vertical_energy_cm, adiabatic_energy_cm, trans_moment
  return ereorg





####################### OPTIONS #########################
if argv[1] == '-h':
  options = ["chk",'opt=verytight freq(savenm) integral(grid=ultrafine) scf(conver=9,novaracc,maxcycle=500,xqc) nosymm']
  options2 = ["chk",'opt=verytight freq(savenm) integral(grid=ultrafine) scf(conver=9,novaracc,maxcycle=500,xqc) nosymm scrf(iefpcm, solvent=methanol)']
  create_hessian_inputs(options2)
elif argv[1] == '-o':
  extract_optgeo()
elif argv[1] == '-g':
  options = ["chk", 'force td(nstates=5,root=1) integral(grid=ultrafine) scf(conver=9,novaracc,maxcycle=500,xqc) nosymm']
  options2 = ["chk", 'force td(nstates=5,root=1) integral(grid=ultrafine) scf(conver=9,novaracc,maxcycle=500,xqc) nosymm scrf(iefpcm, solvent=methanol)']
  create_gradient_inputs(options2)
elif argv[1] == '-go':
  options = ["chk", 'force td(nstates=5,root=1) opt scf(conver=7,maxcycle=1200) nosymm']
  create_gradient_inputs(options)
elif argv[1] == '-a':
  directory_base = read_files()[4]
  hessian_files = listdir("%s/hessian_fchk"%directory_base)
  gradient_files = listdir("%s/gradient_fchk"%directory_base)
#  print "PROCESSED HESSIAN FILES: \n", hessian_files
#  print "PROCESSED GRADIENT FILES: \n", gradient_files
  table(hessian_files, gradient_files)
elif argv[1] == '-v':
  orca_inp(20000, 50000, 20000, 1, 0, 10, 200, "standard")
# orca_inp(start, stop, points, norm, adiab, homo, inhomo)
elif argv[1] == '-e':
  print  ereorg2(argv[2], argv[3])

else:
  print "No", argv[1], "option found!"

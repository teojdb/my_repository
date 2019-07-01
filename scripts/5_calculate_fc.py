#!/usr/bin/python

helptext = """
**************
| FC FACTORS |
|   teojdb   |
************** 
Usage:
./FC.py [option]  FREQ(I)_FREQ(F)_DISPL.FILE   N_STATES
OPTIONS:
-T (switch on temperature effects)
-N (do not include temp. effects) 
"""

### LOAD MODULES ######

from math import sqrt, factorial, pi, e
import numpy
from sys import argv

### DATA INPUT ###

lines = open(argv[2]).readlines()
freq_ini_cm = map(float, [x.strip().split()[0] for x in lines])
freq_fin_cm = map(float, [x.strip().split()[1] for x in lines])
freq_ini_au = [4.55633*10**(-6)*i for i in map(float, [x.strip().split()[0] for x in lines])]
freq_fin_au = [4.55633*10**(-6)*i for i in map(float, [x.strip().split()[1] for x in lines])]
dis_au = map(float, [x.strip().split()[2] for x in lines])
nstates = int(argv[3])
vertical_cm = float(argv[4])
vertical_au = vertical_cm*4.55633*10**(-6)

# adiab_au = vertical_au - corr_au
# corr_au = freq_au**2*dis_au**2
# corr_au = freq_au*dis_dim**2
# dis_au*sqrt(freq_au) = dis_dim

corr_au = 0
for i in range(len(freq_fin_au)):
  corr_au = freq_fin_au[i]**2*dis_au[i]**2 

adiab_au = vertical_au - corr_au

#### CONSTANTS #####

h = 6.626*10**(-34)

### VARAIBLES #####

def variables(omega, omega_prim,d):
  alpha = omega
  alpha_prim = omega_prim
  S = (alpha_prim*alpha*d**2)/(alpha+alpha_prim)
  A = 2*sqrt(alpha*alpha_prim)/(alpha+alpha_prim)
  C = 2*sqrt(alpha)
  D = 2*sqrt(alpha_prim)
  b = - (alpha_prim*sqrt(alpha)*d)/(alpha+alpha_prim)
  b_prim = alpha*sqrt(alpha_prim)*d/(alpha+alpha_prim)
  return alpha, alpha_prim, S, A, C, D, b, b_prim


##### FUNCTIONS ######

def hermite(n, x):
  if n%2 == 0:
    v = 0    
    for l in range(n/2+1):
      v += factorial(n)*((-1)**(0.5*n-l)/(factorial(2*l)*factorial(0.5*n-l)))*(2*x)**(2*l)
    return v
  elif n%2 ==1:
    w = 0
    for l in range(1+(n-1)/2):
      w += factorial(n)*((-1)**((n-1)*0.5 -l)/(factorial(2*l+1)*factorial((n-1)/2 -l)))*(2*x)**(2*l+1)
    return w

def arguments(a=20000.,b=40000.,step=1.):
  arg = []  
  while a <= b:
    arg.append(a)
    a+=step
  return arg

def osc_function(x):
  alpha = variables(omeg, omeg_prim)[0]
  v = overlap[0]
  N_v = (alpha**0.5/(2**v*factorial(v)*pi**0.5))**0.5
  f = N_v*hermite(v,alpha**0.5*x)*e**(-0.5*alpha*x**2)
  return f

def draw_function(file):
  out = open(file, "w")
  for x in arguments(-0.1,0.1,0.01):
    print>>out, x, osc_function(x)
  out.close()


def draw_function2(file):
  for i in range(8):
    out = open("%s%d"%(file,i), "w")
    for x in arguments(-2,2,0.02):
      print>>out, x, Hermite(i,x)*e**(-x**2)
    out.close()


def B(v,v_prim):
  f = (A*e**(-S)/(2**(v+v_prim)*factorial(v)*factorial(v_prim)))**(0.5)
  return f

def Newton(n, k):
  f = factorial(n)/(factorial(k)*factorial(n-k))
  return f

def double(n):
  if n <=0:
    return 1
  else:
    return n*double(n-2)

def I(k, k_prim):
  K = (k+k_prim)/2
  if (k+k_prim)%2 == 1:
    f = 0
  elif (k+k_prim)%2 == 0:
    f = double(2*K -1)/((alpha+alpha_prim)**K)
  return f

alpha, alpha_prim, S, A, C, D, b, b_prim = variables(freq_ini_au[0],freq_fin_au[0],dis_au[0])


def overlap(v,v_prim):
  M = numpy.zeros((v+1,v_prim+1))
  const = B(v,v_prim)
  for k in range(v+1):
    for k_prim in range(v_prim+1):
      M[k][k_prim] += const*Newton(v,k)*Newton(v_prim,k_prim)*hermite(v-k,b)*hermite(v_prim-k_prim,b_prim)*(C**k)*(D**k_prim)*I(k,k_prim) 
#      print "*****Analysing k, k_prim******", k, k_prim 
#      print "B",B(k,k_prim)
#      print "Newton",Newton(v_prim,k_prim)
#      print "hermite", hermite(v-k,b),"hermite",hermite(v_prim-k_prim,b_prim),"C",(C**k),"D",(D**k_prim),"I", I(k,k_prim)
#  print M
  return sum(sum(M))**2
#  fc = sum(sum(M))**2
#  return fc 

#print "00S**7/7!", overlap(0,0)*S**7/(factorial(7))
#print "04", overlap(0,7)
argum = arguments(29000.0,31000.0,5.0)

def gaussian(x,x0,inten,sigma=10.0):
  f=inten*(1/(sigma*2.5066282746))*e**(-0.5*((x-x0)/sigma)*((x-x0)/sigma))
  return f

#for a in argum:
#  print a, gaussian(a,30000.0,2.0,100.0)

##### MAIN MODULE ###########

if argv[1] == "-N":
  out = open("GAUSSIANS", "w")
  out2 = open("SPECTRUM3", "w")
  argum = arguments(25000.0,35000.0,5.0)
  values = numpy.zeros((1,len(argum)))
  print values
  result = {}
  for i in range(len(freq_ini_au)):
    alpha, alpha_prim, S, A, C, D, b, b_prim = variables(freq_ini_au[i],freq_fin_au[i],dis_au[i])
    FC = {}
    for s in range(nstates+1):
      o = overlap(0,s)
      FC["<0|%d>,%f"%(s,freq_ini_au[i]+s*freq_ini_au[i])] = o
      x0 = 219474.63*(freq_ini_au[i]+s*freq_ini_au[i]+adiab_au)
      print>>out, "<0|%d>  %10.3f  %10.5f"%(s,x0,o)
      
      for j in range(len(argum)):
        y = gaussian(argum[j],x0,o,80.65)
        values[0][j] += y
        print>>out, "g  %10.3f  %10.5f"%(argum[j], y)

  for a in range(len(argum)):
    print>>out2, "%10.1f   %10.5f"%(argum[a],values[0][a])
  
#    result["%.2fcm-1"%freq_ini_cm[i]] = FC

elif argv[1] == "-T":
  result = {}
  for i in range(len(freq_ini_au)):
    alpha, alpha_prim, S, A, C, D, b, b_prim = variables(freq_ini_au[i],freq_fin_au[i],dis_au[i])
    FC = numpy.zeros((nstates+1,nstates+1))
    for si in range(nstates+1):
      for sf in range(nstates+1):
        FC[si][sf] += overlap(si,sf)
    result["%.2fcm-1"%freq_ini_au[i]] = FC
  print result




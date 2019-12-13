#!/usr/bin/python
import sys
from os import walk
import os
import pickle
import shutil
import linecache
#import glob
#import fnmatch
#from subprocess import call
#from subprocess import Popen, PIPE
from matplotlib.dates import date2num 
import datetime
import numpy
import matplotlib.pyplot as plt
#import matplotlib.dates as dates
from os.path import basename
from os.path import splitext
from operator import truediv
import scipy.ndimage
from numpy import exp, loadtxt, pi, sqrt
from lmfit import Model
import subprocess


def fit_function(x, A, beta, B, mu, sigma):
    return (A * exp(-x/beta) + B * exp(-1.0 * (x - mu)**2 / (2 * sigma**2)))

def expo(x, scale, decay,offset):
    """exponential"""
    return scale*exp(-x/decay)+offset

def gaussian(x, amp, cen, wid):
    """1-d gaussian: gaussian(x, amp, cen, wid)"""
    return (amp / (sqrt(2*pi) * wid)) * exp(-(x-cen)**2 / (2*wid**2))

def line(x, slope, intercept):
    """a line"""
    return slope*x + intercept

def return_energy(b,a0,a1,a2):
    en=b*b*a2+b*a1+a0
    return en

def find_ind(lst, a): #Finding the indices of matching elements in list 
    return [i for i, x in enumerate(lst) if x==a]

def WriteSPE(h,e,c):
    ofilename=OUTDIR+"/"+filename_data+".spe"
    nh=h.replace(h.split()[6],"%8.6f"%(cf))
    exestr="sed -i "+"\'s/"+h.split()[6]+"/"+"%8.6f"%(cf)+"/g\' "+filename_data
    if bFix:
        os.system(exestr)
    filename  = open(ofilename,'w')
    print >>filename, "%s"%(nh),
    
    for i in range(0,len(c)):
        print >>filename,"%10.6f %10.6f"%(e[i],c[i])
    print "Output file ->",ofilename   

def SplitFileName(fn):
    basename, ext = os.path.splitext(fn)
    path,fname = os.path.split(basename)
    if len(path)==0:
        path="."
    return path,fname,ext

def ReadSP(infilename):
    binv = []
    countv = []
    icount = 0
    f = open(infilename, 'r')
    header = f.readline()
    for line in f:
        line = line.strip()
        columns = line.split()
        binv += [float(columns[0])]
        countv += [float(columns[1])]
    countv[-1] = 1
    return header,binv,countv

def FitSP_Peak(h,b1,c,ext,ppos):
    global cf
    a0=float(h.split()[5])
    a1=float(h.split()[6])
    cf=a1
    a2=float(h.split()[7])
    if "spe" in ext:
        b = b1
    else:
        b = [return_energy(b1[i],a0,a1,a2) for i in range(len(b1))]
    print "pars: ",a0," ",a1," ",a2
    b2=numpy.array(b)
    c2=numpy.array(c)
    cond = (b2>=ppos-100) & (b2<=ppos+100)
    b = b2[ cond ]
    c = c2[ cond ]
    
    ofilename=OUTDIR+"/"+filename_data+".png"
    sh_marker=["+",".","^","*","p","s","x","D","h","o","+",".","^","*","p","s","x","D","h","o","+",".","^","*","p","s"]
    sh_color=['b','g','r','c','m','y','k','orange','darkblue','darkgreen','darkred','darkcyan','darkmagenta','deeppink','firebrick','gold','darkviolet','lavenderblush','dodgerblue','indigo','limegreen']
    #gmodel = Model(gaussian) + Model(line)
    #gmodel = Model(gaussian) + Model(expo)
    gmodel = Model(fit_function)
    gmodel.set_param_hint('A', value=10000)
    gmodel.set_param_hint('beta', value=0.001)
    gmodel.set_param_hint('B', value=500)
    gmodel.set_param_hint('mu', value=ppos)
    gmodel.set_param_hint('sigma', value=20)
    
    pars = gmodel.make_params()
    result = gmodel.fit(c, pars,x=b)
    print(result.fit_report())
    #print result.params
    print "%6.2f %6.4f %8.6f -> %8.6f"%(result.params["mu"].value,result.params["sigma"].value,a1,ppos*a1/result.params["mu"].value)
    if not "spe" in ext:
        cf=ppos*a1/result.params["mu"].value
    energy = [return_energy(b1[i],a0,cf,a2) for i in range(len(b1))]
    
    fig = plt.gcf()
    plt.title(h,fontsize=14,fontweight='bold')
    if "spe" in ext:
        plt.xlabel('Energy ',fontsize=14,fontweight='bold')
    else:
        plt.xlabel('Bins ',fontsize=14,fontweight='bold')
    plt.ylabel('Counts ',fontsize=14,fontweight='bold')
    plt.tick_params(labelsize=18)
    #plt.text(energyBins[binid][0]+10, 0.8*maxv2,strv,fontsize=14,fontweight='bold')

    plt.plot(b,c,marker=sh_marker[0],linestyle='None',color=sh_color[0],label=filename_data)
    plt.plot(b, result.best_fit, 'r-', label="fit %6.2f %6.4f"%(result.params["mu"].value,result.params["sigma"].value))

    legend = plt.legend(loc='upper right', shadow=True)
    frame = legend.get_frame()
    frame.set_facecolor('0.90')
    for label in legend.get_texts():
        label.set_fontsize('large')
    for label in legend.get_lines():
        label.set_linewidth(1.5)  # the legend line width

    plt.grid(True)
    #plt.yscale('log')
    plt.autoscale(enable=True, axis='y')
    fig.set_size_inches(16,12)
    plt.savefig(ofilename,dpi=100)
    print "Output file ->",ofilename
    #plt.show()
    return energy

print ("Usage: %s filename.spe" %str(sys.argv[0]) )
global OUTDIR,filename_data,bFix
filename_data=str(sys.argv[1])
OUTDIR="./outsp"
if (not os.path.isdir(OUTDIR)):
    os.mkdir(OUTDIR,0755)
bFix=True
bFix=False
    
rebin_val=0.5
ppos=662
if "57Co" in filename_data:
    ppos=122.1    
elif "137Cs" in filename_data:
    ppos=661.7
elif "22Na" in filename_data:
    ppos=511
    #ppos=1274
elif "60Co" in filename_data:
    ppos=1173.2
    #ppos=1332.5
else:
    print "Specify the source peak in the filename"
    exit()
    
if os.path.isfile(filename_data):
    [path,fname,ext]= SplitFileName(filename_data)
    print path,fname,ext
else:
    print "file",filename_data," does not exists; exiting"
    exit()

[headvals,binvals,counts]=ReadSP(filename_data)

b2 = numpy.asarray(binvals, dtype = float)
b = scipy.ndimage.interpolation.zoom(b2, rebin_val)

c2 = numpy.asarray(counts, dtype = float)
c = scipy.ndimage.interpolation.zoom(c2, rebin_val)

e=FitSP_Peak(headvals,b,c,ext,ppos)
if not "spe" in ext:
    WriteSPE(headvals,e,c)



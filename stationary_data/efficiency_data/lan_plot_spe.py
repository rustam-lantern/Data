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

def add_lists(a,b):
    c=[a_i + b_i for a_i, b_i in zip(a, b)]
    return c

def ReturnSpList(path,bn,ext):
    rlist=[]
    ind_=find_ind(bn,"_")
    bnc=bn[:ind_[-1]]
    for (root, dirnames, filenames) in walk(path):
        for filename in fnmatch.filter(filenames,bnc+"*"+ext):
            rlist+=[path+"/"+filename]
    return rlist

def find_ind(lst, a): #Finding the indices of matching elements in list 
    return [i for i, x in enumerate(lst) if x==a]

def SplitFileName(fn):
    basename, ext = os.path.splitext(fn)
    path,fname = os.path.split(basename)
    if len(path)==0:
        path="."
    return path,fname,ext

def WriteSP(bin,counts,outfilename):
    filename  = open(outfilename,'w')
    for i in range(0,len(bins)):
        print >>filename,"%10.6f %10.6f"%(bin[i],counts[i])

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

def PlotSPE(h,b,c,ext):
    
    #b1=numpy.array(b1)
    #c1=numpy.array(c1)
    #cond = (b1>=600) & (b1<=800)
    #b = b1[ cond ]
    #c = c1[ cond ]
    
    ofilename=OUTDIR+"/"+bsfile[0]+".png"
    sh_marker=["+",".","^","*","p","s","x","D","h","o","+",".","^","*","p","s","x","D","h","o","+",".","^","*","p","s"]
    sh_color=['b','g','r','c','m','y','k','orange','darkblue','darkgreen','darkred','darkcyan','darkmagenta','deeppink','firebrick','gold','darkviolet','lavenderblush','dodgerblue','indigo','limegreen']
    #gmodel = Model(gaussian) + Model(line)
    #gmodel = Model(gaussian) + Model(expo)
    gmodel = Model(fit_function)
    gmodel.set_param_hint('A', value=10000)
    gmodel.set_param_hint('beta', value=0.001)
    gmodel.set_param_hint('B', value=500)
    gmodel.set_param_hint('mu', value=660)
    gmodel.set_param_hint('sigma', value=30)
    
    pars = gmodel.make_params()
    result = gmodel.fit(c, pars,x=b)
    print(result.fit_report())
    #print result.params
    print "%6.2f %6.4f"%(result.params["mu"].value,result.params["sigma"].value)

    fig = plt.gcf()
    plt.title(h[0],fontsize=14,fontweight='bold')
    if "spe" in ext:
        plt.xlabel('Energy ',fontsize=14,fontweight='bold')
    else:
        plt.xlabel('Bins ',fontsize=14,fontweight='bold')
    plt.ylabel('Counts ',fontsize=14,fontweight='bold')
    plt.tick_params(labelsize=18)
    #plt.text(energyBins[binid][0]+10, 0.8*maxv2,strv,fontsize=14,fontweight='bold')

    plt.plot(b,c,marker=sh_marker[0],linestyle='None',color=sh_color[0],label=bsfile)
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
    plt.show()


    
print ("Usage: %s spfile1 spfile2 ... " %str(sys.argv[0]) )
NARG=(len(sys.argv)-1)
fn=str(sys.argv[1])
global bsfile
spfile=[]
bsfile=[]

if os.path.isfile(fn):
    [path,fname,ext]= SplitFileName(fn)
    print path,fname,ext
else:
    print "file",fn," does not exists; exiting"
    exit()


for i in range(0,NARG):
    j=i+1
    spfile+=[sys.argv[j]]
global OUTDIR
OUTDIR="./outsp"
if (not os.path.isdir(OUTDIR)):
    os.mkdir(OUTDIR,0755)

#rebin_val=0.5
rebin_val=1
headvals = []
binvals = []
counts = []
for i in range(0,len(spfile)):
    [header,bins,cc]=ReadSP(spfile[i])
    bsfile+=[os.path.basename(spfile[i])]
    headvals+=[header]
    binvals+=bins
    if i==1:
        #cc1=[i * 10 for i in cc]
        counts+= cc
    else:
        counts+= cc
print len(counts)
#PlotSP(headvals,binvals,counts,ext)

#b2 = numpy.asarray(binvals, dtype = float)
#b = scipy.ndimage.interpolation.zoom(b2, rebin_val)

#c2 = numpy.asarray(counts, dtype = float)
#c = scipy.ndimage.interpolation.zoom(c2, rebin_val)
b=binvals
c=counts

print len(b),len(c)
PlotSPE(headvals,b,c,ext)

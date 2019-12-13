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

global energyBins

energyBins = [[30,52],
[52,70],
[70,94],
[85,112],
[112,139],
[125,147],
[147,172],
[172,199],
[199,222],
[222,248],
[248,281],
[281,319],
[299,334],
[334,380],
[380,434],
[434,490],
[490,545],
[545,609],
[578,634],
[634,691],
[691,767],
[767,850],
[850,963],
[963,1039],
[1039,1130],
[1130,1215],
[1215,1287],
[1287,1379],
[1379,1510],
[1510,1635],
[1635,1895],
[1895,2537],
[2537,2691]]


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

def GetLiveTime(h):
    NC=int(h.split()[2])
    RT=float(h.split()[3])
    LT=float(h.split()[4])

    return NC,RT,LT

def ConvertSPE_SB(energy,c):
    bincounts = []
    for binid in range(0,len(energyBins)):
        #bincounts += [(sum(float(c[i]) for i, x in enumerate(energy) if (x > energyBins[binid][0]) & (x < energyBins[binid][1])))]
        bincounts += [(sum(c[i] for i, x in enumerate(energy) if (x > energyBins[binid][0]) & (x <= energyBins[binid][1])))]
    #print sum(c),sum(bincounts)
    return bincounts


def ConvertSPE_SBR(energy,c,LT):
    bincounts = []
    for binid in range(len(energyBins)):
        bincounts += [(sum(float(c[i]) for i, x in enumerate(energy) if (x > energyBins[binid][0]) & (x <= energyBins[binid][1])))/LT]
    return bincounts

def WriteSPBCR(c):
    ofilename=OUTDIR+"/"+filename_data+".sbt"
    filename  = open(ofilename,'w')
    for i in range(0,len(c)):
        print >>filename,"%2d %10.6f"%(i+1,c[i])
    print "Output file ->",ofilename   

def SplitFileName(fn):
    basename, ext = os.path.splitext(fn)
    path,fname = os.path.split(basename)
    if len(path)==0:
        path="."
    return path,fname,ext

def ReturnHeaderVars(h):
    
    TS=h.split()[0]
    FN=h.split()[1]
    NC=int(h.split()[2])
    RT=float(h.split()[3])
    LT=float(h.split()[4])
    a0=float(h.split()[5])
    a1=float(h.split()[6])
    a2=float(h.split()[7])
    return TS,FN,NC,RT,LT,a0,a1,a2

def ReadSPE(infilename):
    binv = []
    countv = []
    count_rates = []
    
    f = open(infilename, 'r')
    header = f.readline()
    LT=float(header.split()[4])
    for line in f:
        line = line.strip()
        columns = line.split()
        binv += [float(columns[0])]
        countv += [float(columns[1])]
        count_rates += [float(columns[1])/LT]
    #countv[-1] = 1
    return header,binv,countv,count_rates

def PlotSPE(h,b,c,s):
    ofilename=OUTDIR+"/"+filename_data+"_1.png"
    sh_marker=["+",".","^","*","p","s","x","D","h","o","+",".","^","*","p","s","x","D","h","o","+",".","^","*","p","s"]
    sh_color=['b','g','r','c','m','y','k','orange','darkblue','darkgreen','darkred','darkcyan','darkmagenta','deeppink','firebrick','gold','darkviolet','lavenderblush','dodgerblue','indigo','limegreen']
    c[-1]=1
    bsum=sum(c)
    fig, ax = plt.subplots(2,1)
    ax[0].plot(b,c,marker=sh_marker[1],linestyle='None',color=sh_color[0],label=filename_data)
    ax[0].set_title(filename_data,fontsize=14,fontweight='bold')
    #ax[0].plot(b,c,marker=sh_marker[1],linestyle='None',color=sh_color[0],label=filename_data)
    ax[0].set_ylabel('Counts ',fontsize=14,fontweight='bold')
    ax[0].grid(True)
    ax[0].tick_params(labelsize=18)
    ax[0].set_xlabel('Energy ',fontsize=14,fontweight='bold')
    ax[1].plot(list(range(1, 34)),s,marker=sh_marker[1],linestyle='-',color=sh_color[0],markersize=14,label=filename_data)
    ax[1].set_xlabel('Energy Bins ',fontsize=14,fontweight='bold')
    ax[1].set_ylabel('Counts ',fontsize=14,fontweight='bold')
    ax[1].grid(True)
    #legend = plt.legend(loc='upper right', shadow=True)
    #frame = legend.get_frame()
    #frame.set_facecolor('0.90')
    #for label in legend.get_texts():
    #    label.set_fontsize('large')
    #for label in legend.get_lines():
    #    label.set_linewidth(1.5)  # the legend line width
    plt.tick_params(labelsize=18)
    plt.grid(True)
    plt.autoscale(enable=True, axis='y')
    fig.set_size_inches(16,12)
    plt.savefig(ofilename,dpi=100)
    plt.show()


def PlotSPEC(h,b,c,s):
    ofilename=OUTDIR+"/"+filename_data+"_2.png"
    sh_marker=["+",".","^","*","p","s","x","D","h","o","+",".","^","*","p","s","x","D","h","o","+",".","^","*","p","s"]
    sh_color=['b','g','r','c','m','y','k','orange','darkblue','darkgreen','darkred','darkcyan','darkmagenta','deeppink','firebrick','gold','darkviolet','lavenderblush','dodgerblue','indigo','limegreen']
    c[-1]=1
    bsum=sum(c)
    fig, ax = plt.subplots(2,1)
    ax[0].plot(b,c,marker=sh_marker[1],linestyle='None',color=sh_color[0],label=filename_data)
    ax[0].set_title(filename_data,fontsize=14,fontweight='bold')
    #ax[0].plot(b,c,marker=sh_marker[1],linestyle='None',color=sh_color[0],label=filename_data)
    ax[0].set_ylabel('Count Rates (cps) ',fontsize=14,fontweight='bold')
    ax[0].grid(True)
    ax[0].tick_params(labelsize=18)
    ax[0].set_xlabel('Energy ',fontsize=14,fontweight='bold')
    ax[1].plot(list(range(1, 34)),s,marker=sh_marker[1],linestyle='-',color=sh_color[0],markersize=14,label=filename_data)
    ax[1].set_xlabel('Energy Bins ',fontsize=14,fontweight='bold')
    ax[1].set_ylabel('Count Rates (cps)',fontsize=14,fontweight='bold')
    ax[1].grid(True)
    plt.tick_params(labelsize=18)
    plt.grid(True)
    plt.autoscale(enable=True, axis='y')
    fig.set_size_inches(16,12)
    plt.savefig(ofilename,dpi=100)
    plt.show()

def FitSP_Peak(h,b1,c,ext,ppos):
    global cf
    a0=float(h.split()[5])
    a1=float(h.split()[6])
    #cf=a1
    a2=float(h.split()[7])
    if "spe" in ext:
        b = b1
    else:
        b = [return_energy(b1[i],a0,a1,a2) for i in range(len(b1))]
    print "pars: ",a0," ",a1," ",a2
    b2=numpy.array(b)
    c2=numpy.array(c)
    cond = (b2>=ppos-70) & (b2<=ppos+70)
    #b = b2[ cond ]
    #c = c2[ cond ]
    
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
    cfo=ppos*a1/result.params["mu"].value
    mval=result.params["mu"].value
    sval=result.params["sigma"].value      
    print "%6.2f %6.4f %8.6f -> %8.6f"%(mval,sval,a1,ppos*a1/result.params["mu"].value)
    if bFix:
        cf=cfo
    else:
        cf=a1
        
    print "%6.2f %6.4f %8.6f -> %8.6f : %8.6f (%4.2f %% off)"%(mval,sval,a1,cf,cfo,100*(cf-cfo)/cfo)
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
    plt.plot(b, result.best_fit, 'r-', label="fit %6.2f %6.4f"%(mval,sval))
    plt.xlim([mval-5*sval,mval+5*sval])

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
    plt.show()
    return energy

print ("Usage: %s filename.spe" %str(sys.argv[0]) )
global OUTDIR,filename_data,NC,RT,LT
filename_data=str(sys.argv[1])
OUTDIR="./outsp"
if (not os.path.isdir(OUTDIR)):
    os.mkdir(OUTDIR,0755)
if os.path.isfile(filename_data):
    [path,fname,ext]= SplitFileName(filename_data)
    print path,fname,ext
else:
    print "file",filename_data," does not exists; exiting"
    exit()


[headvals,binvals,counts,count_rates]=ReadSPE(filename_data)
[TS,FN,NC,RT,LT,a0,a1,a2]=ReturnHeaderVars(headvals)
sbins=ConvertSPE_SB(binvals,counts)
sbins_rates=ConvertSPE_SBR(binvals,counts,LT)
print NC,sum(counts),sum(sbins),RT,LT
PlotSPE(headvals,binvals,counts,sbins)
PlotSPEC(headvals,binvals,count_rates,sbins_rates)
WriteSPBCR(sbins_rates)

#b2 = numpy.asarray(binvals, dtype = float)
#b = scipy.ndimage.interpolation.zoom(b2, rebin_val)

#c2 = numpy.asarray(counts, dtype = float)
#c = scipy.ndimage.interpolation.zoom(c2, rebin_val)

#e=FitSP_Peak(headvals,b,c,ext,ppos)
#if not "spe" in ext:
#    WriteSPE(headvals,e,c)



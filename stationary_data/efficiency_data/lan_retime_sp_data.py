#!/usr/bin/python
import sys
from os import walk
import os
import pickle
import subprocess
import shutil
import glob
import fnmatch
from subprocess import call
from subprocess import Popen, PIPE
from matplotlib.dates import date2num 
import datetime as DT
import numpy
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from os.path import basename
from os.path import splitext
from operator import truediv
import scipy.ndimage

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



def add_lists(a,b):
    c=[a_i + b_i for a_i, b_i in zip(a, b)]
    return c

def find_ind(lst, a): #Finding the indices of matching elements in list 
    return [i for i, x in enumerate(lst) if x==a]

def count_freq_occurrence(list):
    d = {x:list.count(x) for x in list}
    return d.keys(),d.values()

#def rebin_data(arr,rfc):
#    b2 = numpy.asarray(arr, dtype = float)
#    b = scipy.ndimage.interpolation.zoom(b2, rdc)
#    return b

def return_energy(b,a0,a1,a2):
    en=b*b*a2+b*a1+a0
    return en

def Fill_Bins(bins,counts):
    counter = 0
    for i in range(0,4096):
        if i in bins:
            print i+1,counts[counter]
            counter +=1
        else:
            print i+1,0
        #for j in range(0,len(bins)):
        #    print i,bins[j],counts[j]

    #print counter
    
def Save_Sp_Data(header,bins,counts):
    outfile=OUTDIR+"/"+filename_data+".sp"
    ofilename  = open(outfile,'w')
    print >>ofilename, "%s"%(header),
    counter = 0
    for i in range(0,4096):
        if i in bins:
            #print i+1,counts[counter]
            if i==0 or i==4095:
                print >>ofilename,"%4d %d"%(i+1,0)
            else:
                print >>ofilename,"%4d %d"%(i+1,counts[counter])
            counter +=1
        else:
            #print i+1,0
            print >>ofilename,"%4d 0"%(i+1)

def Read_Data(infilename,integration_time):
    global header, timestamp,msectime,realtime,livetime,channel,energy,LT,RT,a0,a1,a2
    timestamp=[]
    msectime=[]
    realtime=[]
    livetime=[]
    channel=[]
    energy =[]
    f = open(infilename, 'r')
    header = f.readline()
    for num, line in enumerate(f, 1):
        line = line.strip()
        columns = line.split()
        timestamp += [columns[0]]
        msectime += [int(columns[1])]
        realtime += [int(columns[2])]
        livetime += [int(columns[3])]
        channel += [int(columns[4])]
        energy  += [float(columns[5])]
        if int(columns[2])>int(integration_time*1000):
            break  
        
    RT=float(realtime[-1])/1000.
    LT=float(livetime[-1])/1000.
    #print header
    a0=float(header.split()[2])
    a1=float(header.split()[3])
    a2=float(header.split()[4])
    print "Read data file:",infilename," entries:",num

def timeSampleIndex(time,nsamples,tsampling):
    index_range = [[0]*(2)]*(int(nsamples/tsampling))
    icount =0
    for isample in range(0,nsamples,tsampling):
        index_list = [i for i, x in enumerate(time) if (x > icount*tsampling) & (x <= (icount+1)*tsampling)]
        #print icount,icount*tsampling,(icount+1)*tsampling
        index_range[icount] = [index_list[0],index_list[-1]]
        icount+=1

    return index_range

def bin_data():
    bincounts = []
    for binid in range(len(energyBins)):
        bincounts += [(sum(1 for i, x in enumerate(energy) if (x > energyBins[binid][0]) & (x <= energyBins[binid][1])))]
    return bincounts

def bin_dataC():
    bincounts = []
    for binid in range(len(energyBins)):
        bincounts += [(sum(1 for i, x in enumerate(energy) if (x > energyBins[binid][0]) & (x <= energyBins[binid][1])))/LT]
    return bincounts


def time_dataC(sampling_time,binid):#count rates
    time = [float(msectime[x] - msectime[0])/1000. for x in range(len(channel))]
    nsamples=int(RT/sampling_time)
    index_range=timeSampleIndex(time,int(RT),sampling_time)
    bincounts_time = []
    ntime = []
    for isample in range(nsamples):
        t1= time[index_range[isample][0]]
        t2= time[index_range[isample][1]]
        lt_t = [float(livetime[i])/float(realtime[i]) for i in range(len(time)) if (time[i] >= t1) & (time[i] <= t2) &  (energy[i] > energyBins[binid][0]) &   (energy[i] <= energyBins[binid][1])] #livetime for each time point/bin
        if sum(lt_t)>0:
            for i in range(len(lt_t)):
                if lt_t[i]==0:
                    lt_t[i]=0.001
            bincounts_time += [sum(1./float(lt_t[i]) for i in range(len(lt_t)))/sampling_time]
            #bincounts_time += [sum(1./float(lt_t[i]) for i in range(len(lt_t)))]
        else:
            bincounts_time += [0]
             
        #print isample,t1,t2,t1+(t2-t1)/2,(sum(1 for i in range(len(time)) if (time[i] >= t1) & (time[i] <= t2) &  (energy[i] > energyBins[binid][0]) &   (energy[i] <= energyBins[binid][1]) )),lt_t
        ntime += [t1+(t2-t1)/2]
        
    return ntime,bincounts_time


def time_data(sampling_time,binid):#Just count
    time = [float(msectime[x] - msectime[0])/1000. for x in range(len(channel))]
    nsamples=int(RT/sampling_time)
    index_range=timeSampleIndex(time,int(RT),sampling_time)
    bincounts_time = []
    ntime = []
    for isample in range(nsamples):
        t1= time[index_range[isample][0]]
        t2= time[index_range[isample][1]]
        bincounts_time += [(sum(1 for i in range(len(time)) if (time[i] >= t1) & (time[i] <= t2) &  (energy[i] > energyBins[binid][0]) &   (energy[i] <= energyBins[binid][1]) ))]
        #print isample,t1,t2,t1+(t2-t1)/2
        ntime += [t1+(t2-t1)/2] 
    return ntime,bincounts_time

def time_data_old(sampling_time,binid):
    time = [float(msectime[x] - msectime[0])/1000. for x in range(len(channel))]
    ntime = []
    nsamples=int(RT/sampling_time)
    index_range=timeSampleIndex(time,nsamples,sampling_time)
    bincounts_time= [[0]*(len(energyBins))]*nsamples#2D group indexing
    for binid in range(len(energyBins)):
        bincounts += [(sum(1 for i, x in enumerate(channel) if (x > energyBins[binid][0]) & (x <= energyBins[binid][1])))/LT]
        tmpcounts =[]
        for isample in range(nsamples):
            tmpcounts += [(sum(1 for i in range(len(time)) if (time[i] >= time[index_range[isample][0]]) & (time[i] <= time[index_range[isample][1]]) &  (energy[i] > energyBins[binid][0]) &   (energy[i] <= energyBins[binid][1]) ))]
        bincounts_time[binid] = tmpcounts
    
    return time,bincounts,bincounts_time
        

def PlotSPE(h,b,c,ext):
    ofilename=OUTDIR+"/"+filename_data+".png"
    sh_marker=["+",".","^","*","p","s","x","D","h","o","+",".","^","*","p","s","x","D","h","o","+",".","^","*","p","s"]
    sh_color=['b','g','r','c','m','y','k','orange','darkblue','darkgreen','darkred','darkcyan','darkmagenta','deeppink','firebrick','gold','darkviolet','lavenderblush','dodgerblue','indigo','limegreen']

    bsum=[0]*4096
    c[-1]=1
    bsum=sum(c)
    energy = [return_energy(b[i],a0,a1,a2) for i in range(len(b))]
    fig, ax = plt.subplots(2,1)
    ax[0].plot(b,c,marker=sh_marker[1],linestyle='None',color=sh_color[0],label=filename_data)
    ax[0].set_title(h+" %d"%(bsum)+" %6.1f"%(realtime[-1]/1000.)+ " %6.1f"%(livetime[-1]/1000.),fontsize=14,fontweight='bold')
    ax[0].plot(b,c,marker=sh_marker[1],linestyle='None',color=sh_color[0],label=filename_data)
    ax[0].set_ylabel('Counts ',fontsize=14,fontweight='bold')
    ax[0].grid(True)
    ax[0].tick_params(labelsize=18)
    ax[0].set_xlabel('Bins ',fontsize=14,fontweight='bold')
    ax[1].plot(energy,c,marker=sh_marker[1],linestyle='None',color=sh_color[0],label=filename_data)
    ax[1].set_xlabel('Energy ',fontsize=14,fontweight='bold')
    ax[1].set_ylabel('Counts ',fontsize=14,fontweight='bold')
    ax[1].grid(True)
    legend = plt.legend(loc='upper right', shadow=True)
    frame = legend.get_frame()
    frame.set_facecolor('0.90')
    for label in legend.get_texts():
        label.set_fontsize('large')
    for label in legend.get_lines():
        label.set_linewidth(1.5)  # the legend line width
    plt.tick_params(labelsize=18)
    plt.grid(True)
    plt.autoscale(enable=True, axis='y')
    fig.set_size_inches(16,12)
    plt.savefig(ofilename,dpi=100)
    
    #plt.show()
     
 
def PlotBinData(bindata):
    ofilename=OUTDIR+"/"+filename_data+"_bin.png"
    sh_marker=["+",".","^","*","p","s","x","D","h","o","+",".","^","*","p","s","x","D","h","o","+",".","^","*","p","s"]
    sh_color=['b','g','r','c','m','y','k','orange','darkblue','darkgreen','darkred','darkcyan','darkmagenta','deeppink','firebrick','gold','darkviolet','lavenderblush','dodgerblue','indigo','limegreen']

    fig = plt.figure()
    plt.xlabel('Bins ',fontsize=14,fontweight='bold')
    plt.ylabel('Count Rates (cps) ',fontsize=14,fontweight='bold')
    plt.tick_params(labelsize=18)
    plt.plot(list(range(1, 34)),bindata,marker=sh_marker[1],linestyle='-',color=sh_color[0],label=filename_data,markersize=14)
    #plt.title(h+" %d"%(bsum)+" %6.1f"%(realtime[-1]/1000.)+ " %6.1f"%(livetime[-1]/1000.),fontsize=14,fontweight='bold')
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
    #plt.show()

def PlotTimeData(time,bindata,sampling_time,binid):
    ofilename=OUTDIR+"/"+filename_data+"_tbin.png"
    sh_marker=["+",".","^","*","p","s","x","D","h","o","+",".","^","*","p","s","x","D","h","o","+",".","^","*","p","s"]
    sh_color=['b','g','r','c','m','y','k','orange','darkblue','darkgreen','darkred','darkcyan','darkmagenta','deeppink','firebrick','gold','darkviolet','lavenderblush','dodgerblue','indigo','limegreen']

    fig = plt.figure()
    plt.xlabel('Time ',fontsize=14,fontweight='bold')
    plt.ylabel('Count Rates (cps) ',fontsize=14,fontweight='bold')
    plt.tick_params(labelsize=18)
    plt.plot(time,bindata,marker=sh_marker[1],linestyle='None',color=sh_color[0],label=filename_data)
    plt.title(" Bin=%d ;"%(binid+1)+" SamplingTime=%6.1fs"%(sampling_time),fontsize=14,fontweight='bold')
    legend = plt.legend(loc='upper right', shadow=True)
    frame = legend.get_frame()
    frame.set_facecolor('0.90')
    for label in legend.get_texts():
        label.set_fontsize('large')
    for label in legend.get_lines():
        label.set_linewidth(1.5)  # the legend line width

    plt.grid(True)
    #plt.yscale('log')
    plt.ylim([0,20])
    #plt.autoscale(enable=True, axis='y')
    fig.set_size_inches(16,12)
    plt.savefig(ofilename,dpi=100)
    #plt.show()
    
print ("Usage: %s filename integration_time" %str(sys.argv[0]) )

global OUTDIR,filename_data
OUTDIR="./outsp"
if (not os.path.isdir(OUTDIR)):
    os.mkdir(OUTDIR,0755)


filename_data=str(sys.argv[1])
ext=".sp"
if len (sys.argv) == 3 :
    integration_time=int(sys.argv[2])
else:
    integration_time=-1
print  integration_time
Read_Data(filename_data,integration_time)
[bins,counts]=count_freq_occurrence(channel)
print "Finished Event Counting"
Save_Sp_Data(header,bins,counts)
#PlotSPE(header,bins,counts,ext)





#!/usr/bin/python
import pylab as plb
import matplotlib.pyplot as plt
import math as math
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
#data = plb.loadtxt('part 2.csv', skiprows=2, delimiter=',')


x=ar([-2.5, -2.4,   -2.3,   -2.2,   -2.1,   -2, -1.9,   -1.8,   -1.7,   -1.6,   -1.5,   -1.4,   -1.3,   -1.2,   -1.1,   -1, -0.9,   -0.8,   -0.7,   -0.6,   -0.5,   -0.4,   -0.3,   -0.2,   -0.1,   0,  0.1,    0.2,    0.3,    0.4,    0.5,    0.6,    0.7,    0.8,    0.9,    1,  1.1,    1.2,    1.3,    1.4,    1.5,    1.6,    1.7,    1.8 ,1.9,   2,  2.1,    2.2,    2.3,    2.4,    2.5])
y=ar([2.98317E-6,   2.96168E-6, 2.9548E-6,  2.9535E-6,  2.95508E-6, 2.95839E-6, 2.96301E-6, 2.96861E-6, 2.97448E-6, 2.98083E-6, 2.98819E-6, 2.99527E-6, 3.00295E-6, 3.01103E-6, 3.01969E-6, 3.0286E-6,  3.03798E-6, 3.04744E-6, 3.05716E-6, 3.06697E-6, 3.0769E-6,  3.08716E-6, 3.09771E-6, 3.10858E-6, 3.11936E-6, 3.13004E-6, 3.14179E-6, 0.000003153,    3.16427E-6, 3.17543E-6, 3.18663E-6, 3.19628E-6, 3.2053E-6,  3.21354E-6, 3.21884E-6, 3.22233E-6, 3.22257E-6, 3.21947E-6, 3.21332E-6, 3.20529E-6, 3.19538E-6, 3.18451E-6, 3.17442E-6, 3.16384E-6, 3.15269E-6, 3.14036E-6, 3.12767E-6, 3.11465E-6, 3.10153E-6, 3.08893E-6, 3.07638E-6])


n=sum(y)
mean = sum(x*y)/n                   
sigma=math.sqrt(sum(y*(x-mean)**2)/n)

def gaus(x,a,x0,sigma, offset):
    return a*exp(-(x-x0)**2/(2*sigma**2)) + offset

popt,pcov = curve_fit(gaus,x,y,p0=[1,mean,sigma, 0.0])
print(popt)

plt.plot(x,y,'g+:',label='data')
plt.plot(x,gaus(x,*popt),'b*:',label='fit')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Position (m)')
plt.show()

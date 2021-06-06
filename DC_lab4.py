#! /usr/bin/python3
import sys
import math
import numpy as np
import scipy.stats
import scipy.special as sp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import warnings
warnings.filterwarnings("ignore")
import modulo as mtd

                                ##  DATA    ##

CC=.22 #\muF+-10%
AA=.26**2/4*np.pi #m^2

                                ## PARTE 1 ##
d1=.20 #cm
D1=np.array([[0.5,1.0,1.5,2.0,2.5,3.0],[.61,1.12,1.64,2.17,2.64,3.15]]).T #U_c(kV)//Uair(V)

Uc2=1.5 #kV
D2=np.array([[.10,.15,.20,.25,.30,.35],[2.80,2.28,1.74,1.40,1.24,1.10]]).T	#d(cm)//Uair(V)

                                ## PARTE 2 ##

# plástico
d3=.98 #cm
D3=np.array([[.5,1.0,1.5,2.0,2.5,3.0],[.60,1.00,1.46,1.78,2.23,2.68],[.15,.32,.52,.61,.81,.96]]).T	#Uc(kV)//Upl(V)//Uair(V)

# vidrio
d4=.27 #cm
D4=np.array([[.5,1.0,1.5,2.0],[2.50,3.13,3.57,4.69],[.58,1.02,1.43,1.73]]).T	#Uc(kV)//Uvidrio(V)//Uair(V)
#D4=np.array(([[.5,1.0,1.5,2.0,2.5],[2.75,6.77,10.018,10.685,10.7],[1.232,1.213,2.126,2.665,2.913]])).T
#D4=np.array([[.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0],[.5,.92,1.35,1.8,2.3,2.8,3.1,3.7],[.16,.32,.51,.62,.78,.95,1.12,1.3]]).T

###################################################################################
###################################################################################

def lmsq(M,b):
    qK=sum(M[:,0]*M[:,b])/sum(M[:,0]**2)
    ym=sum(M[:,b])/len(M)
    st=sum((M[:,b]-ym)**2)
    sr=sum((M[:,b]-qK*M[:,0])**2)
    return qK,((st-sr)/st)

#1
P1=np.zeros((len(D1),3))
for i in range(len(D1)):
    P1[i,0]=D1[i,0]
    P1[i,1]=D1[i,1]*CC*1000
    P1[i,2]=P1[i,1]/D1[i,0]/AA*d1
#print(P1)
cLT1=mtd.lingen(P1,mtd.pl1,1)
sLT1=lmsq(P1,1)


#2
P2=np.zeros((len(D2),3))
for i in range(len(D2)):
    P2[i,0]=1/D2[i,0]
    P2[i,1]=D2[i,1]*CC*1000
    P2[i,2]=P2[i,1]/Uc2/AA*D2[i,0]
#print(P2)
cLT2=mtd.lingen(P2,mtd.pl1,1)
sLT2=lmsq(P2,1)


#3
P3=np.zeros((len(D3),5))
for i in range(len(D3)):
    P3[i,0]=D3[i,0]
    P3[i,1]=D3[i,1]*CC*1000
    P3[i,2]=P3[i,1]/D3[i,0]/AA*d3
    P3[i,3]=D3[i,2]*CC*1000
    P3[i,4]=P3[i,3]/D3[i,0]/AA*d3
#print(P3)
sLT3a=lmsq(P3,1)
sLT3b=lmsq(P3,3)
cLT3a=mtd.lingen(P3,mtd.pl1,1)
cLT3b=mtd.lingen(P3,mtd.pl1,3)


#4
P4=np.zeros((len(D4),5))
for i in range(len(D4)):
    P4[i,0]=D4[i,0]
    P4[i,1]=D4[i,1]*CC*1000
    P4[i,2]=P4[i,1]/D4[i,0]/AA*d4
    P4[i,3]=D4[i,2]*CC*1000
    P4[i,4]=P4[i,3]/D4[i,0]/AA*d4
#print(P4)
sLT4a=lmsq(P4,1)
sLT4b=lmsq(P4,3)
cLT4a=mtd.lingen(P4,mtd.pl1,1)
cLT4b=mtd.lingen(P4,mtd.pl1,3)




##  GRAFICAS    ##
plt.figure(1)
plt.plot(P1[:,0],P1[:,1],'o',markersize=2,color='C0')
plt.plot(P1[:,0],sLT1[0]*P1[:,0],'--',lw=.5,color='C0',label='$Q_a=%1.2f\,U_c, \, R^2=%1.4f$'%(sLT1[0],sLT1[1]))
plt.ylabel('$Q \; [nAs]$')
plt.xlabel('$U_c \; [kV]$')
lines1=plt.gca().get_lines()
legend1=plt.legend([lines1[i] for i in [0]],['Aire'],title='Data',loc=2)
plt.legend(loc=4)
plt.gca().add_artist(legend1)


plt.figure(2)
plt.plot(P2[:,0],P2[:,1],'o',markersize=2,color='C0')
plt.plot(P2[:,0],sLT2[0]*P2[:,0],'--',lw=.5,color='C0',label='$Q_a=%1.2f\,d^{-1}, \, R^2=%1.4f$'%(sLT2[0],sLT2[1]))
plt.ylabel('$Q \; [nAs]$')
plt.xlabel('$d^{-1} \; [cm^{-1}]$')
lines2=plt.gca().get_lines()
legend2=plt.legend([lines2[i] for i in [0]],['Aire'],title='Data',loc=2)
plt.legend(loc=4)
plt.gca().add_artist(legend2)


plt.figure(3)
plt.plot(P3[:,0],P3[:,1],'o',markersize=2,color='C1')
plt.plot(P3[:,0],P3[:,3],'o',markersize=2,color='C0')
plt.plot(P3[:,0],sLT3a[0]*P3[:,0],'--',lw=.5,color='C1',label='$Q_{p}=%1.2f\, U_c, \, R^2=%1.4f$'%(sLT3a[0],sLT3a[1]))
plt.plot(P3[:,0],sLT3b[0]*P3[:,0],'--',lw=.5,color='C0',label='$Q_{a}=%1.2f\, U_c, \, R^2=%1.4f$'%(sLT3b[0],sLT3b[1]))
#plt.plot(P3[:,0],cLT3a[0][1]*P3[:,0]+cLT3a[0][0],'--',lw=.5,color='C1',label='$Q_{p}=%1.2f\,U_c+ %1.2f, \, R^2=%1.4f$'%(cLT3a[0][1],cLT3a[0][0],cLT3a[1][0]))
#plt.plot(P3[:,0],cLT3b[0][1]*P3[:,0]+cLT3b[0][0],'--',lw=.5,color='C0',label='$Q_{a}=%1.2f\,U_c+ %1.2f, \, R^2=%1.4f$'%(cLT3b[0][1],cLT3b[0][0],cLT3b[1][0]))
plt.ylabel('$Q \; [nAs]$')
plt.xlabel('$U_c \; [kV]$')
lines3=plt.gca().get_lines()
legend3=plt.legend([lines3[i] for i in [0,1]],['Plástico','Aire'],title='Data',loc=2)
plt.legend(loc=7)
plt.gca().add_artist(legend3)


plt.figure(4)
plt.plot(P4[:,0],P4[:,1],'o',markersize=2,color='C3')
plt.plot(P4[:,0],P4[:,3],'o',markersize=2,color='C0')
#plt.plot(P4[:,0],sLT4a[0]*P4[:,0],'--',lw=.5,color='C3',label='$Q_{v}=%1.2f\, U_c, \, R^2=%1.4f$'%(sLT4a[0],sLT4a[1]))
plt.plot(P4[:,0],sLT4b[0]*P4[:,0],'--',lw=.5,color='C0',label='$Q_{a}=%1.2f\, U_c, \, R^2=%1.4f$'%(sLT4b[0],sLT4b[1]))
plt.plot(P4[:,0],cLT4a[0][1]*P4[:,0]+cLT4a[0][0],'--',lw=.5,color='C3',label='$Q_{v}=%1.2f\,U_c+ %1.2f, \, R^2=%1.4f$'%(cLT4a[0][1],cLT4a[0][0],cLT4a[1][0]))
#plt.plot(P4[:,0],cLT4b[0][1]*P4[:,0]+cLT4b[0][0],'--',lw=.5,color='C0',label='$Q_{a}=%1.2f\,U_c+ %1.2f, \, R^2=%1.4f$'%(cLT4b[0][1],cLT4b[0][0],cLT4b[1][0]))
plt.ylabel('$Q \; [nAs]$')
plt.xlabel('$U_c \; [kV]$')
lines4=plt.gca().get_lines()
legend4=plt.legend([lines4[i] for i in [0,1]],['Vidrio','Aire'],title='Data',loc=2)
plt.legend(loc=7)
plt.gca().add_artist(legend4)


#legend1=plt.legend([lines1[i] for i in [1]],['$\dfrac{Q}{U_c}=\epsilon_0\cdot A \cdot \dfrac{1}{d}$'],loc=5)


plt.show()

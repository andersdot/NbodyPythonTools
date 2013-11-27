import scipy as sp
#from scipy import trapz
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab as pl
from pylab import *
import sys

#ApJ Komatsu et al. 2009
omega_M=.274
omega_L=.726
omega_K=1.-omega_M-omega_L
omega_K=0
h=.705
age=13.72

"""
print ' *** *** *** *** *** *** *** *** *** *** *** *** '
print '                 Welcome to Cosmography          '
print '             Writen by Alex Fry 2009             '
print '    Please visit http://www.astro.washington.edu/users/bastidas/cosmography.html for info                      '

def Variation():
   def Input(min,max,msg,error_msg1,error_msg2):
        while True:
                try:
                        value=max+666
                        while value<min or value>max:
                                value=float(raw_input(msg))
                                if value<min:
                                        print error_msg1
                                if value>max:
                                        print error_msg2
                        return value
                        #break
                except ValueError:
                        print "Not valid input, must be an integer. Try again: "


   h=Input(0,100,'Enter the value of Hubble constant in units of ?','Please, at least 0... ','Please, less than 1000...')
   Omega_M=Input(0,1,'Enter the value of Omega matter: ','Please, at least 0... ','Please, less than unity...')
   Omega_L=Input(0,1-Omega_M,'Enter the value of Omega Lambda: ','Please, at least 0... ','Please, less than unity-Omega_M...')
   Omega_K=1.-Omega_M-Omega_L
   return Omega_M,Omega_L,Omega_K

ask=1
while (ask==1):
   temp=str(raw_input("Do you want a standard LCDM cosmology (yes/no)?"))
   if (temp=="n" or temp=="no"):
   	omega_M,omega_L,omega_K=Variation()
	ask=666
   if (temp=="y" or temp=="yes"):
	ask=666
   else: print "Try again..."
"""
class Cosmology(object):
   "The parameters of your cosmology"
   def __init__(self,omega_M,omega_L,omega_K,h):
	self.omega_M=omega_M
	self.omega_L=omega_L	  
	self.omega_K=omega_K	  
	self.Ho=100.*h #h km s^-1 Mpc^-1
	self.c=3.0*np.power(10.,5) # km/s
	self.Dh=self.c/self.Ho

   def Comoving_Distance(self,z):
	#print self.Ho, "is Ho. \n"
	Dh=self.Dh
	E_z=np.zeros(len(z))
	Dc=np.zeros(len(z))
	for n in range (len(z)):
	   E_z[n]=np.sqrt(omega_M*np.power( (1.+z[n]) ,3.)+omega_K*np.power( (1.+z[n]) ,2.)+omega_L)
	   Dc[n]=Dh*sp.trapz(  (1./E_z[0:n]) ,z[0:n])	
	   #Dc[n]=Dh*sp.integrate.simps( (1.E_z[0:n]) , z[0:n] ,dx=666,axis=-1,even='avg')
	if (omega_K>0):
	   Dm=Dh*( 1./np.sqrt(omega_K))*np.sinh(np.sqrt(omega_K)*(Dc/Dh) )
	   print "Omega k > 0 \n"	
	if (omega_K==0):
	   Dm=Dc
	   #print "Omega k = 0 \n"	
	if (omega_K<0):
	   Dm=Dh*( 1./np.sqrt(abs(omega_K)))*np.sin(np.sqrt(abs(omega_K))*(Dc/Dh) )
	   print "Omega k < 0 \n"	
	return Dm

   def Angular_Diameter_Distance(self,z,Dm):
	Da=Dm/(1.+z)
	#Da=Da/self.Dh
	return Da

   def Luminosity_Distance(self,z,Dm):
	Dl=(1.+z)*Dm
	#Dl=Dl/self.Dh
	return Dl

   def Lookback_Time(self,zt):	
	lenzt=len(zt)
	th=1/self.Ho
	E_zt=np.zeros(lenzt)
	tL=np.zeros(lenzt)
	t=np.zeros(lenzt)
	for n in range (lenzt):
	   E_zt[n]=np.sqrt(omega_M*np.power( (1.+zt[n]) ,3)+omega_K*np.power( (1.+zt[n]) ,2)+omega_L)
	   integrating_quantity=(1.+zt[0:n])*E_zt[0:n]
	   tL[n]=th*sp.trapz(  (1./integrating_quantity) ,zt[0:n])	
	return tL/th
"""
#z meshes to numerically integrating on
#join rough and fine to be more efficent and remain accurate
zt=np.arange(.001,10.0,.001)
z_fine=np.arange(0.0000,.1,.0001)
z_rough=np.arange(.1,10.,.001)
#print z_fine[len(z_fine)-1], "zf \n"
#print z_rough[0], "zf \n"
z=np.concatenate((z_fine,z_rough),axis=1)
cosmo=Cosmology(omega_M,omega_L,omega_K,h)
Dm=cosmo.Comoving_Distance(z)
Dm=Dm[2:len(Dm)]
z=z[2:len(z)]
Da=cosmo.Angular_Diameter_Distance(z,Dm)
Dl=cosmo.Luminosity_Distance(z,Dm)
tL=cosmo.Lookback_Time(zt)

###Luminoisty of an L* galaxy
#zeros in distance modulus creating errors, caused by numerical artifacts
#define a new Dl and remove zeros from it
DlM=Dl
#f=np.nonzero(DlM==0.)
#DlM[f]=DlM[  f[0][len(f[0])-1 ]+1] #make all zeros equal to first nonzero...
#print len(f[0]),"f0 \n"
dM=5.*np.log10(DlM*1000000./10.)

###Angle in arcseconds subtended by 1 kiloparsec
#Da=size/theta    theta is in radians here, 1 radian=57.2957795 deg
# 1 deg = 60 arc min or 3600 arcsec
#Da is in Mpc, multiply by 1000000 to get kpc
za=z
theta=    3600.*57.2957795*(1.)/(Da*1000000.)


###plot###
do_plot=1
plot_E=666
#mpl.rc('text', usetex=True)
#mpl.rc('font', family='serif')
#plt.figure(figsize=(8,10))
#mpl.rc('lines', lw=.1)
plt.subplots_adjust(wspace=0.19)

if (plot_E==1):
   subplot(111)
   print Dm, " dm \n"
   plt.semilogx(z,Dm,linewidth=.5)
   #plot(z,Dm,linewidth=.5)
   #plt.loglog(z,Dm,linewidth=.5)
   xlabel('Z')
   show()

if (do_plot==1):
 
   subplot(321)
   plt.loglog(z,Dl,linewidth=.5)
   #plot(z,Dl,linewidth=.5)
   xlabel('Z')
   ylabel('Luminosity Distance [Mpc]')
   grid(True)
   plt.xlim(xmin=.001)
  
   subplot(322)
   plt.semilogx(za,Da,linewidth=.5)
   #plt.loglog(z,Da,linewidth=.5)
   #plot(z,Da,linewidth=.5)
   xlabel('Z')
   ylabel('Angular Diamter Distance [Mpc]')
   grid(True)
   plt.xlim(xmin=.001)
   #plt.axis([x,x,y,y])
   
   subplot(323)
   #plot(z,dM-22.06,linewidth=.5)
   #plot(z,dM,linewidth=.5)
   plt.semilogx(z,dM-22.06,linewidth=.5)
   xlabel('Z')
   ylabel('Apparent L* Magnitude')
   grid(True)
   plt.xlim(xmin=.001)

   subplot(324)
   #plt.semilogx(za,theta,linewidth=.5)
   plt.loglog(za,theta,linewidth=.5)
   xlabel('Z')
   ylabel('Angle subtended by 1 kpc [arcseconds]')
   plt.xlim(xmin=.001)
   grid(True)
   #plt.axis([x,x,y,y])

   subplot(325)
   plt.semilogx(zt,age*tL,linewidth=.5)
   #plot(zt,tL,linewidth=.5)
   xlabel('Z')
   ylabel('Lookback Time [Gyr]')
   grid(True)
   
   subplot(326)
   #plt.semilogx(zt,np.max(tL)-tL,linewidth=.5)
   plt.semilogx(zt,age*(1-tL),linewidth=.5)
   #plot(zt,np.max(tL)-tL,linewidth=.5)
   xlabel('Z')
   ylabel('Age [Gyr]')
   grid(True)

   #show()


###
"""

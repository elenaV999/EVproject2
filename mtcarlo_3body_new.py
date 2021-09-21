########################################################################################
# mtcarlo_3body: script to generate the initial conditions for binary-single scattering
# following Hut & Bahcall 1983
# co-authored by Marco Dall'Amico and Michela Mapelli in 2020
########################################################################################
# -*- coding: utf-8 -*-
import numpy as np
import random as rnd
import math
import matplotlib.pyplot as plt


###################################### PARAMETERS ######################################

Ntot = 200     #generate inital con dor Ntot binary-single encounters
rnd.seed(a=90) #42

#physical constants
G  = 6.667e-8   #G in cgs
pc = 3.086e18   #pc in cgs
msun = 1.989e33 #msun in cgs
rsun = 6.95e10  #rsun in cgs
AU = 1.49e13    #AU in cgs

#scale quantities for physical --> Nbody units conversion (by multplying) 
L_scale = 1./pc                          #1cm/1pc
M_scale = 1./msun                        #1g/1M_sun
T_scale = (G * L_scale**3/M_scale)**0.5  #tscale in nbody units
V_scale = (L_scale/T_scale)              # vel in nbody units

D = 0.01*pc     #initial distance of the intruder from binary COM
                #0.01 pc in cm

#scale quantities for Nbody --> physical units conversion (by multplying)
print("Tscale= ",1./T_scale," s   Vscale= ",1./V_scale, "   clight= ",2.99792e10*V_scale)
fscale="Nbody_units_scale.txt"
fs=open(fscale,"w")
fs.write("Tscale= "+str(1./T_scale)+" s   Vscale= "+str(1./V_scale)+"   clight= "+str(2.99792e10*V_scale)+"\n") #stores Nbody units scale parameters

                                         
######################### EXTREMES of the DISTRIBUTIONS ################################

#----- mass function -----
Mmin = 3.0*msun  #g
Mmax = 70.0*msun #g

a = 2.3          #Salpeter power-law	
#NOTA: may try to change the slope of the mass function beacuse there are not a lot of constraints on the BH mass function. Can also change the extrema Mmin and Mmax, based on the models seen

#---- semimajor axis-----
amin  = 10.      #R_sun
amax  = 2000.    #R_sun

#---- orbital phase------
fmin=0.0
fmax=2.*math.pi

#--- impact parameters ---
bmin=0.0   #AU
bmax=1.0e1 #AU

#---- psi angle------
psimin=0.0
psimax=2.*math.pi

#--- phi angle------
phimin=0.0
phimax=2.*math.pi

#---- cosinus of theta angle ----
costmin=-1.0
costmax=1.0

#---- mean and standard deviation of velocity ----
vmean=0.0
vstd=5.0*1e5  #cm/s


###################################### FUNCTIONS ######################################## 

def salpeter(y,a,mini,maxi): #mass function of black holes (power law with index alpha between mmin and mmax
    p=1.-a
    x=y*(maxi**p-mini**p)+mini**p
    x=x**(1./p)
    return x

def sign(x):    #true if x>0,     
    if x>=0.0:
        return int(1)
    else:
        return int(0)

def uniformedev(x,minimo,maximo):  # x = x(max-min) + min                  
    x*=(maximo-minimo)
    x+=minimo
    return x

def gauss(x1,x2,med,sigma):	 
    p=sigma*np.sqrt(-2.*np.log(1.-x1))*np.cos(2.*math.pi*x2)
    p=p+med
    return p

def anomaly(x,ecc,E):	# g(x) = 0	function to be solved for E
    g = E - ecc * np.sin(E) - x                    
    return g

def bisec(x,ecc,tol):	
    E1=1.0
    E2=100.
    g2=100.
    while(abs(g2)>tol):
        delta=0.5*(E2+E1)
        g1 = anomaly(x,ecc,E1)
        g2 = anomaly(x,ecc,E2)
        if(sign(g1)== sign(g2)): #if g1 and g2 have same sign search other points
            if(sign(g1)==1):
                if((g1)>(g2)):
                    E1 = E2
                    E2 = E2 + delta
                else:
                    E2 = E1
                    E1 = E1 - delta	
            else:
                if((g1)<(g2)):
                    E1 = E2
                    E2 = E2 + delta
                    
                else:
                    E2 = E1
                    E1 = E1 - delta         
        else:
            E3 = (E2 + E1)*0.5
            g3 = anomaly(x,ecc,E3)
            if(sign(g3)==sign(g1)):
                E1=E3
            else:
                E2=E3               
    return ((E1+E2)*0.5)


#########################################################################################
                                     # MAIN #
#########################################################################################

#------------------------------------  BINARY AND INTRUDER MASS ------------------------------------#
# generate masses distributed according to power law with index a between Mmin and Mmax
i=0

m1=[]
m2=[]
m3=[]

while(i<Ntot):
    x=np.random.rand()
    y=salpeter(x,a,Mmin,Mmax)
    m3.append(y)       

    x=np.random.rand()
    y1=salpeter(x,a,Mmin,Mmax)


    x=np.random.rand()
    y2=salpeter(x,a,Mmin,Mmax)


    if(y1>=y2):
        m1.append(y1)
        m2.append(y2)
    else:
        m1.append(y2)
        m2.append(y1) #all masses here are in g
    i+=1



    
#------------------------------------  SEMI-MAJOR AXIS ------------------------------------ #
# generate semi-major axis according to Sana et al. 2012

u     = -0.325  #exponent of the p-law distribution
u1    = u+1.
r  = np.random.random(size=Ntot)
a  = (amin**u1 + (amax**u1 - amin**u1) * r )**(1./u1)  #R_sun
a = a*rsun      #convert to cm 

#print('a =',len(a),type(a))

#-------------------------------------  ECCENTRICITY -------------------------------------- #

#generate binary orbital eccentricity based on thermal distribution P(e)\propto{}e^2

ecc=[]
i=0
while(i<Ntot):
    p=rnd.random()
    p=uniformedev(p,0.,1.)
    p=p**0.5
    ecc.append(p)
    i+=1
ecc=np.array(ecc,float)    
#print('ecc =',len(ecc),type(ecc))
    
#------------------------------------ ORBITAL PHASE ----------------------------------------#

#generate psi angle and binary orbital phase based on Hut & Bahcall 1983 http://articles.adsabs.harvard.edu/pdf/1983ApJ...268..319H
#(E-esinE), range [0,3.14), to derive the spheric anomaly E
#          f= 2 * atan(((1+e)/(1-e))^0.5 * tan(E/2))

tol=1e-7
Efinale=[]
fase=[]
s=0
i=0

while(i<Ntot):
    p=rnd.random()
    p=uniformedev(p,fmin,fmax)
    Efinale.append(bisec(p, ecc[i],tol))                                                      #call bisection and solves eccentric anomaly
    f=2.0 * np.arctan((((1.+ecc[i])/(1.-ecc[i])))**0.5 * np.tan(Efinale[i]/2.))               #;//convert to phasis
    if(Efinale[i]!=Efinale[-1]):
        print(Efinale[i],Efinale[-1])
    fase.append(f)
    i+=1
  

#print('fase =',len(fase),type(fase))



################################## ENCOUNTER PROPERTIES ########################################

#------------------------------------- IMPACT PARAMETER ---------------------------------------#

#generate impact parameter based on Hut & Bahcall 1983 http://articles.adsabs.harvard.edu/pdf/1983ApJ...268..319H
bmax2=bmax**2
bmin2=bmin**2
    
impact=[]

i=0
while(i<Ntot):
    p=rnd.random()
    p=uniformedev(p,bmin2,bmax2)
    p=p**0.5
    impact.append(p)
    i+=1
impact  = np.array(impact)
impact = impact *AU         #convert to cm 

#print('b =',impact.shape,type(impact))

#---------------------------------------  PSI ANGLE --------------------------------------------#

psi=[]
i=0
while(i<Ntot):
    p=rnd.random()
    p=uniformedev(p,psimin,psimax)
    psi.append(p)
    i+=1


#---------------------------------------  PHI ANGLE --------------------------------------------#  

#generate phi angle based on Hut & Bahcall 1983 http://articles.adsabs.harvard.edu/pdf/1983ApJ...268..319H

phi=[]

i=0
while(i<Ntot):
    p=rnd.random()
    p=uniformedev(p,phimin,phimax)
    phi.append(p)
    i+=1


#print('phi =',len(phi),type(phi))


#---------------------------------------  THETA ANGLE -------------------------------------------#  

#generate theta angle based on Hut & Bahcall 1983 http://articles.adsabs.harvard.edu/pdf/1983ApJ...268..319H
# NOTE: there is a bug in Hut & Bahcall - theta should be between 0 and pi, not between 0 and pi/2, otherwise
# the scattering are asymmetric with respect to orbital angular momentum of the binary. MM fixed this on 2020/09/25

theta=[]

i=0
while(i<Ntot):
    p=rnd.random()
    p=uniformedev(p,costmin,costmax)
    p=np.arccos(p)
    theta.append(p)
    i+=1


#print('theta =',len(theta),type(theta))

#-------------------------------------  RELATIVE VELOCITY ----------------------------------------#  

#generate relative velocity module, based on Maxwellian with vstd=5 km/s (typical of young star clusters)

vel=[]
i=0

while(i<Ntot):
    p=[]
    j=0
    while(j<3):
        p1=rnd.random()
        p2=rnd.random()
        p.append(gauss(p1,p2,vmean,vstd))
        j+=1
    vel.append((p[0]*p[0]+p[1]*p[1]+p[2]*p[2])**0.5)
    i+=1

vel = np.array(vel)                                                                                
print(np.mean(vel/1e5),np.median(vel/1e5))

#vel*= V_scale   
#print('vel =',vel.shape,type(vel))

############################################# CHECK-POINT 1 ######################################### 
# print for check m1 (g), m2 (g), m3 (g), a (cm), ecc, b (cm), phi, theta, psi, fase and vel (cm/s)
stream=open("mtcarlo_ICs.dat","w");
for i in range(len(impact)):
    stream.write(str(m1[i])+' '+str(m2[i])+' '+str(m3[i])+' '+str(a[i])+' '+str(ecc[i])+' '+str(impact[i])+' '+str(phi[i])+' '+str(theta[i])+' '+str(psi[i])+' '+str(fase[i])+' '+str(vel[i])+'\n') 


######################################### CHANGE OF COORDINATES #####################################
# goes from initial parameters to Cartesian
#m1 (m2)=primary (secondary) member of binary, m3=intruder
#pedex 1,2,3 always indicate primary member, secondary member, intruder
    
x1=np.zeros(Ntot,float)
y1=np.zeros(Ntot,float)
z1=np.zeros(Ntot,float)

vx1=np.zeros(Ntot,float)
vy1=np.zeros(Ntot,float)
vz1=np.zeros(Ntot,float)

x2=np.zeros(Ntot,float)
y2=np.zeros(Ntot,float)
z2=np.zeros(Ntot,float)

vx2=np.zeros(Ntot,float)
vy2=np.zeros(Ntot,float)
vz2=np.zeros(Ntot,float)

x3=np.zeros(Ntot,float)
y3=np.zeros(Ntot,float)
z3=np.zeros(Ntot,float)

vx3=np.zeros(Ntot,float)
vy3=np.zeros(Ntot,float)
vz3=np.zeros(Ntot,float)

i=0

while(i<Ntot):
    x1[i] = -m2[i]/(m1[i]+m2[i]) * a[i] * (1.-ecc[i]**2) * np.cos(fase[i])/(1.+(ecc[i] * np.cos(fase[i])))

    y1[i] = -m2[i]/(m1[i]+m2[i]) * a[i] * (1.-ecc[i]**2) * np.sin(fase[i])/(1.+(ecc[i] * np.cos(fase[i])))

    z1[i] = 0.0
    
    vx1[i] = -m2[i]/(m1[i]+m2[i]) * (ecc[i] * np.cos(fase[i]) /(1+ ecc[i] * np.cos(fase[i])) - 1.) * np.sin(fase[i]) * (1+ ecc[i] * np.cos(fase[i])) * (G * (m1[i]+m2[i]) / (a[i]*(1- ecc[i]*ecc[i])))**0.5
    
    vy1[i] = -m2[i]/(m1[i]+m2[i]) * (ecc[i] * (np.sin(fase[i]))**2 /(1+ ecc[i] * np.cos(fase[i])) + np.cos(fase[i])) * (1+ ecc[i] * np.cos(fase[i])) * (G * (m1[i]+m2[i]) / (a[i]*(1- ecc[i]*ecc[i])))**0.5                                      

    vz1[i] = 0.0
    
    x2[i] = m1[i]/(m1[i]+m2[i]) * a[i] * (1.-ecc[i]**2) * np.cos(fase[i])/(1.+(ecc[i] * np.cos(fase[i])))

    y2[i] = m1[i]/(m1[i]+m2[i]) * a[i] * (1.-ecc[i]**2) * np.sin(fase[i])/(1.+(ecc[i] * np.cos(fase[i])))

    z2[i] = 0.0
    
    vx2[i] = m1[i]/(m1[i]+m2[i]) * (ecc[i] * np.cos(fase[i]) /(1+ ecc[i] * np.cos(fase[i])) - 1.) * np.sin(fase[i]) * (1+ ecc[i] * np.cos(fase[i])) * (G * (m1[i]+m2[i]) / (a[i]*(1- ecc[i]*ecc[i])))**0.5   

    vy2[i] = m1[i]/(m1[i]+m2[i]) * (ecc[i] * (np.sin(fase[i]))**2 /(1+ ecc[i] * np.cos(fase[i])) + np.cos(fase[i])) * (1+ ecc[i] * np.cos(fase[i])) * (G * (m1[i]+m2[i]) / (a[i]*(1- ecc[i]*ecc[i])))**0.5                                      

    vz2[i] = 0.0

    
    x3[i] = D * ( np.sin(phi[i]) * (impact[i]/D) * np.cos(psi[i]) - np.cos(phi[i]) * ( (1.0 - (impact[i]/D)**2)**0.5 * np.sin(theta[i]) +  (impact[i]/D) * np.cos(theta[i]) * np.sin(psi[i])))

    y3[i] = -D * ( np.sin(phi[i]) * ((1.0 - (impact[i]/D)**2)**0.5 * np.sin(theta[i]) + (impact[i]/D) * np.cos(theta[i]) * np.sin(psi[i])) + np.cos(phi[i]) * np.cos(psi[i]) * (impact[i]/D))

    z3[i] =  D * ( -np.cos(theta[i]) * (1.0- (impact[i]/D)**2 )**0.5 + (impact[i]/D) * np.sin(theta[i]) * np.sin(psi[i]))     

    vx3[i] = vel[i] * np.sin(theta[i]) * np.cos(phi[i])
    vy3[i] = vel[i] * np.sin(theta[i]) * np.sin(phi[i])
    vz3[i] = vel[i] * np.cos(theta[i]) 
    i+=1

#################################CONVERT FROM CGC TO NBODY UNITS########################
i=0
while(i<Ntot):
    m1[i]*=M_scale
    x1[i]*=L_scale
    y1[i]*=L_scale
    z1[i]*=L_scale
    vx1[i]*=V_scale
    vy1[i]*=V_scale
    vz1[i]*=V_scale
    
    m2[i]*=M_scale
    x2[i]*=L_scale
    y2[i]*=L_scale
    z2[i]*=L_scale
    vx2[i]*=V_scale
    vy2[i]*=V_scale
    vz2[i]*=V_scale
    
    m3[i]*=M_scale
    x3[i]*=L_scale
    y3[i]*=L_scale
    z3[i]*=L_scale
    vx3[i]*=V_scale
    vy3[i]*=V_scale
    vz3[i]*=V_scale
    i+=1
    ############################################# CHECK-POINT 2 ################################################## 
# print file for N-body code
#first row: primary mass, x, y,z, vx,vy,vz in N-body units
#second row: secondary mass, x, y,z, vx,vy,vz in N-body units
#third row: intruder mass, x, y,z, vx,vy,vz in N-body units
#then repeat with the next three lines
for i in range(len(impact)):
    kk=str(i)
    stream2=open("input/threebody_ICs_"+kk+".dat","w");
    stream2.write(str(m1[i])+" "+str(x1[i])+' '+str(y1[i])+' '+str(z1[i])+' '+str(vx1[i])+' '+str(vy1[i])+' '+str(vz1[i])+'\n'+str(m2[i])+' '+str(x2[i])+' '+str(y2[i])+' '+str(z2[i])+' '+str(vx2[i])+' '+str(vy2[i])+' '+str(vz2[i])+'\n'+str(m3[i])+' '+str(x3[i])+' '+str(y3[i])+' '+str(z3[i])+' '+str(vx3[i])+' '+str(vy3[i])+' '+str(vz3[i])+'\n')

########################################## PLOT PARAMETERS ####################################################   

n_bins=50

fig,ax=plt.subplots(4,2)
ax[0][0].hist(impact**2, bins=n_bins)
ax[0][0].set_xlabel('$b^2\,{}(\mathrm{cm}^2)$')

ax[1][0].hist(phi, bins=n_bins)
ax[1][0].set_xlabel('$\phi{}$ (rad)')


ax[2][0].hist(theta, bins=n_bins)
ax[2][0].set_xlabel('$\\theta{}$ (rad)')

ax[3][0].hist(ecc**2,bins=n_bins)
ax[3][0].set_xlabel('$e^2$')

ax[0][1].hist(psi, bins=n_bins)
ax[0][1].set_xlabel('$\psi{}$ (rad)')


ax[1][1].hist(fase, bins=n_bins)
ax[1][1].set_xlabel('phase (rad)')

ax[2][1].hist(vel, bins=n_bins)
ax[2][1].set_xlabel('$v\,{}(\mathrm{cm}\,{}\mathrm{s}^{-1})$ ')

ax[3][1].hist(Efinale,bins=n_bins)
ax[3][1].set_xlabel('Eccentric anomaly')

plt.tight_layout()





fig2,axs=plt.subplots(2,2)
axs[0][0].plot(x1,y1,'r.')
axs[0][0].plot(x2,y2,'b.')
axs[0][0].plot(x3,y3,'g.')
axs[0][0].set_xlabel('$x$ ')
axs[0][0].set_ylabel('$y$ ')


axs[1][0].plot(x1,z1,'r.')
axs[1][0].plot(x2,z2,'b.')
axs[1][0].plot(x3,z3,'g.')
axs[1][0].set_xlabel('$x$ ')
axs[1][0].set_ylabel('$z$ ')


axs[0][1].plot(vx1,vy1,'r.')
axs[0][1].plot(vx2,vy2,'b.')
axs[0][1].plot(vx3,vy3,'g.')
axs[0][1].set_xlabel('$v_x$ ')
axs[0][1].set_ylabel('$v_y$ ')


axs[1][1].plot(vx1,vz1,'r.')
axs[1][1].plot(vx2,vz2,'b.')
axs[1][1].plot(vx3,vz3,'g.')
axs[1][1].set_xlabel('$v_x$ ')
axs[1][1].set_ylabel('$v_z$ ')


plt.tight_layout()


fig3,ax = plt.subplots(2,2)



ax[0][0].hist(a, bins=n_bins)
ax[0][0].set_xlabel('$a\,{}(\mathrm{cm})$')

ax[1][0].hist(m1, bins=n_bins)
ax[1][0].set_xlabel('$M_1\,{}(\mathrm{M}_\odot)$ ')

ax[1][1].hist(m2, bins=n_bins)
ax[1][1].set_xlabel('$M_2\,{}(\mathrm{M}_\odot)$ ')

ax[0][1].hist(m3, bins=n_bins)
ax[0][1].set_xlabel('$M_3\,{}(\mathrm{M}_\odot)$ ')

plt.tight_layout()


plt.show()

import numpy as np
import matplotlib.pyplot as plt
import os
import re


Tscale=470772692629564.4
pcinAU=2.06264806247096e5 #A.U.


n=np.arange(0,303)	#list of all directories
dir=[]
for n in n: 
	dir_='RUN1_200/n'+str(n)+'/'
	dir.append(dir_) 
	
Ntot=len(dir)  #n° of 3-body encounters	

exchange, original , longexchange, mlarger, m3smaller, merge_ini, merge_fin, count=0 , 0,0,0,0, 0, 0, 0  #counters

NSAME=[] #n° of timesteps that the exchange lasts
mdistr , m_single_ini, m_single_fin, m_binary_ini , m_binary_fin = [],[],[],[],[]
q_ini, qex, qor =[],[], []  #min(m1,m2)/m3	or: m_f/m3   (intruder)
sep_ini, sep_fin, ecc_ini, ecc_fin, tgw_in, tgw_fin =[],[],[],[], [],[]



for i in range(len(dir)) : 	#loop over all directories

	#read initial masses
	m=np.genfromtxt(dir[i]+"INPUT.TXT",usecols=(0),unpack=True)   
	m1=m[0] 	#primary
	m2=m[1]   #companion
	m3=m[2]   #intruder
	
	#read binary file
	fname='88_find_binaries.dat'
	f=open(dir[i]+fname , "r") 
	
	#variable 	   1=time   2=ID1  3=ID2   4=a	5=e   
	#a=re.compile("^\s+(\S+)\D+(\d)\D+(\d)\D+(\S+)\D+(\S+)")  #expression to search in file
	a=re.compile("^\s+(\S+)\s+\D+\s+(\d)\D\s+(\d)\D\s+\S+\s+(\S+)\s+\S+\s+(\S+)")
	t,id1,id2,sep,ecc=[],[],[],[],[]
	for line in f: 
		b=a.search(line)
		if(b!=None):
			t.append(eval(b.group(1)))
			id1.append(eval(b.group(2)))
			id2.append(eval(b.group(3)))
			sep.append(eval(b.group(4))) #semimajor axis
			ecc.append(eval(b.group(5))) #eccentricity
	f.close()
	
	t=np.array(t)
	id1=np.array(id1)
	id2=np.array(id2)
	sep=np.array(sep)
	ecc=np.array(ecc)


	#MASSES at t_final
	if    (id1[-1]+id2[-1])==3 :  m3_idx=2	#single star mass at t_final 
	elif  (id1[-1]+id2[-1])==4 : 	m3_idx=1		
	elif  (id1[-1]+id2[-1])==5 :  m3_idx=0		
	else : print('ERROR: FINAL INTRUDER MASS NOT FOUND')
	
	m1_f=m[id1[-1]-1]
	m2_f=m[id2[-1]-1]
	m3_f=m[m3_idx]
	if m1_f==m2_f or m1_f==m3_f or m2_f==m3_f : print('ERROR in final masses') #check that the masses are all different
	
	
	#MASSES:
	mdistr.extend([m1,m2,m3])  #initial 
	m_single_ini.append(m3) 
	m_binary_ini.extend([m1,m2]) 
	
	m_single_fin.append(m3_f)  #final
	m_binary_fin.extend( [m1_f, m2_f]) 
	
	if min(m1, m2)/m3  > 1: m3smaller+=1  #binary less massive satr > intruder mass
	
	#semimajor axis and eccentricity	
	sep_ini.append(sep[0]) 
	ecc_ini.append(ecc[0])
	sep_fin.append(sep[-1]) 
	ecc_fin.append(ecc[-1])
	#if (sep[-1]*pcinAU) > 20 : print('a final large: ', dir[i], sep[0]*pcinAU)


	#GW timescale
	clight= 4573359.917978041
	G=1.
	const=5/256*clight**5/G**3   #in NB units
	tgwin=const*sep[0]**4*(1-ecc[0]**5)**3.5/(  m1*m2*(m1+m2))*Tscale/(60*60*24*365)/1e9
	tgwfin=const*sep[-1]**4*(1-ecc[-1]**5)**3.5/(m1_f*m2_f*(m1_f+m2_f) )*Tscale/(60*60*24*365)/1e9
	
	if tgwin<14: merge_ini+=1
	if tgwfin<14: merge_fin+=1	
	
	tgw_in.append(tgwin)
	tgw_fin.append(tgwfin)
	
	
#####count n° of binaries for each case
	idsum=id1+id2
	Nsame = np.count_nonzero(idsum==idsum[0])	#count how many elemets in idsum are = to the first one
	NNsame=len(t)-Nsame
	texc = NNsame*1e-6*Tscale/(60*60*24*365)
	NSAME.append(texc)
	
	
	#ends with the ORIGINAL binary	
	if idsum[0]==idsum[-1] : 
		qor.append(min(m1, m2)/m3 )	#lower binary mass / intruder (ORIGINAL BINARIES)
		
		if Nsame==len(t): #no exchange at all
			original+=1	
 
		elif  Nsame>970 : #exchange for a very short time
			original+=1 	 

		else : 
			longexchange+=1   #exchange for a long time
			original+=1
			print('LONG EXCHANGE, ends as original: ', dir[i], 't_exc = ', texc)  #77(Nsame=470), 186(Nsame=937)

	#ends with an EXCHANGED binary		
	else : 
		exchange+=1	
		qex.append(m3_f/m3)  #exchanged binary star /  intruder (only for exchanges)
		if m3_f/m3 > 1  : mlarger+=1 	#exchange happens with a more massive star
		if m3_f==m1 : 
			count+=1    #exchange with the most massive star	

sep_ini=np.array(sep_ini)
sep_fin=np.array(sep_fin)

#################  print some info about the encounters:###########################################################
print('\n\nn° of encounters = ', Ntot)
print('exchanged binaries: N = ', exchange , '	', exchange/Ntot*100,'%')
print('original binaries: N = ', original , '	', original/Ntot*100,'%')
print('original, with long time exchanged: N = ', longexchange , '	', longexchange/Ntot*100,'%')

print('\nexchange with more massive star : N = ', mlarger, '	% of exchanges: ', mlarger/exchange*100 , ' total %', mlarger/Ntot*100)
print('nexchange with the most massive star in the binary:', count, '	', count/exchange*100,'%')		

print('\nm_intruder < both m_binary: ',m3smaller,'  ',m3smaller/Ntot*100,' %','\nexchanges: ',mlarger,'	',mlarger/m3smaller*100,' %') 
#calculated on the cases where the intruder has smaller mass than both stars in the binary
print('m_intruder > at least one m_binary: ', Ntot-m3smaller,'  ', (Ntot-m3smaller)/Ntot*100,'\nexchanges', exchange-mlarger,'	',(exchange-mlarger)/(Ntot-m3smaller)*100 ) 
#over the cases in which at least one star from the binary has smaller mass than the intruder (more cases)

print('\nmergers before a Hubble time: ',merge_ini, merge_fin)

##########################################################################################################
##########################################################################################################

########### MASS RATIO: intruder/binary 
lw=2
nbins=30
dens=False

plt.figure(figsize=(8,6), dpi=100)
plt.hist(qex, bins=nbins,  color='tab:red',  histtype='step', linewidth=lw, density=dens,range=(0,4), zorder=1)  # mf3/m3        - exchanged binaries
plt.hist(qor, bins=nbins,  color='tab:blue',  histtype='step', linewidth=lw, density=dens, range=(0,4), zorder=0) # min(m1,m2)/m3 - original binaries
#plt.hist(q_ini, bins=nbins,  color='black',  histtype='step', linewidth=lw, density=dens) #q=min(m_binary) / m_intruder - all systems

plt.yscale('log')
plt.xlabel('q',  fontsize=15)
plt.ylabel('Number' , fontsize=14)
plt.legend(['exchange', 'no exchange'],  fontsize=14)
plt.axvline(x=1., color='black', linestyle='--')
plt.tight_layout()
plt.show()


########### BH BINARY TOTAL MASS, INITIAL vs FINAL
mbins=30
#mbins=np.logspace(np.log10(3),np.log10(100),nbins)
plt.figure(figsize=(8,6), dpi=100)
plt.hist(m_single_fin, bins=mbins,  color='tab:blue',  histtype='step', linewidth=lw, density=dens, range=(0,70)) 
plt.hist(m_binary_fin, bins=mbins,  color='tab:red',  histtype='step', linewidth=lw, density=dens, range=(0,70))   
#plt.hist(m_binary_ini, bins=mbins,  color='tab:red',  histtype='step', linewidth=lw, density=dens, ls='--', alpha=0.6, range=(0,70))    
#plt.hist(m_single_ini, bins=mbins,  color='tab:blue',  histtype='step', linewidth=lw, density=dens, ls='--', alpha=0.6, range=(0,70)) 

#plt.xscale('log')
plt.yscale('log')
plt.xlabel('M [M$_\odot$]', fontsize=15)
plt.ylabel('Number', fontsize=14)
plt.legend(['single BHs', 'BHs in binaries'],  fontsize=14)
#plt.legend( ['binary initial','binary final','single initial','single final'])
plt.tight_layout()
plt.show()


##########################################################################################################
##########################################################################################################

lw=2.
lw2=1.5
nbins=30
dens=False

fig2,ax=plt.subplots(1,3, figsize=(15,4.6))

########### SEMIMAJOR AXIS 
ax[0].hist(sep_ini*pcinAU, bins=nbins,  color='black', histtype='step', range=(0,40) ,linewidth=lw2, density=dens,  ls='--')
ax[0].hist(sep_fin*pcinAU, bins=nbins,  color='orangered',  histtype='step', range=(0,40), linewidth=lw, density=dens,  alpha=0.7)
#ax[1][0].set_yscale('log')
ax[0].set_xlabel('a [A.U.]', fontsize=15)
ax[0].set_ylabel('Number', fontsize=15)
ax[0].set_title('semi-major axis', fontsize=15, fontweight='heavy')


########### ECCENTRICITY 
nbins=15
ax[1].hist(ecc_ini, bins=nbins,color='black', histtype='step', linewidth=lw2,  density=dens, ls='--')
ax[1].hist(ecc_fin , bins=nbins, color='orangered',  histtype='step', linewidth=lw,  density=dens, alpha=0.7)
#ax[1].set_yscale('log')
ax[1].set_xlabel('e', fontsize=15)
ax[1].set_title('eccentricity', fontsize=15, fontweight='heavy')


########## GW EMISSION MERGER TIMESCALE
logbins=np.logspace(-4,11,nbins)
ax[2].hist(tgw_in, bins=logbins,  color='black', histtype='step', linewidth=lw2, density=dens,   ls='--')
ax[2].hist(tgw_fin, bins=logbins,  color='orangered',  histtype='step', linewidth=lw, density=dens, alpha=0.7)
#ax[2].set_yscale('log')
ax[2].set_xscale('log')
ax[2].set_xlabel('t$_{GW}$ [Gyr]', fontsize=15)
ax[2].set_title('GW merger timescale', fontsize=15, fontweight='heavy')
ax[2].legend( ['initial','final'],loc='upper left')

plt.tight_layout()
plt.show()

##########################################################################################################
##########################################################################################################

### time for which the binaries stay exchanged
'''
plt.hist(NSAME, bins=100)
plt.xlabel('Nsame')
plt.title('duration of the exchange')
plt.show()
'''


#!/usr/bin/env python
########################################################################
# Explore the Rm-magnetic conductivity implications.
# 7/25/17 - theory4.py
############

import matplotlib.pyplot as plt
from matplotlib.pyplot import plot,ion,show
from matplotlib import ticker
import numpy as np
from numpy import exp
from math import acos
import scipy.constants as c
import astropy.constants as ac
import argparse
ion()  #interactive mode?
from sys import exit
pi=np.pi
kb=1.38064852e-23  #boltzmann
e=1.6021765e-19  #electron charge
lorenz=pi**2/3.*(kb/e)**2

parser = argparse.ArgumentParser()
parser.add_argument('-s',dest='save',action='store_true')  
parser.add_argument('-m',dest='mode',type=int,
                    default=2,help='mantle temperature mode: 1=T_{LM}=constant, and DT_{LM}=eta_LM*T_cmb; 2=variable T_man, Q_{sec,man}=Q_{sec,core}; 3=T_{LM}=constant')  
parser.add_argument('-lorenz',dest='lorenz',type=float,
                    default=lorenz,help='Lorenz number to use in Wiedemann-Franz Law')
parser.add_argument('-sn',dest='savename',type=str,
                    default='',help='suffix of Figure savename') 
parser.add_argument('-sf',dest='saveformat',type=str,
                    default='eps',help='Figure file format') 
parser.add_argument('-Q_cmb_e',dest='Q_cmb_e',type=float,
                    default=12e12,help='CMB heat flow today in W.') 
args = parser.parse_args()
lorenz = args.lorenz
if args.save:
    fontsize=20
    plt.rcParams['xtick.labelsize']=fontsize
    plt.rcParams['ytick.labelsize']=fontsize
    plt.rcParams['font.size']=fontsize+5
    figsuffix='_m%d'%args.mode
print('mode=%d'%args.mode)
sb=c.Stefan_Boltzmann
R_g=c.R
mu_0=c.mu_0
mass_e=5.972e24  #[kg] mass of E.
mass_core=1.95e24 #[kg] mass of core.
mass_man=mass_e-mass_core
r_e=ac.R_earth.value
area_e=4.*np.pi*r_e**2
r_cmb=3481e3
r_core=r_cmb
r_icb=1221e3
area_cmb=4*pi*pow(r_cmb,2.)
area_core=area_cmb
vol_core=4./3*pi*pow(r_cmb,3.)
vol_man=4./3*pi*(pow(r_e,3.)-pow(r_cmb,3.))
rho_core=mass_core/vol_core
rho_man=mass_man/vol_man
rho_cmb=9903.  #PREM
rho_lm=5500.  #PREM above D''
omega=2*np.pi/(3600.*24)  #[rad/day] E Rotation rate.
cp_man=1260.  #mantle cp.
#cp_core=650.
cp_core=715. # (Davies 15 PEPI)
thermexp=1.35e-5  #thermal expansivity  Davies 2015 PEPI.
thermexp_core=thermexp
thermexp_man=1e-5
k_lm=10.       #LM thermal conductivity
kappa_lm=k_lm/(rho_lm*cp_man)
#kappa_lm=1e-5
#g=10.          #gravity.
g=10.68          #gravity.
g_i=5.
delta_rho_ic=590. #[kg/m3] D16 supp.
latent_i=750e3   #[J/K]
gamma_i=-1.3e-9  #[K/Pa] D16 supp.
#H=cp/(thermexp*g)  #scale height H=cp/(alpha*g)
#lorentz=2.44e-8
Q_m_e=39e12  #mantle Q E
T_m_e=2300.  #mid-mantle E
T_um_e=1630. #upper mantle adiabat.
T_lm_e=2500. # lower mantle E.
T_cmb_e=4e3
delta_T_e=T_cmb_e-T_lm_e
T_icn=T_cmb_e+100.  #T_cmb of ICN.
sigma_core_e = 1e6 #nominal
#k_core_e=100.  #nominal
k_core_e=sigma_core_e*lorenz*T_cmb_e
#Q_ad_e=15e12  #nominal
eta_lm=1.-T_lm_e/T_cmb_e
ra_crit=660.
act_visc=3e5
#Q_cmb_e=15e12
Q_cmb_e=args.Q_cmb_e   #12e12
#visc_man_e=9.855e18  #nominal mantle viscosity
T_tbl_e=(T_cmb_e+T_lm_e)/2.  #viscosity is func of T in TBL.
if args.mode==1: visc_man_e=(thermexp_man*g*rho_lm*cp_man/ra_crit)*(area_cmb*k_lm**(2./3)*(eta_lm*T_cmb_e)**(4./3)/Q_cmb_e)**(3.)
if args.mode==2: visc_man_e=(thermexp_man*g*rho_lm*cp_man/(ra_crit))*(area_cmb*k_lm**(2./3)*(delta_T_e)**(4./3)/Q_cmb_e)**3
if args.mode==3: visc_man_e=(thermexp_man*g*rho_lm*cp_man/(ra_crit))*(area_cmb*k_lm**(2./3)*(delta_T_e)**(4./3)/Q_cmb_e)**3

visc_ref=visc_man_e/exp(act_visc/(R_g*T_tbl_e))
gamma_ad_core=g*thermexp/cp_core  #adidbatic gradient core:dT/dz (Poirier 7.22)
Q_ad_e = area_cmb*k_core_e*gamma_ad_core*T_cmb_e
r_ic_e=1221e3
T_ic_e=5500.
heatcap_core_man = (mass_core*cp_core)/(mass_man*cp_man)
print('T_cmb_e=%d Q_cmb_e=%.2e Q_ad_e=%.2e k_core_e=%.2f visc_man_e=%.2e visc_ref=%.2e gamma_ad_core=%.2e'%(T_cmb_e,Q_cmb_e,Q_ad_e,k_core_e,visc_man_e,visc_ref,gamma_ad_core))

### Rm coefficient.
rm_crit=40.
d=r_cmb-r_icb
c_const=0.5    #Form: Nu=c*Ra^beta  c=?
beta=1./3      #same.
ro_0=0.85      #Rossby coef, OC06 (18)
alpha=2./5     #Rossby exponent, OC06 (18)

### Variables
T_cmb_lo=3.0e3; T_cmb_hi=5e3+1; delta_T_cmb=10
T_cmb_l0=3.5e3; T_cmb_hi=4.5e3+1
T_cmb=np.arange(T_cmb_lo,T_cmb_hi,delta_T_cmb)
nt=len(T_cmb)
i_e_T_cmb = np.where(T_cmb-T_cmb_e>=0)[0].min()

# function to compute T_lm given T_cmb, mode.
def f_T_lm(mode):
    # Set T_lm
    if mode==1: T_lm=T_lm_e+np.zeros(nt)
    if mode==2:
        T_lm_lo=T_lm_e+heatcap_core_man*(T_cmb_lo-T_cmb_e)
        T_lm_hi=T_lm_e+heatcap_core_man*(T_cmb_hi-T_cmb_e)
        delta_T_lm=heatcap_core_man*delta_T_cmb
        T_lm=np.arange(T_lm_lo,T_lm_hi,delta_T_lm)
        if len(T_lm)!=nt:
            print('Error: len(T_lm)=%d not equal to len(T_cmb)=%d '%(len(T_lm),nt))
            exit(0)
    if mode==3: T_lm=T_lm_e+np.zeros(nt)
    return(T_lm)
    
# function to compute delta_tbl given T_cmb,mode.
def f_delta_tbl(mode,T_cmb,visc_tbl,delta_T):
    if mode==1: delta_tbl=pow(ra_crit*kappa_lm*visc_tbl/(thermexp_man*g*eta_lm*T_cmb),1./3)
    if mode==2: delta_tbl=pow(ra_crit*kappa_lm*visc_tbl/(thermexp_man*g*delta_T),1./3)
    if mode==3: delta_tbl=pow(ra_crit*kappa_lm*visc_tbl/(thermexp_man*g*delta_T),1./3)
    return(delta_tbl)
    
# function to compute Q_cmb given T_cmb,mode,delta_T.
def f_Q_cmb(mode,T_cmb,delta_T,visc_tbl,delta_tbl):
    if mode==1: Q_cmb=area_cmb*k_lm**(2./3)*(eta_lm*T_cmb)**(4./3)*(thermexp_man*g*rho_man*cp_man/(ra_crit*visc_tbl))**(1./3)
    if mode==2: Q_cmb=area_cmb*k_lm*(delta_T)/delta_tbl  
    if mode==3: Q_cmb=area_cmb*k_lm*(delta_T)/delta_tbl 
    return(Q_cmb)

# function to compute sigma_crit
def f_sigma_crit(mode,T_cmb,delta_T,visc_tbl):
    if mode==1: sigma_2=k_lm*pow(eta_lm,4./3)/(lorenz*gamma_ad_core)*pow(thermexp_man*g/(ra_crit*kappa_lm),1./3)*pow(visc_tbl*T_cmb**2.,-1./3)
    if mode==2 or mode==3:
        #sigma_2=k_lm/(lorenz*gamma_ad_core)*pow(thermexp_man*g/(ra_crit*kappa_lm*visc_tbl),1./3)*pow(delta_T,4./3)/T_cmb**2.
        sigma_2=k_lm**(2./3)/(lorenz*gamma_ad_core)*pow(thermexp_man*g*rho_lm*cp_man/(ra_crit),1./3)*pow(delta_T,4./3)/(T_cmb**2.*visc_tbl**(1./3))
    return(sigma_2)

def f_visc_tbl(T_tbl):
    return(visc_ref*exp(act_visc/(R_g*T_tbl)))  # use ave TBL temp.
    
# T_lm
T_lm = f_T_lm(args.mode)
# Thermal Boundary Layer
i_e=(np.abs(T_cmb-T_cmb_e)).argmin()  #get index of where T_cmb=T_cmb_e
delta_T=T_cmb-T_lm
T_tbl=(T_cmb+T_lm)/2. # average of LM and CMB.
visc_tbl = f_visc_tbl(T_tbl)
delta_tbl = f_delta_tbl(args.mode,T_cmb,visc_tbl,delta_T)
Q_cmb = f_Q_cmb(args.mode,T_cmb,delta_T,visc_tbl,delta_tbl)
# Critical conductivity: f(T_cmb)
sigma_2 = f_sigma_crit(args.mode,T_cmb,delta_T,visc_tbl)
print('sigma_2[i_e]=%e delta_LM[i_e]=%.4e'%(sigma_2[i_e],delta_tbl[i_e]))

# Compare modes.
T_lm_m2 = f_T_lm(2)
T_lm_m3 = f_T_lm(3)
delta_T_m2 = T_cmb-T_lm_m2
delta_T_m3 = T_cmb-T_lm_m3
T_tbl_m2 = (T_cmb+T_lm_m2)/2.
T_tbl_m3 = (T_cmb+T_lm_m3)/2.
visc_tbl_m2 = f_visc_tbl(T_tbl_m2)
visc_tbl_m3 = f_visc_tbl(T_tbl_m3)
delta_tbl_m2 = f_delta_tbl(2,T_cmb,visc_tbl_m2,delta_T_m2)
delta_tbl_m3 = f_delta_tbl(3,T_cmb,visc_tbl_m3,delta_T_m3)
Q_cmb_m2 = f_Q_cmb(2,T_cmb,delta_T_m2,visc_tbl_m2,delta_tbl_m2)
Q_cmb_m3 = f_Q_cmb(3,T_cmb,delta_T_m3,visc_tbl_m3,delta_tbl_m3)
sigma_2_m2 = f_sigma_crit(2,T_cmb,delta_T_m2,visc_tbl_m2)
sigma_2_m3 = f_sigma_crit(3,T_cmb,delta_T_m3,visc_tbl_m3)

# Sigma Array
sigma_min=1e4; sigma_max=4e6; sigma_step=1e3
sigma_arr=np.arange(sigma_min,sigma_max,sigma_step)  #(1e4,4e6,1e4)
nsig=len(sigma_arr)
i_e_sigma=np.min(np.where(sigma_arr>=sigma_core_e))
#print('sigma_arr=',sigma_arr*1e-6)

# Rm
d=r_cmb
r_star = r_cmb/r_ic_e
rm_0=mu_0*ro_0*pow(d,2*(1.-2.*alpha))*pow(omega,1-3*alpha)*pow(r_cmb,2.*alpha)*pow( (thermexp_core*g*lorenz*gamma_ad_core)/(rho_cmb*cp_core),alpha) 
print('Rm_0=%.2e'%rm_0)

rm=np.zeros((nt,nsig))
for i in range(nsig):
    ind_real=np.argwhere(sigma_2 >=sigma_arr[i])   #only compute fractional powers of positive nums
    rm[ind_real,i]=rm_0*pow(sigma_arr[i],1.+alpha)*pow(T_cmb[ind_real],2.*alpha)*pow(sigma_2[ind_real]/sigma_arr[i]-1.,alpha)

# Plotting specifics
colors=['k','b','g','r','c','m','y']
ncol=len(colors)
nsig_short=ncol-1
i_short=np.linspace(0,nsig-1,nsig_short)  #create short index array.
sig2plot=np.array([0.6e6,0.8e6,1.0e6,1.2e6,1.4e6])
n_sig2plot=len(sig2plot)
i_2plot=np.zeros(n_sig2plot)
for i in range(n_sig2plot): i_2plot[i]=np.argwhere(sig2plot[i]==sigma_arr)
i_2plot=i_2plot.astype(int)

### Special Temps to plot
T2plot=np.array([3.8e3,3.9e3,4e3,4.1e3,4.2e3])
n_T2plot=len(T2plot)
i_T2plot=np.zeros(n_T2plot)
for i in range(n_T2plot): i_T2plot[i]=np.argwhere(T2plot[i]==T_cmb)
i_T2plot=i_T2plot.astype(int)

##############################################    
### Compute individual terms
d_e=r_cmb-r_ic_e
k_core=np.zeros([nt,nsig])
Q_ad=np.zeros([nt,nsig])
q_conv=np.zeros([nt,nsig])
bflux_thermal=np.zeros([nt,nsig])
bflux_comp=np.zeros([nt,nsig])
bflux_total=np.zeros([nt,nsig])
Ra_Q=np.zeros([nt,nsig])+1e-10
u=np.zeros([nt,nsig])+1e-10
rm2=np.zeros([nt,nsig])+1e-10
sigma_crit=sigma_2   #where Rm=0
sigma_dynamo1=(rm_crit/rm_0)/pow(T_cmb**2*sigma_crit,alpha)
sigma_dynamo2=sigma_crit/(1.+pow(sigma_dynamo1/sigma_crit,1./alpha))   #where Rm=Rm_crit
sigma_peak=(1.-1./(1.+1./alpha))*sigma_crit
rm_peak=np.zeros(nt)

# No IC for T_cmb<T_icn
icn=np.zeros([nt])  # boolean switch for IC: y=1,n=0
for j in range(nt):
    if T_cmb[j] <= T_icn: icn[j]=1.  #yes IC for T_cmb<T_icn
i_icn = np.min(np.where(T_cmb>=T_icn)[0])  #T_cmb[i_icn]=T_icn
# Loop through sigma_arr
for i in range(nsig):
    # Approximation with g_i fixed.
    k_core[:,i]=sigma_arr[i]*lorenz*T_cmb[:]
    Q_ad[:,i]=area_cmb*k_core[:,i]*gamma_ad_core*T_cmb
    q_conv[:,i]=(Q_cmb[:]-Q_ad[:,i])/area_cmb
    bflux_thermal[:,i]=thermexp*g/(rho_core*cp_core)*q_conv[:,i]
    bflux_comp[:,i]=icn[:]*g_i*(delta_rho_ic/rho_core+thermexp*latent_i/cp_core)*Q_cmb[:]/(-gamma_i*rho_core*g*mass_core*cp_core)
#    bflux_total[:,i]=bflux_thermal[:,i]+bflux_comp[:,i]
    bflux_total[:,i]=bflux_thermal[:,i]
    Ra_Q[:,i]=bflux_total[:,i]*r_cmb**2./(pow(d_e,4.)*pow(omega,3.)) #theory7
    ind_pos=np.argwhere(Ra_Q[:,i] >= 0.)  #keep only positive nums for fractional exponent.
    u[ind_pos,i]=ro_0*omega*d_e*pow(Ra_Q[ind_pos,i],alpha)
    rm2[:,i]=mu_0*sigma_arr[i]*u[:,i]*d_e
for j in range(nt):  #loop T_cmb to get each rm_peak value
    rm_peak[j]=rm2[j,np.argmin(np.abs(sigma_arr-sigma_peak[j]))]  #rm[j,index where sigma=peak]
    
# Compute sigma(T) from Ohta 2016
resist_1520=17e-8
resist_2450=31e-8
resist_3750=40.4e-8
coef_b=(resist_3750-resist_1520)/(3750.-1520.)  #slope
coef_a=resist_1520
T_ref=1520.
resist_fe=coef_a+coef_b*(T_cmb-T_ref)
f_melt=0.2    #is this 20% resistivity increase upon melting?
resist_le=1./((1+f_melt)/86.9e-8-1./40.4e-8)  #effect of light elements on resistivity
resist=(1+f_melt)/(1./resist_fe+1./resist_le)  #resistivitiy of Fe and Le effect.
sigma_T=1./resist  #T-dependent electrical conductivity
# Compute Rm(sigma(T),T), note same T array is used in sigma_T so only 1D array needed.
rm3=rm_0*pow(sigma_T[:],1.+alpha)*pow(T_cmb[:],2.*alpha)*pow(sigma_2[:]/sigma_T[:]-1.,alpha)
rm3=np.nan_to_num(rm3)  #replace nan with zero.

### Print E-values.
print('Rm(i_e)=%f sigma_crit(i_e)=%.2e'%(rm2[i_e,i_e_sigma],sigma_crit[i_e]))
### Print ICN values.
i_icn=i_e+10
i_icn_sigma=np.min(np.where(sigma_T[i_icn]<=sigma_arr))
print('at T=T_ICN=%d: sigma_T=%.4e Rm=%.2f'%(T_cmb[i_icn],sigma_arr[i_icn_sigma],rm2[i_icn,i_icn_sigma]))

##############################################
##############################################
##############################################

### Plot Q_cmb,Rm,log(Rm),Buoyancy Flux vs T_cmb
print('Figure 2')
nfig=2
plt.figure(nfig,figsize=(15,10))
# T_cmb vs Q_cmb
xlabel=r'$T_{cmb}$ [K]'
plt.subplot(221)
plot(T_cmb,Q_cmb*1e-12,colors[0],label=r'$Q_{cmb}$',linewidth=2)
for i in range(n_sig2plot):
    plot(T_cmb,Q_ad[:,i_2plot[i]]*1e-12,colors[i+1],label=r'$Q_{a,c}$, $\sigma=$'+str(sigma_arr[i_2plot[i]]*1e-6)+r'$\times 10^6$',linewidth=2)
plt.xlabel(xlabel)
plt.xlim([3.5e3,4.5e3])
plt.ylim([0,30])
plt.ylabel(r'$Q_{cmb}$ [TW]')
plt.legend(loc='best',fontsize=15,ncol=2,frameon=False)
# Tcmb vs Rm
plt.subplot(223) #bottom left.
for i in range(n_sig2plot):
    plot(T_cmb,rm2[:,i_2plot[i]]+1e-20,label='sigma='+str(sigma_arr[i_2plot[i]]*1e-6)+'x10^6',linewidth=2.)
plt.xlabel(xlabel)
plt.xlim([3.5e3,4.5e3]); plt.ylim([0,8e3])
plt.ylabel('Rm')
# Tcmb vs Rm (ylog)

plt.subplot(222) #top right
for i in range(n_sig2plot):
    plot(T_cmb,rm2[:,i_2plot[i]]+1e-20,label='sigma='+str(sigma_arr[i_2plot[i]]*1e-6)+'x10^6',linewidth=2.)
plt.xlabel(xlabel)
plt.ylabel('Rm')
plt.yscale('log')
plt.ylim(1.,1e4)
# Buoyancy flux
plt.subplot(224)
plot(T_cmb,bflux_thermal[:,0])
for i in range(n_sig2plot):
    plot(T_cmb,bflux_thermal[:,i_2plot[i]],label='sigma='+str(sigma_arr[i_2plot[i]]*1e-6)+'x10^6',linewidth=2.)
plot(T_cmb,bflux_comp[:,0],label='bflux_comp',linewidth=2.)
plt.xlabel(xlabel)
plt.ylabel('Buoynacy Flux')
plt.yscale('log')
plt.tight_layout()

if args.save:
    figpref='Figure%d'%nfig
    plt.savefig(figpref+figsuffix+args.savename+'.'+args.saveformat,format=args.saveformat)



# Plot sigma vs Tcmb alone
print('Figure 4')
nfig=4
plt.figure(nfig,figsize=(15,10))
# T_cmb vs sigma_2
yscale=1e-6
alpha_color=0.6
# Fill No Dynamo
plt.fill_between(T_cmb[i_icn:],sigma_dynamo2[i_icn:]*yscale,10.,color='grey',alpha=alpha_color)
# Fill Thermal Dynamo
plt.fill_between(T_cmb[i_icn:],0,sigma_dynamo2[i_icn:]*yscale,color='red',alpha=alpha_color)
# Fill Thermo-compositional Dynamo
plt.fill_between(T_cmb[:i_icn+1],0,sigma_dynamo2[:i_icn+1]*yscale,color='green',alpha=alpha_color)
# Fill Compositional Dynamo
plt.fill_between(T_cmb[:i_icn+1],sigma_dynamo2[:i_icn+1]*yscale,10.,color='blue',alpha=alpha_color)
# Plot curves.
plot(T_cmb,sigma_dynamo2*yscale,linewidth=2)
plot(T_cmb,sigma_crit*yscale,linewidth=2)
plot(T_cmb,sigma_T*yscale,color='g',linewidth=2)
plt.xlabel(r'$T_{cmb}$ [K]')
plt.ylabel(r'$\sigma$ [$10^{6}$ $\Omega^{-1}$m$^{-1}$]')
plt.xlim([3.5e3,4.5e3]); plt.ylim([0.5,1.5])

if args.save:
    figpref='Figure%d'%nfig
    plt.savefig(figpref+figsuffix+args.savename+'.'+args.saveformat,format=args.saveformat)

# Plot Rm vs sigma
print('Figure 5')
nfig=5
plt.figure(nfig,figsize=(30,10))
# Use loglog scales.
plt.subplot(121)
for i in range(n_T2plot):
    plt.loglog(sigma_arr,rm2[i_T2plot[i],:],label=str(T_cmb[i_T2plot[i]].astype(int))+' K',
                   linewidth=2.,color=colors[i])
    plt.axvline(sigma_crit[i_T2plot[i]],color=colors[i])
    plt.axvline(sigma_dynamo2[i_T2plot[i]],linestyle='-',color=colors[i])
plot((1,1e10),(rm_crit,rm_crit),color='0.5')
plt.xlabel(r'Electrical Conductivity $\sigma$ [$\Omega^{-1}m^{-1}$]'); plt.ylabel('Rm')
plt.xlim(1e4,1e7); plt.ylim(1e1,1e4)
plt.axhspan(5e2,2e3,color='grey',alpha=0.1)
plt.legend(loc='best')
# Same but with linear scales.
plt.subplot(122)
sigma_scale = 1e-6
for i in range(n_T2plot):
    plot(sigma_arr*sigma_scale,rm2[i_T2plot[i],:],label=str(T_cmb[i_T2plot[i]].astype(int))+' K',
         linewidth=2.,color=colors[i])
    plot((sigma_dynamo2[i_T2plot[i]]*sigma_scale),(rm_crit),'o',color=colors[i]) #sigma_D2
    plot((sigma_dynamo1[i_T2plot[i]]*sigma_scale),(rm_crit),'o',color=colors[i]) #sigma_D1
    plot((sigma_peak[i_T2plot[i]]*sigma_scale),(rm_peak[i_T2plot[i]]),'o',color=colors[i]) #sigma_peak
plt.axhline(rm_crit,color='0.5')
plt.xlabel(r'Electrical Conductivity $\sigma$ [$\times10^6$ $\Omega^{-1}m^{-1}$]'); plt.ylabel('Rm')
plt.xlim(0,1.2)
plt.ylim(0,5e3)
plt.axhspan(5e2,2e3,color='grey',alpha=0.1)
plt.legend(loc='best',frameon=False,fontsize=20)
if args.save:
    figpref='Figure%d'%nfig
    plt.savefig(figpref+figsuffix+args.savename+'.'+args.saveformat,format=args.saveformat)

### Find where thermal stratification equals compositional buoyancy.
# Quadratic adiabat and liquidus.  THIS IS DIFFERENT THAN THE ABOVE VERSION!!!!!!!!!!!!!!!!!!
T_ad1=(T_ic_e-T_cmb_e)/(r_cmb**2-r_ic_e**2)
T_liq1=3527.  #Du2017 table S4
T_liq2=(T_ic_e-T_liq1)/(r_cmb**2-r_ic_e**2)
q_cmb=Q_cmb/area_cmb
q_cmb2=k_lm**(2./3)*(eta_lm*T_cmb)**(4./3)*(thermexp_man*rho_man*g*cp_man/(ra_crit*visc_tbl))**(1./3)
energy_ic=rho_core*(latent_i)   #+gravitational??
# IC
r_ic_2=np.zeros(nt)
for i in range(nt): r_ic_2[i]=max(0,(T_liq1-T_cmb[i])/(T_liq2-T_ad1)+r_cmb**2)  #>=0. to avoid sqrt of a neg
r_ic=np.sqrt(r_ic_2)   #Du17 Supp eq 10
r_ic=np.clip(r_ic,0,r_cmb)
T_icn=T_cmb[np.min(np.where(r_ic==0))]  #T_cmb(ICN)
d_ic=r_cmb-r_ic  #shell thickness with variable r_ic
i_solid=1+np.max(np.where(d_ic<=0.)[0])  #1+index of T_cmb where r_ic=r_cmb, fully solid.
# Compute r_ic, g_i
r_ic_dot=np.zeros(nt)
bflux_comp=np.zeros(nt)
Gamma_ic=np.zeros(nt)
for i in range(nt):
    if r_ic[i]>0 and r_ic[i]<r_cmb:
        r_ic_dot[i]=-q_cmb[i]/(2./3*cp_core*r_cmb*rho_core*r_ic[i]*(T_ad1-T_liq2)-energy_ic*(r_ic[i]/r_cmb)**2)  #10
        bflux_comp[i]=g*r_ic[i]/r_cmb*(delta_rho_ic/rho_core+thermexp*latent_i/cp_core)*r_ic_dot[i]
        Gamma_ic[i]=(cp_core/thermexp*delta_rho_ic/rho_core+latent_i)/(2./3*cp_core*r_cmb**2*(T_ad1-T_liq2)-energy_ic/rho_core*r_ic[i]/r_cmb)
bflux_total=np.zeros((nt,nsig))
# Thermo-compositional combined.
Ra_Q_tc=np.zeros((nt,nsig))
u_tc=np.zeros((nt,nsig))
rm_tc=np.zeros((nt,nsig))
for isig in range(nsig):
    bflux_total[i_solid:,isig]=bflux_thermal[i_solid:,isig]+bflux_comp[i_solid:]
    #Ra_Q_tc[i_solid:,isig]=bflux_total[i_solid:,isig]*(r_cmb**2/d_ic[i_solid:]**4)/pow(omega,3.) #Aubert2009 (13).
    Ra_Q_tc[:,isig]=bflux_total[:,isig]*(r_cmb/r_ic_e)/d_e**2/pow(omega,3.)
    ind_pos=range(nt)
    if np.min(Ra_Q_tc[:,isig])<0:
        ind_pos=np.argwhere(Ra_Q_tc[:,isig] >= 0.)  #keep only positive nums for fractional exponent.
    #print('isig=%d'%isig)
    u_tc[ind_pos,isig]=ro_0*omega*d_ic[ind_pos]*pow(Ra_Q_tc[ind_pos,isig],alpha)
    rm_tc[ind_pos,isig]=mu_0*sigma_arr[isig]*u_tc[ind_pos,isig]*d_ic[ind_pos]
    
# Compute sigma(T) for "thermal stratification=compositional convection".
sigma_buoybal=1./(2*lorenz*T_cmb*T_ad1*r_cmb)*(q_cmb+(rho_core*cp_core)/(thermexp*g)*bflux_comp)

#### Print some detailed numbers.
i_T_cold = 55
print('The coldest temperature before the core completely solidifies is T_cmb=%d, which corresponds to Q_cmb=%.2e W.'%(T_cmb[i_T_cold],Q_cmb[i_T_cold]))
print('For sigma=%.2e , F_th(T_cmb=%d)=%.2e, F_comp(T_cmb=%d)=%.2e'%(sigma_arr[i_e_sigma],T_cmb[i_e],bflux_thermal[i_e,i_e_sigma],T_cmb[i_e],bflux_comp[i_e]))
print('For sigma=%.2e , F_th(T_cmb=%d)=%.2e, F_comp(T_cmb=%d)=%.2e'%(sigma_arr[i_e_sigma],T_cmb[i_T_cold],bflux_thermal[i_T_cold,i_e_sigma],T_cmb[i_T_cold],bflux_comp[i_T_cold]))
if args.mode==1: i_bal=373.
if args.mode==2 or args.mode==3: i_bal=np.max(np.where(bflux_thermal[i_T_cold,:]+bflux_comp[i_T_cold]>=0.)[0])
print('For sigma_balance=%.2e, F_th(T_cmb=%d)=%.2e, F_comp(T_cmb=%d)=%.2e'%(sigma_arr[i_bal],T_cmb[i_T_cold],bflux_thermal[i_T_cold,i_bal],T_cmb[i_T_cold],bflux_comp[i_T_cold]))


# Plot scaled general form: f(x)
# From Appendix of paper.
print('Figure 9')
nfig=9
plt.figure(nfig,figsize=(15,10))
num=1e2
x=np.linspace(1e-9,1,num)
betai=[1./3,1,1]
alphai=[2./5,2./5,4./9]
labels=[r'$\alpha=2/5$, $\beta=1/3$',r'$\alpha=2/5$, $\beta=1$',r'$\alpha=4/9$, $\beta=1$']
colors=['blue','green','orange']
npar=len(betai)
f=np.zeros((npar,num))
x1=0.; x2=1.
x3=np.zeros(npar)  #peak x
f3=np.zeros(npar)  #peak f
for i in range(3):
    f[i]=x**(alphai[i]+1.)*(x**(-betai[i])-1)**alphai[i]
    x3[i]=(1-betai[i]/(1+1./alphai[i]))**(1./betai[i])
    f3[i]=x3[i]**(1+alphai[i])*(x3[i]**(-betai[i])-1.)**(alphai[i])
    plot(x,f[i],label=labels[i],color=colors[i])
    plot(x3[i],f3[i],'o',color=colors[i])
plt.ylabel('f')
plt.xlabel('x')
plt.legend(loc='best')
if args.save:
    figpref='Figure%d'%nfig
    plt.savefig(figpref+figsuffix+args.savename+'.'+args.saveformat,format=args.saveformat)


### Compare modes 2 and 3.
print('Figure 11')
nfig=11
plt.figure(nfig,figsize=(15,10))
plt.subplot(221)
plot(T_cmb,T_lm_m2,label=r'$Q_{sec,man}=Q_{sec,core}$')
plot(T_cmb,T_lm_m3,label=r'$T_{LM}=2500$ K')
plt.xlim([3.5e3,4.5e3])
plt.legend(loc='best',fontsize=12)
plt.xlabel(r'$T_{cmb}$'); plt.ylabel(r'$T_{LM}$')
plt.subplot(222)
plot(T_cmb,delta_T_m2,label='m2')
plot(T_cmb,delta_T_m3,label='m3')
plt.xlim([3.5e3,4.5e3])
plt.xlabel(r'$T_{cmb}$'); plt.ylabel(r'$\Delta T_{LM}$')
plt.subplot(223)
plot(T_cmb,Q_cmb_m2*1e-12,label='m2')
plot(T_cmb,Q_cmb_m3*1e-12,label='m3')
plt.xlim([3.5e3,4.5e3])
plt.xlabel(r'$T_{cmb}$'); plt.ylabel(r'$Q_{cmb}$ [TW]')
plt.subplot(224)
plot(T_cmb,sigma_2_m2*1e-6,label='m2')
plot(T_cmb,sigma_2_m3*1e-6,label='m3')
plt.axvline(T_icn,linestyle='--',color='k',linewidth=2.)  #vertical at T_icn
plot(T_cmb,sigma_T*yscale,color='g',linewidth=2)
plt.xlabel(r'$T_{cmb}$')
plt.ylabel(r'$\sigma$ [$10^{6}$ $\Omega^{-1}$m$^{-1}$]')
plt.xlim([3.5e3,4.5e3]); plt.ylim([0.5,1.5])
plt.tight_layout()
if args.save:
    figpref='Figure%d'%nfig
    plt.savefig(figpref+figsuffix+args.savename+'.'+args.saveformat,format=args.saveformat)





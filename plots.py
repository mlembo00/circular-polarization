import sys
import numpy as np
import pylab as plt
from matplotlib import gridspec
import camb
from camb import model, initialpower

twopi = 2.0*np.pi
fourpi = 4.0*np.pi
lmax_v = 1999
lmax_e_b = 1999
tau=0.05905

betaV_squared = 3.15e-2 #4.0e-2
betaE_squared = 37.5 #55.0

betaE_squared2 = 0.1

pars = camb.read_ini('../../Documents/Work-University/Codes/camb-puliti/CAMB-0.1.7/params.ini')
# pars = camb.CAMBparams(WantTensors = True,  max_l = 5500, max_eta_k = 10000.0,  max_l_tensor = 4000,  max_eta_k_tensor = 8000.0, YHe = 0.241, num_nu_massless = 2.046, share_delta_neff = True); 
# pars.set_cosmology(H0=67.32, ombh2=0.02236, omch2=0.11201, mnu=0.06, omk=0, tau=0.05905);  
# pars.InitPower.set_params(As=2.10100312e-9, ns=0.96508, r=0.07); 
pars.Reion.Reionization=True 
# #print(pars)   
# pars.set_for_lmax(3500, lens_potential_accuracy=0); 
results = camb.get_results(pars)
powers =results.get_cmb_power_spectra(pars, CMB_unit='muK')
cl_tilde_R=powers['unlensed_total']

tildedle_R = cl_tilde_R[2:, 1]
tildedlb_R = cl_tilde_R[2:, 2]

pars.Reion.Reionization=False 
results = camb.get_results(pars)
powers =results.get_cmb_power_spectra(pars, CMB_unit='muK')
cl_tilde=powers['unlensed_total']

tildedle = cl_tilde[2:, 1]*np.exp(-2.0*tau)
tildedlb = cl_tilde[2:, 2]*np.exp(-2.0*tau)
lmax = len(tildedle)

sR_e = tildedle_R - tildedle
sR_b = tildedlb_R - tildedlb

prefac = twopi/(np.arange(2, lmax+1)*np.arange(3, lmax+2))
prefac_p1 = twopi/(np.arange(3, lmax+2)*np.arange(4, lmax+3))
prefac_m1 = twopi/(np.arange(1, lmax)*np.arange(2, lmax+1))
prefac_p2 = twopi/(np.arange(4, lmax+3)*np.arange(5, lmax+4))
prefac_m2 = twopi/(np.arange(2, lmax)*np.arange(1, lmax-1))
prefac_m2 = np.insert(prefac_m2, 0, 0)

tildedle_p1 = tildedle[1:]
tildedle_m1 = np.insert(tildedle, 0, 0)
tildedlb_p1 = tildedlb[1:]
tildedlb_m1 = np.insert(tildedlb, 0, 0)
tildedlb_p2 = tildedlb[2:]
tildedlb_m2 = np.insert(tildedlb, 0, 0)
tildedlb_m2 = np.insert(tildedlb_m2, 0, 0)

fw3j_ee_bb = np.loadtxt('./code_cl-camb2modified-cl/fw3j_ee_bb.txt', unpack=True)
K11 = fw3j_ee_bb[0]
K22p1 = fw3j_ee_bb[1]
K22p1 = K22p1[1:]
K22m1 = fw3j_ee_bb[2]

fw3j_vv = np.loadtxt('./code_cl-camb2modified-cl/fw3j_vv.txt', unpack=True)
K44 = fw3j_vv[0]
K33p1 = fw3j_vv[1]
K33p1 = K33p1[1:]
K33m1 = fw3j_vv[2]
K44p2 = fw3j_vv[3]
K44p2 = K44p2[2:]
K44m2 = fw3j_vv[4]


cle = prefac[:lmax_e_b]*tildedle[:lmax_e_b] + 1/fourpi * betaV_squared * (K11[:lmax_e_b]*prefac[:lmax_e_b]*tildedle[:lmax_e_b] + \
			K22p1[:lmax_e_b]*prefac_p1[:lmax_e_b]*tildedlb_p1[:lmax_e_b] + K22m1[:lmax_e_b]*prefac_m1[:lmax_e_b]*tildedlb_m1[:lmax_e_b]) 


clb = prefac[:lmax_e_b]*tildedlb[:lmax_e_b] + 1/fourpi * betaV_squared * (K11[:lmax_e_b]*prefac[:lmax_e_b]*tildedlb[:lmax_e_b] + \
			K22p1[:lmax_e_b]*prefac_p1[:lmax_e_b]*tildedle_p1[:lmax_e_b] + K22m1[:lmax_e_b]*prefac_m1[:lmax_e_b]*tildedle_m1[:lmax_e_b])


clv = 1/np.pi * betaE_squared * (K44[:lmax_v]*prefac[:lmax_v]*tildedlb[:lmax_v] + K44p2[:lmax_v]*prefac_p2[:lmax_v]*tildedlb_p2[:lmax_v] \
										+ K44m2[:lmax_v]*prefac_m2[:lmax_v]*tildedlb_m2[:lmax_v] + \
											K33p1[:lmax_v]*prefac_p1[:lmax_v]*tildedle_p1[:lmax_v] + K33m1[:lmax_v]*prefac_m1[:lmax_v]*tildedle_m1[:lmax_v])

clv2 = 1/np.pi * betaE_squared2 * (K44[:lmax_v]*prefac[:lmax_v]*tildedlb[:lmax_v] + K44p2[:lmax_v]*prefac_p2[:lmax_v]*tildedlb_p2[:lmax_v] \
										+ K44m2[:lmax_v]*prefac_m2[:lmax_v]*tildedlb_m2[:lmax_v] + \
											K33p1[:lmax_v]*prefac_p1[:lmax_v]*tildedle_p1[:lmax_v] + K33m1[:lmax_v]*prefac_m1[:lmax_v]*tildedle_m1[:lmax_v])


dle = cle*(np.arange(2, lmax_e_b+2)*np.arange(3, lmax_e_b+3)/twopi) + sR_e[:len(cle)]
dlb = clb*(np.arange(2, lmax_e_b+2)*np.arange(3, lmax_e_b+3)/twopi) + sR_b[:len(clb)]
dlv = clv*(np.arange(2, lmax_v+2)*np.arange(3, lmax_v+3)/twopi)
dlv2 = clv2*(np.arange(2, lmax_v+2)*np.arange(3, lmax_v+3)/twopi)
ell = np.arange(2,lmax_e_b+2)

#lclass,bpclass,errclass=np.loadtxt('./data/Vmodes_CLASS.txt',unpack=True)



fig = plt.figure()
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
gs = gridspec.GridSpec(2, 1, height_ratios=[2.15, 0.95]) 
spectrum = plt.subplot(gs[0])
spectrum.loglog(ell, dle, c = 'tab:blue', lw='2', label=r'$EE$') #c = ''mediumblue''
spectrum.loglog(ell, tildedle_R[:lmax_e_b], c = 'tab:blue', lw='2', ls='--', label=r'$\widetilde{EE}$')
spectrum.loglog(ell, dlb, c = 'tab:orange', lw='2',label=r'$BB$')  #firebrick
spectrum.loglog(ell, tildedlb_R[:lmax_e_b], c = 'tab:orange',lw='2', ls='--', label=r'$\widetilde{BB}$') #forestgreen
spectrum.loglog(ell, dlv, c = 'tab:green', lw='2',label=r'$VV$') #darkgreen
#spectrum.loglog(ell, dlv2, c = 'tab:red', lw='2',label=r'$VV$ ($\beta_E^2$ = 0.1)') #darkgreen
#spectrum.errorbar(lclass,bpclass+errclass*2, c='tab:red',fmt='o',xerr=6,yerr=10**(np.log10(bpclass+errclass*2)+0.15)-(bpclass+errclass*2),uplims=np.ones(len(lclass),dtype=bool),label='CLASS 40GHz')
residual = plt.subplot(gs[1], sharex = spectrum)
residual.loglog(ell, dle - tildedle_R[:lmax_e_b], c = 'tab:blue', lw='2', label=r'$EE$')
plt.setp(spectrum.get_xticklabels(), visible=False)
leg1=spectrum.legend(loc='upper left', ncol=3,fontsize=12)
leg2=residual.legend(loc='upper left',fontsize=12)
leg1.get_frame().set_linewidth(0.0)
leg2.get_frame().set_linewidth(0.0)
spectrum.set_ylabel(r'$\mathcal{D}_\ell\;[\mu K^2]$', fontsize=18, labelpad=10)
residual.set_ylabel(r'$\Delta  \mathcal{D}_\ell\;[\mu K^2]$', fontsize=18, labelpad=10)
residual.set_xlabel(r'Multipole - $\ell$', fontsize=18)
spectrum.tick_params('y', labelrotation=90)
#spectrum.set_ylim(1.5e-7,2e3)
spectrum.set_ylim(2e-6,2e3)
residual.set_ylim(1e-8,5e-5)
#yaxis.set_tick_params(rotation=90)
#set_ticklabels(rotation=90)
#tick_params('y', labelrotation=90)
residual.tick_params('y', labelrotation=90)
#set_xticklabels(rotation=90)
#tick_params('y', labelrotation=90)
#plt.yticks(rotation=90, va='center')
plt.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.12, hspace=.0)
plt.yticks(va="center")
plt.show()

#sys.exit(0)

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

lclass,bpclass,errclass=np.loadtxt('./data/Vmodes_CLASS.txt',unpack=True)

lspider95,bpspider95,errspider95ul=np.loadtxt('./data/Spider_Vmodes_95GHz_spectrum.txt',unpack=True,usecols=(0,1,4))
lspider150,bpspider150,errspider150ul=np.loadtxt('./data/Spider_Vmodes_150GHz_spectrum.txt',unpack=True,usecols=(0,1,4))

plt.errorbar(lspider95,errspider95ul,fmt='o',xerr=12.5,yerr=10**(np.log10(errspider95ul)+0.15)-errspider95ul,uplims=np.ones(len(lspider95),dtype=bool),label='SPIDER 95GHz')
plt.errorbar(lspider150,errspider150ul,fmt='o',xerr=12.5,yerr=10**(np.log10(errspider150ul)+0.15)-errspider150ul,uplims=np.ones(len(lspider150),dtype=bool),label='SPIDER 150GHz')
plt.errorbar(lclass,bpclass+errclass*2,fmt='o',xerr=6,yerr=10**(np.log10(bpclass+errclass*2)+0.15)-(bpclass+errclass*2),uplims=np.ones(len(lclass),dtype=bool),label='CLASS 40GHz')

leg1=plt.legend(loc=4,ncol=1,fontsize=12)
leg1.get_frame().set_linewidth(0.0)

plt.plot(ell,dlv,'k',lw=2, label=r'$\beta_E^2$ = %s' % betaE_squared)
plt.plot(ell, dlv2, 'k', ls='--', lw='2',label=r'$\beta_E^2$ = 0.1')

leg2=plt.legend(loc=4,ncol=2,fontsize=12)
leg2.get_frame().set_linewidth(0.0)

plt.yscale('log')



plt.xlim(1.5,310)
plt.ylim(1e-4,5e3)
plt.ylabel(r'$\mathcal{D}^{VV}_\ell\;[\mu K^2]$',fontsize=18, labelpad=10)
plt.yticks(rotation='vertical',va='center')
plt.xlabel(r'Multipole - $\ell$',fontsize=18)
plt.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.12, hspace=.0)
#plt.subplots_adjust(top=0.998,left=0.105,bottom=0.125,right=0.998)
plt.show()



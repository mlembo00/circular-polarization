import sys
import numpy as np
import pylab as plt
from astropy.table import Table
from astropy.io import fits
from matplotlib import gridspec
import camb
from camb import model, initialpower

# author M.Lembo (using pythoncamb)

# argv1=betaV_squared, argv2=betaE_squared

# Values used in the paper
# betaV_squared = 3.15e-2
# betaE_squared = 0.140

betaV_squared = float(sys.argv[1])
betaE_squared = float(sys.argv[2])

directory =  '/Users/mac/Desktop/CircularPolarisation-project/plots/GFEmodified-spectra'# 'path_of_your_directory'
Want_fig = False  # if False no fig. are produced.

twopi = 2.0 * np.pi
fourpi = 4.0 * np.pi
# this correspond to lmax=2000, this value can be change. The max value is 3987.
lmax_v = 1999
lmax_e_b = 1999
tau = 0.05905

pars = camb.read_ini('/Users/mac/Documents/Work-University/Codes/camb-puliti/CAMB-0.1.7/params.ini')
pars.Reion.Reionization = True
results = camb.get_results(pars)
powers = results.get_cmb_power_spectra(pars, CMB_unit='muK')
cl_tilde_R = powers['unlensed_total']

tildedle_R = cl_tilde_R[2:, 1]
tildedlb_R = cl_tilde_R[2:, 2]
tildedlte_R = cl_tilde_R[2:, 3]

pars.Reion.Reionization = False
results = camb.get_results(pars)
powers = results.get_cmb_power_spectra(pars, CMB_unit='muK')
cl_tilde = powers['unlensed_total']

tildedle = cl_tilde[2:, 1] * np.exp(-2.0 * tau)
tildedlb = cl_tilde[2:, 2] * np.exp(-2.0 * tau)
tildedlte = cl_tilde[2:, 3] * np.exp(-2.0 * tau)
lmax = len(tildedle)

sR_e = tildedle_R - tildedle
sR_b = tildedlb_R - tildedlb
sR_te = tildedlte_R - tildedlte

prefac = twopi / (np.arange(2, lmax + 1) * np.arange(3, lmax + 2))
prefac_p1 = twopi / (np.arange(3, lmax + 2) * np.arange(4, lmax + 3))
prefac_m1 = twopi / (np.arange(1, lmax) * np.arange(2, lmax + 1))
prefac_p2 = twopi / (np.arange(4, lmax + 3) * np.arange(5, lmax + 4))
prefac_m2 = twopi / (np.arange(2, lmax) * np.arange(1, lmax - 1))
prefac_m2 = np.insert(prefac_m2, 0, 0)

tildedle_p1 = tildedle[1:]
tildedle_m1 = np.insert(tildedle, 0, 0)
tildedlb_p1 = tildedlb[1:]
tildedlb_m1 = np.insert(tildedlb, 0, 0)
tildedlb_p2 = tildedlb[2:]
tildedlb_m2 = np.insert(tildedlb, 0, 0)
tildedlb_m2 = np.insert(tildedlb_m2, 0, 0)

fw3j_ee_bb = np.loadtxt('%s/fw3j_ee_bb.txt' % directory, unpack=True)
K11 = fw3j_ee_bb[0]
K22p1 = fw3j_ee_bb[1]
K22p1 = K22p1[1:]
K22m1 = fw3j_ee_bb[2]

fw3j_vv = np.loadtxt('%s/fw3j_vv.txt' % directory, unpack=True)
K44 = fw3j_vv[0]
K33p1 = fw3j_vv[1]
K33p1 = K33p1[1:]
K33m1 = fw3j_vv[2]
K44p2 = fw3j_vv[3]
K44p2 = K44p2[2:]
K44m2 = fw3j_vv[4]

cle = prefac[:lmax_e_b] * tildedle[:lmax_e_b] - \
      1.0 / fourpi * (betaV_squared + betaE_squared) * prefac[:lmax_e_b] * tildedle[:lmax_e_b] + \
      1.0 / fourpi * betaV_squared * (K11[:lmax_e_b] * prefac[:lmax_e_b] * tildedle[:lmax_e_b] +
                                      K22p1[:lmax_e_b] * prefac_p1[:lmax_e_b] * tildedlb_p1[:lmax_e_b] +
                                      K22m1[:lmax_e_b] * prefac_m1[:lmax_e_b] * tildedlb_m1[:lmax_e_b])

clb = prefac[:lmax_e_b] * tildedlb[:lmax_e_b] - \
      1.0 / fourpi * (betaV_squared + betaE_squared) * prefac[:lmax_e_b] * tildedlb[:lmax_e_b] + \
      1.0 / fourpi * betaV_squared * (K11[:lmax_e_b] * prefac[:lmax_e_b] * tildedlb[:lmax_e_b] +
                                      K22p1[:lmax_e_b] * prefac_p1[:lmax_e_b] * tildedle_p1[:lmax_e_b] +
                                      K22m1[:lmax_e_b] * prefac_m1[:lmax_e_b] * tildedle_m1[:lmax_e_b])

clte = prefac[:lmax_e_b] * tildedlte[:lmax_e_b] - \
       0.5 / fourpi * (betaV_squared + betaE_squared) * prefac[:lmax_e_b] * tildedlte[:lmax_e_b]

clv = 1 / np.pi * betaE_squared * (
            K44[:lmax_v] * prefac[:lmax_v] * tildedlb[:lmax_v] +
            K44p2[:lmax_v] * prefac_p2[:lmax_v] * tildedlb_p2[:lmax_v]
            + K44m2[:lmax_v] * prefac_m2[:lmax_v] * tildedlb_m2[:lmax_v] +
            K33p1[:lmax_v] * prefac_p1[:lmax_v] * tildedle_p1[:lmax_v] +
            K33m1[:lmax_v] * prefac_m1[:lmax_v] * tildedle_m1[:lmax_v])

dle = cle * (np.arange(2, lmax_e_b + 2) * np.arange(3, lmax_e_b + 3) / twopi) + sR_e[:len(cle)]
dlb = clb * (np.arange(2, lmax_e_b + 2) * np.arange(3, lmax_e_b + 3) / twopi) + sR_b[:len(clb)]
dlte = clte * (np.arange(2, lmax_e_b + 2) * np.arange(3, lmax_e_b + 3) / twopi) + sR_te[:len(clte)]
dlv = clv * (np.arange(2, lmax_v + 2) * np.arange(3, lmax_v + 3) / twopi)
ell = np.arange(2, lmax_e_b + 2)

# SR_X: We only rotate the polarization at the last scattering.
# So we are adding the reionization part of the spectra unchanged "a posteriori".

spectra = Table([ell, dle, dlb, dlte, dlv], names=("ell", "EE", "BB", "TE", "VV"))
spectra.write("%s/spectra_lmax%s.fits" % (directory, len(ell)+1), format="fits")

if Want_fig == True:
    fig = plt.figure()
    plt.rc('xtick', labelsize=14)
    plt.rc('ytick', labelsize=14)
    gs = gridspec.GridSpec(2, 1, height_ratios=[2.15, 0.95])
    spectrum = plt.subplot(gs[0])
    spectrum.loglog(ell, dlv, c='none', lw='2', label=' ')
    spectrum.loglog(ell, dlv, c='tab:red', lw='2', label=r'$VV$')
    spectrum.loglog(ell, abs(dlte), c='tab:green', lw='2', label=r'$TE$')
    spectrum.loglog(ell, abs(tildedlte_R[:lmax_e_b]), c='tab:green', lw='2', ls='--', label=r'$\widetilde{TE}$')
    spectrum.loglog(ell, dlb, c='tab:orange', lw='2', label=r'$BB$')
    spectrum.loglog(ell, tildedlb_R[:lmax_e_b], c='tab:orange', lw='2', ls='--', label=r'$\widetilde{BB}$')
    spectrum.loglog(ell, dle, c='tab:blue', lw='2', label=r'$EE$')
    spectrum.loglog(ell, tildedle_R[:lmax_e_b], c='tab:blue', lw='2', ls='--', label=r'$\widetilde{EE}$')
    residual = plt.subplot(gs[1], sharex=spectrum)
    residual.loglog(ell, abs(dlte - tildedlte_R[:lmax_e_b]), c='tab:green', lw='2', label=r'$TE$')
    residual.loglog(ell, abs(dle - tildedle_R[:lmax_e_b]), c='tab:blue', lw='2', label=r'$EE$')
    plt.setp(spectrum.get_xticklabels(), visible=False)
    leg1 = spectrum.legend(loc='lower right', ncol=4, fontsize=12)
    leg2 = residual.legend(loc='lower right', ncol=2, fontsize=12)
    leg1.get_frame().set_linewidth(0.0)
    leg2.get_frame().set_linewidth(0.0)
    spectrum.set_ylabel(r'$\mathcal{D}_\ell\;[\mu K^2]$', fontsize=18, labelpad=10)
    residual.set_ylabel(r'$\Delta  \mathcal{D}_\ell\;[\mu K^2]$', fontsize=18, labelpad=10)
    residual.set_xlabel(r'Multipole - $\ell$', fontsize=18)
    spectrum.tick_params('y', labelrotation=90)
    spectrum.set_ylim(2e-8, 3e2)
    residual.set_ylim(2e-8, 7e1)
    residual.tick_params('y', labelrotation=90)
    plt.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.12, hspace=.0)
    plt.yticks(va="center")
    plt.savefig('%s/spettri.pdf' % directory, bbox_inches='tight', quality=100, format='pdf')
    plt.show()

import sys
import numpy as np
import camb
from camb import model, initialpower
import pylab as plt
from scipy import interpolate


twopi = 2.0*np.pi
lmax_v = 119
lik_exist = True
### choosen values for \beta^2_E ###
start = 0.0
end = 200.0
nstep = 1000000
betaEsquared_vec = np.linspace(start, end, nstep)

if lik_exist == False:
	### reading data ###
	data = np.loadtxt('Vmodes_CLASS.txt', unpack=True)

	bin_center = data[0]
	print(bin_center)
	dlv_data = data[1]
	print(dlv_data)
	sigma = data[2]
	deltabin = 12
	nbin = len(bin_center)
	print('num of bin:', nbin)

	### reading tilde cl ###

	pars = camb.read_ini('../../Documents/Work-University/Codes/camb-puliti/CAMB-0.1.7/params.ini')
	#pars.Reion.Reionization=True 
	results = camb.get_results(pars)
	powers =results.get_cmb_power_spectra(pars, CMB_unit='muK')
	cl_tilde=powers['unlensed_total']

	tildedle = cl_tilde[2:, 1]
	tildedlb = cl_tilde[2:, 2]
	lmax = len(tildedle)
	prefac = twopi/(np.arange(2, lmax+1)*np.arange(3, lmax+2))
	prefac_p1 = twopi/(np.arange(3, lmax+2)*np.arange(4, lmax+3))
	prefac_m1 = twopi/(np.arange(1, lmax)*np.arange(2, lmax+1))
	prefac_p2 = twopi/(np.arange(4, lmax+3)*np.arange(5, lmax+4))
	prefac_m2 = twopi/(np.arange(2, lmax)*np.arange(1, lmax-1))
	prefac_m2 = np.insert(prefac_m2, 0, 0)

	tildedle_p1 = tildedle[1:]
	tildedle_m1 = np.insert(tildedle, 0, 0)
	tildedlb_p2 = tildedlb[2:]
	tildedlb_m2 = np.insert(tildedlb, 0, 0)
	tildedlb_m2 = np.insert(tildedlb_m2, 0, 0)

	### reading fw3j_vv file ###

	fw3j_vv = np.loadtxt('./code_cl-camb2modified-cl/fw3j_vv.txt', unpack=True)
	K44 = fw3j_vv[0]
	K33p1 = fw3j_vv[1]
	K33p1 = K33p1[1:]
	K33m1 = fw3j_vv[2]
	K44p2 = fw3j_vv[3]
	K44p2 = K44p2[2:]
	K44m2 = fw3j_vv[4]


	LnLike_vec = []

	for betaE_squared in betaEsquared_vec:
		print(r'$\beta^2_E$ = ', betaE_squared)

		clv = 1/np.pi * betaE_squared * (K44[:lmax_v]*prefac[:lmax_v]*tildedlb[:lmax_v] + K44p2[:lmax_v]*prefac_p2[:lmax_v]*tildedlb_p2[:lmax_v] \
										+ K44m2[:lmax_v]*prefac_m2[:lmax_v]*tildedlb_m2[:lmax_v] + \
											K33p1[:lmax_v]*prefac_p1[:lmax_v]*tildedle_p1[:lmax_v] + K33m1[:lmax_v]*prefac_m1[:lmax_v]*tildedle_m1[:lmax_v])

		dlv = clv*(np.arange(2, lmax_v+2)*np.arange(3, lmax_v+3)/twopi)
		dlv = np.insert(dlv, 0, 0)
		#np.savetxt('test_v.txt',dlv)

		bin_min = np.arange(1,120,12, dtype='int')
		bin_max = np.arange(12,121,12, dtype='int')
	
		LnLike = 0
		for i in range(nbin):
			dlv_bin = np.sum(dlv[bin_min[i]:bin_max[i]])/deltabin
			#print(dlv_bin)
			#print(dlv_data[i])
	
			LnLike = LnLike + ((dlv_bin - dlv_data[i])**2 /(sigma[i])**2)
		
	
		LnLike_vec = np.append(LnLike_vec, LnLike) 
	
	LnLike_vec = -0.5*LnLike_vec
	LnLike_norm_vec = LnLike_vec - max(LnLike_vec)

	#print(LnLike_vec)
	np.savetxt('VLnlik_%s_%s_%s.txt' %(start, end, nstep), LnLike_norm_vec)
	#np.savetxt('Vlik_%s_%s_%s.txt' %(start, end, nstep), np.exp(-0.5*LnLike_vec))
	
	plt.plot(betaEsquared_vec, np.exp(LnLike_norm_vec))
	plt.show()



	sys.exit(0)

else:

	lim1 = 1 - 0.6827 #0.6827/2.0
	lim2 = 1 - 0.9545 #0.9545/2.0
	lim3 = 1 - 0.9973 #0.9973/2.0
	set_lim1 = False
	set_lim2 = False
	set_lim3 = False
	
	Lnl = np.loadtxt('VLnlik_%s_%s_%s.txt' %(start, end, nstep))
	l = np.exp(Lnl)
	p = np.linspace(start, end, nstep) 
	index = np.where(l == l.max())[0][0]
	f= open("constraints.txt","w+")
	f.write("best fit value: %5.3f\r\n" % p[index])  
	f.close()
	print('best fit value:', p[index])
	ll = l[:] #l[index:]
	pp = p[:] #p[index:]
	#plt.plot(pp,ll)
	#plt.show()
	max_steps = len(ll)-1
	step = 10

	total_area = np.trapz(ll,pp) 

	for i in range(max_steps, 0, -step):
		#print(i)
		if (np.trapz(ll[i:],pp[i:]) >= total_area*lim1) and (set_lim1 != True) :
			print('area:', np.trapz(ll[i:],pp[i:]), 'to be compared with:', total_area*lim1)
			#print('area:', np.trapz(ll[:i],pp[:i]), 'to be compared with:', total_area*lim1)
			print('1 sigma constraint:', pp[i])
			f=open("constraints.txt", "a+")
			f.write("1 sigma constraint: %5.3f\r\n" % pp[i])
			f.close()
			set_lim1 = True

		if (np.trapz(ll[i:],pp[i:]) >= total_area*lim2) and (set_lim2 != True) :
			print('area:', np.trapz(ll[i:],pp[i:]), 'to be compared with:', total_area*lim2)
			print('2 sigma constraint:', pp[i])
			f=open("constraints.txt", "a+")
			f.write("2 sigma constraint: %5.3f\r\n" % pp[i])
			f.close()
			set_lim2 = True

		if (np.trapz(ll[i:],pp[i:]) >= total_area*lim3) and (set_lim3 != True) :
			print('area:', np.trapz(ll[i:],pp[i:]), 'to be compared with:', total_area*lim3)
			print('3 sigma constraint:', pp[i])
			f=open("constraints.txt", "a+")
			f.write("3 sigma constraint: %5.3f\r\n" % pp[i])
			f.close()
			set_lim3 = True

		if (set_lim1== True) and (set_lim2 ==True) and (set_lim3 == True):
			sys.exit(0)

	
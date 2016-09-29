import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats, integrate
from scipy.stats import norm
from sklearn.neighbors import KernelDensity
from matplotlib  import pyplot,mpl

mpl.rcParams['xtick.major.pad'] = '10'
mpl.rcParams['ytick.major.pad'] = '10'

def tick_func_deg2kpc(X):
	kpc = 20.3e6*np.arctan(X*np.pi/180.)/1e3
	return ["%.1f" % z for z in kpc]

def tick_func_arcmin2kpc(X):
	kpc = 20.3e6*np.arctan(X/60.*np.pi/180.)/1e3
	return ["%.1f" % z for z in kpc]

for_ra = 54.6208#54.55229 image coords #54.6208 #ngc1399 ra
for_dec = -35.450742#-35.536447 image coords #-35.05056 #ngc1399 dec

ra      = np.genfromtxt('NGFS_FCC_cat_non_nucleated.dat', dtype=float, comments='#', delimiter='|', missing_values='-', skip_header=1, usecols=(2))
dec     = np.genfromtxt('NGFS_FCC_cat_non_nucleated.dat', dtype=float, comments='#', delimiter='|', missing_values='-', skip_header=1, usecols=(3))
ra_nuc  = np.genfromtxt('NGFS_FCC_cat_nucleated.dat', dtype=float, comments='#', delimiter='|', missing_values='-', skip_header=1, usecols=(2))
dec_nuc = np.genfromtxt('NGFS_FCC_cat_nucleated.dat', dtype=float, comments='#', delimiter='|', missing_values='-', skip_header=1, usecols=(3))

#ra, dec = np.loadtxt('finallist_new_sm_fcc.txt',unpack=True,comments='#')
#ra_nuc, dec_nuc = np.loadtxt('finallist_new_nuc_fcc.txt',unpack=True,comments='#')

dra = ra-for_ra# - ra # - for_ra
ddec = dec - for_dec #for_dec - dec# - for_dec

dra_nuc = ra_nuc-for_ra#-ra_nuc# - for_ra
ddec_nuc = dec_nuc - for_dec# - dec_nuc# - for_dec

dra_tot = np.concatenate([dra,dra_nuc])
ddec_tot = np.concatenate([ddec,ddec_nuc])

# plt.figure()
# plt.plot([np.max(dra_tot)+0.1,np.min(dra_tot)-0.1],[0,0],'k:')
# plt.plot([0,0],[np.min(ddec_tot)-0.1,np.max(ddec_tot)+0.1],'k:')
# #plt.plot(dra_tot,ddec_tot,'ko')
# plt.plot(dra,ddec,'bo')
# plt.plot(dra_nuc,ddec_nuc,'ro')
# plt.xlim(np.max(dra_tot)+0.1,np.min(dra_tot)-0.1)
# plt.ylim(np.min(ddec_tot)-0.1,np.max(ddec_tot)+0.1)
# plt.xlabel('$\Delta$RA [deg]')
# plt.ylabel('$\Delta$Dec [deg]')
# plt.savefig('dwarf_spatial_dist.pdf')

ang = []
for ii in range(len(dra)):
	if ddec[ii] == 0 and dra[ii] < 0:
		ang.append(270.)
	if ddec[ii] == 0 and dra[ii] > 0:
		ang.append(90.)
	if ddec[ii] > 0 and dra[ii] == 0:
		ang.append(0.)
	if ddec[ii] < 0 and dra[ii] == 0:
		ang.append(180.)
	
	if dra[ii] < 0 and ddec[ii] > 0:
		ang.append(270.+np.arctan(np.abs(ddec[ii]/dra[ii]))*(180./np.pi))
	if dra[ii] > 0 and ddec[ii] > 0:
		ang.append(np.arctan(dra[ii]/ddec[ii])*(180./np.pi))
	if dra[ii] > 0 and ddec[ii] < 0:
		ang.append(90.+np.arctan(np.abs(ddec[ii]/dra[ii]))*(180./np.pi))
	if dra[ii] < 0 and ddec[ii] < 0:
		ang.append(180.+np.arctan(np.abs(dra[ii]/ddec[ii]))*(180./np.pi))

ang_nuc = []
for ii in range(len(dra_nuc)):
	print ii, ddec_nuc[ii], dra_nuc[ii]
	if ddec_nuc[ii] == 0 and dra_nuc[ii] < 0:
		ang_nuc.append(270.)
	if ddec_nuc[ii] == 0 and dra_nuc[ii] > 0:
		ang_nuc.append(90.)
	if ddec_nuc[ii] > 0 and dra_nuc[ii] == 0:
		ang_nuc.append(0.)
	if ddec_nuc[ii] < 0 and dra_nuc[ii] == 0:
		ang_nuc.append(180.)
	
	if dra_nuc[ii] < 0 and ddec_nuc[ii] > 0:
		ang_nuc.append(270.+np.arctan(np.abs(ddec_nuc[ii]/dra_nuc[ii]))*(180./np.pi))
	if dra_nuc[ii] > 0 and ddec_nuc[ii] > 0:
		ang_nuc.append(np.arctan(dra_nuc[ii]/ddec_nuc[ii])*(180./np.pi))
	if dra_nuc[ii] > 0 and ddec_nuc[ii] < 0:
		ang_nuc.append(90.+np.arctan(np.abs(ddec_nuc[ii]/dra_nuc[ii]))*(180./np.pi))
	if dra_nuc[ii] < 0 and ddec_nuc[ii] < 0:
		ang_nuc.append(180.+np.arctan(np.abs(dra_nuc[ii]/ddec_nuc[ii]))*(180./np.pi))

ang = np.array(ang)
ang_nuc = np.array(ang_nuc)
ang_tot = np.concatenate([ang,ang_nuc])

rgc = np.sqrt((for_ra-ra)**2+(for_dec-dec)**2)*60. #radial distance from ngc1399 in arcmin
rgc_nuc = np.sqrt((for_ra-ra_nuc)**2+(for_dec-dec_nuc)**2)*60.
rgc_tot = np.concatenate([rgc,rgc_nuc])

# print len(ang), np.min(ang), np.max(ang)
# print len(rgc), np.min(rgc), np.max(rgc)
# 
# print len(ang_tot), np.min(ang_tot), np.max(ang_tot)
# print len(rgc_tot), np.min(rgc_tot), np.max(rgc_tot)
# 
# print len(ang_nuc), np.min(ang_nuc), np.max(ang_nuc)
# print len(rgc_nuc), np.min(rgc_nuc), np.max(rgc_nuc)
# #exit()

rad_lim = 60.

mask = np.where(rgc <= rad_lim,1,0)
rgc = np.compress(mask,rgc)
ang = np.compress(mask,ang)
mask = np.where(rgc_nuc <= rad_lim,1,0)
rgc_nuc = np.compress(mask,rgc_nuc)
ang_nuc = np.compress(mask,ang_nuc)
mask = np.where(rgc_tot <= rad_lim,1,0)
rgc_tot = np.compress(mask,rgc_tot)
ang_tot = np.compress(mask,ang_tot)

print np.min(rgc_tot), np.max(rgc_tot)

#Starting at 0.1, increase radius until 20 objects are included
#Compute surface density (Rgc < r), record radius
#Increase radius until another 20, and repeat until all dwarfs are found
drad = 0.1
num_win_tot = 20
num_win = num_win_tot
print num_win_tot

rad = drad
dens_tot = []
rad_tot = []
num_tot = []
num_chk = 0
rad_d = 0.
while rad <= rad_lim:
	num = 0
	mask = np.where(rgc_tot > rad_d,1,0)
	rgc_tot_cnt = np.compress(mask,rgc_tot)
	mask = np.where(rgc_tot_cnt <= rad,1,0)
	rgc_tot_cnt = np.compress(mask,rgc_tot_cnt)
	num = len(rgc_tot_cnt)
#	print rad_d, rad, rad-rad_d, num
	if num >= num_win:
		num_chk+=num
		num_tot.append(num_chk)
# 		print rad_d, rad, rad-rad_d, num
		dens_tot.append(num/(np.pi*(rad**2-rad_d**2)))
		rad_tot.append((rad_d+rad)/2.)
		rad_d=rad
		rad+=drad
	if num < num_win and rad+drad<=rad_lim:
		rad+=drad
	if rad+drad > rad_lim:
		num_chk+=num
		num_tot.append(num_chk)
		dens_tot.append(num/(np.pi*(rad**2-rad_d**2)))
		rad_tot.append((rad_d+rad)/2.)
		rad_d=rad
		rad+=drad

num_win = float(len(rgc))/len(rgc_tot)*num_win_tot
print num_win

rad = drad
dens_sm = []
rad_sm = []
num_sm = []
num_chk = 0
rad_d = 0.
while rad <= rad_lim:
	num = 0
	mask = np.where(rgc > rad_d,1,0)
	rgc_sm_cnt = np.compress(mask,rgc)
	mask = np.where(rgc_sm_cnt <= rad,1,0)
	rgc_sm_cnt = np.compress(mask,rgc_sm_cnt)
	num = len(rgc_sm_cnt)
	if num >= num_win:
		num_chk+=num
		num_sm.append(num_chk)
# 		print rad_d, rad, rad-rad_d, num
		dens_sm.append(num/(np.pi*(rad**2-rad_d**2)))
		rad_sm.append((rad_d+rad)/2.)
		rad_d=rad
		rad+=drad
	if num < num_win and rad+drad<=rad_lim:
		rad+=drad
	if rad+drad > rad_lim:
		num_chk+=num
		num_sm.append(num_chk)
# 		print rad_d, rad, rad-rad_d, num
		dens_sm.append(num/(np.pi*(rad**2-rad_d**2)))
		rad_sm.append((rad_d+rad)/2.)
		rad_d=rad
		rad+=drad

num_win = float(len(rgc_nuc))/len(rgc_tot)*num_win_tot
print num_win

rad = drad
dens_nuc = []
rad_nuc = []
num_nuc = []
num_chk = 0
rad_d = 0.
while rad <= rad_lim:
	num = 0
	mask = np.where(rgc_nuc > rad_d,1,0)
	rgc_nuc_cnt = np.compress(mask,rgc_nuc)
	mask = np.where(rgc_nuc_cnt <= rad,1,0)
	rgc_nuc_cnt = np.compress(mask,rgc_nuc_cnt)
	num = len(rgc_nuc_cnt)
	if num >= num_win:
		num_chk+=num
		num_nuc.append(num_chk)
# 		print rad_d, rad, rad-rad_d, num
		dens_nuc.append(num/(np.pi*(rad**2-rad_d**2)))
		rad_nuc.append((rad_d+rad)/2.)
		rad_d=rad
		rad+=drad
	if num < num_win and rad+drad<=rad_lim:
		rad+=drad
	if rad+drad > rad_lim:
		num_chk+=num
		num_nuc.append(num_chk)
# 		print rad_d, rad, rad-rad_d, num
		dens_nuc.append(num/(np.pi*(rad**2-rad_d**2)))
		rad_nuc.append((rad_d+rad)/2.)
		rad_d=rad
		rad+=drad

radii = np.array(rad_sm)
dens = np.array(dens_sm)

radii_nuc = np.array(rad_nuc)
dens_nuc = np.array(dens_nuc)

radii_tot = np.array(rad_tot)
dens_tot = np.array(dens_tot)


print "Average density of total sample = %.4f" % np.mean(dens_tot)
print "Average density of non-nucleated sample = %.4f" % np.mean(dens)
print "Average density of nucleated sample = %.4f" % np.mean(dens_nuc)



fig = plt.figure(figsize=(12,10))
#plt.plot([0.1,100],[0,0],'k:')

# #==========================================================================================================================================#                                                               
# theta,dist,w,w_err = np.genfromtxt('NGFS_FCC_all_dwarfs_two-point_correlation.dat',dtype=float,comments='#',delimiter='|',unpack=True,skip_header=1,usecols=(1,2,3,4))
# ax1 = fig.add_subplot(311)
# plt.xlim(0,1.03)
# plt.ylim(-0.175,0.5)
# ax1.plot(theta,w,linestyle='-',color='#6b4423',linewidth=3,ms=12,zorder=10,alpha=0.9)
# ax1.fill_between(theta, w+w_err, w-w_err, facecolor='#e0cda7',zorder=1,alpha=0.7)
# #ax1.errorbar(dist, w, yerr=w_err, fmt='o', facecolor='#e0cda7')
# ax1.set_xlabel(r'$\theta$ [deg]',fontsize=18)
# ax1.set_ylabel(r'$w(\theta)$',fontsize=18)
# ax1.legend(numpoints=1,frameon=False,loc='upper right',prop={'size':12})
# 
# axis = plt.gca()
# xmin, xmax = axis.get_xlim()
# ymin, ymax = axis.get_ylim()
# ax2 = ax1.twiny()
# #new_tick_locs = np.array([10,20,30,40,50,60])#ax1.get_xticks()
# new_tick_locs = np.linspace(xmin, xmax, 10)
# ax2.set_xticks(new_tick_locs)
# ax2.set_xticklabels(tick_func_deg2kpc(new_tick_locs))
# ax2.set_xlabel(r'$R_{gc}$ [kpc]',fontsize=12)
# ax2.tick_params('both',length=6,width=2,which='major',labelsize=10)
# ax2.tick_params('both',length=4,width=1,which='minor')
# 
# #==========================================================================================================================================#                                                               
# ax3 = fig.add_subplot(312)
# plt.xlim(0,62)
# ax3.plot(radii_tot,dens_tot,linestyle='-',color='#6b4423',linewidth=3,ms=12,zorder=-1,alpha=0.5)
# ax3.plot(radii,dens        ,linestyle='-',color='#bb1515',linewidth=3,ms=12,zorder=-1,alpha=0.5)
# ax3.plot(radii_nuc,dens_nuc,linestyle='-',color='#2a334f',linewidth=3,ms=12,zorder=-1,alpha=0.5)
# ax3.scatter(radii_tot,dens_tot,marker='o',color='#6b4423',s=120,zorder=10,alpha=1.0)
# ax3.scatter(radii,dens        ,marker='o',color='#bb1515',s=120,zorder=10,alpha=1.0)
# ax3.scatter(radii_nuc,dens_nuc,marker='o',color='#2a334f',s=120,zorder=10,alpha=1.0)
# # tot, = plt.plot(radii_lim,fit_func(radii_lim,*popt_tot),'k-',lw=4)
# ax3.set_xlabel(r'$R_{gc}$ [arcmin]',fontsize=18)
# ax3.set_ylabel(r'$\Sigma$ [arcmin$^{-2}$]',fontsize=18)
# ax3.set_yscale('log')
# ax3.set_ylim(1e-3,2e-1)
# ax3.legend((r'Total',r'Non-nucleated',r'Nucleated'),numpoints=1,frameon=False,loc='upper left',prop={'size':12})
# ax3.tick_params('both',length=6,width=2,which='major',labelsize=12)
# ax3.tick_params('both',length=4,width=1,which='minor')
# 
# axis = plt.gca()
# xmin, xmax = axis.get_xlim()
# ymin, ymax = axis.get_ylim()
# ax4 = ax3.twiny()
# #new_tick_locs = np.array([10,20,30,40,50,60])#ax1.get_xticks()
# new_tick_locs = np.linspace(xmin, xmax, 10)
# ax4.set_xticks(new_tick_locs)
# ax4.set_xticklabels(tick_func_arcmin2kpc(new_tick_locs))
# ax4.set_xlabel(r'$R_{gc}$ [kpc]',fontsize=12)
# ax4.tick_params('both',length=6,width=2,which='major',labelsize=10)
# ax4.tick_params('both',length=4,width=1,which='minor')
# 
# #==========================================================================================================================================#                                                               
ang_bs = 15.
#plt.subplot(313)
plt.xlim(0,360)
plt.ylim(-1.5,22)
axis = plt.gca()
xmin, xmax = axis.get_xlim()
ymin, ymax = axis.get_ylim()

plt.hist(ang_tot,bins=np.arange(np.min(ang_tot),np.max(ang_tot)+ang_bs,ang_bs),histtype='stepfilled',color='#6b4423',linestyle='solid',linewidth=2,zorder='1',alpha=0.25)
plt.hist(ang,bins=np.arange(np.min(ang_tot),np.max(ang_tot)+ang_bs,ang_bs),histtype='stepfilled',color='#bb1515',linewidth=0.5,zorder='2',alpha=0.4)
plt.hist(ang_nuc,bins=np.arange(np.min(ang_tot),np.max(ang_tot)+ang_bs,ang_bs),histtype='stepfilled',color='#2a334f',linewidth=0.5,zorder='3',alpha=0.4)
#plt.hist(ang,color='b',histtype='step',lw=2)
#plt.hist(ang_nuc,color='r',histtype='step',lw=2)
#plt.hist(ang_tot,color='k',histtype='step',lw=2)
plt.legend((r'Total',r'Non-nucleated',r'Nucleated'),frameon=False,loc='upper left',prop={'size':20})

bandw = 25.
samplingsteps = 3600
### ALL
mirr_ang_tot = np.concatenate((ang_tot-360,ang_tot,ang_tot+360)) # wrap vector before and after
normdat = len(mirr_ang_tot)*ang_bs
x_in = mirr_ang_tot[:, np.newaxis]
x_smpl = np.linspace(xmin-bandw*2, xmax+bandw*2, samplingsteps)[:, np.newaxis]
#for kernel in ['gaussian', 'tophat', 'epanechnikov']:
kde = KernelDensity(kernel='epanechnikov', bandwidth=bandw).fit(x_in)
log_dens = kde.score_samples(x_smpl)
plt.plot(x_smpl[:, 0],np.exp(log_dens)*normdat,linestyle='-',color='#6b4423',linewidth=5,alpha=0.5)
#plt.scatter(ang_tot, 21.4 - 0.25 * np.random.random(ang_tot.shape[0]),marker='o',s=50,edgecolor='none',zorder='5',facecolor='#6b4423',alpha=0.4)

### NON-NUCLEATED
mirr_ang = np.concatenate((ang-360,ang,ang+360)) # wrap vector before and after
normdat = len(mirr_ang)*ang_bs
x_in = mirr_ang[:, np.newaxis]
x_smpl = np.linspace(xmin-bandw*2, xmax+bandw*2, samplingsteps)[:, np.newaxis]
kde = KernelDensity(kernel='epanechnikov', bandwidth=bandw).fit(x_in)
log_dens = kde.score_samples(x_smpl)
plt.plot(x_smpl[:, 0],np.exp(log_dens)*normdat,linestyle='-',color='#bb1515',linewidth=4,alpha=0.5)
plt.scatter(ang, -0.4 - 0.25 * np.random.random(ang.shape[0]),marker='o',s=60,edgecolor='none',zorder='5',facecolor='#bb1515',alpha=0.4)

### NUCLEATED
mirr_ang_nuc = np.concatenate((ang_nuc-360,ang_nuc,ang_nuc+360)) # wrap vector before and after
normdat = len(mirr_ang_nuc)*ang_bs
x_in = mirr_ang_nuc[:, np.newaxis]
x_smpl = np.linspace(xmin-bandw*2, xmax+bandw*2, samplingsteps)[:, np.newaxis]
kde = KernelDensity(kernel='epanechnikov', bandwidth=bandw).fit(x_in)
log_dens = kde.score_samples(x_smpl)
plt.plot(x_smpl[:, 0],np.exp(log_dens)*normdat,linestyle='-',color='#2a334f',linewidth=4,alpha=0.5)
plt.scatter(ang_nuc, -0.9 - 0.25 * np.random.random(ang_nuc.shape[0]),marker='o',s=60,edgecolor='none',zorder='5',facecolor='#2a334f',alpha=0.4)

#plt.legend((r'KDE(Total)'),loc='top right',frameon=False)

plt.tick_params('both',length=6,width=2,which='major',labelsize=22)
plt.tick_params('both',length=4,width=1,which='minor')
plt.xlabel(r'$\Phi$ [deg E of N]',fontdict={'fontsize':22})
plt.ylabel(r'Number',fontdict={'fontsize':22})
plt.tight_layout()

plt.savefig('dwarf_azimuthal_distrib.pdf')
plt.show()
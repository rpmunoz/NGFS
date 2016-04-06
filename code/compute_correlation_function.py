import numpy as np
from matplotlib import pyplot as plt
from astropy.table import Table, Column, join
from astroML.correlation import bootstrap_two_point_angular,two_point_angular
from astropy.io import ascii
from tabulate import tabulate
import astropy.units as u
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=8, usetex=True)

decam_scale=0.263*u.arcsec
fornax_distance=20.0*u.Mpc
cat_dir='../catalogs'

cat_ngfs_nuc=ascii.read(cat_dir+'/NGFS_FCC_cat_nucleated.dat', Reader=ascii.FixedWidth)
cat_ngfs_non=ascii.read(cat_dir+'/NGFS_FCC_cat_non_nucleated.dat', Reader=ascii.FixedWidth)

#bins=np.linspace(0.005,1,11)
bins=np.arange(0.02,1.3,0.08)

corr=two_point_angular(cat_ngfs_nuc['RA'], cat_ngfs_nuc['DEC'], bins, method='landy-szalay')
print '\nNucleated dwarf galaxies'
print 'Angle(deg)   w(theta)'
for j in np.arange(len(bins)-1):
	print "{:.2f}   {:.2f}".format(bins[j], corr[j])
print 'Uniform distribution? ', np.allclose(corr, 0, atol=0.02)

#bins=np.linspace(0.005,1,11)
#corr=two_point_angular(cat_ngfs_non['RA'], cat_ngfs_non['DEC'], bins, method='landy-szalay')
#print '\nNon nucleated dwarf galaxies'
#print 'Angle(deg)   w(theta)'
#for j in np.arange(len(bins)-1):
#	print "{:.2f}   {:.2f}".format(bins[j], corr[j])
#print 'Uniform distribution? ', np.allclose(corr, 0, atol=0.02)

# Now the non-nucleated dwarfs

result = bootstrap_two_point_angular( cat_ngfs_non['RA'], cat_ngfs_non['DEC'], bins, method='landy-szalay', Nbootstraps=5000)
bin_centers = 0.5 * (bins[1:] + bins[:-1])
(corr, corr_err, bootstraps)=result
print 'Non-nucleated dwarfs - Uniform distribution?', np.allclose(corr, 0, atol=0.02)

fig = plt.figure(figsize=(8, 6))
plt.xlabel(r'$\theta\ (deg)$')
plt.ylabel(r'$\hat{w}(\theta)$')
plt.errorbar(bin_centers, corr, corr_err, fmt='.k', ecolor='gray', lw=1)
plt.savefig('../catalogs/NGFS_FCC_non_nucleated_dwarfs_two-point_correlation.pdf')

correlation_function=Table()
correlation_function['theta']=bin_centers
correlation_function['theta_kpc']=2*fornax_distance.to('kpc')*np.tan( (bin_centers*u.deg).to('rad')/2. )
correlation_function['w']=corr
correlation_function['w_err']=corr_err

correlation_function['theta'].format = '6.2f'
correlation_function['theta'].unit = 'degree'
correlation_function['theta_kpc'].format = '6.2f'
correlation_function['theta_kpc'].unit = 'kpc'
correlation_function['w'].format = '8.4f'
correlation_function['w_err'].format = '8.4f'
ascii.write(correlation_function, '../catalogs/NGFS_FCC_non_nucleated_dwarfs_two-point_correlation.dat', format='fixed_width')  


# Now all the dwarfs

result = bootstrap_two_point_angular( np.concatenate( (np.asarray(cat_ngfs_non['RA']), np.asarray(cat_ngfs_nuc['RA'])), axis=0), np.concatenate( (np.asarray(cat_ngfs_non['DEC']), np.asarray(cat_ngfs_nuc['DEC'])), axis=0), bins, method='landy-szalay', Nbootstraps=5000)

bin_centers = 0.5 * (bins[1:] + bins[:-1])
(corr, corr_err, bootstraps)=result
print 'All dwarfs - Uniform distribution?', np.allclose(corr, 0, atol=0.02)

fig = plt.figure(figsize=(8, 6))
plt.xlabel(r'$\theta\ (deg)$')
plt.ylabel(r'$\hat{w}(\theta)$')
plt.errorbar(bin_centers, corr, corr_err, fmt='.k', ecolor='gray', lw=1)
plt.savefig('../catalogs/NGFS_FCC_all_dwarfs_two-point_correlation.pdf')

correlation_function=Table()
correlation_function['theta']=bin_centers
correlation_function['theta_kpc']=2*fornax_distance.to('kpc')*np.tan( (bin_centers*u.deg).to('rad')/2. )
correlation_function['w']=corr
correlation_function['w_err']=corr_err

correlation_function['theta'].format = '6.2f'
correlation_function['theta'].unit = 'degree'
correlation_function['theta_kpc'].format = '6.2f'
correlation_function['theta_kpc'].unit = 'kpc'
correlation_function['w'].format = '8.4f'
correlation_function['w_err'].format = '8.4f'
ascii.write(correlation_function, '../catalogs/NGFS_FCC_all_dwarfs_two-point_correlation.dat', format='fixed_width')  

#print '\nAll dwarf galaxies'
#print 'Angle(deg)\t  Distance(kpc)\t  w(theta)'
#for j in np.arange(len(bins)-1):
#	print "{:6.2f}\t  {:6.2f}\t  {:5.2f}".format(bins[j], (20.3*u.Mpc).to('kpc')*(bins[j]*u.deg).to('rad'), corr[j])
#print 'Uniform distribution?', np.allclose(corr, 0, atol=0.02)

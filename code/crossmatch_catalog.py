import sys,os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy.table import Table, Column, join, vstack, hstack
import aplpy
import skimage
import skimage.io
from pyavm import AVM
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

im_rgb_sd_file= '../catalogs/ngfs_tile1_rgb_asinh_v2_sd.png'
cat_ferguson=ascii.read('../catalogs/Ferguson_1989_cat_all.dat')
cat_paul_nuc=ascii.read('../catalogs/Paul_nucleated.dat')
cat_paul_non=ascii.read('../catalogs/Paul_non_nucleated.dat')

cat_paul_galfit_nuc=ascii.read('../catalogs/Paul_GALFIT_nucleated.dat', fill_values=('-', '0.0'))
cat_paul_galfit_non=ascii.read('../catalogs/Paul_GALFIT_non_nucleated.dat', fill_values=('-', '0.0'))

cat_mieske=ascii.read('../catalogs/Mieske_2007_cat_dwarf.dat')
cat_richtler=ascii.read('../catalogs/Richtler_2015_private.dat')

decam_scale=0.263*u.arcsec
fornax_distance=20.0*u.Mpc

cat_ferguson.add_column( Column(['FCC'+str(j) for j in cat_ferguson['FCC']], name='ID'), index=0)

#coo=SkyCoord(cat_richtler['RA_hms'], cat_richtler['DEC_dms'], unit=(u.hourangle, u.deg))
#cat_richtler.add_column( Column(coo.ra.deg, name='RA'), index=1)
#cat_richtler.add_column( Column(coo.dec.deg, name='DEC'), index=2)
#ascii.write(cat_richtler, '../catalogs/Richtler_2015_private.dat', format='fixed_width')

# We revise the Mieske dwarfs
mieske_revised=Table(names=('ID', 'MType', 'Memb'), dtype=('S20', 'S10', 'f8'))
mieske_revised.add_row(('WFLSB4-3', 'd', 0.))
mieske_revised.add_row(('WFLSB11-8', 'd', 0.))
mieske_revised.add_row(('WFLSB1-9', 'd', 0.))
mieske_revised.add_row(('WFLSB10-5', 'dN?', 0.))

ferguson_revised=Table(names=('ID', 'MType', 'RA', 'DEC'), dtype=('S20', 'S10', 'f8', 'f8'))
# We revise the FCC non-nucleated
ferguson_revised.add_row(('FCC131','dN?', 0, 0))
ferguson_revised.add_row(('FCC187','d', 0, 0))
ferguson_revised.add_row(('FCC246','dN', 0, 0))
ferguson_revised.add_row(('FCC191','dN', 0, 0))
ferguson_revised.add_row(('FCC192','dN', 0, 0))
ferguson_revised.add_row(('FCC141','dE', 53.738396, -35.206796))
ferguson_revised.add_row(('FCC196','dN', 0, 0))
ferguson_revised.add_row(('FCC259','dN', 0, 0))
ferguson_revised.add_row(('FCC244','dN', 0, 0))
ferguson_revised.add_row(('FCC218','dN', 0, 0))
ferguson_revised.add_row(('FCC194','dN', 0, 0))
ferguson_revised.add_row(('FCC162','E?', 0, 0))
ferguson_revised.add_row(('FCC165','dN', 0, 0))
ferguson_revised.add_row(('FCC168','dN', 0, 0))
ferguson_revised.add_row(('FCC229','dN', 0, 0))
ferguson_revised.add_row(('FCC236','dN', 0, 0))

# We revise the FCC nucleated
ferguson_revised.add_row(('FCC238','dN', 55.071428, -36.534886 ))
ferguson_revised.add_row(('FCC243','dN', 55.1124, -36.499083 ))
ferguson_revised.add_row(('FCC245','dN', 55.140862, -35.022865 ))
ferguson_revised.add_row(('FCC188','dN', 54.26876, -35.590182 ))
ferguson_revised.add_row(('FCC136','dN?', 0, 0))
ferguson_revised.add_row(('FCC252','dN', 55.209934, -35.74846 ))
ferguson_revised.add_row(('FCC195','d', 0, 0))
ferguson_revised.add_row(('FCC257','d', 0, 0))
ferguson_revised.add_row(('FCC266','dN', 55.421939, -35.170375 ))
ferguson_revised.add_row(('FCC202','dN', 54.527266, -35.439882 ))
ferguson_revised.add_row(('FCC150','dN?', 53.850193, -36.363752 ))
ferguson_revised.add_row(('FCC264','dS0', 55.381966, -35.589611 ))
ferguson_revised.add_row(('FCC211','E?', 54.589284, -35.259921 ))
ferguson_revised.add_row(('FCC274','dN', 55.571741, -35.540775 ))
ferguson_revised.add_row(('FCC159','d', 0, 0))
ferguson_revised.add_row(('FCC222','dN', 54.805246, -35.371481 ))
ferguson_revised.add_row(('FCC164','d', 54.053605, -36.166272 ))
ferguson_revised.add_row(('FCC207','d', 54.580167, -35.12905 ))
ferguson_revised.add_row(('FCC230','d', 55.005457, -34.758155 ))


# We made visual inspection of Mieske 2007 dwarf galaxies and we replace the main catalog
for j in np.arange(len(mieske_revised)):
	gv=np.where(cat_mieske['ID'] == mieske_revised['ID'][j])
	cat_mieske['Memb'][gv]=mieske_revised['Memb'][j]
	print gv, cat_mieske['Memb'][gv]

# Now we do the crossmatch
im_data = np.flipud(skimage.io.imread(im_rgb_sd_file))
im_size= im_data.shape
avm=AVM.from_image(im_rgb_sd_file)
w = avm.to_wcs()
w.naxis1=im_size[1]
w.naxis2=im_size[0]

fig = plt.figure(figsize=(14.,8.05), dpi=300)
f1 = aplpy.FITSFigure(w, figure=fig)

gv_tile=[]
cat_coo = f1.world2pixel(np.asarray(cat_ferguson['RA']), np.asarray(cat_ferguson['DEC']))
for j in range(len(cat_ferguson)):
	if (0. < cat_coo[0][j] < im_size[1]) and (0. < cat_coo[1][j] < im_size[0]):
		temp=np.median(im_data[cat_coo[1][j]-2:cat_coo[1][j]+2,cat_coo[0][j]-2:cat_coo[0][j]+2,:])
		if temp > 0 : gv_tile.append(j)

gv_dwarf=[]
for j in range(len(cat_ferguson)):
	if cat_ferguson['MType'][j].startswith('d'): gv_dwarf.append(j)

gv_dwarf_tile=list( set(gv_tile) & set(gv_dwarf) )
cat_ferguson_dwarf=cat_ferguson[gv_dwarf_tile]
print 'Number of Ferguson dwarfs inside the Tile ', len(gv_dwarf_tile)

gv_nuc=[]
gv_non=[]
for j in range(len(cat_ferguson_dwarf)):
	if cat_ferguson_dwarf['MType'][j].startswith('d') and (cat_ferguson_dwarf['MType'][j].endswith('N') or cat_ferguson_dwarf['MType'][j].endswith('N?')): gv_nuc.append(j)
	else : gv_non.append(j)

print '\nTotal number of non-nucleated FCC dwarf galaxies in Tile1: ', len(gv_non)
print 'Total number of nucleated FCC dwarf galaxies in Tile1: ', len(gv_nuc)

gv_mieske_memb=np.where( cat_mieske['Memb'] <= 2.5)
gv_mieske_non=np.where( cat_mieske['Memb'] > 2.5)

cat_ferguson=cat_ferguson[gv_tile]
cat_paul=join(cat_paul_non, cat_paul_nuc, join_type='outer')
cat_mieske_memb=cat_mieske[gv_mieske_memb]
cat_mieske_non=cat_mieske[gv_mieske_non]

cat_ferguson_coo = SkyCoord(cat_ferguson['RA'], cat_ferguson['DEC'], unit="deg")
cat_paul_coo = SkyCoord(cat_paul['RA'], cat_paul['DEC'], unit="deg") #SkyCoord( np.concatenate( ([np.asarray(cat_paul_non['RA']), np.asarray(cat_paul_nuc['RA'])]) ), np.concatenate( ([np.asarray(cat_paul_non['DEC']), np.asarray(cat_paul_nuc['DEC'])]) ), unir="deg")
cat_mieske_memb_coo = SkyCoord(cat_mieske_memb['RA'], cat_mieske_memb['DEC'], unit="deg")
cat_mieske_non_coo = SkyCoord(cat_mieske_non['RA'], cat_mieske_non['DEC'], unit="deg")

cat_ferguson_nuc_coo = SkyCoord(cat_ferguson['RA'][gv_nuc], cat_ferguson['DEC'][gv_nuc], unit="deg")
cat_ferguson_non_coo = SkyCoord(cat_ferguson['RA'][gv_non], cat_ferguson['DEC'][gv_non], unit="deg")
cat_paul_non_coo = SkyCoord(cat_paul_non['RA'], cat_paul_non['DEC'], unit="deg")
cat_paul_nuc_coo = SkyCoord(cat_paul_nuc['RA'], cat_paul_nuc['DEC'], unit="deg")

cat_paul_galfit_non_coo = SkyCoord(cat_paul_galfit_non['RA'], cat_paul_galfit_non['DEC'], unit="deg")
cat_paul_galfit_nuc_coo = SkyCoord(cat_paul_galfit_nuc['RA'], cat_paul_galfit_nuc['DEC'], unit="deg")

print '\nCrossmatching Paul against FCC dwarf galaxy catalog'
print '---------------------------------'

idx, sep2d, sep3d = cat_paul_non_coo.match_to_catalog_sky(cat_ferguson_coo)
print '\nPaul non-nucleated already in FCC dwarf catalog: ', len(np.where(sep2d.arcsecond <= 12.)[0])
print 'Paul non-nucleated galaxies not present in the FCC dwarf catalog: ', len(np.where(sep2d.arcsecond > 12.)[0])

idx, sep2d, sep3d = cat_paul_nuc_coo.match_to_catalog_sky(cat_ferguson_coo)
print '\nPaul nucleated already in FCC dwarf catalog: ', len(np.where(sep2d.arcsecond <= 12.)[0])
print 'Paul nucleated galaxies not present in the FCC dwarf catalog: ', len(np.where(sep2d.arcsecond > 12.)[0])

print '\nCrossmatching FCC agaings Paul non-nucleated and nucleated dwarf galaxy catalog'
print '---------------------------------'

idx, sep2d, sep3d = cat_ferguson_non_coo.match_to_catalog_sky(cat_paul_non_coo)
gv_ferguson=np.where(sep2d.arcsecond <= 12.)[0]
bv_ferguson=np.where(sep2d.arcsecond > 12.)[0]
print '\nFerguson non-nucleated galaxies already in Paul non-nucleated catalog: ', len(gv_ferguson)
print '\nFerguson non-nucleated galaxies not present in Paul non-nucleated catalog: ', len(bv_ferguson)
print cat_ferguson_dwarf[gv_non][bv_ferguson]

idx, sep2d, sep3d = cat_ferguson_nuc_coo.match_to_catalog_sky(cat_paul_nuc_coo)
gv_ferguson=np.where(sep2d.arcsecond <= 12.)[0]
bv_ferguson=np.where(sep2d.arcsecond > 12.)[0]
print '\nFerguson nucleated galaxies already in Paul nucleated catalog: ', len(gv_ferguson)
print '\nFerguson nucleated galaxies not present in Paul nucleated catalog: ', len(bv_ferguson)
print cat_ferguson_dwarf[gv_nuc][bv_ferguson]

# We crossmatch the Ferguson dwarf catalog agains Paul nucleated and non-nucleated

idx, sep2d, sep3d = cat_ferguson_coo.match_to_catalog_sky(cat_paul_coo)
gv_ferguson=np.where(sep2d.arcsecond <= 12.)[0]
bv_ferguson=np.where(sep2d.arcsecond > 12.)[0]
print '\nFerguson dwarf galaxies already in Paul dwarf catalog: ', len(gv_ferguson)
print '\nFerguson dwarf galaxies not present in Paul dwarf catalog: ', len(bv_ferguson)
print '\nDo visual inspection of the following galaxies'
print cat_ferguson[bv_ferguson]

plt.figure()
plt.hist(sep2d*3600., bins=30, range=[0.,100.])
plt.xlabel('Separation angle (arcsec)')
plt.ylabel('Number')
plt.savefig('Crossmatch FCC original dwarf galaxies with Paul dwarf catalog.pdf')

# Now let's create the Tables

cat_fornax_nuc=Table()
cat_fornax_non=Table()

# Now the nucleated

coo=SkyCoord( cat_paul_nuc['RA'], cat_paul_nuc['DEC'], unit="deg")
cat_fornax_nuc['ID'] = ['NGFS'+ coo[j].ra.to_string(unit=u.hour, sep='',pad=1,precision=0, alwayssign=0, fields=3)+ coo[j].dec.to_string(unit=u.degree, sep='',pad=1,precision=0, alwayssign=1, fields=3) for j in np.arange(len(coo))]
cat_fornax_nuc['RA'] = np.asarray(coo.ra.deg)
cat_fornax_nuc['DEC'] = np.asarray(coo.dec.deg)
cat_fornax_nuc['RA_hms'] = np.asarray( coo.ra.to_string(unit=u.hour, sep=':', pad=1,precision=2) )
cat_fornax_nuc['DEC_dms'] = np.asarray( coo.dec.to_string(unit=u.deg, sep=':', pad=1,precision=2) )

coo=SkyCoord( cat_fornax_nuc['RA'], cat_fornax_nuc['DEC'], unit="deg")
idx, sep2d, sep3d = coo.match_to_catalog_sky(cat_ferguson_coo)
gv_fornax=np.where(sep2d.arcsecond < 12.)[0]
gv_ferguson=idx[np.where(sep2d.arcsecond < 12.)[0]]
cat_fornax_nuc['ID'][gv_fornax]=cat_ferguson['ID'][gv_ferguson]

for j in np.arange(len(ferguson_revised)):
	if np.sum(cat_fornax_nuc['ID']==ferguson_revised['ID'][j])==0 and ferguson_revised['MType'][j].startswith('d') and (ferguson_revised['MType'][j].endswith('N') or ferguson_revised['MType'][j].endswith('N?')):
		gv=np.where(cat_ferguson['ID'] == ferguson_revised['ID'][j])
		coo = ([cat_ferguson['RA'][gv][0],cat_ferguson['DEC'][gv][0]] if (ferguson_revised['RA'][j]==0 and ferguson_revised['DEC'][j]==0) else [ferguson_revised['RA'][j], ferguson_revised['DEC'][j]])
		cat_fornax_nuc.add_row( (cat_ferguson['ID'][gv][0], coo[0], coo[1], Angle(coo[0]*u.deg).to_string(unit=u.hour, sep=':', pad=1,precision=2), Angle(coo[1]*u.deg).to_string(unit=u.deg, sep=':', pad=1,precision=2)) )

for j in np.arange(len(mieske_revised)):
	if mieske_revised['MType'][j].startswith('d') and (mieske_revised['MType'][j].endswith('N') or mieske_revised['MType'][j].endswith('N?')):
		gv=np.where(cat_mieske['ID'] == mieske_revised['ID'][j])
		cat_fornax_nuc.add_row( (cat_mieske['ID'][gv][0], cat_mieske['RA'][gv][0], cat_mieske['DEC'][gv][0], Angle(cat_mieske['RA'][gv][0]*u.deg).to_string(unit=u.hour, sep=':', pad=1,precision=2), Angle(cat_mieske['DEC'][gv][0]*u.deg).to_string(unit=u.deg, sep=':', pad=1,precision=2)) )

coo=SkyCoord( cat_fornax_nuc['RA'], cat_fornax_nuc['DEC'], unit="deg")
cat_fornax_nuc['ID'] = ['NGFS'+ coo[j].ra.to_string(unit=u.hour, sep='',pad=1,precision=0, alwayssign=0, fields=3)+ coo[j].dec.to_string(unit=u.degree, sep='',pad=1,precision=0, alwayssign=1, fields=3) for j in np.arange(len(coo))]

cat_fornax_nuc.add_column( Column( np.full(len(cat_fornax_nuc),np.nan), name='m_i' ) )
cat_fornax_nuc.add_column( Column( np.full(len(cat_fornax_nuc),np.nan), name='M_i' ) )
cat_fornax_nuc.add_column( Column( np.full(len(cat_fornax_nuc),np.nan), name='sersic_n' ) )
cat_fornax_nuc.add_column( Column( np.full(len(cat_fornax_nuc),np.nan), name='reff_arcsec' ) )
cat_fornax_nuc.add_column( Column( np.full(len(cat_fornax_nuc),np.nan), name='reff_kpc' ) )
cat_fornax_nuc.add_column( Column( ['             ' for j in np.arange(len(cat_fornax_nuc)) ], name='Reference' ) )
cat_fornax_nuc.add_column( Column( ['             ' for j in np.arange(len(cat_fornax_nuc)) ], name='MType' ) )

cat_fornax_nuc_coo=SkyCoord( cat_fornax_nuc['RA'], cat_fornax_nuc['DEC'], unit="deg")
idx, sep2d, sep3d = cat_fornax_nuc_coo.match_to_catalog_sky(cat_paul_galfit_nuc_coo)
gv_fornax=np.where(sep2d.arcsecond <= 12.)[0]
gv_paul=idx[np.where(sep2d.arcsecond <= 12.)[0]]

cat_fornax_nuc['m_i'][gv_fornax] = cat_paul_galfit_nuc['mag'][gv_paul]
cat_fornax_nuc['M_i'][gv_fornax] = cat_paul_galfit_nuc['mag'][gv_paul] - 31.51
cat_fornax_nuc['sersic_n'][gv_fornax] = cat_paul_galfit_nuc['sersic_n'][gv_paul]
cat_fornax_nuc['reff_arcsec'][gv_fornax] = cat_paul_galfit_nuc['reff'][gv_paul]*decam_scale
cat_fornax_nuc['reff_kpc'][gv_fornax] = cat_paul_galfit_nuc['reff'][gv_paul]* fornax_distance.to('kpc') * decam_scale.to('rad')

idx, sep2d, sep3d = cat_fornax_nuc_coo.match_to_catalog_sky(cat_ferguson_coo)
gv_fornax=np.where(sep2d.arcsecond <= 12.)[0]
gv_ferguson=idx[np.where(sep2d.arcsecond <= 12.)[0]]
cat_fornax_nuc['Reference'][gv_fornax] = cat_ferguson['ID'][gv_ferguson]
cat_fornax_nuc['MType'][gv_fornax] = cat_ferguson['MType'][gv_ferguson]

idx, sep2d, sep3d = cat_fornax_nuc_coo.match_to_catalog_sky(cat_mieske_memb_coo)
gv_fornax=np.where(sep2d.arcsecond <= 12.)[0]
gv_mieske=idx[np.where(sep2d.arcsecond <= 12.)[0]]
for j in np.arange(len(gv_fornax)):
	print not bool(cat_fornax_nuc['Reference'][gv_fornax][j].strip())
	if not bool(cat_fornax_nuc['Reference'][gv_fornax][j].strip()):
		cat_fornax_nuc['Reference'][gv_fornax[j]] = cat_mieske_memb['ID'][gv_mieske[j]]

plt.figure()
plt.hist(sep2d*3600., bins=30, range=[0.,100.])
plt.xlabel('Separation angle (arcsec)')
plt.ylabel('Number')
plt.savefig('Crossmatch Fornax nucleated against Mieske 2007 dwarf galaxies.pdf')


# Now the non-nucleated
coo=SkyCoord( cat_paul_non['RA'], cat_paul_non['DEC'], unit="deg")
cat_fornax_non['ID'] = ['NGFS'+ coo[j].ra.to_string(unit=u.hour, sep='',pad=1,precision=0, alwayssign=0, fields=3)+ coo[j].dec.to_string(unit=u.degree, sep='',pad=1,precision=0, alwayssign=1, fields=3) for j in np.arange(len(coo))]
cat_fornax_non['RA'] = np.asarray(coo.ra.deg)
cat_fornax_non['DEC'] = np.asarray(coo.dec.deg)
cat_fornax_non['RA_hms'] = np.asarray( coo.ra.to_string(unit=u.hour, sep=':', pad=1,precision=2) )
cat_fornax_non['DEC_dms'] = np.asarray( coo.dec.to_string(unit=u.deg, sep=':', pad=1,precision=2) )

coo=SkyCoord( cat_fornax_non['RA'], cat_fornax_non['DEC'], unit="deg")
idx, sep2d, sep3d = coo.match_to_catalog_sky(cat_ferguson_coo)
gv_fornax=np.where(sep2d.arcsecond < 12.)[0]
gv_ferguson=idx[np.where(sep2d.arcsecond < 12.)[0]]
cat_fornax_non['ID'][gv_fornax]=cat_ferguson['ID'][gv_ferguson]

for j in np.arange(len(ferguson_revised)):
	if np.sum(cat_fornax_non['ID']==ferguson_revised['ID'][j])==0 and ferguson_revised['MType'][j].startswith('d') and not (ferguson_revised['MType'][j].endswith('N') or ferguson_revised['MType'][j].endswith('N?')):
		gv=np.where(cat_ferguson['ID'] == ferguson_revised['ID'][j])
		coo = ([cat_ferguson['RA'][gv][0],cat_ferguson['DEC'][gv][0]] if (ferguson_revised['RA'][j]==0 and ferguson_revised['DEC'][j]==0) else [ferguson_revised['RA'][j], ferguson_revised['DEC'][j]])
		cat_fornax_non.add_row( (cat_ferguson['ID'][gv][0], coo[0], coo[1], Angle(coo[0]*u.deg).to_string(unit=u.hour, sep=':', pad=1,precision=2), Angle(coo[1]*u.deg).to_string(unit=u.deg, sep=':', pad=1,precision=2)) )

for j in np.arange(len(mieske_revised)):
	if mieske_revised['MType'][j].startswith('d') and not (mieske_revised['MType'][j].endswith('N') or mieske_revised['MType'][j].endswith('N?')):
		gv=np.where(cat_mieske['ID'] == mieske_revised['ID'][j])
		cat_fornax_non.add_row( (cat_mieske['ID'][gv][0], cat_mieske['RA'][gv][0], cat_mieske['DEC'][gv][0], Angle(cat_mieske['RA'][gv][0]*u.deg).to_string(unit=u.hour, sep=':', pad=1,precision=2), Angle(cat_mieske['DEC'][gv][0]*u.deg).to_string(unit=u.deg, sep=':', pad=1,precision=2)) )


coo=SkyCoord( cat_fornax_non['RA'], cat_fornax_non['DEC'], unit="deg")
cat_fornax_non['ID'] = ['NGFS'+ coo[j].ra.to_string(unit=u.hour, sep='',pad=1,precision=0, alwayssign=0, fields=3)+ coo[j].dec.to_string(unit=u.degree, sep='',pad=1,precision=0, alwayssign=1, fields=3) for j in np.arange(len(coo))]

cat_fornax_non.add_column( Column( np.full(len(cat_fornax_non),np.nan), name='m_i' ) )
cat_fornax_non.add_column( Column( np.full(len(cat_fornax_non),np.nan), name='M_i' ) )
cat_fornax_non.add_column( Column( np.full(len(cat_fornax_non),np.nan), name='sersic_n' ) )
cat_fornax_non.add_column( Column( np.full(len(cat_fornax_non),np.nan), name='reff_arcsec' ) )
cat_fornax_non.add_column( Column( np.full(len(cat_fornax_non),np.nan), name='reff_kpc' ) )
cat_fornax_non.add_column( Column( ['             ' for j in np.arange(len(cat_fornax_non)) ], name='Reference' ) )
cat_fornax_non.add_column( Column( ['             ' for j in np.arange(len(cat_fornax_non)) ], name='MType' ) )

cat_fornax_non_coo=SkyCoord( cat_fornax_non['RA'], cat_fornax_non['DEC'], unit="deg")
idx, sep2d, sep3d = cat_fornax_non_coo.match_to_catalog_sky(cat_paul_galfit_non_coo)
gv_fornax=np.where(sep2d.arcsecond <= 12.)[0]
gv_paul=idx[np.where(sep2d.arcsecond <= 12.)[0]]

cat_fornax_non['m_i'][gv_fornax] = cat_paul_galfit_non['mag'][gv_paul]
cat_fornax_non['M_i'][gv_fornax] = cat_paul_galfit_non['mag'][gv_paul] - 31.51
cat_fornax_non['sersic_n'][gv_fornax] = cat_paul_galfit_non['sersic_n'][gv_paul]
cat_fornax_non['reff_arcsec'][gv_fornax] = cat_paul_galfit_non['reff'][gv_paul]*decam_scale
cat_fornax_non['reff_kpc'][gv_fornax] = cat_paul_galfit_non['reff'][gv_paul]* fornax_distance.to('kpc') * decam_scale.to('rad')

idx, sep2d, sep3d = cat_fornax_non_coo.match_to_catalog_sky(cat_ferguson_coo)
gv_fornax=np.where(sep2d.arcsecond <= 15.)[0]
gv_ferguson=idx[np.where(sep2d.arcsecond <= 15.)[0]]
cat_fornax_non['Reference'][gv_fornax] = cat_ferguson['ID'][gv_ferguson]
cat_fornax_non['MType'][gv_fornax] = cat_ferguson['MType'][gv_ferguson]

idx, sep2d, sep3d = cat_fornax_non_coo.match_to_catalog_sky(cat_mieske_memb_coo)
gv_fornax=np.where(sep2d.arcsecond <= 15.)[0]
gv_mieske=idx[np.where(sep2d.arcsecond <= 15.)[0]]
for j in np.arange(len(gv_fornax)):
	print not bool(cat_fornax_non['Reference'][gv_fornax][j].strip())
	if not bool(cat_fornax_non['Reference'][gv_fornax][j].strip()):
		cat_fornax_non['Reference'][gv_fornax[j]] = cat_mieske_memb['ID'][gv_mieske[j]]

plt.figure()
plt.hist(sep2d*3600., bins=30, range=[0.,100.])
plt.xlabel('Separation angle (arcsec)')
plt.ylabel('Number')
plt.savefig('Crossmatch Fornax non-nucleated againts Mieske 2007 dwarf galaxies.pdf')

cat_fornax_nuc['RA'].format = '12.6f'
cat_fornax_nuc['DEC'].format = '12.6f'
cat_fornax_nuc['m_i'].format = '6.2f'
cat_fornax_nuc['M_i'].format = '6.2f'
cat_fornax_nuc['sersic_n'].format = '6.2f'
cat_fornax_nuc['reff_arcsec'].format = '8.3f'
cat_fornax_nuc['reff_kpc'].format = '8.3f'

cat_fornax_non['RA'].format = '12.6f'
cat_fornax_non['DEC'].format = '12.6f'
cat_fornax_non['m_i'].format = '6.2f'
cat_fornax_non['M_i'].format = '6.2f'
cat_fornax_non['sersic_n'].format = '6.2f'
cat_fornax_non['reff_arcsec'].format = '8.3f'
cat_fornax_non['reff_kpc'].format = '8.3f'

cat_fornax_nuc.sort('reff_arcsec')
cat_fornax_non.sort('reff_arcsec')

ascii.write(cat_fornax_nuc, '../catalogs/NGFS_FCC_cat_nucleated.dat', format='fixed_width')  
ascii.write(cat_fornax_non, '../catalogs/NGFS_FCC_cat_non_nucleated.dat', format='fixed_width')  

cat_fornax_nuc.add_column( Column( ['$\odot$' for j in np.arange(len(cat_fornax_nuc)) ], name='Type' ), index=len(cat_fornax_nuc.keys())-2 )
cat_fornax_non.add_column( Column( ['$\medcircle$' for j in np.arange(len(cat_fornax_non)) ], name='Type' ), index=len(cat_fornax_non.keys())-2 )
cat_fornax=join(cat_fornax_nuc, cat_fornax_non, join_type='outer')

cat_fornax.sort('RA')

ascii.write(cat_fornax, '../catalogs/NGFS_FCC_cat_dwarf_deluxe.tex', exclude_names=['RA','DEC','MType'], format='aastex', latexdict = {'tabletype': 'deluxetable*', 'units': {'ID':'', 'RA_hms': 'h:m:s', 'DEC_dms':'d:m:s', 'm_i':'mag','M_i':'mag','sersic_n':'','reff_arcsec':'arcsec','reff_kpc':'kpc','Type':''}}, caption='Dwarf galaxies in the core of the Fornax Cluster\label{tab:dwarftable}')  

sys.exit()

cat_fornax_coo=SkyCoord( cat_fornax['RA'], cat_fornax['DEC'], unit="deg")
cat_richtler_coo=SkyCoord( cat_richtler['RA'], cat_richtler['DEC'], unit="deg")
idx, sep2d, sep3d = cat_richtler_coo.match_to_catalog_sky(cat_fornax_coo)
gv_richtler=np.where(sep2d.arcsecond <= 12.)[0]
gv_fornax=idx[np.where(sep2d.arcsecond <= 12.)[0]]
bv_richtler=np.where(sep2d.arcsecond > 12.)[0]

print 'Number of Richtler 2015 private communication galaxies AND NGFS ', len(gv)
hstack( [cat_richtler['ID','RA','DEC'][gv_richtler], cat_fornax['ID','M_i','Reference'][gv_fornax]] )

print 'Number of Richtler 2015 private communication galaxies NOT NGFS ', len(bv)
cat_richtler['ID','RA','DEC'][bv_richtler]

# Now Mieske

cat_fornax_coo=SkyCoord( cat_fornax['RA'], cat_fornax['DEC'], unit="deg")
cat_mieske_coo=SkyCoord( cat_mieske['RA'], cat_mieske['DEC'], unit="deg")
idx, sep2d, sep3d = cat_richtler_coo.match_to_catalog_sky(cat_fornax_coo)
gv_richtler=np.where(sep2d.arcsecond <= 12.)[0]
gv_fornax=idx[np.where(sep2d.arcsecond <= 12.)[0]]
bv_richtler=np.where(sep2d.arcsecond > 12.)[0]

print 'Number of Richtler 2015 private communication galaxies AND NGFS ', len(gv)
hstack( [cat_richtler['ID','RA','DEC'][gv_richtler], cat_fornax['ID','M_i','Reference'][gv_fornax]] )

print 'Number of Richtler 2015 private communication galaxies NOT NGFS ', len(bv)
cat_richtler['ID','RA','DEC'][bv_richtler]



#idx, sep2d, sep3d = cat_mieske_memb_coo.match_to_catalog_sky(cat_fornax_coo)
#gv_mieske=np.where(sep2d.arcsecond <= 15.)[0]
#gv_fornax=idx[np.where(sep2d.arcsecond <= 15.)[0]]
#print 'Number of Mieske 2007 dwarfs detected in the NGFS survey ', len(gv_fornax)
#print tabulate(cat_fornax[gv_fornax], headers="keys")

# Now comapring against Mieske et al 2007

gv=[]
for j in np.arange(len(cat_fornax_nuc)):
 if cat_fornax_nuc['ID'][j].startswith('FCC'): gv.append(j)
cat_fcc_nuc=cat_fornax_nuc[gv]

gv=[]
for j in np.arange(len(cat_fornax_nuc)):
 if cat_fornax_nuc['ID'][j].startswith('NGFS'): gv.append(j)
cat_ngfs_nuc=cat_fornax_nuc[gv]

gv=[]
for j in np.arange(len(cat_fornax_non)):
 if cat_fornax_non['ID'][j].startswith('FCC'): gv.append(j)
cat_fcc_non=cat_fornax_non[gv]

gv=[]
for j in np.arange(len(cat_fornax_non)):
 if cat_fornax_non['ID'][j].startswith('NGFS'): gv.append(j)
cat_ngfs_non=cat_fornax_non[gv]

gv_mieske_fcc=[]
for j in np.arange(len(cat_mieske)):
	if cat_mieske['ID'][j].startswith('FCC'): gv_mieske_fcc.append(j)
cat_mieske_fcc=cat_mieske[gv_mieske_fcc]

cat_fcc_nuc = SkyCoord(cat_fcc_nuc['RA'], cat_fcc_nuc['DEC'], unit="deg")
cat_fcc_non = SkyCoord(cat_fcc_non['RA'], cat_fcc_non['DEC'], unit="deg")

idx, sep2d, sep3d = cat_mieske_coo.match_to_catalog_sky(cat_ferguson_coo)
gv_mieske=np.where(sep2d.arcsecond <= 15.)[0]
print 'Mieske dwarf galaxies already in FCC catalog: ', len(gv_mieske)

idx, sep2d, sep3d = cat_mieske_coo.match_to_catalog_sky(cat_fcc_non)
gv_mieske=np.where(sep2d.arcsecond <= 15.)[0]
print 'Mieske non-nucleated dwarf galaxies already in FCC catalog: ', len(gv_mieske)

idx, sep2d, sep3d = cat_mieske_coo.match_to_catalog_sky(cat_fcc_nuc)
gv_mieske=np.where(sep2d.arcsecond <= 15.)[0]
print 'Mieske nucleated dwarf galaxies already in FCC catalog: ', len(gv_mieske)


sys.exit()

plt.figure()
plt.hist(sep2d*3600., bins=30, range=[0.,100.])
plt.xlabel('Separation angle (arcsec)')
plt.ylabel('Number')
plt.savefig('Crossmatch Mieske 2007 with FCC original dwarf galaxies.pdf')

gv_mieske=np.where(sep2d.arcsecond <= 15.)[0]
bv_mieske=np.where(sep2d.arcsecond > 15.)[0]
print 'Mieske dwarf galaxies already in FCC catalog: ', len(gv_mieske)

ascii.write(cat_mieske[gv_mieske], '../catalogs/Mieske_2007_FCC.dat', format='fixed_width')

cat_mieske_nonfcc = cat_mieske[bv_mieske]
cat_mieske_nonfcc_coo = SkyCoord(cat_mieske_nonfcc['RA'], cat_mieske_nonfcc['DEC'], unit="deg")
cat_fornax_coo=SkyCoord( cat_fornax['RA'], cat_fornax['DEC'], unit="deg")

idx, sep2d, sep3d = cat_mieske_nonfcc_coo.match_to_catalog_sky(cat_fornax_coo)
plt.figure()
plt.hist(sep2d*3600., bins=30, range=[0.,100.])
plt.xlabel('Separation angle (arcsec)')
plt.ylabel('Number')
plt.savefig('Crossmatch Mieske 2007 non-FCC with Fornax master dwarf catalog.pdf')

cat_fornax_nuc_coo = SkyCoord(cat_fornax_nuc['RA'], cat_fornax_nuc['DEC'], unit="deg")
cat_fornax_non_coo = SkyCoord(cat_fornax_non['RA'], cat_fornax_non['DEC'], unit="deg")

gv_mieske=np.where(sep2d.arcsecond <= 15.)[0]
bv_mieske=np.where(sep2d.arcsecond > 15.)[0]


sys.exit()

idx, sep2d, sep3d = cat_paul_non_coo.match_to_catalog_sky(cat_ferguson_coo)
print '\nNon-nucleated FCC galaxies: ', len(np.where(sep2d.arcsecond <= 12.)[0])
print 'Non-nucleated new NGFS Paul galaxies: ', len(np.where(sep2d.arcsecond > 12.)[0])

plt.figure()
plt.hist(sep2d*3600., bins=30, range=[0.,100.])
plt.xlabel('Separation angle (arcsec)')
plt.ylabel('Number')
plt.savefig('Crossmatch non-nucleated NGFS_Paul_FCC_orig.pdf')

idx, sep2d, sep3d = cat_paul_nuc_coo.match_to_catalog_sky(cat_ferguson_coo)
print '\nNucleated FCC galaxies: ', len(np.where(sep2d.arcsecond <= 12.)[0])
print 'Nucleated new NGFS Paul galaxies: ', len(np.where(sep2d.arcsecond > 12.)[0])

plt.figure()
plt.hist(sep2d*3600., bins=30, range=[0.,100.])
plt.xlabel('Separation angle (arcsec)')
plt.ylabel('Number')
plt.savefig('Crossmatch nucleated NGFS_Paul_FCC_orig.pdf')

# Now crossmatching Paul and FCC galaxies

sys.exit()

cat_ngfs_non=ascii.read('../catalogs/NGFS_FCC_cat_non_nucleated.dat', Reader=ascii.FixedWidth)
cat_ngfs_nuc=ascii.read('../catalogs/NGFS_FCC_cat_nucleated.dat', Reader=ascii.FixedWidth)

cat_ferguson_coo = SkyCoord(cat_ferguson['RA'], cat_ferguson['DEC'], unit="deg")
cat_ngfs_non_coo = SkyCoord(cat_ngfs_non['RA'], cat_ngfs_non['DEC'], unit="deg")
cat_ngfs_nuc_coo = SkyCoord(cat_ngfs_nuc['RA'], cat_ngfs_nuc['DEC'], unit="deg")

idx, sep2d, sep3d = cat_ngfs_non_coo.match_to_catalog_sky(cat_ferguson_coo)
print '\nNon-nucleated FCC galaxies: ', len(np.where(sep2d.arcsecond <= 12.)[0])
print 'Non-nucleated new NGFS galaxies: ', len(np.where(sep2d.arcsecond > 12.)[0])

plt.figure()
plt.hist(sep2d*3600., bins=30, range=[0.,100.])
plt.xlabel('Separation angle (arcsec)')
plt.ylabel('Number')
plt.savefig('Crossmatch non-nucleated NGFS_FCC.pdf')

idx, sep2d, sep3d = cat_ngfs_nuc_coo.match_to_catalog_sky(cat_ferguson_coo)
print '\nNucleated FCC galaxies: ', len(np.where(sep2d.arcsecond <= 12.)[0])
print 'Nucleated new NGFS galaxies: ', len(np.where(sep2d.arcsecond > 12.)[0])
print np.max(sep2d < 20.)

plt.figure()
plt.hist(sep2d*3600., bins=30, range=[0.,100.])
plt.xlabel('Separation angle (arcsec)')
plt.ylabel('Number')
plt.savefig('Crossmatch nucleated NGFS_FCC.pdf')

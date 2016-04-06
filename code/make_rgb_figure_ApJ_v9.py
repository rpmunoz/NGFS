import sys,os
import pyfits
import astropy
from astropy.io import fits
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy.table import Table, Column
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import aplpy
#import skimage.io as io
import skimage
import skimage.io
from skimage import exposure
from skimage import transform
from PIL import Image
from PIL import ImageEnhance
from tabulate import tabulate
from astroquery.simbad import Simbad
from astroquery.ned import Ned
import astropy.units as u
import matplotlib.gridspec as gridspec
from pyavm import AVM
u.set_enabled_equivalencies(u.dimensionless_angles())

stack_dir='/Volumes/Q6/NGFS/DECam/stacks'
cat_dir='/Volumes/Q6/NGFS/DECam/catalogs'

stack_tile=['1'] #['1','2','3','4','5','6','7','10','13']
stack_tile_ref='1'
stack_filter=['i','g','u']
stack_filter_ref='i'
stack_version='2'
do_rgb_quality='low'

decam_scale=0.263*u.arcsec
fornax_distance=20.0*u.Mpc
cat_dir='../catalogs'

cat_ferguson=ascii.read(cat_dir+'/Ferguson_1989_cat_all.dat')
cat_mieske=ascii.read(cat_dir+'/Mieske_2007_cat_dwarf.dat')
cat_ngfs_nuc=ascii.read(cat_dir+'/NGFS_FCC_cat_nucleated.dat', Reader=ascii.FixedWidth)
cat_ngfs_non=ascii.read(cat_dir+'/NGFS_FCC_cat_non_nucleated.dat', Reader=ascii.FixedWidth)

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
ferguson_revised.add_row(('FCC141','dN', 53.735265, -35.190913))
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

dwarf_group_coo=SkyCoord('3h39m53s','-35d30m0s')
dwarf_group_radius=Angle(2.6*u.arcmin)
#dwarf_zoom_coo=SkyCoord('3h40m03s','-35d27m53s')
dwarf_zoom_coo=SkyCoord('3h39m42s','-35d28m06s')
dwarf_zoom_radius=Angle(25*u.arcsec)

ngc_offset = {'NGC1399': Angle([7,-3]*u.arcmin),
	'NGC1404': Angle([-7,-3]*u.arcmin),
	'NGC1427': Angle([6,-3]*u.arcmin),
	'NGC1380': Angle([7,-3]*u.arcmin),
	'NGC1374': Angle([6,2.5]*u.arcmin),
	'NGC1379': Angle([3,-3]*u.arcmin),
	'NGC1387': Angle([5,3]*u.arcmin),
	'NGC1386': Angle([6,-2]*u.arcmin),
	}

ngfs_offset = {'NGFS034003-352754': Angle([-60,-35]*u.arcsec),
	'NGFS034003-352920': Angle([0,25]*u.arcsec),
	'NGFS034002-352930': Angle([-60,-25]*u.arcsec),
	'NGFS033942-352806': Angle([60,-35]*u.arcsec),
	'NGFS033942-353155': Angle([60,-25]*u.arcsec),
	}

mieske_offset = {'WFLSB1-5': Angle([0,-35]*u.arcsec),
	'WFLSB1-3': Angle([0,-25]*u.arcsec),
	}

stack_im_file=['ss_fornax_tileX_i_long.003.fits','ss_fornax_tileX_g_long_ALIGNi.003.fits','ss_fornax_tileX_u_long_ALIGNi.003.fits']
stack_weight_file=['ss_fornax_tileX_i_long.003.WEIGHT.fits','ss_fornax_tileX_g_long_ALIGNi.003.WEIGHT.fits','ss_fornax_tileX_u_long_ALIGNi.003.WEIGHT.fits']
stack_mask_file=['ss_fornax_tileX_i_long.003.MASK.fits','ss_fornax_tileX_g_long_ALIGNi.003.MASK.fits','ss_fornax_tileX_u_long_ALIGNi.003.MASK.fits']


# We made visual inspection of Ferguson 1989 dwarf galaxies and we update the catalog
cat_ferguson.add_column( Column(['FCC'+str(j) for j in cat_ferguson['FCC']], name='ID'), index=0)
for j in np.arange(len(ferguson_revised)):
	gv=np.where(cat_ferguson['ID']==ferguson_revised['ID'][j])
	cat_ferguson['MType'][gv]=ferguson_revised['MType'][j]
	if not (ferguson_revised['RA'][j]==0. and ferguson_revised['DEC'][j]==0.):
		print 'Updating the coordinates of galaxy ', ferguson_revised['ID'][j]
		cat_ferguson['RA'][gv]=ferguson_revised['RA'][j]
		cat_ferguson['DEC'][gv]=ferguson_revised['DEC'][j]

# We made visual inspection of Mieske 2007 dwarf galaxies and we update the catalog
for j in np.arange(len(mieske_revised)):
	gv=np.where(cat_mieske['ID'] == mieske_revised['ID'][j])
	print 'Updating the membership of galaxy ', cat_mieske['ID'][gv]
	cat_mieske['Memb'][gv]=mieske_revised['Memb'][j]


stack_im_file_full = [[stack_dir+'/'+im_file.replace('tileX','tile'+im_tile) for im_file in stack_im_file] for im_tile in stack_tile]
stack_weight_file_full = [[stack_dir+'/'+im_file.replace('tileX','tile'+im_tile) for im_file in stack_weight_file] for im_tile in stack_tile]
stack_mask_file_full = [[stack_dir+'/'+im_file.replace('tileX','tile'+im_tile) for im_file in stack_mask_file] for im_tile in stack_tile]
stack_rgb_file_full=[stack_dir+'/ngfs_tile'+im_tile+'_rgb.fits' for im_tile in stack_tile]

#stack_cat_file=['Paul_tile1.dat']
#stack_cat_file_full=[cat_dir+'/'+cat_file for cat_file in stack_cat_file]

for i in range(len(stack_im_file_full)):

	stack_im_file=stack_im_file_full[i]
	stack_weight_file=stack_weight_file_full[i]
	stack_mask_file=stack_mask_file_full[i]
	stack_rgb_file=stack_rgb_file_full[i]
	stack_rgb_limit=np.zeros((3,2), dtype=np.float32)

	im_rgb_file=stack_rgb_file.replace('.fits','_asinh_v'+stack_version+'.jpg')

	if do_rgb_quality == 'medium':
		nim_scale_factor=10.
		nim_rgb_file='../figures/NGFS_tile1_rgb_asinh_v2_medium.png'
		im_rgb_labels_file='../figures/NGFS_tile1_rgb_labels_medium.pdf'
	elif do_rgb_quality == 'low':
		nim_scale_factor=30.
		nim_rgb_file='../figures/NGFS_tile1_rgb_asinh_v2_low.jpg'
		im_rgb_labels_file='../figures/NGFS_tile1_rgb_labels_low.pdf'

	im_dwarf_group_file='../figures/ngfs_tile1_rgb_asinh_v2_dwarf_group.jpg' #stack_rgb_file.replace('.fits','_asinh_v'+stack_version+'_dwarf_group.png')
	im_dwarf_zoom_file='../figures/ngfs_tile1_rgb_asinh_v2_dwarf_zoom.png' #stack_rgb_file.replace('.fits','_asinh_v'+stack_version+'_dwarf_zoom.png')

	if not os.path.exists(im_rgb_file):
	
		hist_nbin=200
		hist_percentile=[0.25,99.5] #[0.25,99.8] #[0.25,99.5]  #([0.25,99.5],[0.25,99.55],[0.22,99.8])
		
		hdulist = pyfits.open(stack_im_file[0])
		im_h=hdulist[0].header
		hdulist.close()
		
		nx = int(im_h['NAXIS1'])
		ny = int(im_h['NAXIS2'])
		im_data_cube = np.zeros((3, ny, nx), dtype=np.float32)
		
	#	for f in stack_im_file:
	#		if not os.path.exists(f):
	#			raise Exception("File does not exist : " + f)
		
		if not os.path.exists(stack_rgb_file):
		
			print '\nCreating rgb cube ', stack_rgb_file
			for j in range(len(stack_filter)):
		
				print '\nProcessing image ', stack_im_file[j], ' - filter ', stack_filter[j]
				im_file=stack_im_file[j]
				weight_file=stack_weight_file[j]
				mask_file=stack_mask_file[j]
			
				if os.path.exists(im_file):
					print '\nReading image file ', im_file
					hdulist = pyfits.open(im_file)
					im_data=hdulist[0].data
					im_h=hdulist[0].header
					hdulist.close()
			
					print 'Reading weight file ', weight_file
					hdulist = pyfits.open(weight_file)
					weight_data=hdulist[0].data
					hdulist.close()
			
					bv_weight = (weight_data==0.)
					im_data[bv_weight]=np.nan
			
					print 'Image size along X-axis ', im_h['NAXIS1']
					print 'Image size along Y-axis ', im_h['NAXIS2']
			
					ypad = ny - im_data.shape[0]
					xpad = nx - im_data.shape[1]
					pad_width = ((0, ypad), (0, xpad))
	
					im_data_pad = np.pad(im_data, pad_width, mode='constant', constant_values=[np.nan])
					im_data_cube[j,:,:]=im_data_pad
				else:
					print 'The image file does not exists ', im_file
		
					nim_data = np.empty((ny, nx), dtype=np.float32)
					nim_data[:] = np.nan
					try:
						gv = stack_tile.index(stack_tile_ref)
						im_file=stack_im_file_full[gv][j]
						weight_file=stack_weight_file_full[gv][j]
						mask_file=stack_mask_file_full[gv][j]
	
						print '\nReading image file ', im_file
						hdulist = pyfits.open(im_file)
						im_data=hdulist[0].data
						hdulist.close()
			
						print 'Reading weight file ', weight_file
						hdulist = pyfits.open(weight_file)
						weight_data=hdulist[0].data
						hdulist.close()
			
						print 'Reading mask file ', mask_file
						hdulist = pyfits.open(mask_file)
						mask_data=hdulist[0].data
						hdulist.close()
	
						gv = np.logical_and(weight_data>0,mask_data==0)
				
						if np.count_nonzero(gv)>1:
							sky_data=im_data[gv]
							hist_limit= np.abs(np.percentile(sky_data, 0.001))
							hist_xmin=-1.*hist_limit
							hist_xmax=1.*hist_limit
							hist, bin_edges = np.histogram( sky_data, bins=np.linspace(hist_xmin, hist_xmax, hist_nbin), range=(hist_xmin, hist_xmax))
							plt.bar(bin_edges[:-1], hist, width = bin_edges[1]-bin_edges[0], linewidth=0., log=True)
							plt.xlim(hist_xmin, hist_xmax)
							plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
							plt.xlabel('Flux')
							plt.ylabel('Number count')
							plt.savefig(stack_dir+'/ngfs_tile'+stack_tile_ref+'_histogram_sky.pdf')
	
							weight_file=stack_weight_file[0]
							print 'Reading weight file ', weight_file
							hdulist = pyfits.open(weight_file)
							weight_data=hdulist[0].data
							hdulist.close()
	
							print 'Random sample from the reference image'
							gv=(weight_data>0)
							nim_data[gv]=np.random.choice(sky_data, np.count_nonzero(gv))
							im_data_cube[j,:,:]=nim_data
	
					except ValueError:
						print 'Adding gaussian model'
						im_data_cube[j,:,:]= 6. * np.random.randn(ny, nx)
						im_data_cube[j,:,:][bv_weight]=np.nan
	
			
			print 'Saving the rgb cube file ', stack_rgb_file
			pyfits.writeto(stack_rgb_file, im_data_cube, header=im_h)
		
	
		if not os.path.exists(im_rgb_file):
	
			print "Reading the RGB fits cube file: ", stack_rgb_file		
			hdulist = pyfits.open(stack_rgb_file)
			im_data_cube=hdulist[0].data
			im_h_cube=hdulist[0].header
			hdulist.close()
		
			w, h = 2*matplotlib.figure.figaspect(1.2)
			fig = plt.figure(figsize=(w,h))
			for j in range(len(stack_filter)):
				print '\nProcessing filter ', stack_filter[j]
				im_data=im_data_cube[j,:,:]
				gv=np.isfinite(im_data)
	
				if np.count_nonzero(gv)>1:
	
					if os.path.exists(stack_im_file[j]):
						print 'Computing the rgb limit using the actual data'
	
						stack_rgb_limit[j,:]= np.percentile(im_data[gv], hist_percentile)  #[plot_xmin, plot_xmax]
						if stack_rgb_limit[j,0]==stack_rgb_limit[j,1]:
							print 'The data seems to contain a constant value. Filter '+stack_filter[j]
							continue
	
						hist_xmin=np.amin(im_data[gv])
						hist_xmax=np.amax(im_data[gv])
		
						print 'Percentiles for computing RGB min and max ', hist_percentile
						print 'Histogram min and max ', hist_xmin, hist_xmax
						print 'RGB image min and max ', stack_rgb_limit[j,:]
		
		#hist_min=np.nanmin(im_data)
		#hist_max=np.nanmax(im_data)
		#fig = plt.figure()
		#ax = fig.add_subplot(1,1,1,)
		#n, bins, patches = ax.hist( np.ravel(im_data), bins=2., range=(hist_min, hist_max), histtype='bar')
		
						hist, bin_edges = np.histogram( im_data, bins=np.linspace(hist_xmin, hist_xmax, hist_nbin), range=(hist_xmin, hist_xmax))
			
		#	ax = fig.add_subplot(2,3,1+2*i,)
						plt.subplot(3,2,1+2*j)
						plt.bar(bin_edges[:-1], hist, width = bin_edges[1]-bin_edges[0], linewidth=0., log=True)
						plt.xlim(hist_xmin, hist_xmax)
						plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
		#	ax.set_yscale('log')
						plt.title('Tile '+stack_tile[i]+' - '+stack_filter[j])
						plt.xlabel('Flux')
						plt.ylabel('Number count')
		
						if stack_filter[j]=='i':
							hist, bin_edges = np.histogram( im_data, bins=np.linspace(stack_rgb_limit[j,0], stack_rgb_limit[j,1], hist_nbin), range=stack_rgb_limit[j,:] )
							stack_rgb_ratio_ref= [ hist[0]*1./np.amax(hist), hist[-1]*1./np.amax(hist) ]
							gv_max=np.argmax(hist)
							stack_rgb_limit_ratio_ref= (bin_edges[gv_max]-stack_rgb_limit[j,0])*1./(stack_rgb_limit[j,1]-stack_rgb_limit[j,0])
							print 'stack_rgb_ratio_ref ', stack_rgb_ratio_ref
							print 'stack_rgb_limit_ratio_ref ', stack_rgb_limit_ratio_ref
						else:
		#		if stack_filter[i]=='u': stack_rgb_ratio=[ stack_rgb_ratio_ref[0]**0.25, stack_rgb_ratio_ref[1] ]
		#		else: stack_rgb_ratio=stack_rgb_ratio_ref
							stack_rgb_ratio=stack_rgb_ratio_ref
							print 'stack_rgb_ratio ', stack_rgb_ratio
		
							delta=stack_rgb_limit[j,1]-stack_rgb_limit[j,0]
							if stack_filter[j]=='u':
								stack_rgb_limit[j,0] -= delta/8.
								stack_rgb_limit[j,1] += 2*delta
							else:
								stack_rgb_limit[j,0] -= delta/8.
								stack_rgb_limit[j,1] += delta/4.
		
							print "Guess RGB image min and max ", stack_rgb_limit[j,:]
					
							hist, bin_edges = np.histogram( im_data, bins=np.linspace(stack_rgb_limit[j,0], stack_rgb_limit[j,1], 2*hist_nbin), range=stack_rgb_limit[j,:] )
							gv_max=np.argmax(hist)
		#		stack_rgb_limit[i,:]= [ bin_edges[ np.argmin( np.abs(hist[0:gv_max]*1./np.amax(hist) - stack_rgb_ratio[0]) )], bin_edges[ gv_max + np.argmin( np.abs(hist[gv_max:-1]*1./np.amax(hist) - stack_rgb_ratio[1]) ) ] ]
							gv_limit_max=gv_max + np.argmin( np.abs(hist[gv_max:-1]*1./np.amax(hist) - stack_rgb_ratio[1]) )
							gv_limit_min=np.argmin( np.abs( (bin_edges[gv_max]- bin_edges[0:gv_max])*1./(bin_edges[gv_limit_max]-bin_edges[0:gv_max]) - stack_rgb_limit_ratio_ref) )
							stack_rgb_limit[j,:]=bin_edges[[gv_limit_min,gv_limit_max]]
							print "Computed RGB limit ", stack_rgb_limit[j,:]
	#						print "New stack_rgb_ratio ", hist[gv_limit_min]*1./hist[gv_max], hist[gv_limit_max]*1./hist[gv_max]
	
					else:
						print 'Computing the rgb limit using the reference tile'
	
						gv_tile = stack_tile.index(stack_tile_ref)
						gv_filter = stack_filter.index(stack_filter_ref)
						hdulist = pyfits.open(stack_rgb_file_full[gv_tile])
						im_h=hdulist[0].header
						hdulist.close()
	
						stack_rgb_limit[j,:]=np.asarray(im_h['RGB_'+stack_filter[j]].split(',')).astype(np.float)
						delta=stack_rgb_limit[j,1]-stack_rgb_limit[j,0]
						stack_rgb_limit[j,0] -= delta/4.
						stack_rgb_limit[j,1] += delta/4.
	
						hist, bin_edges = np.histogram( im_data, bins=np.linspace(stack_rgb_limit[j,0], stack_rgb_limit[j,1], hist_nbin), range=stack_rgb_limit[j,:] )
						gv_max=np.argmax(hist)
	
						stack_rgb_limit[j,1]= np.asarray(im_h['RGB_'+stack_filter[j]].split(',')).astype(np.float)[1]/np.asarray(im_h['RGB_'+stack_filter_ref].split(',')).astype(np.float)[1]*stack_rgb_limit[gv_filter,1]
						stack_rgb_limit[j,0]= bin_edges[ np.argmin( np.abs( (bin_edges[gv_max]- bin_edges[0:gv_max])*1./(stack_rgb_limit[j,1]-bin_edges[0:gv_max]) - stack_rgb_limit_ratio_ref) ) ]
						print "Computed RGB limit ", stack_rgb_limit[j,:]
	#					stack_rgb_limit[j,:]=np.asarray(im_h['RGB_'+stack_filter[j]].split(',')).astype(np.float)
						
					hist, bin_edges = np.histogram( im_data, bins=np.linspace(stack_rgb_limit[j,0], stack_rgb_limit[j,1], hist_nbin), range=stack_rgb_limit[j,:] )
		
		#	ax = fig.add_subplot(2,3,1+2*i+1,)
					print 'Adding RGB limits to header'
					print stack_rgb_limit[j,:]
					im_h_cube['RGB_'+stack_filter[j]]= ','.join(np.char.mod('%.2f', stack_rgb_limit[j,:]))
	
					plt.subplot(3,2,1+2*j+1)
					plt.bar(bin_edges[:-1], hist, width = bin_edges[1]-bin_edges[0], linewidth=0., log=True)
					plt.xlim(stack_rgb_limit[j,0], stack_rgb_limit[j,1])
					plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
		#	ax.set_yscale('log')
					plt.title('Tile '+stack_tile[i]+' - '+stack_filter[j])
					plt.xlabel('Flux')
					plt.ylabel('Number count')
	
				else:
					print '\nNo valid pixels for filter '+stack_filter[j]
	
			print 'Updating the rgb cube file ', stack_rgb_file
			pyfits.writeto(stack_rgb_file, im_data_cube, header=im_h_cube, clobber=True)
		
			plt.tight_layout()
			plt.savefig(stack_dir+'/ngfs_tile'+stack_tile[i]+'_histogram_v'+stack_version+'.pdf')
			plt.close(fig)
		
			print "stack_rgb_limit: \n", stack_rgb_limit
			
			im_rgb_file=stack_rgb_file.replace('.fits','_asinh_v'+stack_version+'.jpg')
			print 'Creating RGB image file ', im_rgb_file
			aplpy.make_rgb_image(stack_rgb_file, im_rgb_file, stretch_r='arcsinh', stretch_g='arcsinh', stretch_b='arcsinh', vmin_r=stack_rgb_limit[0,0], vmin_g=stack_rgb_limit[1,0], vmin_b=stack_rgb_limit[2,0], vmax_r=stack_rgb_limit[0,1], vmax_g=stack_rgb_limit[1,1], vmax_b=stack_rgb_limit[2,1], vmid_r=-0.07, vmid_g=-0.07, vmid_b=-0.07, make_nans_transparent=True, embed_avm_tags=True)


	if not os.path.exists(nim_rgb_file):
		print 'Reading full resolution image and downsampling it for Aplpy'
		im_data = skimage.io.imread(im_rgb_file)
		im_data = skimage.transform.rescale(im_data, 1/nim_scale_factor, order=1)
		skimage.io.imsave(nim_rgb_file, im_data)
		im_data=0

		hdu=fits.open(stack_im_file_full[i][0])
		hdu[0].header['NAXIS1'] = np.int(hdu[0].header['NAXIS1']/nim_scale_factor)
		hdu[0].header['NAXIS2'] = np.int(hdu[0].header['NAXIS2']/nim_scale_factor)
		hdu[0].header['CD1_1'] *= nim_scale_factor
		hdu[0].header['CD2_2'] *= nim_scale_factor
		hdu[0].header['CRPIX1'] = (hdu[0].header['NAXIS1']+1)/2.
		hdu[0].header['CRPIX2'] = (hdu[0].header['NAXIS2']+1)/2.
		w=astropy.wcs.WCS(hdu[0].header)
		w.naxis1=hdu[0].header['NAXIS1']
		w.naxis2=hdu[0].header['NAXIS2']

		avm = AVM.from_wcs(w)
		avm.Spatial.Equinox='J2000'
		avm.embed(nim_rgb_file, nim_rgb_file)

	if not (os.path.exists(im_dwarf_group_file) and os.path.exists(im_dwarf_zoom_file)) :
		print 'Creating dwarf group and dwarf zoom images'
		im_data = np.flipud(skimage.io.imread(im_rgb_file))

		hdu=fits.open(stack_im_file_full[i][0])
		w=astropy.wcs.WCS(hdu[0].header)
		im_crop_coo=w.wcs_world2pix([[dwarf_group_coo.ra.deg,(dwarf_group_coo.dec+dwarf_group_radius).deg],[dwarf_group_coo.ra.deg,(dwarf_group_coo.dec-dwarf_group_radius).deg]], 1)
		im_crop_size=(np.abs(im_crop_coo[0,1]-im_crop_coo[1,1])*np.asarray([1.2,1.])).astype(int)
		im_crop_coo=(w.wcs_world2pix([[dwarf_group_coo.ra.deg,dwarf_group_coo.dec.deg]], 1)[0]).astype(int)
		im_crop_data=im_data[im_crop_coo[1]-im_crop_size[1]/2:im_crop_coo[1]+im_crop_size[1]/2,im_crop_coo[0]-im_crop_size[0]/2:im_crop_coo[0]+im_crop_size[0]/2]
		skimage.io.imsave(im_dwarf_group_file, np.flipud(im_crop_data))

		hdu=fits.open(stack_im_file_full[i][0])
		hdu[0].header['NAXIS1'] = im_crop_data.shape[1]
		hdu[0].header['NAXIS2'] = im_crop_data.shape[0]
		hdu[0].header['CRPIX1'] -= (im_crop_coo[0]-im_crop_size[0]/2)
		hdu[0].header['CRPIX2'] -= (im_crop_coo[1]-im_crop_size[1]/2)
		w=astropy.wcs.WCS(hdu[0].header)
		w.naxis1=hdu[0].header['NAXIS1']
		w.naxis2=hdu[0].header['NAXIS2']
		avm = AVM.from_wcs(w)
		avm.Spatial.Equinox='J2000'
		avm.embed(im_dwarf_group_file, im_dwarf_group_file)

		hdu=fits.open(stack_im_file_full[i][0])
		w=astropy.wcs.WCS(hdu[0].header)
		im_crop_coo=w.wcs_world2pix([[dwarf_zoom_coo.ra.deg,(dwarf_zoom_coo.dec+dwarf_zoom_radius).deg],[dwarf_zoom_coo.ra.deg,(dwarf_zoom_coo.dec-dwarf_zoom_radius).deg]], 1)
		im_crop_size=(np.abs(im_crop_coo[0,1]-im_crop_coo[1,1])*np.asarray([1.2,1.])).astype(int)
		im_crop_coo=(w.wcs_world2pix([[dwarf_zoom_coo.ra.deg,dwarf_zoom_coo.dec.deg]], 1)[0]).astype(int)
		im_crop_data=im_data[im_crop_coo[1]-im_crop_size[1]/2:im_crop_coo[1]+im_crop_size[1]/2,im_crop_coo[0]-im_crop_size[0]/2:im_crop_coo[0]+im_crop_size[0]/2]
		skimage.io.imsave(im_dwarf_zoom_file, np.flipud(im_crop_data))

		hdu=fits.open(stack_im_file_full[i][0])
		hdu[0].header['NAXIS1'] = im_crop_data.shape[1]
		hdu[0].header['NAXIS2'] = im_crop_data.shape[0]
		hdu[0].header['CRPIX1'] -= (im_crop_coo[0]-im_crop_size[0]/2)
		hdu[0].header['CRPIX2'] -= (im_crop_coo[1]-im_crop_size[1]/2)
		w=astropy.wcs.WCS(hdu[0].header)
		w.naxis1=hdu[0].header['NAXIS1']
		w.naxis2=hdu[0].header['NAXIS2']
		avm = AVM.from_wcs(w)
		avm.Spatial.Equinox='J2000'
		avm.embed(im_dwarf_zoom_file, im_dwarf_zoom_file)

	if os.path.exists(nim_rgb_file):
		print 'Displaying and overplotting the position of dwarf galaxies'
		im_data = np.flipud(skimage.io.imread(nim_rgb_file))
		im_size= im_data.shape
		avm=AVM.from_image(nim_rgb_file)
		w = avm.to_wcs()
		w.naxis1=im_size[1]
		w.naxis2=im_size[0]

		fig = plt.figure(figsize=(14.,18.))
#		fig = plt.figure(figsize=(14.,16.5), dpi=300)
		gs=gridspec.GridSpec(3, 2) #, width_ratios=[1.2,1.2], height_ratios=[1,1,1])
#		gs=gridspec.GridSpec(8, 6)
		gs.update(wspace=0.2, hspace=0.2)
		ax1 = plt.subplot(gs[0:2,:])
		ax2 = plt.subplot(gs[2,0])
		ax3 = plt.subplot(gs[2,1])
#		ax1 = plt.subplot(gs[0:6,:])
#		ax2 = plt.subplot(gs[6:8,0:3])
#		ax3 = plt.subplot(gs[6:8,3:6])
		gs.tight_layout(fig, rect=[0.06, 0.01, 1., 1.])

		ax1.set_aspect('equal')
		ax2.axis('off')
		ax3.axis('off')
		ax1.axes.get_xaxis().set_ticks([])
		ax1.axes.get_yaxis().set_ticks([])
		ax2.axes.get_xaxis().set_ticks([])
		ax2.axes.get_yaxis().set_ticks([])
		ax3.axes.get_xaxis().set_ticks([])
		ax3.axes.get_yaxis().set_ticks([])


		f1 = aplpy.FITSFigure(w, figure=fig, subplot=list(ax1.get_position().bounds) ) # [0.11,0.07,0.58,0.92]) # wcs, hdu[0])
		f1.tick_labels.set_xformat('hh:mm')
		f1.tick_labels.set_yformat('dd:mm')
		f1.ticks.set_minor_frequency(6)
		f1.ticks.set_length(8) 
		f1.ticks.set_linewidth(1.5)   
		f1.show_rgb(nim_rgb_file)

		f1.show_rectangles([dwarf_group_coo.ra.deg], [dwarf_group_coo.dec.deg], 2*1.2*dwarf_group_radius.deg, 2*dwarf_group_radius.deg, edgecolor='orange', linewidth=0.5, alpha=0.5)

		# Now reading the original Ferguson catalog and generating the dwarf catalog

		gv_dwarf=[]
		for j in range(len(cat_ferguson)):
			if cat_ferguson['MType'][j].startswith('d'): gv_dwarf.append(j)

		cat_ferguson_dwarf=cat_ferguson[gv_dwarf]
		ascii.write(cat_ferguson_dwarf, cat_dir+'/ferguson_1989_cat_dwarfs.dat', format='commented_header')

		cat_coo = f1.world2pixel(np.asarray(cat_ferguson['RA'][gv_dwarf]), np.asarray(cat_ferguson['DEC'][gv_dwarf]))
		gv_dwarf_tile=[]
		for j in range(len(gv_dwarf)):
			if (0. < cat_coo[0][j] < im_size[1]) and (0. < cat_coo[1][j] < im_size[0]):
				temp=np.median(im_data[cat_coo[1][j]-2:cat_coo[1][j]+2,cat_coo[0][j]-2:cat_coo[0][j]+2,:])
				if temp > 0 : gv_dwarf_tile.append(j)

		print 'Number of Ferguson dwarfs inside the Tile ', len(gv_dwarf_tile)
		cat_ferguson_dwarf_tile=cat_ferguson[gv_dwarf][gv_dwarf_tile]
		ascii.write(cat_ferguson_dwarf_tile, cat_dir+'/ferguson_1989_cat_dwarfs_tile.dat', format='commented_header')	

		gv_nuc=[]
		gv_non=[]
		for j in range(len(cat_ferguson_dwarf)):
			if cat_ferguson_dwarf['MType'][j].endswith('N') or cat_ferguson_dwarf['MType'][j].endswith('N?'): gv_nuc.append(j)
			else : gv_non.append(j)

		ascii.write(cat_ferguson_dwarf[gv_nuc], cat_dir+'/ferguson_1989_cat_nucleated.dat', format='commented_header')
		ascii.write(cat_ferguson_dwarf[gv_non], cat_dir+'/ferguson_1989_cat_non_nucleated.dat', format='commented_header')

		gv_nuc=[]
		gv_non=[]
		for j in range(len(cat_ferguson_dwarf_tile)):
			if cat_ferguson_dwarf_tile['MType'][j].endswith('N') or cat_ferguson_dwarf_tile['MType'][j].endswith('N?'): gv_nuc.append(j)
			else : gv_non.append(j)

		ascii.write(cat_ferguson_dwarf_tile[gv_nuc], cat_dir+'/ferguson_1989_cat_nucleated_tile.dat', format='commented_header')
		ascii.write(cat_ferguson_dwarf_tile[gv_non], cat_dir+'/ferguson_1989_cat_non_nucleated_tile.dat', format='commented_header')

		# Now crossmatching Ferguson dwarfs with the NGFS catalogs
		
		gv_ngfs_nuc_ferguson=[]
		gv_ngfs_nuc_mieske=[]
		gv_ngfs_nuc_new=[]
		for j in range(len(cat_ngfs_nuc)):
			if not bool(cat_ngfs_nuc['Reference'][j]): gv_ngfs_nuc_new.append(j)
			elif cat_ngfs_nuc['Reference'][j].startswith('FCC'): gv_ngfs_nuc_ferguson.append(j)
			else: gv_ngfs_nuc_mieske.append(j)

		gv_ngfs_non_ferguson=[]
		gv_ngfs_non_mieske=[]
		gv_ngfs_non_new=[]
		for j in range(len(cat_ngfs_non)):
			if not bool(cat_ngfs_non['Reference'][j]): gv_ngfs_non_new.append(j)
			elif cat_ngfs_non['Reference'][j].startswith('FCC'): gv_ngfs_non_ferguson.append(j)
			else: gv_ngfs_non_mieske.append(j)


		f1.show_circles(np.asarray(cat_ngfs_non['RA'][gv_ngfs_non_ferguson]), np.asarray(cat_ngfs_non['DEC'][gv_ngfs_non_ferguson]), 40./3600, edgecolor='white', linewidth=1., alpha=0.4)
		f1.show_circles(np.asarray(cat_ngfs_non['RA'][gv_ngfs_non_mieske]), np.asarray(cat_ngfs_non['DEC'][gv_ngfs_non_mieske]), 60./3600, edgecolor='white', linewidth=1., alpha=0.4)
		f1.show_circles(np.asarray(cat_ngfs_non['RA'][gv_ngfs_non_new]), np.asarray(cat_ngfs_non['DEC'][gv_ngfs_non_new]), 60./3600, edgecolor='red', linewidth=1., alpha=0.6)

		f1.show_circles(np.asarray(cat_ngfs_nuc['RA'][gv_ngfs_nuc_ferguson]), np.asarray(cat_ngfs_nuc['DEC'][gv_ngfs_nuc_ferguson]), 40./3600, edgecolor='white', linewidth=1., alpha=0.4)
		f1.show_circles(np.asarray(cat_ngfs_nuc['RA'][gv_ngfs_nuc_ferguson]), np.asarray(cat_ngfs_nuc['DEC'][gv_ngfs_nuc_ferguson]), 15./3600, edgecolor='none', facecolor='white', alpha=0.4)
		f1.show_circles(np.asarray(cat_ngfs_nuc['RA'][gv_ngfs_nuc_mieske]), np.asarray(cat_ngfs_nuc['DEC'][gv_ngfs_nuc_mieske]), 60./3600, edgecolor='white', linewidth=1., alpha=0.4)
		f1.show_circles(np.asarray(cat_ngfs_nuc['RA'][gv_ngfs_nuc_mieske]), np.asarray(cat_ngfs_nuc['DEC'][gv_ngfs_nuc_mieske]), 20./3600, edgecolor='none', facecolor='white', alpha=0.4)
		f1.show_circles(np.asarray(cat_ngfs_nuc['RA'][gv_ngfs_nuc_new]), np.asarray(cat_ngfs_nuc['DEC'][gv_ngfs_nuc_new]), 60./3600, edgecolor='red', linewidth=1., alpha=0.6)
		f1.show_circles(np.asarray(cat_ngfs_nuc['RA'][gv_ngfs_nuc_new]), np.asarray(cat_ngfs_nuc['DEC'][gv_ngfs_nuc_new]), 20./3600, edgecolor='none', facecolor='red', alpha=0.4)

		customSimbad = Simbad()
		customSimbad.add_votable_fields('ra(d)','dec(d)','flux(V)','flux(K)')
		result = customSimbad.query_criteria('region(circle, NGC 1399, 1.2d) & cat="ngc"')
		result.sort('FLUX_K')
		for j in range(8):
			f1.add_label(result['RA_d'][j]+ngc_offset[result['MAIN_ID'][j].replace(' ','')][0].deg, result['DEC_d'][j]+ngc_offset[result['MAIN_ID'][j].replace(' ','')][1].deg,  result['MAIN_ID'][j].replace(' ',''), size=7, color='silver', horizontalalignment='center', alpha=0.5, family='sans-serif', style='italic')

		label='NGC1427-A'
		label_offset=Angle([6,2.5]*u.arcmin)
		result = customSimbad.query_object(label)
		f1.add_label(result['RA_d'][0]+label_offset[0].deg, result['DEC_d'][0]+label_offset[1].deg, label, size=7, color='silver', horizontalalignment='center', alpha=0.5, family='sans-serif', style='italic')

		f1.add_scalebar((100*u.kpc/fornax_distance.to('kpc')*u.rad).to('arcmin'), '', color='white')
		f1.scalebar.set_corner('bottom left')
		f1.scalebar.set_label('')
#		f1.scalebar.set_alpha(0.9)
#		f1.scalebar.set_font(size='small', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')

		label_coo=SkyCoord('3h43m5s','-36d28m0s')
		f1.add_label( label_coo.ra.deg, label_coo.dec.deg, '100 kpc', size=12, color='white', horizontalalignment='center', alpha=0.9, family='sans-serif', style='normal')

		im_data = skimage.io.imread(im_dwarf_group_file)
		im_size = im_data.shape
		avm=AVM.from_image(im_dwarf_group_file)
		w = avm.to_wcs()
		w.naxis1=im_size[1]
		w.naxis2=im_size[0]
		f2 = aplpy.FITSFigure(w, figure=fig, subplot=list(ax2.get_position().bounds) ) # [0.11,0.07,0.58,0.92]) # wcs, hdu[0])
		f2.show_rgb(im_dwarf_group_file)
		f2.frame.set_linewidth(1)
		f2.frame.set_color('white')
		f2.axis_labels.hide()
		f2.tick_labels.hide()
		f2.ticks.hide()

		f2.show_rectangles([dwarf_zoom_coo.ra.deg], [dwarf_zoom_coo.dec.deg], 2*1.2*dwarf_zoom_radius.deg, 2*dwarf_zoom_radius.deg, edgecolor='orange', linewidth=1, alpha=0.5)

		f2.show_circles(np.asarray(cat_ngfs_non['RA'][gv_ngfs_non_ferguson]), np.asarray(cat_ngfs_non['DEC'][gv_ngfs_non_ferguson]), 12./3600, edgecolor='white', linewidth=2., alpha=0.4)
		f2.show_circles(np.asarray(cat_ngfs_non['RA'][gv_ngfs_non_mieske]), np.asarray(cat_ngfs_non['DEC'][gv_ngfs_non_mieske]), 15./3600, edgecolor='white', linewidth=1.8, alpha=0.4)
		f2.show_circles(np.asarray(cat_ngfs_non['RA'][gv_ngfs_non_new]), np.asarray(cat_ngfs_non['DEC'][gv_ngfs_non_new]), 15./3600, edgecolor='red', linewidth=1.8, alpha=0.4)

		f2.show_circles(np.asarray(cat_ngfs_nuc['RA'][gv_ngfs_nuc_ferguson]), np.asarray(cat_ngfs_nuc['DEC'][gv_ngfs_nuc_ferguson]), 12./3600, edgecolor='white', linewidth=2., alpha=0.4)
		f2.show_circles(np.asarray(cat_ngfs_nuc['RA'][gv_ngfs_nuc_ferguson]), np.asarray(cat_ngfs_nuc['DEC'][gv_ngfs_nuc_ferguson]), 5./3600, edgecolor='none', facecolor='white', alpha=0.08)
		f2.show_circles(np.asarray(cat_ngfs_nuc['RA'][gv_ngfs_nuc_mieske]), np.asarray(cat_ngfs_nuc['DEC'][gv_ngfs_nuc_mieske]), 15./3600, edgecolor='white', linewidth=1.8, alpha=0.4)
		f2.show_circles(np.asarray(cat_ngfs_nuc['RA'][gv_ngfs_nuc_mieske]), np.asarray(cat_ngfs_nuc['DEC'][gv_ngfs_nuc_mieske]), 5./3600, edgecolor='none', facecolor='white', alpha=0.08)
		f2.show_circles(np.asarray(cat_ngfs_nuc['RA'][gv_ngfs_nuc_new]), np.asarray(cat_ngfs_nuc['DEC'][gv_ngfs_nuc_new]), 15./3600, edgecolor='red', linewidth=1.8, alpha=0.4)
		f2.show_circles(np.asarray(cat_ngfs_nuc['RA'][gv_ngfs_nuc_new]), np.asarray(cat_ngfs_nuc['DEC'][gv_ngfs_nuc_new]), 5./3600, edgecolor='none', facecolor='red', alpha=0.08)

		print 'FCC inside dwarf region box'
		label=['FCC'+str(x) for x in cat_ferguson['FCC']]
		coo=SkyCoord( cat_ferguson['RA'], cat_ferguson['DEC'], unit="deg")
		coo_x,coo_y =f2.world2pixel(coo.ra.value, coo.dec.value)
		for j in range(len(cat_ferguson)):
			if (0. < coo_x[j] < im_size[1]) and (0. < coo_y[j] < im_size[0]):
				print coo_x[j], coo_y[j]
				f2.add_label(coo[j].ra.value, coo[j].dec.value-30./3600, label[j], size=7, color='silver', horizontalalignment='center', alpha=0.5, family='sans-serif', style='italic')

		print 'Mieske inside dwarf region box'
		coo=SkyCoord( cat_ngfs_non['RA'][gv_ngfs_non_mieske], cat_ngfs_non['DEC'][gv_ngfs_non_mieske], unit="deg")
		coo_x,coo_y =f2.world2pixel(coo.ra.value, coo.dec.value)
		for j in range(len(coo)):
			if (0. < coo_x[j] < im_size[1]) and (0. < coo_y[j] < im_size[0]):
				print coo_x[j], coo_y[j]
				label=cat_ngfs_non['Reference'][gv_ngfs_non_mieske[j]]
				f2.add_label(coo[j].ra.value+mieske_offset[label][0].deg, coo[j].dec.value+mieske_offset[label][1].deg, label, size=7, color='silver', horizontalalignment='center', alpha=0.5, family='sans-serif', style='italic')

		print 'NGFS non-nucleated inside dwarf region box'
		coo=SkyCoord( cat_ngfs_non['RA'][gv_ngfs_non_new], cat_ngfs_non['DEC'][gv_ngfs_non_new], unit="deg")
		coo_x,coo_y =f2.world2pixel(coo.ra.value, coo.dec.value)
		for j in range(len(gv_ngfs_non_new)):
			if (0. < coo_x[j] < im_size[1]) and (0. < coo_y[j] < im_size[0]):
				print coo_x[j], coo_y[j]
				label='NGFS'+ coo[j].ra.to_string(unit=u.hour, sep='',pad=1,precision=0, alwayssign=0, fields=3)+ coo[j].dec.to_string(unit=u.degree, sep='',pad=1,precision=0, alwayssign=1, fields=3)
				f2.add_label(coo[j].ra.value+ngfs_offset[label][0].deg, coo[j].dec.value+ngfs_offset[label][1].deg, label, size=7, color='silver', horizontalalignment='center', alpha=0.5, family='sans-serif', style='italic')

		length=np.abs(dwarf_group_radius.to('deg')*2*1.2/np.cos(dwarf_group_coo.dec.rad))  #np.abs(im_size[1]*w.wcs.cdelt[0]/np.cos(dwarf_group_coo.dec.rad))*u.deg
		length_phy=fornax_distance.to('kpc')*length.to(u.rad).to('')

		f2.show_arrows([dwarf_group_coo.ra.deg+length.value/20.], [dwarf_group_coo.dec.deg-dwarf_group_radius.deg*1.05], [length.value/2.*9/10.], [0], edgecolor='none', facecolor='black', length_includes_head=1, width=7, head_width=30, head_length=30)
		f2.show_arrows([dwarf_group_coo.ra.deg-length.value/20.], [dwarf_group_coo.dec.deg-dwarf_group_radius.deg*1.05], [-1*length.value/2.*9/10.*0.99], [0], edgecolor='none', facecolor='black', length_includes_head=1, width=7, head_width=30, head_length=30)
		f2.add_label( dwarf_group_coo.ra.deg, dwarf_group_coo.dec.deg-dwarf_group_radius.deg*1.05, "{:.1f} '".format(length.value*60.), size=9, color='black', horizontalalignment='center', alpha=1.)

		f2.show_arrows([dwarf_group_coo.ra.deg+length.value*2/20.], [dwarf_group_coo.dec.deg+dwarf_group_radius.deg*1.05], [length.value/2.*8/10.*0.99], [0], edgecolor='none', facecolor='black', length_includes_head=1, width=7, head_width=30, head_length=30)
		f2.show_arrows([dwarf_group_coo.ra.deg-length.value*2/20.], [dwarf_group_coo.dec.deg+dwarf_group_radius.deg*1.05], [-1*length.value/2.*8/10.], [0], edgecolor='none', facecolor='black', length_includes_head=1, width=7, head_width=30, head_length=30)
		f2.add_label( dwarf_group_coo.ra.deg, dwarf_group_coo.dec.deg+dwarf_group_radius.deg*1.05, "{:.1f} kpc".format(length_phy.value), size=9, color='black', horizontalalignment='center', alpha=1.)

		f2.add_scalebar((5*u.kpc/fornax_distance.to('kpc')*u.rad).to('arcmin'), '', color='white')
		f2.scalebar.set_corner('bottom left')
		f2.scalebar.set_label('')
#		f2.scalebar.set_alpha(0.7)
#		f2.scalebar.set_font(size='small', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')

		f2.add_label( dwarf_group_coo.ra.deg+dwarf_group_radius.deg*1.18, dwarf_group_coo.dec.deg-dwarf_group_radius.deg*0.85, '5 kpc', size=12, color='white', horizontalalignment='center', alpha=0.9, family='sans-serif', style='normal')

		im_data = skimage.io.imread(im_dwarf_zoom_file)
		im_size=im_data.shape
		avm=AVM.from_image(im_dwarf_zoom_file)
		w = avm.to_wcs()
		w.naxis1=im_size[1]
		w.naxis2=im_size[0]
		f3 = aplpy.FITSFigure(w, figure=fig, subplot=list(ax3.get_position().bounds))
		f3.show_rgb(im_dwarf_zoom_file)
		f3.frame.set_linewidth(1)
		f3.frame.set_color('white')
		f3.axis_labels.hide()
		f3.tick_labels.hide()
		f3.ticks.hide()

		print 'NGFS non-nucleated inside dwarf zoom box'
		coo=SkyCoord( cat_ngfs_non['RA'][gv_ngfs_non_new], cat_ngfs_non['DEC'][gv_ngfs_non_new], unit="deg")
		coo_x,coo_y =f3.world2pixel(coo.ra.value, coo.dec.value)
		for j in range(len(gv_ngfs_non_new)):
			if (0. < coo_x[j] < im_size[1]) and (0. < coo_y[j] < im_size[0]):
				print coo_x[j], coo_y[j]
				label='NGFS'+ coo[j].ra.to_string(unit=u.hour, sep='',pad=1,precision=0, alwayssign=0, fields=3)+ coo[j].dec.to_string(unit=u.degree, sep='',pad=1,precision=0, alwayssign=1, fields=3)
				f3.add_label(coo[j].ra.value-20./3600, coo[j].dec.value-23./3600, label, size=9, color='silver', horizontalalignment='center', alpha=0.5, family='sans-serif', style='italic')

		length=np.abs(dwarf_zoom_radius.to('deg')*2*1.2/np.cos(dwarf_zoom_coo.dec.rad))  #np.abs(im_size[1]*w.wcs.cdelt[0]/np.cos(dwarf_zoom_coo.dec.rad))*u.deg
		length_phy=fornax_distance.to('kpc')*length.to(u.rad).to('')

		f3.show_arrows([dwarf_zoom_coo.ra.deg+length.value/20.], [dwarf_zoom_coo.dec.deg-dwarf_zoom_radius.deg*1.05], [length.value/2.*9/10.], [0], edgecolor='none', facecolor='black', length_includes_head=1, width=1, head_width=5, head_length=5)
		f3.show_arrows([dwarf_zoom_coo.ra.deg-length.value/20.], [dwarf_zoom_coo.dec.deg-dwarf_zoom_radius.deg*1.05], [-1*length.value/2.*9/10.*0.99], [0], edgecolor='none', facecolor='black', length_includes_head=1, width=1, head_width=5, head_length=5)
		f3.add_label( dwarf_zoom_coo.ra.deg, dwarf_zoom_coo.dec.deg-dwarf_zoom_radius.deg*1.05, "{:.1f} '".format(length.to(u.arcmin).value), size=9, color='black', horizontalalignment='center', alpha=1.)

		f3.show_arrows([dwarf_zoom_coo.ra.deg+length.value*2/20.], [dwarf_zoom_coo.dec.deg+dwarf_zoom_radius.deg*1.05], [length.value/2.*8/10.*0.99], [0], edgecolor='none', facecolor='black', length_includes_head=1, width=1, head_width=5, head_length=5)
		f3.show_arrows([dwarf_zoom_coo.ra.deg-length.value*2/20.], [dwarf_zoom_coo.dec.deg+dwarf_zoom_radius.deg*1.05], [-1*length.value/2.*8/10.], [0], edgecolor='none', facecolor='black', length_includes_head=1, width=1, head_width=5, head_length=5)
		f3.add_label( dwarf_zoom_coo.ra.deg, dwarf_zoom_coo.dec.deg+dwarf_zoom_radius.deg*1.05, "{:.1f} kpc".format(length_phy.value), size=9, color='black', horizontalalignment='center', alpha=1.)

		f3.add_scalebar((1*u.kpc/fornax_distance.to('kpc')*u.rad).to('arcmin'), '', color='white')
		f3.scalebar.set_corner('bottom left')
		f3.scalebar.set_label('')
#		f3.scalebar.set_alpha(0.7)
#		f3.scalebar.set_font(size='small', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')

		f3.add_label( dwarf_zoom_coo.ra.deg+dwarf_zoom_radius.deg*1.15, dwarf_zoom_coo.dec.deg-dwarf_zoom_radius.deg*0.85, '1 kpc', size=12, color='white', horizontalalignment='center', alpha=0.9, family='sans-serif', style='normal')

		fig.savefig(im_rgb_labels_file, dpi=70)

#		cat_nucleated.add_column( Table.Column(name='Number', data=np.arange(1,len(cat_coo)+1)), 0)
#		cat_nucleated.add_column( Table.Column(name='Separation', unit='arcseconds', data=sep2d.arcsecond) )
#		print tabulate(cat_nucleated, headers="keys")

		sys.exit()

sys.exit()

aplpy.make_rgb_image(stack_rgb_file, stack_rgb_file.replace('.fits','_v1.jpg'), stretch_r='arcsinh', stretch_g='arcsinh', stretch_b='arcsinh', pmin_r=0.25, pmin_g=0.25, pmin_b=0.25, pmax_r=99.75, pmax_g=99.75, pmax_b=99.75, make_nans_transparent=True, embed_avm_tags=True)
aplpy.make_rgb_image(stack_rgb_file, stack_rgb_file.replace('.fits','_v2.jpg'), stretch_r='arcsinh', stretch_g='arcsinh', stretch_b='arcsinh', pmin_r=0.25, pmin_g=0.25, pmin_b=0.25, pmax_r=85., pmax_g=85., pmax_b=85., make_nans_transparent=True, embed_avm_tags=True)
aplpy.make_rgb_image(stack_rgb_file, stack_rgb_file.replace('.fits','_v3.jpg'), stretch_r='arcsinh', stretch_g='arcsinh', stretch_b='arcsinh', pmin_r=10., pmin_g=10., pmin_b=10., pmax_r=99.75, pmax_g=99.75, pmax_b=99.75, make_nans_transparent=True, embed_avm_tags=True)
aplpy.make_rgb_image(stack_rgb_file, stack_rgb_file.replace('.fits','_v4.jpg'), stretch_r='arcsinh', stretch_g='arcsinh', stretch_b='arcsinh', pmin_r=10., pmin_g=10., pmin_b=10., pmax_r=85., pmax_g=85., pmax_b=85., make_nans_transparent=True, embed_avm_tags=True)
aplpy.make_rgb_image(stack_rgb_file, stack_rgb_file.replace('.fits','_v5.jpg'), stretch_r='arcsinh', stretch_g='arcsinh', stretch_b='arcsinh', pmin_r=5., pmin_g=5., pmin_b=5., pmax_r=85., pmax_g=85., pmax_b=85., make_nans_transparent=True, embed_avm_tags=True)
aplpy.make_rgb_image(stack_rgb_file, stack_rgb_file.replace('.fits','_v6.jpg'), stretch_r='arcsinh', stretch_g='arcsinh', stretch_b='arcsinh', pmin_r=0.25, pmin_g=0.25, pmin_b=0.25, pmax_r=99.75, pmax_g=99.75, pmax_b=99.75, make_nans_transparent=True, embed_avm_tags=True)

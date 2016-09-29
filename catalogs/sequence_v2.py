import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
from matplotlib.patches import Polygon
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


rotated_labels = []
def text_slope_match_line(text, x, y, line):
    global rotated_labels

    # find the slope
    xdata, ydata = line.get_data()

    x1 = xdata[0]
    x2 = xdata[-1]
    y1 = ydata[0]
    y2 = ydata[-1]

    rotated_labels.append({"text":text, "line":line, "p1":np.array((x1, y1)), "p2":np.array((x2, y2))})

def update_text_slopes():
    global rotated_labels

    for label in rotated_labels:
        # slope_degrees is in data coordinates, the text() and annotate() functions need it in screen coordinates
        text, line = label["text"], label["line"]
        p1, p2 = label["p1"], label["p2"]

        # get the line's data transform
        ax = line.get_axes()

        sp1 = ax.transData.transform_point(p1)
        sp2 = ax.transData.transform_point(p2)

        rise = (sp2[1] - sp1[1])
        run = (sp2[0] - sp1[0])

        slope_degrees = math.degrees(math.atan(rise/run))

        text.set_rotation(slope_degrees)

p1 = plt.subplot(1,1,1)
plt.axis([-5.,-25.,0.5,50000.])
axis = plt.gca()
xmin, xmax = axis.get_xlim()
ymin, ymax = axis.get_ylim()
plt.xticks(np.arange(-5.,-25.01,-5.))

majorLocator   = MultipleLocator(5)
majorFormatter = FormatStrFormatter('%d')
minorLocator   = MultipleLocator(1)
p1.xaxis.set_minor_locator(minorLocator)
p1.tick_params(axis='y', which='major', pad=1)

#==========================================================================================================================================#
mag_smooth           = np.genfromtxt('NGFS_FCC_cat_non_nucleated.dat',dtype=float, comments='#', delimiter='|', missing_values='-', skip_header=1, usecols = (6))    #
mag_nucleated        = np.genfromtxt('NGFS_FCC_cat_nucleated.dat' ,   dtype=float, comments='#', delimiter='|', missing_values='-', skip_header=1, usecols = (6))    #
                                                                                                                                                                    #
reff_smooth          = np.genfromtxt('NGFS_FCC_cat_non_nucleated.dat',dtype=float, comments='#', delimiter='|', missing_values='-', skip_header=1, usecols = (9))   #
reff_nucleated       = np.genfromtxt('NGFS_FCC_cat_nucleated.dat' ,   dtype=float, comments='#', delimiter='|', missing_values='-', skip_header=1, usecols = (9))   #
                                                                                                                                                                    #
#reff_NGFSnuc         = np.genfromtxt('NGFS_nuclei.dat'            ,   dtype=float, comments='#', delimiter='', missing_values='-', skip_header=0, usecols = (1))   #
#i_NGFSnuc            = np.genfromtxt('NGFS_nuclei.dat'            ,   dtype=float, comments='#', delimiter='', missing_values='-', skip_header=0, usecols = (0))   #
                                                                                                                                                                    #
reff_vandokkum15     = np.genfromtxt('vandokkum2015_47dwarfs.txt' ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (1))    #
M_g_vandokkum15      = np.genfromtxt('vandokkum2015_47dwarfs.txt' ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (2))    #
                                                                                                                                                                    #
janglee_Mv           = np.genfromtxt('janglee.dat'                ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=0, usecols = (0))    #
janglee_reff         = np.genfromtxt('janglee.dat'                ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=0, usecols = (1))    #
janglee_vi           = np.genfromtxt('janglee.dat'                ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=0, usecols = (2))    #
                                                                                                                                                                    #
globs_reff           = np.genfromtxt('table4.dat'                 ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=0, usecols = (14))   #
globs_g              = np.genfromtxt('table4.dat'                 ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=0, usecols = (8))    #

mwgc_distkpc         = np.genfromtxt('mwgc10_onetab.dat'          ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=0, usecols = (9))   #
mwgc_reff_arcmin     = np.genfromtxt('mwgc10_onetab.dat'          ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=0, usecols = (35))   #
mwgc_Mv              = np.genfromtxt('mwgc10_onetab.dat'          ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=0, usecols = (20))    #
mwgc_VI              = np.genfromtxt('mwgc10_onetab.dat'          ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=0, usecols = (24))    #

n1399gc_reff         = np.genfromtxt('NGC1399_lumsize.dat'        ,   dtype=float, comments='#', delimiter='', missing_values='*', skip_header=0, usecols = (2))   #
n1399gc_v            = np.genfromtxt('NGC1399_lumsize.dat'        ,   dtype=float, comments='#', delimiter='', missing_values='*', skip_header=0, usecols = (0))    #
                                                                                                                                                                    #
peacock_reff         = np.genfromtxt('peacock.dat'                ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (3))    #
peacock_g            = np.genfromtxt('peacock.dat'                ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (0))    #
peacock_gr           = np.genfromtxt('peacock.dat'                ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (1))    #
peacock_ri           = np.genfromtxt('peacock.dat'                ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (2))    #
                                                                                                                                                                    #
ferrarese06_reff     = np.genfromtxt('ferrarese06.dat'            ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (0))    #
ferrarese06_g        = np.genfromtxt('ferrarese06.dat'            ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (2))    #
                                                                                                                                                                    #
nuclei_reff          = np.genfromtxt('nuclei.dat'                 ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (2))    #
nuclei_g             = np.genfromtxt('nuclei.dat'                 ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (0))    #

ucds_reff            = np.genfromtxt('zhang_ucd.dat'              ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (18))    #
ucds_i             	 = np.genfromtxt('zhang_ucd.dat'              ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (11))    #

mucds_reff           = np.genfromtxt('liu_massiveucds.dat'        ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (7))    #
mucds_g            	 = np.genfromtxt('liu_massiveucds.dat'        ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (3))    #
mucds_gi           	 = np.genfromtxt('liu_massiveucds.dat'        ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (5))    #

fielddwarfs_reff     = np.genfromtxt('hopp_tab6.dat'              ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=0, usecols = (13))    #
fielddwarfs_Mg     	 = np.genfromtxt('hopp_tab6.dat'              ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=0, usecols = (4))    #
fielddwarfs_gr     	 = np.genfromtxt('hopp_tab6.dat'              ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=0, usecols = (5))    #
fielddwarfs_ri     	 = np.genfromtxt('hopp_tab6.dat'              ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=0, usecols = (7))    #

fielddwarfs_reff *= 1e3
fielddwarfs_Mi = fielddwarfs_Mg - fielddwarfs_gr - fielddwarfs_ri

janglee_Mi = janglee_Mv - janglee_vi

mwgc_MI 			= mwgc_Mv - 0.8 #mwgc_VI
mwgc_reff_new 		= (mwgc_reff_arcmin/60)*(math.pi/180)*mwgc_distkpc*1e3

n1399gc_reff_new	 = n1399gc_reff[(n1399gc_reff<70) & (n1399gc_reff>0)]
n1399gc_Mi 			 = n1399gc_v[(n1399gc_reff<70) & (n1399gc_reff>0)] - 0.9 - 31.51

mucds_Mi = mucds_g - mucds_gi - 31.51
ucds_Mi              = ucds_i - 31.51
#Mi_NGFSnuc 			 = i_NGFSnuc - 31.51
#reff_NGFSnuc_new	 = (reff_NGFSnuc * 98.417)

#nuclei_reff_new      = np.log10(nuclei_reff * 98.417)
nuclei_reff_new      = (nuclei_reff * 98.417)
nuclei_Mi            = nuclei_g - 31.51 - 1.0

#ferrarese06_reff_new = np.log10(ferrarese06_reff * 79.9942573831108)
ferrarese06_reff_new = (ferrarese06_reff * 79.9942573831108)
ferrarese06_Mi       = ferrarese06_g - 31.0874 - 1.0

peacock_Mi           = peacock_g - peacock_gr - peacock_ri - 24.471
#peacock_reff_new     = np.log10(peacock_reff)
peacock_reff_new     = (peacock_reff)

LSB_reff             = np.genfromtxt('LSB.dat'                   ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=0, usecols = (1))   #
LSB_Mi               = np.genfromtxt('LSB.dat'                   ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=0, usecols = (0))   #
                                                                                                                                                                  #
irvine_DL            = np.genfromtxt('irvine_Es_new.dat'         ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (0))   #
irvine_reff          = np.genfromtxt('irvine_Es_new.dat'         ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (1))   #
irvine_B             = np.genfromtxt('irvine_Es_new.dat'         ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (2))   #
irvine_V             = np.genfromtxt('irvine_Es_new.dat'         ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (3))   #
irvine_R             = np.genfromtxt('irvine_Es_new.dat'         ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (4))   #
irvine_I             = np.genfromtxt('irvine_Es_new.dat'         ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (5))   #

localdSphs_reff             = np.genfromtxt('dSphs.dat'         ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (1))   #
localdSphs_MV               = np.genfromtxt('dSphs.dat'         ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (0))   #

#localdSphs_reff_new = np.log10(localdSphs_reff)
localdSphs_reff_new = (localdSphs_reff)
localdSphs_Mi = localdSphs_MV + 0.308 -1.

modulus = -5 + 5 * np.log10(irvine_DL * 1000000)
#irvine_reff_new = np.log10(irvine_reff * 60 * 1/206264.8 * irvine_DL * 1000000)
irvine_reff_new = (irvine_reff * 60 * 1/206264.8 * irvine_DL * 1000000)
RI = irvine_R - irvine_I
ri = RI - 0.21
r  = irvine_V - 0.44 * (irvine_B - irvine_V) + 0.12
irvine_i  = r - ri
irvine_Mi = irvine_i - modulus

# DECam pixel size is 0.262968 arcsec
#79.9942573831108    pc/arcsec Virgo assuming D=16.5 Mpc

#globs_reff_new        = np.log10(globs_reff * 79.9942573831108)
globs_reff_new        = (globs_reff * 79.9942573831108)
globs_Mi              = globs_g - 1.0 - 31.0874

#reff_vandokkum15_new  = np.log10(reff_vandokkum15 * 1000)
reff_vandokkum15_new  = (reff_vandokkum15 * 1000)
M_i_vandokkum15_new   = M_g_vandokkum15 - 0.9

mag_smooth_new        = mag_smooth - 31.51
mag_nucleated_new     = mag_nucleated - 31.51
#reff_smooth_new       = np.log10(reff_smooth    * 96.96274)
#reff_nucleated_new    = np.log10(reff_nucleated * 96.96274)
reff_smooth_new       = (reff_smooth    * 96.96274)
reff_nucleated_new    = (reff_nucleated * 96.96274)
#==========================================================================================================================================#                                                               

#plt.scatter(Mi_NGFSnuc, reff_NGFSnuc_new,          s=20,  marker='o', c='green', linewidths=0, facecolor='green',zorder='5', edgecolor='green',alpha=0.5)
### Field dwarfs
plt.scatter(fielddwarfs_Mi, fielddwarfs_reff, s=40, marker='D', c='#ff9900', linewidths=0.9,  zorder='3', facecolor='#ff9900',alpha=0.4)
### Coma
plt.scatter(M_i_vandokkum15_new, reff_vandokkum15_new, s=40, marker='D', c='#CD5C5C', linewidths=0.9,  zorder='3', facecolor='#CD5C5C',alpha=0.25)
### GCs
plt.scatter(n1399gc_Mi, n1399gc_reff_new,         s=10, marker='.', c='#5e0200', linewidths=0.0,  zorder='5', facecolor='#5e0200', alpha=0.5)
plt.scatter(globs_Mi, globs_reff_new,             s=3, marker='.', c='#8b8878', linewidths=0.0,  zorder='5', facecolor='#8b8878', alpha=0.25)
plt.scatter(peacock_Mi, peacock_reff_new,         s=3, marker='s', c='#191919', linewidths=1.5,  zorder='5', facecolor='#191919', alpha=0.5)
plt.scatter(mwgc_MI, mwgc_reff_new,               s=5, marker='^', c='#191919', linewidths=1.5,  zorder='5', facecolor='#191919', alpha=0.5)
### other galaxies
plt.scatter(LSB_Mi,np.power(10,LSB_reff),         s=40, marker='D', c='#cc0099', linewidths=0.9,  zorder='5', facecolor='#cc0099', edgecolor='#cc0099',alpha=0.4)
plt.scatter(janglee_Mi,janglee_reff,              s=40, marker='D', c='red', linewidths=0.5, edgecolor='#330000',zorder='5', facecolor='red',alpha=0.4)
plt.scatter(irvine_Mi, irvine_reff_new,           s=45, marker='s', c='#191919', linewidths=0.9,  zorder='5', facecolor='#191919', edgecolor='#191919',alpha=0.5)
plt.scatter(ferrarese06_Mi, ferrarese06_reff_new, s=45, marker='s', c='#5e4b4b', linewidths=0.9,  zorder='5', facecolor='#5e4b4b', edgecolor='#5e4b4b',alpha=0.5)
### nuclei and UCDs
plt.scatter(nuclei_Mi, nuclei_reff_new,           s=40, marker='p', c='#e39e9c', linewidths=0.9,  zorder='5', facecolor='#e39e9c', edgecolor='#5e0200',alpha=0.5)
plt.scatter(ucds_Mi, ucds_reff,                   s=40, marker='H', c='#e39e9c', linewidths=0.9,  zorder='5', facecolor='#e39e9c', edgecolor='#5e0200',alpha=0.5)
plt.scatter(mucds_Mi, mucds_reff,                 s=40, marker='H', c='#e39e9c', linewidths=0.9,  zorder='5', facecolor='#e39e9c', edgecolor='#5e0200',alpha=0.75)
### NGFS
plt.scatter(mag_smooth_new, reff_smooth_new      , s=80, marker='o', c='red', linewidths=0.5, edgecolor='#330000',zorder='5', facecolor='red',alpha=0.4)
plt.scatter(mag_nucleated_new, reff_nucleated_new, s=80, marker='o', c='red', linewidths=0.5, edgecolor='#330000',zorder='5', facecolor='red',alpha=0.4)
plt.scatter(mag_nucleated_new, reff_nucleated_new, s=20,  marker='o', c='black', linewidths=0, facecolor='black',zorder='5', edgecolor='black',alpha=0.5)
### Local Group dwarfs
plt.scatter(localdSphs_Mi, localdSphs_reff_new,   s=40, marker='h', c='#ff9900', linewidths=0.9,  zorder='5', facecolor='#ff9900', edgecolor='#660000',alpha=0.5)
#==========================================================================================================================================#                                                               

### Making own legend
yalign = np.power(10,np.arange(1.75,-0.2,-0.1375))
xalign = -18.6
xtextalign = xalign-0.5
ngfstxt = 7
txt = 5
plt.scatter(xalign,yalign[0], s=72, marker='o', c='red', linewidths=0.5, edgecolor='#330000',zorder='5', facecolor='red',alpha=0.4)
plt.scatter(xalign,yalign[1], s=72, marker='o', c='red', linewidths=0.5, edgecolor='#330000',zorder='5', facecolor='red',alpha=0.4)
plt.scatter(xalign,yalign[1], s=18, marker='o', c='black', linewidths=0, facecolor='black',zorder='5', edgecolor='black',alpha=0.5)
plt.text(xtextalign,yalign[0], r'NGFS dwarfs (non-nucleated)', fontsize=ngfstxt, weight='heavy', verticalalignment='center', alpha=0.5)
plt.text(xtextalign,yalign[1], r'NGFS dwarfs (nucleated)', fontsize=ngfstxt, weight='heavy', verticalalignment='center', alpha=0.5)
### Coma
plt.scatter(xalign,yalign[2], s=36, marker='D', c='#cc0099', linewidths=0.9,  zorder='5', facecolor='#cc0099', edgecolor='#cc0099',alpha=0.4)
plt.text(xtextalign,yalign[2], r'Virgo UDGs (Mihos et al. 2015)', fontsize=txt, verticalalignment='center', alpha=0.5)
plt.scatter(xalign,yalign[3], s=36, marker='D', c='#CD5C5C', linewidths=0.9,  zorder='3', facecolor='#CD5C5C',alpha=0.25)
plt.text(xtextalign,yalign[3], r'Coma UDGs (van Dokkum et al. 2015)', fontsize=txt, verticalalignment='center', alpha=0.5)
### other galaxies
plt.scatter(xalign,yalign[4], s=36, marker='D', c='#ff9900', linewidths=0.9,  zorder='3', facecolor='#ff9900',alpha=0.4)
plt.text(xtextalign,yalign[4], r'Group dwarf galaxies (Hopp & Vennik 2014)', fontsize=txt, verticalalignment='center', alpha=0.5)
plt.scatter(xalign,yalign[5], s=40, marker='s', c='#191919', linewidths=0.9,  zorder='5', facecolor='#191919', edgecolor='#191919',alpha=0.5)
plt.text(xtextalign,yalign[5], r'Carnegie-Irvine galaxies (Ho et al. 2011)', fontsize=txt, verticalalignment='center', alpha=0.5)
plt.scatter(xalign,yalign[6], s=40, marker='s', c='#5e4b4b', linewidths=0.9,  zorder='5', facecolor='#5e4b4b', edgecolor='#5e4b4b',alpha=0.5)
plt.text(xtextalign,yalign[6], r'Virgo/Fornax ETGs (Ferrarese et al. 2012)', fontsize=txt, verticalalignment='center', alpha=0.5)
### Local Group dwarfs
plt.scatter(xalign,yalign[7], s=40, marker='h', c='#ff9900', linewidths=0.9,  zorder='5', facecolor='#ff9900', edgecolor='#660000',alpha=0.5)
plt.text(xtextalign,yalign[7], r'Local Group dSphs (McConnachie et al. 2012)', fontsize=txt, verticalalignment='center', alpha=0.5)
plt.scatter(xalign,yalign[8], s=36, marker='D', c='red', linewidths=0.5, edgecolor='#330000',zorder='5', facecolor='red',alpha=0.4)
plt.text(xtextalign,yalign[8], r'Virgo UFD (Jang & Lee 2014)', fontsize=txt, verticalalignment='center', alpha=0.5)
### nuclei and UCDs
plt.scatter(xalign,yalign[9], s=40, marker='p', c='#e39e9c', linewidths=0.9,  zorder='5', facecolor='#e39e9c', edgecolor='#5e0200',alpha=0.5)
plt.text(xtextalign,yalign[9], r'Fornax NSCs (Turner et al. 2012)', fontsize=txt, verticalalignment='center', alpha=0.5)
plt.scatter(xalign,yalign[10], s=40, marker='H', c='#e39e9c', linewidths=0.9,  zorder='5', facecolor='#e39e9c', edgecolor='#5e0200',alpha=0.5)
plt.text(xtextalign,yalign[10], r'Virgo UCDs (Zhang et al. 2015)', fontsize=txt, verticalalignment='center', alpha=0.5)
plt.scatter(xalign,yalign[11], s=40, marker='H', c='#e39e9c', linewidths=0.9,  zorder='5', facecolor='#e39e9c', edgecolor='#5e0200',alpha=0.75)
plt.text(xtextalign,yalign[11], r'M59/M60 UCDs (Liu et al. 2015)', fontsize=txt, verticalalignment='center', alpha=0.5)
### GCs
plt.scatter(xalign,yalign[12], s=10, marker='.', c='#5e0200', linewidths=0.0,  zorder='5', facecolor='#5e0200', alpha=0.5)
plt.text(xtextalign,yalign[12], r'Fornax GCs (Puzia et al. 2014)', fontsize=txt, verticalalignment='center', alpha=0.5)
plt.scatter(xalign,yalign[13], s=5, marker='.', c='#8b8878', linewidths=0.0,  zorder='5', facecolor='#8b8878', alpha=0.25)
plt.text(xtextalign,yalign[13], r'Virgo GCs (Jordan et al. 2009)', fontsize=txt, verticalalignment='center', alpha=0.5)
plt.scatter(xalign-0.1,yalign[14], s=3, marker='s', c='#191919', linewidths=1.5,  zorder='5', facecolor='#191919', alpha=0.5)
plt.scatter(xalign+0.1,yalign[14], s=5, marker='^', c='#191919', linewidths=1.5,  zorder='5', facecolor='#191919', alpha=0.5)
plt.text(xtextalign,yalign[14], r'MW+M31 GCs (Harris et al. 2010, Peacock et al. 2010)', fontsize=txt, verticalalignment='center', alpha=0.5)
# plt.scatter(xalign,yalign[14], s=5, marker='^', c='#191919', linewidths=1.5,  zorder='5', facecolor='#191919', alpha=0.5)
# plt.text(xtextalign,yalign[14], r'Milky Way GCs (Harris et al. 2010)', fontsize=txt, verticalalignment='center', alpha=0.5)

# verts = [(-18.5, 1.)] + [(-24.5, 1.)] + [(-24.5, 105.)] + [(-18.5, 105.)] + [(-18.5, 1.)] 
# poly = Polygon(verts, facecolor='none', edgecolor='#800000', linewidth=5, zorder='5', alpha=0.2)
# p1.add_patch(poly)

### Draw lines of const. surface-brightness and label them
x=np.linspace(-28,-5, num=100)       	# set range in x
lims = np.array([22.,24.,26.,28.])		# set SB values
for sbval in lims:
	line, = plt.plot(x,np.power(10,(-0.2*(x+31.51)-0.5*np.log10(2*3.14)+np.log10(98.417)+sbval/5.)),color='0.5',linewidth=2,alpha=0.15,zorder='2')
	xlocal = -5.5
	ylocal = np.power(10,(-0.2*(xlocal+31.51)-0.5*np.log10(2*3.14)+np.log10(98.417)+sbval/5.))
	t = plt.annotate('${:2.0f}$'.format(sbval), xy=(xlocal,ylocal), xytext=(-10, -2.5),fontsize=10,weight='light',textcoords='offset points',horizontalalignment='left',verticalalignment='bottom',color='k',zorder='10',alpha=0.5)
	text_slope_match_line(t, xlocal, ylocal, line)

### Align the text along the mu=28 mag/arcsec^2 limit line
line1, = plt.plot(x,np.power(10,(-0.2*(x+31.51)-0.5*np.log10(2*3.14)+np.log10(98.417)+28./5.)),color='0.5',linewidth=0,alpha=0.0,zorder='2')
x1 = -11.
y1 = np.power(10,(-0.2*(x1+31.51)-0.5*np.log10(2*3.14)+np.log10(98.417)+28./5.))
string1 = '$\left\langle \mu_i \\right\\rangle_{\\rm{e}}\,= \,28\,{\\rm{mag/arcsec}}^2$'
t = plt.annotate(string1,xy=(x1,y1),xytext=(-10,0),fontsize=11,weight='light',textcoords='offset points',horizontalalignment='left',verticalalignment='bottom',color='k',zorder='10',alpha=0.5)
text_slope_match_line(t, x1, y1, line1)

### Add polygons
xsmpl=np.linspace(-8,xmax,num=1000)       	# set range in x
verts = [(xmin, 0)] + [(-8., 0)] + list(zip(xsmpl,np.power(10,(-0.2*(xsmpl+31.51)-0.5*np.log10(2*3.14)+np.log10(98.417)+28./5.)))) + [(xmin, ymax)]
poly = Polygon(verts, facecolor='#efe8d2', edgecolor='none', alpha=0.5)
p1.add_patch(poly)

verts = [(-8., 0)] + [(xmax, 0)] + [(xmax, 75.)] + [(-8.0, 75.)] 
poly = Polygon(verts, facecolor='#efe8d2', edgecolor='none', alpha=0.2)
p1.add_patch(poly)

plt.plot((-8.04,-25.), (75.,75.), linewidth=2, color='#efe8d2', label='test', alpha=0.5)

plt.text(-6., 10000.,'NGFS insensitive',rotation='0',fontsize=12,weight='heavy',color='k',zorder='10',alpha=0.15)
plt.text(-10.2, 80.,'NGFS resolution limit',rotation='0',fontsize=10,weight='heavy',color='k',zorder='10',alpha=0.15)
plt.text(-13., 2.25,'NGFS point sources', rotation='0',fontsize=10,weight='heavy',color='k',zorder='10',alpha=0.15)

#axis.xaxis.set_ticks(np.arange(xmin, xmax, 0.712123))
#axis.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.yscale('log')
plt.xlabel(r'$M_i\, [{\rm mag}]$',fontdict={'fontsize':20})
plt.ylabel(r'$r_e\, [{\rm pc}]$',fontdict={'fontsize':20})
#p1.xaxis.set_label_coords(1.05, -0.025)
#plt.title(r'Stellar Systems in the nearby Universe')

plt.tight_layout()
update_text_slopes()

plt.savefig('sequence_v2.pdf',format='pdf')
plt.show()
exit()

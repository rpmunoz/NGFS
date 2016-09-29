#! /Applications/Ureka/variants/common/bin/python                                                                                                                   
                                                                                                                                                                    
#=========================================== PACKAGES ============================================#                                                                 
from      matplotlib          import    pyplot, mpl                                               #                                                                 
import    numpy               as        np                                                        #                                                                 
from      numpy               import    *                                                         #                                                                 
from      matplotlib          import    rc                                                        #                                                                 
import    matplotlib.image    as        mpimg                                                     #                                                                 
from      matplotlib.ticker   import    MaxNLocator, MultipleLocator, FormatStrFormatter          #                                                                 
#=================================================================================================#                                                                 
                                                                                                                                                                    
#====================================== STYLE DEFINITIONS ========================================#                                                                 
                                                                                                  #                                                                 
#FONT                                                                                             #                                                                 
mpl.rcParams['font.size']              = 14                                                       #                                                                 
mpl.rcParams['font.family']            = 'Arial Narrow'                                           #                                                                 
mpl.rcParams['mathtext.fontset']       = 'stixsans'                                               #                                                                 
mpl.rcParams['mathtext.default']       = 'regular'                                                #                                                                 
mpl.rcParams['text.usetex']            = 'False'                                                  #                                                                 
mpl.rcParams['text.latex.unicode']=True                                                           #                                                                 
                                                                                                  #                                                                 
                                                                                                  #                                                                 
#TICKS GENERAL                                                                                    #                                                                 
mpl.rcParams['xtick.labelsize']        = 'medium'                                                 #                                                                 
mpl.rcParams['ytick.labelsize']        = 'medium'                                                 #                                                                 
mpl.rcParams['xtick.color']            = 'k'                                                      #                                                                 
mpl.rcParams['ytick.color']            = 'k'                                                      #                                                                 
#MAJOR TICKS                                                                                      #                                                                 
mpl.rcParams['xtick.major.size']       = '5'                                                      #                                                                 
mpl.rcParams['ytick.major.size']       = '5'                                                      #                                                                 
mpl.rcParams['xtick.major.width']      = '0.75'                                                   #                                                                 
mpl.rcParams['ytick.major.width']      = '0.75'                                                   #                                                                 
mpl.rcParams['xtick.major.pad']        = '5'                                                      #                                                                 
mpl.rcParams['ytick.major.pad']        = '5'                                                      #                                                                 
#MINOR TICKS                                                                                      #                                                                 
mpl.rcParams['xtick.minor.size']       = '3'                                                      #                                                                 
mpl.rcParams['ytick.minor.size']       = '3'                                                      #                                                                 
mpl.rcParams['xtick.minor.width']      = '0.75'                                                   #                                                                 
mpl.rcParams['ytick.minor.width']      = '0.75'                                                   #                                                                 
mpl.rcParams['xtick.minor.pad']        = '5'                                                      #                                                                 
mpl.rcParams['ytick.minor.pad']        = '5'                                                      #                                                                 
#AXES                                                                                             #                                                                 
mpl.rcParams['axes.facecolor']         = 'white'                                                  #                                                                 
mpl.rcParams['axes.edgecolor']         = 'black'                                                  #                                                                 
mpl.rcParams['axes.linewidth']         = '1'                                                      #                                                                 
mpl.rcParams['axes.labelsize']         = 'medium'                                                 #                                                                 
mpl.rcParams['axes.labelcolor']        = 'black'                                                  #                                                                 
mpl.rcParams['axes.axisbelow']         = 'False'                                                  #                                                                 
                                                                                                  #                                                                 
#=================================================================================================#                                                                 
                                                                                                                                                                    
                                                                                                                                                                    
#============================================================ IMPORT DATA =================================================================#                        
                                                                                                                                                                    
mag_smooth           = np.genfromtxt('finallist_smooth.txt'       ,   dtype=float, comments='#', delimiter='', missing_values='-', skip_header=0, usecols = (9))    #
mag_nucleated        = np.genfromtxt('finallist_nucleated.txt'    ,   dtype=float, comments='#', delimiter='', missing_values='-', skip_header=0, usecols = (9))    #
                                                                                                                                                                    #
reff_smooth          = np.genfromtxt('finallist_smooth.txt'       ,   dtype=float, comments='#', delimiter='', missing_values='-', skip_header=0, usecols = (11))   #
reff_nucleated       = np.genfromtxt('finallist_nucleated.txt'    ,   dtype=float, comments='#', delimiter='', missing_values='-', skip_header=0, usecols = (11))   #
                                                                                                                                                                    #
reff_vandokkum15     = np.genfromtxt('vandokkum2015_47dwarfs.txt' ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (1))    #
M_g_vandokkum15      = np.genfromtxt('vandokkum2015_47dwarfs.txt' ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (2))    #
                                                                                                                                                                    #
globs_reff           = np.genfromtxt('table4.dat'                 ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=0, usecols = (14))   #
globs_g              = np.genfromtxt('table4.dat'                 ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=0, usecols = (8))    #
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


nuclei_reff_new      = np.log10(nuclei_reff * 98.417)
nuclei_Mi            = nuclei_g - 31.54 - 1.0


ferrarese06_reff_new = np.log10(ferrarese06_reff * 79.9942573831108)
ferrarese06_Mi       = ferrarese06_g - 31.0874 - 1.0

peacock_Mi           = peacock_g - peacock_gr - peacock_ri - 24.471
peacock_reff_new     = np.log10(peacock_reff)

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

localdSphs_reff_new = np.log10(localdSphs_reff)
localdSphs_Mi = localdSphs_MV + 0.308 -1.

modulus = -5 + 5 * np.log10(irvine_DL * 1000000)
irvine_reff_new = np.log10(irvine_reff * 60 * 1/206264.8 * irvine_DL * 1000000)
RI = irvine_R - irvine_I
ri = RI - 0.21
r  = irvine_V - 0.44 * (irvine_B - irvine_V) + 0.12
irvine_i  = r - ri
irvine_Mi = irvine_i - modulus

#79.9942573831108    pc/arcsec Virgo assuming D=16.5 Mpc

globs_reff_new        = np.log10(globs_reff * 79.9942573831108)
globs_Mi              = globs_g - 1.0 - 31.0874

reff_vandokkum15_new  = np.log10(reff_vandokkum15 * 1000)
M_i_vandokkum15_new   = M_g_vandokkum15 - 0.9

mag_smooth_new        = mag_smooth - 31.54
mag_nucleated_new     = mag_nucleated - 31.54
reff_smooth_new       = np.log10(reff_smooth    * 25.9526)
reff_nucleated_new    = np.log10(reff_nucleated * 25.9526)
#==========================================================================================================================================#                                                               

#============================================================ CALCULATONS =================================================================#                                                               
#==========================================================================================================================================#                                                               
                                                                                                                                                                                                           
                                                                                                                                                                                                           
                                                                                                                                                                                                           
#===================================== SET FIGURE AND AXES =====================================#                                                                                                          
figure     =  pyplot.figure(figsize=(10,10))                                                    #                                                                                                          
                                                                                                #                                                                                                          
diagram1   = figure.add_axes([0.01,0.01,0.50,0.5],     autoscale_on=True)                      #                                                                                                          

#===============================================================================================#                                                                                                          
                                                                                                                                                                                                           
                                                                                                                                                                                                           
#=============================================================== CONFIGURE AXES ===============================================================#                                                           
                                                                                                                                               #                                                           
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ DIAGRAM 1 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#       #                                                           
                                                                                                                                       #       #                                                           
#------------------------------------------------------ SET AX SCALE --------------------------------------------------------------#   #       #                                                           
#diagram1.set_xscale("log")                                                                                                        #   #       #                                                           
#diagram1.set_yscale("log")                                                                                                        #   #       #                                                           
diagram1.set_xlim(-6,-25)                                                                                                       #   #       #                                                           
diagram1.set_ylim(-0.2,4.7)                                                                                                           #   #       #                                                           
##---------------------------------------------------- SET TICK DIRECTIONS ----------------------------------------------------------#   #       #                                                           
diagram1.get_xaxis().set_tick_params(which='both', direction='in')                                                                 #   #       #                                                           
diagram1.get_yaxis().set_tick_params(which='both', direction='in')                                                                 #   #       #                                                           
#---------------------------------------------------- SET TICK LOCATIONS ----------------------------------------------------------#   #       #                                                           
major_xticks_diagram1 = MultipleLocator(2)                                                                                       #   #       #                                                           
minor_xticks_diagram1 = MultipleLocator(1)                                                                                       #   #       #                                                           
major_yticks_diagram1 = MultipleLocator(0.5)                                                                                       #   #       #                                                           
minor_yticks_diagram1 = MultipleLocator(0.1)                                                                                        #   #       #                                                           
diagram1.xaxis.set_major_locator(major_xticks_diagram1)                                                                            #   #       #                                                           
diagram1.xaxis.set_minor_locator(minor_xticks_diagram1)                                                                            #   #       #                                                           
diagram1.yaxis.set_major_locator(major_yticks_diagram1)                                                                           #   #       #                                                           
diagram1.yaxis.set_minor_locator(minor_yticks_diagram1)                                                                           #   #       #                                                           
#---------------------------------------------------- CUSTOM TICK LABELS ----------------------------------------------------------#   #       #                                                           
#diagram1.xaxis.set_ticks_position("top")                                                                                        #   #       #                                                           
#diagram1.xaxis.set_ticks_position("both")                                                                                          #   #       #                                                           
##diagram1.yaxis.set_ticks_position("right")                                                                                        #   #       #                                                           
##diagram1.yaxis.set_ticks_position("both")                                                                                         #   #       #                                                           
#xlabels_diagram1 = diagram1.get_xticklabels()                                                                                      #   #       #                                                           
#ylabels_diagram1 = diagram1.get_yticklabels()                                                                                      #   #       #                                                           
##xlabels_diagram1[0].set_visible(False)                                                                                            #   #       #                                                           
##xlabels_diagram1[6].set_visible(False)                                                                                            #   #       #                                                           
#ylabels_diagram1[1].set_visible(False)                                                                                             #   #       #                                                           
#ylabels_diagram1[6].set_visible(False)                                                                                             #   #       #                                                           
##------------------------------------- HIDE ALL TICK LABELS OR SET LABELS INDIVIDUALLY --------------------------------------------#   #       #                                                           
##diagram1.set_xticklabels([])                                                                                                      #   #       #                                                           
#diagram1.set_yticklabels([])                                                                                                      #   #       #                                                           
#---------------------------------------------------- FORMAT TICK LABELS ----------------------------------------------------------#   #       #                                                           
#format_major_xticks_diagram1 = FormatStrFormatter('%1.3f')                                                                        #   #       #                                                           
#format_major_yticks_diagram1 = FormatStrFormatter('%5.2f')                                                                        #   #       #                                                           
#diagram1.xaxis.set_major_formatter(format_major_xticks_diagram1)                                                                  #   #       #                                                           
#diagram1.yaxis.set_major_formatter(format_major_yticks_diagram1)                                                                  #   #       #                                                           
#----------------------------------------------------------- GRID -----------------------------------------------------------------#   #       #                                                           
#diagram3.xaxis.grid(which='major', dash_capstyle = 'round' , dashes=[0,2,0,2] , color='red' , linewidth=0.5, zorder='2')          #   #       #                                                           
#diagram3.xaxis.grid(which='minor', dash_capstyle = 'round' , dashes=[0,2,0,2] , color='red' , linewidth=0.5, zorder='2')          #   #       #                                                           
#diagram3.yaxis.grid(which='minor', dash_capstyle = 'round' , dashes=[0,2,0,2] , color='red' , linewidth=0.5, zorder='2')          #   #       #                                                           
#diagram3.yaxis.grid(which='major', dash_capstyle = 'round' , dashes=[0,2,0,2] , color='red' , linewidth=0.5, zorder='2')          #   #       #                                                           
#----------------------------------------------------- SET AXIS LABELS ------------------------------------------------------------#   #       #                                                           
#diagram1.xaxis.set_label_position("top")                                                                                       #   #       #                                                           
#diagram1.yaxis.set_label_position("right")                                                                                        #   #       #                                                           
diagram1.set_xlabel('M$_i$ [mag]')                                                                                          #   #       #                                                           
diagram1.set_ylabel('log r$_e$ [pc] ')                                                                                            #   #       #                                                           
#-------------------------------------------------- SET AXIS LABEL PADDING --------------------------------------------------------#   #       #                                                           
diagram1.xaxis.labelpad        = 2                                                                                                 #   #       #                                                           
diagram1.yaxis.labelpad        = 2                                                                                                 #   #       #                                                           
#----------------------------------------------------------------------------------------------------------------------------------#   #       #                                                           
                                                                                                                                       #       #                                                           
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#       #                                                           
                                                                                                                                               #                                                           
#==============================================================================================================================================#                                                           
                                                                                                                                                                                                           
                                                                                                                                                                                                           
#================================================================================= PLOT THINGS ==================================================================================#                         

diagram1.scatter(mag_smooth_new, reff_smooth_new, s=10, marker='^', c='k', linewidths=0.75,  zorder='50', facecolor='black')                                                  #                                                  #                         
diagram1.scatter(mag_nucleated_new, reff_nucleated_new, s=50, marker='.', c='k', linewidths=0.75, facecolor='white', zorder='50', edgecolor='k')                                                  #                                                  #                         
diagram1.scatter(mag_nucleated_new, reff_nucleated_new, s=1.5, marker='.', c='k', linewidths=0.5, facecolor='k', zorder='50', edgecolor='k')                                                  #                                                  #                         

diagram1.scatter(M_i_vandokkum15_new, reff_vandokkum15_new, s=12, marker='x', c=[0., 0.5, 1.], linewidths=1.1,  zorder='5', facecolor=[0., 0.5, 1.])                                                  #                                                  #                         

diagram1.scatter(globs_Mi, globs_reff_new, s=5, marker='.', c=[0.62, 0.02, 0.03], linewidths=0,  zorder='5', facecolor=[0.62, 0.02, 0.03])                                                  #                                                  #                         


diagram1.scatter(peacock_Mi, peacock_reff_new, s=2, marker='s', c=[1, 0.43, 0.48], linewidths=0.0,  zorder='5', facecolor=[1, 0.43, 0.48])                                                  #                                                  #                         

diagram1.scatter(LSB_Mi, LSB_reff, s=40, marker='1', c=[0., 0., 1.], linewidths=1.6,  zorder='5', facecolor=[0., 0., 1.], edgecolor=[0., 0., 1.])                                                  #                                                  #                         

diagram1.scatter(irvine_Mi, irvine_reff_new, s=5, marker='v', c=[0.93, 0.25, 0], linewidths=1.25,  zorder='5', facecolor=[0.93, 0.25, 0], edgecolor=[0.93, 0.25, 0])                                                  #                                                  #                         

diagram1.scatter(ferrarese06_Mi, ferrarese06_reff_new, s=5, marker='<', c=[1., 0.60, 0.07], linewidths=1.25,  zorder='5', facecolor=[1., 0.60, 0.07], edgecolor=[1., 0.60, 0.07])                                                  #                                                  #                         

diagram1.scatter(nuclei_Mi, nuclei_reff_new, s=3, marker='D', c=[1., 0., 0.4], linewidths=0.9,  zorder='5', facecolor=[1., 0., 0.4], edgecolor=[1., 0., 0.4])                                                  #                                                  #                         

diagram1.scatter(localdSphs_Mi, localdSphs_reff_new, s=6, marker='s', c=[0.05, 0.75, 0.91], linewidths=0.9,  zorder='5', facecolor=[0.05, 0.75, 0.91], edgecolor=[0.05, 0.75, 0.91])                                                  #                                                  #                         


x=np.linspace(-28,-5, num=100)       # set range in x                                                                   #                               #                                                       #                         

f22=-0.2*(x+31.54)-0.5*np.log10(2*3.14)+np.log10(98.417)+22./5.
#f23=-0.2*(x+31.54)-0.5*np.log10(2*3.14)+np.log10(98.417)+23./5.
f24=-0.2*(x+31.54)-0.5*np.log10(2*3.14)+np.log10(98.417)+24./5.
#f25=-0.2*(x+31.54)-0.5*np.log10(2*3.14)+np.log10(98.417)+25./5.
f26=-0.2*(x+31.54)-0.5*np.log10(2*3.14)+np.log10(98.417)+26./5.
#f27=-0.2*(x+31.54)-0.5*np.log10(2*3.14)+np.log10(98.417)+27./5.
f28=-0.2*(x+31.54)-0.5*np.log10(2*3.14)+np.log10(98.417)+28./5.
#f29=-0.2*(x+31.54)-0.5*np.log10(2*3.14)+np.log10(98.417)+29./5.
#f30=-0.2*(x+31.54)-0.5*np.log10(2*3.14)+np.log10(98.417)+30./5.

diagram1.plot( x,  f22, dash_capstyle  = 'round' , dashes=[10,3,10,3]     , color=[0.39, 0.39, 0.39]    , linewidth=1., zorder='2')          #                               #                                                       #                         
#diagram1.plot( x,  f23, dash_capstyle  = 'round' , dashes=[10,3,10,3]     , color='k'    , linewidth=0.75, zorder='2')          #                               #                                                       #                         
diagram1.plot( x,  f24, dash_capstyle  = 'round' , dashes=[10,3,10,3]     , color=[0.39, 0.39, 0.39]    , linewidth=1., zorder='2')          #                               #                                                       #                         
#diagram1.plot( x,  f25, dash_capstyle  = 'round' , dashes=[10,3,10,3]     , color='k'    , linewidth=0.75, zorder='2')          #                               #                                                       #                         
diagram1.plot( x,  f26, dash_capstyle  = 'round' , dashes=[10,3,10,3]     , color=[0.39, 0.39, 0.39]    , linewidth=1., zorder='2')          #                               #                                                       #                         
#diagram1.plot( x,  f27, dash_capstyle  = 'round' , dashes=[10,3,10,3]     , color='k'    , linewidth=0.75, zorder='2')          #                               #                                                       #                         
diagram1.plot( x,  f28, dash_capstyle  = 'round' , dashes=[10,3,10,3]     , color=[0.39, 0.39, 0.39]    , linewidth=1., zorder='2')          #                               #                                                       #                         
#diagram1.plot( x,  f29, dash_capstyle  = 'round' , dashes=[10,3,10,3]     , color='k'    , linewidth=0.75, zorder='2')          #                               #                                                       #                         
#diagram1.plot( x,  f30, dash_capstyle  = 'round' , dashes=[10,3,10,3]     , color='k'    , linewidth=0.75, zorder='2')          #                               #                                                       #                         


diagram1.text(0.03,  0.323-0.054,'22'   ,fontsize=7, weight='heavy', color=[0.39,0.39,0.39], zorder='10', transform=diagram1.transAxes, rotation='39')                                    #                                                 #                                              
diagram1.text(0.03,  0.404-0.054,'24'   ,fontsize=7, weight='heavy', color=[0.39,0.39,0.39], zorder='10', transform=diagram1.transAxes, rotation='39')                                    #                                                 #                                              
diagram1.text(0.03,  0.485-0.054,'26'   ,fontsize=7, weight='heavy', color=[0.39,0.39,0.39], zorder='10', transform=diagram1.transAxes, rotation='39')                                    #                                                 #                                              
diagram1.text(0.03,  0.567-0.054,'28'   ,fontsize=7, weight='heavy', color=[0.39,0.39,0.39], zorder='10', transform=diagram1.transAxes, rotation='39')                                    #                                                 #                                              
#
diagram1.text(0.24, 0.86-0.052,'$\left\langle \mu  \\right\\rangle_{\\rm{e}} \, = \, m_{\\rm{tot}}+2.5\,log(2\pi r_{\\rm{e}}^2)$'       ,fontsize=8, weight='heavy', color=[0.39, 0.39, 0.39], zorder='10', transform=diagram1.transAxes, rotation='39')                                    #                                                 #                                              
#================================================================================================================================================================================#                         
                                                                                                                                                                                                           
                                                                                                                                                                                                           
#================================================================ SAVE TO EPS =================================================================#                                                           
#pyplot.show()                                                                                                                                 #                                                           
pyplot.savefig('sequence.pdf', bbox_inches='tight', dpi=300, format='pdf')
#==============================================================================================================================================#                                                           

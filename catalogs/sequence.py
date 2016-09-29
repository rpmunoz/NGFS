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
mpl.rcParams['font.size']              = 12                                                       #                                                                                                        
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
mpl.rcParams['xtick.major.width']      = '0.75'                                                      #                                                                                                        
mpl.rcParams['ytick.major.width']      = '0.75'                                                      #                                                                                                        
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
mpl.rcParams['axes.linewidth']         = '0.75'                                                      #                                                                                                        
mpl.rcParams['axes.labelsize']         = 'medium'                                                 #                                                                                                        
mpl.rcParams['axes.labelcolor']        = 'black'                                                  #                                                                                                        
mpl.rcParams['axes.axisbelow']         = 'False'                                                  #                                                                                                        
                                                                                                  #                                                                                                        
#=================================================================================================#                                                                                                        
                                                                                                                                                                                                           
                                                                                                                                                                                                           
#============================================================ IMPORT DATA =================================================================#                                                               

mag_smooth      = np.genfromtxt('finallist_smooth.txt'     ,   dtype=float, comments='#', delimiter='', missing_values='-', skip_header=0, usecols = (9))   #                                                               
mag_nucleated   = np.genfromtxt('finallist_nucleated.txt'  ,   dtype=float, comments='#', delimiter='', missing_values='-', skip_header=0, usecols = (9))   #                                                               

reff_smooth      = np.genfromtxt('finallist_smooth.txt'     ,   dtype=float, comments='#', delimiter='', missing_values='-', skip_header=0, usecols = (11))   #                                                               
reff_nucleated   = np.genfromtxt('finallist_nucleated.txt'  ,   dtype=float, comments='#', delimiter='', missing_values='-', skip_header=0, usecols = (11))   #                                                               

reff_vandokkum15      = np.genfromtxt('vandokkum2015_47dwarfs.txt'     ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (1))   #                                                               
M_g_vandokkum15       = np.genfromtxt('vandokkum2015_47dwarfs.txt'     ,   dtype=float, comments='#', delimiter='', missing_values='_', skip_header=1, usecols = (2))   #                                                               

reff_vandokkum15_new  =    np.log10(reff_vandokkum15*1000)
M_g_vandokkum15_new   =  M_g_vandokkum15-1

mag_smooth_new       =    mag_smooth    - 31.54
mag_nucleated_new    =    mag_nucleated - 31.54
reff_smooth_new      =    np.log10(reff_smooth    * 25.9526)
reff_nucleated_new   =    np.log10(reff_nucleated * 25.9526)
#==========================================================================================================================================#                                                               



#============================================================ CALCULATONS =================================================================#                                                               
#==========================================================================================================================================#                                                               
                                                                                                                                                                                                           
                                                                                                                                                                                                           
                                                                                                                                                                                                           
#===================================== SET FIGURE AND AXES =====================================#                                                                                                          
figure     =  pyplot.figure(figsize=(10,10))                                                    #                                                                                                          
                                                                                                #                                                                                                          
diagram1   = figure.add_axes([0.01,0.01,0.50,0.5],     autoscale_on=True)                      #                                                                                                          
#diagram2   = figure.add_axes([0.01,0.01,0.90,0.48],     autoscale_on=True)                      #                                                                                                          

#===============================================================================================#                                                                                                          
                                                                                                                                                                                                           
                                                                                                                                                                                                           
#=============================================================== CONFIGURE AXES ===============================================================#                                                           
                                                                                                                                               #                                                           
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ DIAGRAM 1 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#       #                                                           
                                                                                                                                       #       #                                                           
#------------------------------------------------------ SET AX SCALE --------------------------------------------------------------#   #       #                                                           
#diagram1.set_xscale("log")                                                                                                        #   #       #                                                           
#diagram1.set_yscale("log")                                                                                                        #   #       #                                                           
#diagram1.set_xlim(-17,-7)                                                                                                       #   #       #                                                           
#diagram1.set_ylim(1.5,4)                                                                                                           #   #       #                                                           
##---------------------------------------------------- SET TICK DIRECTIONS ----------------------------------------------------------#   #       #                                                           
diagram1.get_xaxis().set_tick_params(which='both', direction='in')                                                                 #   #       #                                                           
diagram1.get_yaxis().set_tick_params(which='both', direction='in')                                                                 #   #       #                                                           
#---------------------------------------------------- SET TICK LOCATIONS ----------------------------------------------------------#   #       #                                                           
major_xticks_diagram1 = MultipleLocator(1)                                                                                       #   #       #                                                           
minor_xticks_diagram1 = MultipleLocator(0.2)                                                                                       #   #       #                                                           
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

diagram1.scatter(mag_smooth_new, reff_smooth_new, s=20, marker='.', c='k', linewidths=0.5,  zorder='5', facecolor='black')                                                  #                                                  #                         
diagram1.scatter(mag_nucleated_new, reff_nucleated_new, s=50, marker='.', c='red', linewidths=0.5, facecolor='white', zorder='5')                                                  #                                                  #                         
diagram1.scatter(mag_nucleated_new, reff_nucleated_new, s=1, marker='.', c='red', linewidths=0.5, facecolor='black', zorder='5')                                                  #                                                  #                         

diagram1.scatter(M_g_vandokkum15_new, reff_vandokkum15_new, s=20, marker='x', c='red', linewidths=1.5,  zorder='5', facecolor='red')                                                  #                                                  #                         

#                         
x=np.linspace(-20,-5, num=100)       # set range in x                                                                   #                               #                                                       #                         
f25=-0.2*(x+31.54)-0.5*np.log10(2*3.14)+np.log10(98.417)+25./5.
              # define function                                                                     #                               #                                                      #                         
f24=-0.2*(x+31.54)-0.5*np.log10(2*3.14)+np.log10(98.417)+24./5.
f26=-0.2*(x+31.54)-0.5*np.log10(2*3.14)+np.log10(98.417)+26./5.
f27=-0.2*(x+31.54)-0.5*np.log10(2*3.14)+np.log10(98.417)+27./5.

diagram1.plot( x,  f24, dash_capstyle  = 'round' , dashes=[3,1,3,1]     , color='k'    , linewidth=0.5, zorder='2')          #                               #                                                       #                         
diagram1.plot( x,  f25, dash_capstyle  = 'round' , dashes=[3,1,3,1]     , color='k'    , linewidth=0.5, zorder='2')          #                               #                                                       #                         
diagram1.plot( x,  f26, dash_capstyle  = 'round' , dashes=[3,1,3,1]     , color='k'    , linewidth=0.5, zorder='2')          #                               #                                                       #                         
diagram1.plot( x,  f27, dash_capstyle  = 'round' , dashes=[3,1,3,1]     , color='k'    , linewidth=0.5, zorder='2')          #                               #                                                       #                         


diagram1.text(0.030,  0.738-0.005,'24'        ,fontsize=8.5, weight='heavy', color=[0.36, 0.56, 0.96], zorder='10', transform=diagram1.transAxes, rotation='-39')                                    #                                                 #                                              
diagram1.text(0.065,  0.790-0.005,'25'        ,fontsize=8.5, weight='heavy', color=[0.36, 0.56, 0.96], zorder='10', transform=diagram1.transAxes, rotation='-39')                                    #                                                 #                                              
diagram1.text(0.100,  0.840-0.005,'26'         ,fontsize=8.5, weight='heavy', color=[0.36, 0.56, 0.96], zorder='10', transform=diagram1.transAxes, rotation='-39')                                    #                                                 #                                              
diagram1.text(0.140,  0.890-0.005,'27'         ,fontsize=8.5, weight='heavy', color=[0.36, 0.56, 0.96], zorder='10', transform=diagram1.transAxes, rotation='-39')                                    #                                                 #                                              

diagram1.text(0.24, 0.846+0.006,'$\left\langle \mu  \\right\\rangle_{\\rm{e}} \, = \, m_{\\rm{tot}}+2.5\,log(2\pi r_{\\rm{e}}^2)$'       ,fontsize=10, weight='heavy', color=[0.36, 0.56, 0.96], zorder='10', transform=diagram1.transAxes, rotation='-39')                                    #                                                 #                                              
#================================================================================================================================================================================#                         
                                                                                                                                                                                                           
                                                                                                                                                                                                           
#================================================================ SAVE TO EPS =================================================================#                                                           
#pyplot.show()                                                                                                                                 #                                                           
pyplot.savefig('sequence.eps', bbox_inches='tight', dpi=300)                                                                                     #                                                           
#==============================================================================================================================================#                                                           

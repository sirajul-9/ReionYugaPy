#Initialise required variables
#**********************************************************#
nbody_path="../N-body-threads"                 #write the n-body output directory path here, don't put slash at the end
halo_cat_path='../FoF-Halo-finder-master'     #write the halo_catalogue directory path here, don't put slash at the end
sfac = 2                                                   #factor for rescaling the  grid
Nthreads=4                                                 #set the number of openmp threads, see for which value you are getting maximum speedup

#parameters of the reionization model
nion=23.21                                                 #parameter Nion, fid. value=23.21
rmfp=20.00                                                 #parameter Rmfp, fid. value=20.00
mmin=10                                                    #parameter Mmin(in units of DM_mass), fid. value=1.09 * 10**9 M_sun


redshifts=[7,8,9]                                 #values of redshifts                           
#**********************************************************#


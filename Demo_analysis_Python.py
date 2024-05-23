# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 22:18:34 2023

@author: Peter Rupprecht, p.t.r.rupprecht+centripetal+propagation@gmail.com
"""

"""

Compute delay maps from astrocytic in vivo calcium imaging data

Code written by Peter Rupprecht (rupprecht@hifo.uzh.ch) in 2023

Code is deposited in this repository: https://github.com/HelmchenLabSoftware/Centripetal_propagation_astrocytes
See this repository for an explanation of the code

Please cite this paper when using the code:
Rupprecht P, Duss SN, Becker D, Lewis CM, Bohacek J, and Helmchen F.
"Centripetal integration of past events by hippocampal astrocytes regulated by locus coeruleus."
https://www.nature.com/articles/s41593-024-01612-8


"""


# import required packages

import numpy as np
import glob
import matplotlib.pyplot as plt
import tifffile
from skimage.registration import phase_cross_correlation
from scipy.ndimage import fourier_shift

# generic pattern of input files; if multiple files are recognized by this
# pattern, separate delay maps will be computed for each file, but the
# final result will be obtained by averaging across these delay maps
filenames = glob.glob('Raw calcium imaging data/FOV_excerpt_recording*.tif');

# frame rate specific for this recording; change to make the scaling of the delay map correct
framerate = 4.419;



# Go through each recording segment; delay maps will be computed for each segment and averaged afterwards
Corr_maps_all = [0]*len(filenames)
for kkk,filename in enumerate(filenames):
    
    print('Computing delay map for: ' + filename)

    
    # Read raw data from hard disk
    movie = tifffile.imread(filename).astype('float')
    
    # Compute reference (mean across FOV)
    mean_activity = np.squeeze(np.mean(np.mean(movie,axis=2),axis=1))

    # Drift correction (align segments with respect to each other, if necessary)
    
    if kkk == 0:
        
        #compute alignment template from first segment
        ref_excerpt = np.mean(movie,axis=0)
        
    else:
        
        #compute shift to alignment template using very simple movement correction
        this_excerpt = np.mean(movie,axis=0)
        
        shift =  phase_cross_correlation(this_excerpt, ref_excerpt)[0]
        
        movie = np.roll(movie,[0,np.int32(-shift[0]),np.int32(-shift[1])],[0,1,2])
        
        # print(shift)
    
    # Compute delay map based on correlation function peaks
    max_delay = 30; # in #frames; do increase value if the frame rate is higher
    
    Corr_map = np.zeros([movie.shape[1],movie.shape[2]]);
    for j in np.arange(movie.shape[1]):
        for k in np.arange(movie.shape[2]):
            # extract time trace of this pixel
            trace =  movie[:,j,k]
            # initialize cross_correlation vector
            cross_correlation = np.zeros([2*max_delay+1,]);
            for kk in np.arange(-max_delay,max_delay):
                # cross-correlate X and Y
                X = mean_activity
                Y = np.roll(trace,kk,axis=0);
                if kk >= 0:
                    cross_correlation[kk+max_delay] = np.corrcoef(X[kk:],Y[kk:],rowvar=False)[0,1]
                else:
                    cross_correlation[kk+max_delay] = np.corrcoef(X[:kk],Y[:kk],rowvar=False)[0,1]

            # find peak of the cross-correlation
            xi = np.argmax(cross_correlation);
            delay = xi - max_delay
            # assign peak delay to the delay map pixel
            Corr_map[j,k] = -delay

    Corr_maps_all[kkk] = Corr_map
    
print('All delay maps are computed.')


# Concatenate delay maps
Corr_maps_all_concatenated = np.expand_dims(Corr_maps_all[0],axis=2)
for Corr_maps_single in Corr_maps_all[1:]:
    Corr_maps_all_concatenated = np.concatenate( (Corr_maps_all_concatenated,np.expand_dims(Corr_maps_single,axis=2)),axis=2)

# Average concatenated delay maps
delay_map = np.mean(Corr_maps_all_concatenated,axis=2)/framerate;


# Visualize results
plt.figure(77), plt.imshow(delay_map,vmin=-2,vmax=2,cmap='viridis_r')
plt.colorbar()
plt.grid(False)
plt.title('Map of delays (s)')


plt.figure(78), plt.imshow(ref_excerpt,vmin=np.quantile(ref_excerpt,0.0),vmax=np.quantile(ref_excerpt,0.95),cmap='Greys_r')
plt.colorbar()
plt.grid(False)
plt.title('Anatomical reference')



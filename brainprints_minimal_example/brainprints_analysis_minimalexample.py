import os
import numpy as np
import pandas as pd
from helpers_brainprints import brainprints_from_FS
from helpers_similaritycalc import compute_similarity

"""
    @author: p23-markushs

    Brainprint similarity
    Calculate brainprint similarity between twins at a given visit
    Minimal example, april 2025
"""

# Paths
path_FS = '/ess/p23/cluster/data/bids_fmri/s2c_gridtask_prisma/bids/derivatives/freesurfer_7.3.2'
pipeline = 'long' # at LCBC either 'cross', 'long', 'long_alt1' or 'long_alt2'
parcellation = 'aparc' # or 'aparc.a2009s' for Destrieux, but then you'll have to calculate some additional features during the FreeSurfer steps

# Load main dataframe
df = pd.read_csv('path_to_main_dataframe.csv', sep=',')

# df needs to contain at least the columns shown in toy dataframe below,
# with sub and ses being the folder names used in Freesurfer structure for a given subject-session combo
# e.g. recon-all output is found under /.../bids/derivatives/freesurfer_7.3.2/sub-1700021/ses-01

# sub	        ses	    twin_pair	twin_number	    zygo
# sub-1700021	ses-01	2	        1	            1
# sub-1700022	ses-01	2	        2	            1
# sub-1700041	ses-01	4	        1	            1
# sub-1700042	ses-01	4	        2	            1
# sub-1700061	ses-01	6	        1	            1
# sub-1700062	ses-01	6	        2 	            1
# sub-1700071	ses-01	7	        1	            1
# sub-1700072	ses-01	7	        2	            1
# sub-1700081	ses-01	8	        1	            2
# sub-1700082	ses-01	8	        2	            2
# sub-1700091	ses-01	9	        1	            1
# sub-1700092	ses-01	9	        2	            1
# sub-1700101	ses-01	10	        1	            2
# sub-1700102	ses-01	10	        2	            2

# ensure sorting of dataframe so that twin number 1s are shown first, then twin number 2s. Allows for simple vertical splitting below
df = df.reset_index(drop=True).sort_values(['twin_number', 'twin_pair', 'ses'])

# feature extraction routines (see helpers_brainprints.py)
# the subfunction "_get_fs_stats", called by "brainprints_from_FS" is where the program looks for feature details. Missing files (routines not being run)
# will lead to errors and warnings, so this function will have to be updated to fit your needs (and we're of course happy to assist here :-) )
# note: the pipeline flag let's you extract features from a longitudinal pipeline, assuming that Freesurfer-processed data is available
# under the subject-session path discussed above but with a long-pipeline-suffix added, 
# e.g. /.../bids/derivatives/freesurfer_7.3.2/sub-1700021/ses-01.long.base
fps_twins, regions = brainprints_from_FS(df, path_FS, parcellation, pipeline) 

# split complete feature sets (extracted over all subjects in df_main) set into twin 1s and twin 2s
fp_twin1, fp_twin2 = np.vsplit(fps_twins, 2)

# Keep certain feature categories based on their prefixes (check region variable to get full overview)
features_to_keep = ['aparcAREA_', 'aparcCT_', 'aparcVOL_', 'aparcMC_'] # underscore important to not pick up aparcAREApial_ etc (if not intended)

# index of features to keep
whichfeatures = np.array([idx for idx, region in enumerate(regions) if any(substring in region for substring in features_to_keep)])

# Brainprint similarity routines (see helpers_similaritycalc.py)
# The output variable mtx represent a similarity matrix estimated using the chosen similarity metric (default is measurement error). Diagonal represent twin pair similarity 
print("Calculating overall brainprint similarity")
mtx, _ = compute_similarity(fp_twin1, fp_twin2, 'error', features=whichfeatures)
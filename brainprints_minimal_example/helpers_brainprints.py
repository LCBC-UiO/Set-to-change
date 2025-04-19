#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: p23-markushs
"""

# Consider adding "sulcal depth" to brainprints_from_FS output (_get_fs_stats) using alternative approach
# See Madan 2019 "Robust estimation of sulcal morphology", https://cmadan.github.io/calcSulc/
#%%

import os
import glob
import numpy as np
import pandas as pd

def brainprints_from_FC(df, path_FC, parcellation='MMP', fisherstd=False, whichdata='edges'):

    fps = []
    
    for sub,ses,tedana_run1,tedana_run2 in zip(df['sub'], df['ses'], df['tedana_run1'],df['tedana_run2']):
        if all((int(tedana_run1)<3, int(tedana_run2)<3)):
            run='allcombined'
        elif (int(tedana_run1)<3):
            run='1'
        elif (int(tedana_run2)<3):
            run='2'
        else:
            print('WARNING, No valid tedana ratings for {0}, {1}'.format(sub,ses))
        
        if whichdata == 'weights':
            file_fconn = os.path.join(path_FC,sub,ses,'task-grid_run-{0}_parcellation-{1}_fconn-matrix.tsv'.format(run,parcellation))
            if parcellation=='MMP':
                default_shape = (379,)
            elif parcellation=='CA':
                default_shape = (718,)
            
            if os.path.isfile(file_fconn):
                fp = np.loadtxt(file_fconn, skiprows=1)
                if fisherstd==True:
                    fp = np.arctanh(fp)                    
                np.fill_diagonal(fp, 0)
                fp = fp.sum(axis=1)
            else:
                print("{0} --> {1} --> run-{2} is missing. Path = {3}".format(sub,ses,run,file_fconn))
                fp = np.full(default_shape, np.nan)                
                
        else:
            file_fconn = os.path.join(path_FC,sub,ses,'task-grid_run-{0}_parcellation-{1}_fconn-flat.txt'.format(run,parcellation))
            if parcellation=='MMP':
                default_shape = (71631,)
            elif parcellation=='CA':
                default_shape = (257403,)
            
            if os.path.isfile(file_fconn):            
                fp = np.loadtxt(file_fconn)
                if fisherstd==True:
                    fp = np.arctanh(fp)                    
            else:
                print("{0} --> {1} --> run-{2} is missing. Path = {3}".format(sub,ses,run,file_fconn))
                fp = np.full(default_shape, np.nan)
                        
        fps.append(fp)

    return np.vstack(fps)

def brainprints_from_FS(df, path_FS, parcellation='aparc', pipeline='cross'):

    fps = []
    for idx, sub in enumerate(df['sub']):
        ses = df.iloc[idx]['ses']
        y_sub = _get_fs_stats(path_FS,sub,ses,parcellation,pipeline,samseg_drop_all=True)
        fps.append(y_sub)
        if idx==0:
            regions = y_sub.index

    return np.vstack(fps), regions

def brainprints_from_MIND(df, path_MIND, pipeline='long', parcellation='aparc'):

    fps = []
    for idx, sub in enumerate(df['sub']):
        ses = df.iloc[idx]['ses']
        file_in_flat = os.path.join(path_MIND,sub,pipeline,'MIND_'+ses+'_parcellation-'+parcellation+'_flat.txt')
        y_sub = np.loadtxt(file_in_flat)
        fps.append(y_sub)

        if idx == 0:
            file_in_mtx = file_in_flat.replace('_flat.txt', '_matrix.tsv')
            regionnames = pd.read_csv(file_in_mtx, delimiter='\t', index_col=0).index.tolist()
            rows, cols = np.triu_indices(len(regionnames),k=1)
            edgenames = ['-'.join((regionnames[row], regionnames[col])) for row,col in zip(rows,cols)]
            region_in_edges_indices = []
            for region_idx in np.arange(len(regionnames)):
                region_in_edges_indices.append(np.where((rows == region_idx) | (cols == region_idx))[0])

    return np.vstack(fps), regionnames, edgenames, region_in_edges_indices

def brainprints_from_shapeDNA(df, path_shapeDNA, pipeline='long', analysis='norm1reweight1'):

    fps = []
    for idx, sub in enumerate(df['sub']):
        ses = df.iloc[idx]['ses']
        file_in = glob.glob(os.path.join(path_shapeDNA,sub,pipeline,analysis,ses+'*.brainprint.csv'))[0]
        df_tmp = pd.read_csv(file_in, sep=',')
        # Drop the first three rows where 'Unnamed: 0' is 'area', 'volume', or 'ev0'
        df_tmp = df_tmp[df_tmp['Unnamed: 0'].isin(['area', 'volume', 'ev0']) == False]
        # Drop the 'Unnamed: 0' column
        df_tmp.drop(columns=['Unnamed: 0'], inplace=True)
        fps.append(df_tmp.values)

    return np.stack(fps,axis=2), df_tmp.columns.tolist()


def DNA_fix_nans(fp):
    # assuming structure (eigenvalues, structures, subjects)
    # Find indices where any measurement is NaN
    structure_subject_nan = np.any(np.isnan(fp), axis=0)
    # np.where produces tuple of indices, unpacking the tuple for zip with *
    for struct, sub in zip(*np.where(structure_subject_nan)):
        if struct >= 30:
            fp[:,struct,sub] = np.nanmean(fp[:,30:,sub], axis=1)
        else:
            print('WARNING, DNA NaN found in non-cortical brainprint')

    return fp


def FC_add_fd(df, path_FC):
    fds_filtered_mean = []

    for sub,ses in zip(df['sub'], df['ses']):
        file_fd_filtered_mean= os.path.join(path_FC,sub,ses,'fd_filtered_mean_run-allcombined.txt')
        if os.path.isfile(file_fd_filtered_mean):
            fd_f_mean = np.loadtxt(file_fd_filtered_mean)
            fds_filtered_mean.append(fd_f_mean)
        else:
            print("{0} --> {1} is missing FD-data. Path = {2}".format(sub,ses,file_fd_filtered_mean))
            fds_filtered_mean.append(np.nan)

    return np.array(fds_filtered_mean)


def _get_fs_stats(path_FS, sub, ses, parcellation, pipeline='cross', samseg_drop_all=False):

    # Define helper function to get the samseg timepoint
    def get_samseg_tp(sub, ses):
        # Handle special cases
        if sub in ['sub-1700101', 'sub-1780951']:
            return 'tp00' + str(int(ses[-1]) - 1)
        else:
            return 'tp00' + ses[-1]

    # Determine suffixes and base paths based on pipeline
    pipeline_suffix = {
        'cross': {'stats_suffix': os.path.join(sub, ses), 'samseg_folder': 'samseg-cross', 'base_suffix': ''},
        'long': {'stats_suffix': os.path.join(sub, ses + '.long.base'), 'samseg_folder': 'samseg-long', 'base_suffix': '.long.base'},
        'long_alt1': {'stats_suffix': os.path.join(sub, ses + '.long.base_alt1'), 'samseg_folder': 'samseg-long_alt1', 'base_suffix': '.long.base_alt1'},
        'long_alt2': {'stats_suffix': os.path.join(sub, ses + '.long.base_alt2'), 'samseg_folder': 'samseg-long_alt2', 'base_suffix': '.long.base_alt2'},
    }

    if pipeline not in pipeline_suffix:
        raise ValueError(f"Unknown pipeline: {pipeline}")

    suffix_info = pipeline_suffix[pipeline]
    stats_dir = os.path.join(path_FS, suffix_info['stats_suffix'], 'stats')
    mri_dir = os.path.join(path_FS, suffix_info['stats_suffix'], 'mri')

    # Determine samseg paths
    if pipeline == 'cross':
        samseg_stats = os.path.join(path_FS, sub, 'samseg-cross', ses, 'samseg.stats')
    else:
        samseg_tp = get_samseg_tp(sub, ses)
        samseg_stats = os.path.join(path_FS, sub, suffix_info['samseg_folder'], samseg_tp, 'samseg.stats')
        if ses == 'ses-01' and not os.path.isfile(samseg_stats):
            print(f'Using cross for {sub}, {ses} as no long-processed samseg data are found')
            samseg_stats = os.path.join(path_FS, sub, 'samseg-cross', ses, 'samseg.stats')

    sbtiv_stats = samseg_stats.replace('samseg.stats', 'sbtiv.stats')

    # Paths for hemisphere-specific stats
    hemis = ['lh', 'rh']
    paths = {}
    for hemi in hemis:
        paths[f'global_white_{hemi}'] = os.path.join(stats_dir, f'{hemi}.global.white.stats')
        paths[f'global_pial_{hemi}'] = os.path.join(stats_dir, f'{hemi}.global.pial.stats')
        paths[f'global_gwr_{hemi}'] = os.path.join(stats_dir, f'{hemi}.global.w-g.pct.stats')
        paths[f'aparc_{hemi}'] = os.path.join(stats_dir, f'{hemi}.{parcellation}.stats')
        paths[f'aparc_pial_{hemi}'] = os.path.join(stats_dir, f'{hemi}.{parcellation}.pial.stats')
        paths[f'aparc_sd_{hemi}'] = os.path.join(stats_dir, f'{hemi}.sulcaldepth.{parcellation}.stats')
        paths[f'gwr_{hemi}'] = os.path.join(stats_dir, f'{hemi}.w-g.pct.stats')
        paths[f'subfields_{hemi}'] = os.path.join(mri_dir, f'{hemi}.hippoSfVolumes-T1-T2multispectral.v22.txt')

    paths['wmparc'] = os.path.join(stats_dir, 'wmparc.stats')
    paths['samseg'] = samseg_stats
    paths['sbtiv'] = sbtiv_stats

    # Check if all required files exist
    if not all(os.path.isfile(p) for p in paths.values()):
        print(f"Missing data for {sub}, {ses}")
        for key, path in paths.items():
            print(f"{key}: {path} exists: {os.path.isfile(path)}")
        return None

    # Initialize list to collect data
    data_series = []

    # Map hemisphere abbreviations to 'left' and 'right'
    hemi_map = {'lh': 'left', 'rh': 'right'}

    # Process hemisphere-specific files
    for hemi in hemis:
        hemi_side = hemi_map[hemi]

        # Read features
        global_area = _get_feature_from_fsstats(paths[f'global_white_{hemi}'], 2)
        global_area.index = [f'globalAREA_{hemi_side}']
        global_area_pial = _get_feature_from_fsstats(paths[f'global_pial_{hemi}'], 2)
        global_area_pial.index = [f'globalAREApial_{hemi_side}']
        global_vol = _get_feature_from_fsstats(paths[f'global_white_{hemi}'], 3)
        global_vol.index = [f'globalVOL_{hemi_side}']
        global_ct = _get_feature_from_fsstats(paths[f'global_white_{hemi}'], 4)
        global_ct.index = [f'globalCT_{hemi_side}']
        global_mc = _get_feature_from_fsstats(paths[f'global_white_{hemi}'], 6)
        global_mc.index = [f'globalMC_{hemi_side}']
        global_mc_pial = _get_feature_from_fsstats(paths[f'global_pial_{hemi}'], 6)
        global_mc_pial.index = [f'globalMCpial_{hemi_side}'] 
        global_fold = _get_feature_from_fsstats(paths[f'global_white_{hemi}'], 8)
        global_fold.index = [f'globalFOLD_{hemi_side}']
        global_fold_pial = _get_feature_from_fsstats(paths[f'global_pial_{hemi}'], 8)
        global_fold_pial.index = [f'globalFOLDpial_{hemi_side}']     
        global_gwr = _get_feature_from_fsstats(paths[f'global_gwr_{hemi}'], 5, skiprows=53, index_col=4)
        global_gwr.index = [f'globalGWR_{hemi_side}']
                
        area = _get_feature_from_fsstats(paths[f'aparc_{hemi}'], 2)
        area.index = [f'aparcAREA_{hemi_side}_{name}' for name in area.index]
        area_pial = _get_feature_from_fsstats(paths[f'aparc_pial_{hemi}'], 2)
        area_pial.index = [f'aparcAREApial_{hemi_side}_{name}' for name in area_pial.index]
        vol = _get_feature_from_fsstats(paths[f'aparc_{hemi}'], 3)
        vol.index = [f'aparcVOL_{hemi_side}_{name}' for name in vol.index]
        ct = _get_feature_from_fsstats(paths[f'aparc_{hemi}'], 4)
        ct.index = [f'aparcCT_{hemi_side}_{name}' for name in ct.index]
        mc = _get_feature_from_fsstats(paths[f'aparc_{hemi}'], 6)
        mc.index = [f'aparcMC_{hemi_side}_{name}' for name in mc.index]
        mc_pial = _get_feature_from_fsstats(paths[f'aparc_pial_{hemi}'], 6)
        mc_pial.index = [f'aparcMCpial_{hemi_side}_{name}' for name in mc_pial.index]
        fold = _get_feature_from_fsstats(paths[f'aparc_{hemi}'], 8)
        fold.index = [f'aparcFOLD_{hemi_side}_{name}' for name in fold.index]
        fold_pial = _get_feature_from_fsstats(paths[f'aparc_pial_{hemi}'], 8)
        fold_pial.index = [f'aparcFOLDpial_{hemi_side}_{name}' for name in fold_pial.index]        
        sd = _get_feature_from_fsstats(paths[f'aparc_sd_{hemi}'], 4, skiprows=60)
        sd.index = [f'aparcSD_{hemi_side}_{name}' for name in sd.index]
        gwr = _get_feature_from_fsstats(paths[f'gwr_{hemi}'], 5, skiprows=57, index_col=4)
        gwr = gwr[~gwr.index.str.startswith('unknown')]
        gwr.index = [f'aparcGWR_{hemi_side}_{name}' for name in gwr.index]
        subfields = pd.read_csv(paths[f'subfields_{hemi}'], header=None, sep=' ', index_col=0)[1].round(3)
        subfields.index = [f'subfields_{hemi_side}_{name}' for name in subfields.index]

        data_series.extend([global_area,
                            global_area_pial,
                            global_vol,
                            global_ct,
                            global_mc,
                            global_mc_pial,
                            global_fold,
                            global_fold_pial,
                            global_gwr,
                            area, 
                            area_pial,
                            vol,
                            ct,
                            mc,
                            mc_pial,
                            fold,
                            fold_pial,
                            sd,
                            gwr,
                            subfields])

    # Process wmparc
    wmparc = _get_feature_from_fsstats(paths['wmparc'], 5, skiprows=65, index_col=4)
    wmparc = wmparc[wmparc.index.str.startswith('wm')]
    wmparc.index = wmparc.index.str.replace('wm-lh-', 'aparcWM_left_').str.replace('wm-rh-', 'aparcWM_right_')
    data_series.append(wmparc)

    # Process samseg and sbtiv
    samseg = pd.read_csv(paths['samseg'], header=None, skiprows=1, sep=',', index_col=0)[1].round(3)
    samseg.index = ['samseg_' + name.split(' ')[-1] for name in samseg.index]
    sbtiv = pd.read_csv(paths['sbtiv'], header=None, sep=',', index_col=0)[1].round(3)
    sbtiv.index = ['sbtiv_' + name.split(' ')[-1] for name in sbtiv.index]

    # Drop unwanted samseg indices
    samseg_drop = [
        'samseg_Soft_Nonbrain_Tissue', 'samseg_Fluid_Inside_Eyes', 'samseg_non-WM-hypointensities',
        'samseg_WM-hypointensities', 'samseg_Optic-Chiasm', 'samseg_Right-choroid-plexus',
        'samseg_Left-choroid-plexus', 'samseg_Left-vessel', 'samseg_Right-vessel'
    ]
    if samseg_drop_all:
        samseg_drop.extend([
            'samseg_Skull', 'samseg_Brain-Stem', 'samseg_CSF', 'samseg_Left-Cerebellum-Cortex',
            'samseg_Right-Cerebellum-Cortex', 'samseg_4th-Ventricle', 'samseg_Left-Cerebral-Cortex',
            'samseg_Right-Cerebral-Cortex', 'samseg_Left-Cerebellum-White-Matter',
            'samseg_Right-Cerebellum-White-Matter', 'samseg_Right-Cerebral-White-Matter',
            'samseg_Left-Cerebral-White-Matter', 'samseg_Right-Inf-Lat-Vent',
            'samseg_Left-Inf-Lat-Vent', 'samseg_3rd-Ventricle', 'samseg_Left-Lateral-Ventricle',
            'samseg_Right-Lateral-Ventricle', 'samseg_5th-Ventricle'
        ])
    samseg = samseg[~samseg.index.isin(samseg_drop)]
    data_series.extend([samseg, sbtiv])

    # Combine all data into a single Series
    combined_series = pd.concat(data_series)

    # Sort index to ensure alignment
    combined_series.sort_index(inplace=True)

    return combined_series

def _get_feature_from_fsstats(file_path,feature_idx,skiprows=61,header=None,sep='\s+',index_col=0):
    df_feature = pd.read_csv(
        file_path,
        skiprows=skiprows,
        header=header,
        sep=sep,
        index_col=index_col)[feature_idx]

    return df_feature

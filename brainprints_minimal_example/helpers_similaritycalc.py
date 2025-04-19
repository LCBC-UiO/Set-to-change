#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: p23-markushs
"""
#%%

import numpy as np
from scipy.stats import zscore, multivariate_normal 
from sklearn.metrics import pairwise_distances
from scipy.linalg import inv


# Function to compute similarity
def compute_similarity(matrix1, matrix2, metric='error', features=None, zscaling=False):

    if features is None or features.size == 0:
        print("No features provided, estimating from all")
    else:
        #print("Features provided, length {}".format(features.shape))        
        matrix1 = matrix1[:, features]
        matrix2 = matrix2[:, features]
        
    if matrix1.ndim == 1:
        if metric == 'error':
            matrix1 =  matrix1.reshape(-1, 1)
            matrix2 =  matrix2.reshape(-1, 1)            
        else:
            matrix1 =  matrix1.reshape(1, -1)
            matrix2 =  matrix2.reshape(1, -1)      

    if zscaling:
        matrix1 = zscore(matrix1)
        matrix2 = zscore(matrix2)

    if metric == 'correlation':
        #print('similarity via correlation')
        # calculate pearson correlation
        sim_matrix = _generate_correlation_map(matrix1, matrix2)
    elif metric == 'error':
        #print('similarity via measurement error')
        sim_matrix = 1-(_generate_measurement_error_map(matrix1, matrix2))
    elif metric == 'mahalanobis':
        #print('similarity via mahalanobis')        
        cov_matrix, inv_cov_matrix = _compute_covariance(matrix1, matrix2)
        dist_matrix = pairwise_distances(matrix1, matrix2, metric=metric, VI=inv_cov_matrix)
        min_distance = np.min(dist_matrix[np.nonzero(dist_matrix)])
        max_distance = np.max(dist_matrix)
        dist_matrix = (dist_matrix - min_distance) / (max_distance - min_distance)
        sim_matrix = 1 - dist_matrix
    elif metric == 'cosine':
        #print('similarity via cosine')
        sim_matrix = 1-(pairwise_distances(matrix1, matrix2, metric=metric))
    else:
        #print('similarity via {}'.format(metric))        
        dist_matrix = pairwise_distances(matrix1, matrix2, metric=metric)
        min_distance = np.min(dist_matrix[np.nonzero(dist_matrix)])
        max_distance = np.max(dist_matrix)
        dist_matrix = (dist_matrix - min_distance) / (max_distance - min_distance)
        sim_matrix = 1 - dist_matrix

    # calculate individual differentiability (within-subject similarity - between subjects similarity)
    within_subj_sim = np.mean(np.diag(sim_matrix))
    between_subj_sim = np.mean(sim_matrix[np.triu_indices(sim_matrix.shape[0], k = 1)])
    idiff = within_subj_sim - between_subj_sim
 
    return sim_matrix, idiff

def calculate_twin_idiff(df, sim_matrix, zygo):
    twin_indices = _get_twin_indices(df, zygo)
    twin_similarities = np.array([sim_matrix[i, j] for i, j in twin_indices])
    avg_twin_similarity = np.mean(twin_similarities)

    all_indices = np.arange(sim_matrix.shape[0])
    other_similarities = []
    for i, j in twin_indices:
        non_twin_indices = np.delete(all_indices, [i, j])
        other_similarities.extend(sim_matrix[i, non_twin_indices])
        other_similarities.extend(sim_matrix[j, non_twin_indices])

    avg_other_similarity = np.mean(other_similarities)
    idiff = avg_twin_similarity - avg_other_similarity
    return idiff, avg_twin_similarity, avg_other_similarity, twin_similarities, other_similarities

def DNA_compute_similarity(training_data, test_data, vectorized=True):
    
    num_eigenvalues, num_structures, num_participants = training_data.shape    
    
    if vectorized == False:        
        overall_similarity_matrix = np.zeros((num_participants, num_participants))
        per_structure_similarity_matrix = np.zeros((num_participants, num_participants, num_structures))
    
        # Compute the global covariance for each structure across all subjects in the training data
        global_covariances = np.var(training_data, axis=(0, 2)) + 1e-6 # Adding a small value for numerical stability
    
        # Iterate over participants in the test dataset
        for i in range(num_participants):
            # Iterate over participants in the training dataset
            for j in range(num_participants):
                log_prob_sum = 0
                # Iterate over structures
                for structure in range(num_structures):
                    x_diff = test_data[:, structure, i] - training_data[:, structure, j]
                    # Calculate the log probability density for the current structure
                    log_prob = multivariate_normal.logpdf(
                        x_diff,
                        mean=np.zeros(num_eigenvalues), # Mean is zero because we're looking at the difference
                        cov=global_covariances[structure] * np.eye(num_eigenvalues)
                    )
                    log_prob_sum += log_prob
                    # Store per-structure similarity
                    per_structure_similarity_matrix[i, j, structure] = -log_prob  # Negating to convert to similarity score
    
                # Store overall similarity
                overall_similarity_matrix[i, j] = -log_prob_sum # Negating to convert to similarity score
                
    else:
        # Expand dimensions for vectorized computation of differences
        test_expanded = test_data[:, :, np.newaxis, :]  # Shape: [num_eigenvalues, num_structures, 1, num_participants]
        train_expanded = training_data[:, :, :, np.newaxis]  # Shape: [num_eigenvalues, num_structures, num_participants, 1]
        
        # Difference matrix: shape [num_eigenvalues, num_structures, num_participants, num_participants]
        diff_matrix = test_expanded - train_expanded
        
        # Compute variances and add a small value for numerical stability
        variances = np.var(training_data, axis=(0, 2)) + 1e-6  # Shape: [num_structures]
        
        # Compute log probabilities
        log_var = np.log(variances)
        const_term = num_eigenvalues * np.log(2 * np.pi)
        
        # Correct variance reshaping to match broadcasting with the difference matrix
        variances_reshaped = variances[:, np.newaxis, np.newaxis]
        
        # Sum across eigenvalues for each structure, resulting in shape [num_structures, num_participants, num_participants]
        log_probs = -0.5 * (const_term + np.sum(log_var) + 
                            np.sum((diff_matrix ** 2) / variances_reshaped, axis=0))
        
        # Sum log probabilities across structures
        overall_similarity_matrix = np.sum(log_probs, axis=0).T
        per_structure_similarity_matrix = log_probs.T
                
    return overall_similarity_matrix, per_structure_similarity_matrix

def calculate_similarity_score(similarity_matrix, score_type='idiff'):
    self_similarity = np.diag(similarity_matrix)
    if score_type == 'max_binary':
        max_other_similarity_row = np.max(np.where(
            np.eye(len(self_similarity)) == 1, -np.inf, similarity_matrix), axis=1)
        max_other_similarity_col = np.max(np.where(
            np.eye(len(self_similarity)) == 1, -np.inf, similarity_matrix), axis=0)
        score_row = np.mean(self_similarity > max_other_similarity_row)
        score_col = np.mean(self_similarity > max_other_similarity_col)
    elif score_type == 'idiff':
        mean_other_similarity_row = np.nanmean(
            np.where(np.eye(len(self_similarity)) == 1, np.nan, similarity_matrix), axis=1)
        mean_other_similarity_col = np.nanmean(
            np.where(np.eye(len(self_similarity)) == 1, np.nan, similarity_matrix), axis=0)
        score_row = np.mean(self_similarity - mean_other_similarity_row)
        score_col = np.mean(self_similarity - mean_other_similarity_col)

    score = (score_row + score_col) / 2  # Average or another combination logic
    return score

def _generate_correlation_map(x, y):
    """Correlate each n with each m.

    Parameters
    ----------
    x : np.array
      Shape N X T.

    y : np.array
      Shape M X T.

    Returns
    -------
    np.array
      N X M array in which each element is a correlation coefficient.

    """
    mu_x = x.mean(1)
    mu_y = y.mean(1)
    n = x.shape[1]
    if n != y.shape[1]:
        raise ValueError('x and y must ' +
                         'have the same number of timepoints.')
    s_x = x.std(1, ddof=n - 1)
    s_y = y.std(1, ddof=n - 1)
    cov = np.dot(x,
                 y.T) - n * np.dot(mu_x[:, np.newaxis],
                                  mu_y[np.newaxis, :])
    return cov / np.dot(s_x[:, np.newaxis], s_y[np.newaxis, :])


def _generate_measurement_error_map(x, y):
    """
    Calculate the mean measurement error between subjects over all features,
    automatically handling negative or zero values by shifting the data.

    Parameters
    ----------
    x : np.array
        Shape N x F.
    y : np.array
        Shape N x F.

    Returns
    -------
    np.array
        N x N array where each element is the mean measurement error between subjects.

    Notes
    -----
    If the input arrays contain negative or zero values, the function shifts the data
    to ensure all values are positive before computing measurement errors.
    """
    if x.shape != y.shape:
        raise ValueError('x and y must have the same shape.')

    # Ensure both x and y are at least 2-dimensional
    if x.ndim == 1:
        x = x[:, np.newaxis]
    if y.ndim == 1:
        y = y[:, np.newaxis]

    # Shift the data if necessary to ensure all values are positive
    min_value = min(np.min(x), np.min(y))
    if min_value <= 0:
        shift_constant = -min_value + 1e-6  # Small epsilon to avoid zero
        x_shifted = x + shift_constant
        y_shifted = y + shift_constant
    else:
        x_shifted = x
        y_shifted = y

    # Initialize the NxN matrix to store mean measurement errors
    N = x_shifted.shape[0]
    mean_measurement_errors = np.zeros((N, N))

    # Calculate the mean measurement error for each subject pair
    for i in range(N):
        for j in range(N):
            # Compute the measurement errors for all features between subjects i and j
            numerator = np.abs(x_shifted[i, :] - y_shifted[j, :])
            denominator = (x_shifted[i, :] + y_shifted[j, :]) / 2

            # Handle cases where denominator is zero
            with np.errstate(divide='ignore', invalid='ignore'):
                measurement_errors = np.divide(numerator, denominator)
                measurement_errors[~np.isfinite(measurement_errors)] = np.nan  # Set inf and -inf to NaN

            # Use nanmean to avoid NaN issues if any denominator was zero
            mean_measurement_errors[i, j] = np.nanmean(measurement_errors)

    return mean_measurement_errors


# Function to compute covariance matrix and its inverse (with Tikhonov regularization to avoid nans)
def _compute_covariance(matrix1, matrix2):
    epsilon = 1e-5  # small regularization factor
    cov_matrix = np.cov(np.vstack([matrix1, matrix2]), rowvar=False)
    cov_matrix += epsilon * np.eye(cov_matrix.shape[0])  # add small value to the diagonal
    inv_cov_matrix = inv(cov_matrix)
    return cov_matrix, inv_cov_matrix

def _get_twin_indices(df, zygo):
    # Find indices for twin pairs based on zygosity
    df['Index'] = df.index
    twins = df[df['zygo'] == zygo]
    twin_pairs = twins.groupby('twin_pair')['Index'].apply(list).values
    return [pair for pair in twin_pairs if len(pair) == 2]
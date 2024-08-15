# -*- coding: utf-8 -*-
"""This module defines multiple functions to compute different dynamics-based
attributes using ProDy."""

import prody as pr
import numpy as np
import sys
import math
import os
import re
from utilities.paths import ANM_folder, GNM_folder, MSF_folder, PRS_folder, CORR_folder, HINGE_folder

__author__ = "Anupam Banerjee"
__date__ = "March 2024"
__maintainer__ = "Anupam Banerjee"
__email__ = "anupam.banerjee@stonybrook.edu"
__status__ = "Alpha"


def cal_MSF(pdb_path, chain, mut_res, n_modes, enm, **kwargs):
    """
    Computes the absolute value, Z-Score, and normalized Z-Score of the MSF of 
    a residue based on the chosen ENM and chosen number of modes.

    Parameters:
    pdb_path (str): Path to the PDB file.
    chain (str): Chain identifier.
    mut_res (int or str): Residue number or 'all' for all residues.
    enm (str): Elastic Network Model to use ('ANM' or 'GNM').
    n_modes (int or str): Number of modes to consider ('all', 'hot', 'soft', or a fraction of total modes).
    kwargs: Additional arguments like 'frac_var'.

    Returns:
    numpy.ndarray or tuple: MSF values for all residues or a specific residue.
    """
    supported_arguments = {'pdb_path', 'chain', 'mut_res', 'enm', 'n_modes', 'frac_var'}
    for key in kwargs:
        if key not in supported_arguments:
            raise ValueError(f"Unsupported argument: {key}")

    pdb_code = re.search(r'([^/]+)\.pdb$', pdb_path)
    if not pdb_code:
        raise ValueError(f"Invalid PDB file path: {pdb_path}")

    pdb_name = f"{pdb_code.group(1)}_{chain}"
    msf_file_path = os.path.join(MSF_folder, f"{pdb_name}_{enm}_{n_modes}_all.npy")

    # Check if the MSF file exists
    if os.path.exists(msf_file_path):
        result = np.load(msf_file_path)
        if mut_res == 'all':
            return result
        else:
            res_index = np.where(result[:, 0] == mut_res)[0]
            if res_index.size == 0:
                raise ValueError(f"Residue number {mut_res} not found in the MSF file.")
            res_index = res_index[0]
            res_fluc = result[res_index, 1]
            res_zscore = result[res_index, 2]
            res_norm_zscore = result[res_index, 3]
            return res_fluc, res_zscore, res_norm_zscore

    sel_enm = None

    if enm == 'ANM':
        anm_file_path = os.path.join(ANM_folder, f"{pdb_name}.anm.npz")
        if os.path.exists(anm_file_path):
            sel_enm = pr.loadModel(anm_file_path)
        else:
            pdb = pr.parsePDB(pdb_path, chain=chain, subset='ca')
            sel_enm = pr.ANM()
            sel_enm.buildHessian(pdb)
            sel_enm.calcModes(n_modes='all')
            pr.saveModel(sel_enm, filename=anm_file_path)

    elif enm == 'GNM':
        gnm_file_path = os.path.join(GNM_folder, f"{pdb_name}.gnm.npz")
        if os.path.exists(gnm_file_path):
            sel_enm = pr.loadModel(gnm_file_path)
        else:
            pdb = pr.parsePDB(pdb_path, chain=chain, subset='ca')
            sel_enm = pr.GNM()
            sel_enm.buildKirchhoff(pdb)
            sel_enm.calcModes(n_modes='all')
            pr.saveModel(sel_enm, filename=gnm_file_path)

    else:
        raise ValueError(f"Unsupported ENM type: {enm}")

    tot_modes = sel_enm.numModes()
    if isinstance(n_modes, str):
        if n_modes == 'all':
            sq_fluc = pr.calcSqFlucts(sel_enm)
        elif n_modes == 'hot':
            num_modes = 1 + math.ceil(0.02 * tot_modes)
            sq_fluc = pr.calcSqFlucts(sel_enm[-num_modes:])
        elif n_modes == 'soft':
            num_modes = math.ceil(0.02 * tot_modes)
            sq_fluc = pr.calcSqFlucts(sel_enm[:num_modes])
        else:
            raise ValueError(f"Invalid n_modes argument: {n_modes}")
    elif isinstance(n_modes, (int, float)):
        if n_modes > 0 and n_modes < 1:
            num_modes = math.ceil(n_modes * tot_modes)
            sq_fluc = pr.calcSqFlucts(sel_enm[:num_modes])
        elif n_modes >= 1 and n_modes <= tot_modes:
            sq_fluc = pr.calcSqFlucts(sel_enm[:int(n_modes)])
        else:
            raise ValueError(f"Invalid n_modes argument: {n_modes}")
    elif 'frac_var' in kwargs:
        frac_var = kwargs['frac_var']
        if 0 < frac_var <= 1:
            fra_var = pr.calcFractVariance(sel_enm)
            cumulative_sum = np.cumsum(fra_var)
            min_mode = 1 + np.argmax(cumulative_sum >= frac_var)
            sq_fluc = pr.calcSqFlucts(sel_enm[:min_mode])
        else:
            raise ValueError(f"Invalid frac_var argument: {frac_var}")
    else:
        sq_fluc = pr.calcSqFlucts(sel_enm)

    pdb = pr.parsePDB(pdb_path, chain=chain, subset='ca')
    resnums = pdb.getResnums()

    mean_fluc = sq_fluc.mean()
    std_fluc = sq_fluc.std()
    max_fluc = sq_fluc.max()
    min_fluc = sq_fluc.min()
    zscore_min = (min_fluc - mean_fluc) / std_fluc
    zscore_max = (max_fluc - mean_fluc) / std_fluc

    res_zscore = (sq_fluc - mean_fluc) / std_fluc
    res_norm_zscore = (res_zscore - zscore_min) / (zscore_max - zscore_min)

    result = np.vstack((resnums, sq_fluc, res_zscore, res_norm_zscore)).T
    np.save(msf_file_path, result)

    if mut_res == 'all':
        return result

    res_index = np.where(result[:, 0] == mut_res)[0]
    if res_index.size == 0:
        raise ValueError(f"Residue number {mut_res} not found in the computed MSF data.")
    res_index = res_index[0]
    res_fluc = result[res_index, 1]
    res_zscore = result[res_index, 2]
    res_norm_zscore = result[res_index, 3]

    return res_fluc, res_zscore, res_norm_zscore


def cal_effecsens(pdb_path, chain, mut_res, n_modes, enm, **kwargs):
    """
    Computes the absolute value, Z-Score, and normalized Z-Score of the effectiveness and 
    sensitivity of a residue based on PRS.

    Parameters:
    pdb_path (str): Path to the PDB file.
    chain (str): Chain identifier.
    mut_res (int or str): Residue number or 'all' for all residues.
    enm (str): Elastic Network Model to use ('ANM' or 'GNM').
    n_modes (int or str): Number of modes to consider ('all', 'soft', or a fraction of total modes).
    kwargs: Additional arguments like 'frac_var'.

    Returns:
    numpy.ndarray or tuple: PRS values for all residues or a specific residue.
    """
    supported_arguments = {'pdb_path', 'chain', 'mut_res', 'enm', 'n_modes', 'frac_var'}
    for key in kwargs:
        if key not in supported_arguments:
            raise ValueError(f"Unsupported argument: {key}")

    pdb_code = re.search(r'([^/]+)\.pdb$', pdb_path)
    if not pdb_code:
        raise ValueError(f"Invalid PDB file path: {pdb_path}")

    pdb_name = f"{pdb_code.group(1)}_{chain}"
    prs_file_path = os.path.join(PRS_folder, f"{pdb_name}_{enm}_{n_modes}_all.npy")

    # Check if the PRS file exists
    if os.path.exists(prs_file_path):
        result = np.load(prs_file_path)
        if mut_res == 'all':
            return result
        else:
            res_index = np.where(result[:, 0] == mut_res)[0]
            if res_index.size == 0:
                raise ValueError(f"Residue number {mut_res} not found in the PRS file.")
            res_index = res_index[0]
            res_effec = result[res_index, 1]
            res_zscore_effec = result[res_index, 2]
            res_norm_zscore_effec = result[res_index, 3]
            res_sens = result[res_index, 4]
            res_zscore_sens = result[res_index, 5]
            res_norm_zscore_sens = result[res_index, 6]
            return res_effec, res_zscore_effec, res_norm_zscore_effec, res_sens, res_zscore_sens, res_norm_zscore_sens

    sel_enm = None

    if enm == 'ANM':
        anm_file_path = os.path.join(ANM_folder, f"{pdb_name}.anm.npz")
        if os.path.exists(anm_file_path):
            sel_enm = pr.loadModel(anm_file_path)
        else:
            pdb = pr.parsePDB(pdb_path, chain=chain, subset='ca')
            sel_enm = pr.ANM()
            sel_enm.buildHessian(pdb)
            sel_enm.calcModes(n_modes='all')
            pr.saveModel(sel_enm, filename=anm_file_path)

    elif enm == 'GNM':
        gnm_file_path = os.path.join(GNM_folder, f"{pdb_name}.gnm.npz")
        if os.path.exists(gnm_file_path):
            sel_enm = pr.loadModel(gnm_file_path)
        else:
            pdb = pr.parsePDB(pdb_path, chain=chain, subset='ca')
            sel_enm = pr.GNM()
            sel_enm.buildKirchhoff(pdb)
            sel_enm.calcModes(n_modes='all')
            pr.saveModel(sel_enm, filename=gnm_file_path)

    else:
        raise ValueError(f"Unsupported ENM type: {enm}")

    tot_modes = sel_enm.numModes()
    if isinstance(n_modes, str):
        if n_modes == 'all':
            prs_mtrx, effec, sens = pr.calcPerturbResponse(sel_enm)
        elif n_modes == 'soft':
            num_modes = math.ceil(0.02 * tot_modes)
            prs_mtrx, effec, sens = pr.calcPerturbResponse(sel_enm[:num_modes])
        else:
            raise ValueError(f"Invalid n_modes argument: {n_modes}")
    elif isinstance(n_modes, (int, float)):
        if n_modes > 0 and n_modes < 1:
            num_modes = math.ceil(n_modes * tot_modes)
            prs_mtrx, effec, sens = pr.calcPerturbResponse(sel_enm[:num_modes])
        elif n_modes >= 1 and n_modes <= tot_modes:
            prs_mtrx, effec, sens = pr.calcPerturbResponse(sel_enm[:int(n_modes)])
        else:
            raise ValueError(f"Invalid n_modes argument: {n_modes}")
    elif 'frac_var' in kwargs:
        frac_var = kwargs['frac_var']
        if 0 < frac_var <= 1:
            fra_var = pr.calcFractVariance(sel_enm)
            cumulative_sum = np.cumsum(fra_var)
            min_mode = 1 + np.argmax(cumulative_sum >= frac_var)
            prs_mtrx, effec, sens = pr.calcPerturbResponse(sel_enm[:min_mode])
        else:
            raise ValueError(f"Invalid frac_var argument: {frac_var}")
    else:
        prs_mtrx, effec, sens = pr.calcPerturbResponse(sel_enm)

    pdb = pr.parsePDB(pdb_path, chain=chain, subset='ca')
    resnums = pdb.getResnums()

    max_effec, max_sens = effec.max(), sens.max()
    min_effec, min_sens = effec.min(), sens.min()
    mean_effec, mean_sens = effec.mean(), sens.mean()
    std_effec, std_sens = effec.std(), sens.std()
    zscore_min_effec, zscore_min_sens = (min_effec - mean_effec) / std_effec, (min_sens - mean_sens) / std_sens
    zscore_max_effec, zscore_max_sens = (max_effec - mean_effec) / std_effec, (max_sens - mean_sens) / std_sens

    res_zscore_effec = (effec - mean_effec) / std_effec
    res_norm_zscore_effec = (res_zscore_effec - zscore_min_effec) / (zscore_max_effec - zscore_min_effec)
    res_zscore_sens = (sens - mean_sens) / std_sens
    res_norm_zscore_sens = (res_zscore_sens - zscore_min_sens) / (zscore_max_sens - zscore_min_sens)

    result = np.vstack((resnums, effec, res_zscore_effec, res_norm_zscore_effec, sens, res_zscore_sens, res_norm_zscore_sens)).T
    np.save(prs_file_path, result)

    if mut_res == 'all':
        return result

    res_index = np.where(result[:, 0] == mut_res)[0]
    if res_index.size == 0:
        raise ValueError(f"Residue number {mut_res} not found in the computed PRS data.")
    res_index = res_index[0]
    res_effec = result[res_index, 1]
    res_zscore_effec = result[res_index, 2]
    res_norm_zscore_effec = result[res_index, 3]
    res_sens = result[res_index, 4]
    res_zscore_sens = result[res_index, 5]
    res_norm_zscore_sens = result[res_index, 6]

    return res_effec, res_zscore_effec, res_norm_zscore_effec, res_sens, res_zscore_sens, res_norm_zscore_sens


def cal_hingedist(pdb_path, chain, mut_res, n_modes, enm='GNM', **kwargs): 
    """
    Computes the sequence distance of the mutated residue from the nearest hinge site across n_modes.
    The default value of n_modes is 5, and the default ENM is GNM.

    Parameters:
    pdb_path (str): Path to the PDB file.
    chain (str): Chain identifier.
    mut_res (int or str): Residue number or 'all' for all residues.
    enm (str): Elastic Network Model to use ('ANM' or 'GNM').
    n_modes (int or str): Number of modes to consider ('all', 'soft', or a fraction of total modes).
    kwargs: Additional arguments like 'frac_var'.

    Returns:
    int or np.ndarray: If mut_res is a residue number, returns 1 if seq_dist_min < 2, otherwise 0.
                       If mut_res is 'all', returns a numpy array with residue numbers and binary values.
    """
    supported_arguments = {'pdb_path', 'chain', 'mut_res', 'enm', 'n_modes'}
    for key in kwargs:
        if key not in supported_arguments:
            raise ValueError(f"Unsupported argument: {key}")

    pdb_code = re.search(r'([^/]+)\.pdb$', pdb_path)
    if not pdb_code:
        raise ValueError(f"Invalid PDB file path: {pdb_path}")

    pdb_name = f"{pdb_code.group(1)}_{chain}"
    hinge_file_path = os.path.join(HINGE_folder, f"{pdb_name}_{enm}_{n_modes}_all.npy")

    # Check if the hinge file exists
    if os.path.exists(hinge_file_path):
        result = np.load(hinge_file_path)
        if mut_res == 'all':
            return result
        else:
            res_index = np.where(result[:, 0] == mut_res)[0]
            if res_index.size == 0:
                raise ValueError(f"Residue number {mut_res} not found in the hinge file.")
            return int(result[res_index[0], 1])

    sel_enm = None

    pdb = pr.parsePDB(pdb_path, chain=chain, subset='ca')

    if enm == 'ANM':
        anm_file_path = os.path.join(ANM_folder, f"{pdb_name}.anm.npz")
        if os.path.exists(anm_file_path):
            sel_enm = pr.loadModel(anm_file_path)
        else:
            sel_enm = pr.ANM()
            sel_enm.buildHessian(pdb)
            sel_enm.calcModes(n_modes='all')
            pr.saveModel(sel_enm, filename=anm_file_path)

    elif enm == 'GNM':
        gnm_file_path = os.path.join(GNM_folder, f"{pdb_name}.gnm.npz")
        if os.path.exists(gnm_file_path):
            sel_enm = pr.loadModel(gnm_file_path)
        else:
            sel_enm = pr.GNM()
            sel_enm.buildKirchhoff(pdb)
            sel_enm.calcModes(n_modes='all')
            pr.saveModel(sel_enm, filename=gnm_file_path)

    else:
        raise ValueError(f"Unsupported ENM type: {enm}")

    tot_modes = sel_enm.numModes()
    if 'n_modes' in kwargs:
        n_modes = kwargs['n_modes']
        if isinstance(n_modes, str):
            if n_modes == 'all':
                num_modes = tot_modes
            elif n_modes == 'soft':
                num_modes = math.ceil(0.02 * tot_modes)
            else:
                raise ValueError(f"Invalid n_modes argument: {n_modes}")
        elif isinstance(n_modes, (int, float)):
            if n_modes > 0 and n_modes < 1:
                num_modes = math.ceil(n_modes * tot_modes)
            elif n_modes >= 1 and n_modes <= tot_modes:
                num_modes = int(n_modes)
            else:
                raise ValueError(f"Invalid n_modes argument: {n_modes}")
        else:
            raise ValueError(f"Unsupported n_modes type: {type(n_modes)}")
    else:
        num_modes = 5

    hinges = pr.calcHinges(sel_enm[:num_modes])
    hinge_res = pdb.getResnums()[hinges]

    result = []
    for res in pdb.getResnums():
        seq_dist = np.abs(res - hinge_res).min()
        result.append([res, 1 if seq_dist < 2 else 0])

    result = np.array(result)
    np.save(hinge_file_path, result)

    if mut_res == 'all':
        return result
    else:
        res_index = np.where(result[:, 0] == mut_res)[0]
        if res_index.size == 0:
            raise ValueError(f"Residue number {mut_res} not found.")
        return int(result[res_index[0], 1])


def cal_meancorr(pdb_path, chain, mut_res, n_modes, enm, **kwargs):
    """
    Computes the absolute value, Z-Score, and normalized Z-Score of the mean correlation of 
    a residue with all other residues based on the chosen ENM and chosen number of modes.

    Parameters:
    pdb_path (str): Path to the PDB file.
    chain (str): Chain identifier.
    mut_res (int or str): Residue number or 'all' for all residues.
    enm (str): Elastic Network Model to use ('ANM' or 'GNM').
    n_modes (int or str): Number of modes to consider ('all', 'soft', or a fraction of total modes).
    kwargs: Additional arguments like 'frac_var'.

    Returns:
    numpy.ndarray or tuple: Correlation values for all residues or a specific residue.
    """
    supported_arguments = {'pdb_path', 'chain', 'mut_res', 'enm', 'n_modes', 'frac_var'}
    for key in kwargs:
        if key not in supported_arguments:
            raise ValueError(f"Unsupported argument: {key}")

    pdb_code = re.search(r'([^/]+)\.pdb$', pdb_path)
    if not pdb_code:
        raise ValueError(f"Invalid PDB file path: {pdb_path}")

    pdb_name = f"{pdb_code.group(1)}_{chain}"
    corr_file_path = os.path.join(CORR_folder, f"{pdb_name}_{enm}_{n_modes}_all.npy")

    # Check if the correlation file exists
    if os.path.exists(corr_file_path):
        result = np.load(corr_file_path)
        if mut_res == 'all':
            return result
        else:
            res_index = np.where(result[:, 0] == mut_res)[0]
            if res_index.size == 0:
                raise ValueError(f"Residue number {mut_res} not found in the correlation file.")
            res_index = res_index[0]
            res_avg_corr = result[res_index, 1]
            res_zscore = result[res_index, 2]
            res_norm_zscore = result[res_index, 3]
            return res_avg_corr, res_zscore, res_norm_zscore

    sel_enm = None

    if enm == 'ANM':
        anm_file_path = os.path.join(ANM_folder, f"{pdb_name}.anm.npz")
        if os.path.exists(anm_file_path):
            sel_enm = pr.loadModel(anm_file_path)
        else:
            pdb = pr.parsePDB(pdb_path, chain=chain, subset='ca')
            sel_enm = pr.ANM()
            sel_enm.buildHessian(pdb)
            sel_enm.calcModes(n_modes='all')
            pr.saveModel(sel_enm, filename=anm_file_path)

    elif enm == 'GNM':
        gnm_file_path = os.path.join(GNM_folder, f"{pdb_name}.gnm.npz")
        if os.path.exists(gnm_file_path):
            sel_enm = pr.loadModel(gnm_file_path)
        else:
            pdb = pr.parsePDB(pdb_path, chain=chain, subset='ca')
            sel_enm = pr.GNM()
            sel_enm.buildKirchhoff(pdb)
            sel_enm.calcModes(n_modes='all')
            pr.saveModel(sel_enm, filename=gnm_file_path)

    else:
        raise ValueError(f"Unsupported ENM type: {enm}")

    tot_modes = sel_enm.numModes()
    if isinstance(n_modes, str):
        if n_modes == 'all':
            cross_corr = pr.calcCrossCorr(sel_enm)
        elif n_modes == 'soft':
            num_modes = math.ceil(0.02 * tot_modes)
            cross_corr = pr.calcCrossCorr(sel_enm[:num_modes])
        else:
            raise ValueError(f"Invalid n_modes argument: {n_modes}")
    elif isinstance(n_modes, (int, float)):
        if n_modes > 0 and n_modes < 1:
            num_modes = math.ceil(n_modes * tot_modes)
            cross_corr = pr.calcCrossCorr(sel_enm[:num_modes])
        elif n_modes >= 1 and n_modes <= tot_modes:
            cross_corr = pr.calcCrossCorr(sel_enm[:int(n_modes)])
        else:
            raise ValueError(f"Invalid n_modes argument: {n_modes}")
    elif 'frac_var' in kwargs:
        frac_var = kwargs['frac_var']
        if 0 < frac_var <= 1:
            fra_var = pr.calcFractVariance(sel_enm)
            cumulative_sum = np.cumsum(fra_var)
            min_mode = 1 + np.argmax(cumulative_sum >= frac_var)
            cross_corr = pr.calcCrossCorr(sel_enm[:min_mode])
        else:
            raise ValueError(f"Invalid frac_var argument: {frac_var}")
    else:
        cross_corr = pr.calcCrossCorr(sel_enm)

    pdb = pr.parsePDB(pdb_path, chain=chain, subset='ca')
    resnums = pdb.getResnums()

    avg_corr = np.abs(cross_corr).mean(axis=1)

    max_avg_corr = avg_corr.max()
    min_avg_corr = avg_corr.min()
    mean_avg_corr = avg_corr.mean()
    std_avg_corr = avg_corr.std()
    zscore_min = (min_avg_corr - mean_avg_corr) / std_avg_corr
    zscore_max = (max_avg_corr - mean_avg_corr) / std_avg_corr

    res_zscore = (avg_corr - mean_avg_corr) / std_avg_corr
    res_norm_zscore = (res_zscore - zscore_min) / (zscore_max - zscore_min)

    result = np.vstack((resnums, avg_corr, res_zscore, res_norm_zscore)).T
    np.save(corr_file_path, result)

    if mut_res == 'all':
        return result

    res_index = np.where(result[:, 0] == mut_res)[0]
    if res_index.size == 0:
        raise ValueError(f"Residue number {mut_res} not found in the computed correlation data.")
    res_index = res_index[0]
    res_avg_corr = result[res_index, 1]
    res_zscore = result[res_index, 2]
    res_norm_zscore = result[res_index, 3]

    return res_avg_corr, res_zscore, res_norm_zscore


#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison

import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import pickle
from numba import jit
import os
import sys


@jit(nopython=True)
def p_generate(val):
    '''
    Iteratively grabbing scaffold pairs from a list
    '''
    p1 = val*2
    p2 = p1+1
    return p1,p2


@jit(nopython=True)
def prob_compare(prob, rate_cutoff):
    '''
    Filter machine learning probabilities
    '''
    prob_check = False
    if prob >= rate_cutoff:
        prob_check = True
    return prob_check


@jit(nopython=True)
def prob_generate_hybrid(ET, NN):
    '''
    Combine ET and NN probabilities
    '''
    prob = (ET + NN)/2
    return prob


@jit(nopython=True)
def cd_compare(p, c, network):
    '''
    Generate the probability cutoff and compare to maximum threshold
    '''
    r = c/p
    r_check = False
    if r <= network:
        r_check = True
    return r_check, r


def machine_stuff(metrics, presets, model_method, new_pairs, cohen_list, iterations):
    '''
    Main wrapper for machine learning and cutoff checks
    '''
    anno = pd.read_csv(metrics, sep='\t', header=0)

    base_scripts = os.path.dirname(os.path.abspath(__file__)).replace('/scripts','')
    base_bin = os.path.dirname(os.path.abspath(__file__)).replace('/bin','')
    if os.path.exists(f'{base_scripts}/models/vRhyme_machine_model_NN.sav'):
        ml_path = f'{base_scripts}/models/'
    else: # pip
        ml_path = f'{base_bin}/lib/python{sys.version_info.major}.{sys.version_info.minor}/site-packages/vRhyme/models/'

    if model_method == 'hybrid' or model_method == 'ET':
        with open(f'{ml_path}vRhyme_machine_model_ET.sav', 'rb') as read_model_ET:
            model_ET = pickle.load(read_model_ET)

        probs_ET = model_ET.predict_proba(anno)
        model_ET = None # close
        probs_ET = [p[1] for p in probs_ET] # just the probability of "same" genome according to model
        length = len(probs_ET)

    if model_method == 'hybrid' or model_method == 'NN':
        with open(f'{ml_path}vRhyme_machine_model_NN.sav', 'rb') as read_model_NN:
            model_NN = pickle.load(read_model_NN)

        probs_NN = model_NN.predict_proba(anno)
        model_NN = None # close
        probs_NN = [p[1] for p in probs_NN] # just the probability of "same" genome according to model
        length = len(probs_NN)

    anno = None  # close
    net_full = {}

    net_full = {}
    ml_ratelist = []
    cohen_ratelist = []
    net_ratelist = []
    for i in range(0,iterations):
        ml_ratelist.append(presets[i][0])
        cohen_ratelist.append(presets[i][1])
        net_ratelist.append(presets[i][2])
        net_full[i] = []

    for n in range(0,length):

        checklist = []

        if model_method == 'hybrid':
            prob_ET = probs_ET[n]
            prob_NN = probs_NN[n]
            prob = prob_generate_hybrid(prob_ET, prob_NN)
            for i in range(0,iterations):
                rate_cutoff = ml_ratelist[i]
                prob_check = prob_compare(prob, rate_cutoff)
                checklist.append(prob_check)

        elif model_method == 'ET':
            prob = probs_ET[n]
            for i in range(0,iterations):
                rate_cutoff = ml_ratelist[i]
                prob_check = prob_compare(prob, rate_cutoff)
                checklist.append(prob_check)

        elif model_method == 'NN':
            prob = probs_NN[n]
            for i in range(0,iterations):
                rate_cutoff = ml_ratelist[i]
                prob_check = prob_compare(prob, rate_cutoff)
                checklist.append(prob_check)

        if True in checklist:
            p1,p2 = p_generate(n)
            p1 = new_pairs[p1]
            p2 = new_pairs[p2]
            cd = cohen_list[n]
            for i in range(0,iterations):
                if checklist[i] == True:
                    rate_cutoff = cohen_ratelist[i]
                    if cd <= rate_cutoff:
                        rate_cutoff = net_ratelist[i]
                        ratio_check, ratio = cd_compare(prob, cd, rate_cutoff)
                        if ratio_check == True:
                            net_full[i].append((p1,p2,ratio))

    new_pairs = None
    cohen_list = None
    probs_ET = None
    probs_NN = None

    return net_full



#
#
#

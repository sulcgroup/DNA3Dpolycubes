from nupack import *
import random as rnd 
import sys
import os
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from copy import deepcopy
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from scipy.special import comb
from itertools import combinations
from numba import njit


complement_translation = str.maketrans('ACGT', 'TGCA')


def get_compseq(seq):
    return seq.translate(complement_translation)[::-1]


# load seq input and add its complementary pair 
# return a list with all the seq
def load_ensemble(fname):
    seqs = []
    with open(fname) as inf:
        for line in inf:
            if len(line.strip().split()) >= 1:
                seq = line.strip().split()[0]
                #seq = seq[::-1] #uncomment to make it take 3-5 input 
                comp = get_compseq(seq)

                if comp == seq: # don't want homopolymers 
                    continue

                if seq not in seqs:
                    #seqs.append(seq)
                    seqs.insert(0,seq)
                if comp not in seqs:
                    seqs.append(get_compseq(seq))

    # should filter the ensemble here if we want
    rnd.shuffle(seqs)
    
    return seqs


def add_spacers(seqs, spacer):
    seqs_spaced = []
    for seq in seqs:
        seqs_spaced.append(spacer + seq)
    return seqs_spaced


def filter_out_seqs_with_secondary_struct(seqs, my_model, spacer, filter_struct, constraint_seqs):
    mfe_structures = [mfe(strands=[seq], model=my_model) for seq in seqs]
    non_seqs = [seq for seq,energy in zip(seqs, mfe_structures) if str(energy[0].structure).count('(') > filter_struct]
    
    if spacer is not None:
        non_seqs_no_spacer = [seq.replace(spacer, '') for seq in non_seqs]
        comp_non_seqs_no_spacer = [get_compseq(seq) for seq in non_seqs_no_spacer]
        comp_non_seqs = [f'{spacer}{seq}' for seq in comp_non_seqs_no_spacer]
    else:
        comp_non_seqs = [get_compseq(seq) for seq in non_seqs]
        
    has_sec = set(non_seqs + comp_non_seqs)
    print(f'{len(has_sec)} sequences removed out of {len(seqs)} after secondary structure filter')
    
    seqs_new = [seq for seq in seqs if seq not in has_sec]
    
    if len(constraint_seqs) > 0:
        for seq in constraint_seqs:
            if spacer is not None:
                comp_seq = f"{spacer}{get_compseq(seq)}"
                seq = f"{spacer}{seq}"
            else:
                comp_seq = get_compseq(seq)
                
            if seq not in seqs_new:
                seqs_new.append(seq)
            if comp_seq not in seqs_new:
                seqs_new.append(comp_seq)

    
    print(f'{len(seqs_new)} sequecnes remaining')
    
    return seqs_new


# compute the dGs for ant seq pair, complementary or not
def fill_dGs(seqs, my_model, dG_file, spacer, filter_struct, constraint_seqs):
    if not os.path.exists(dG_file):

        if spacer is not None:
            print(f'Adding {spacer} 5` spacer')
            old_seqs = deepcopy(seqs)
            seqs = add_spacers(seqs, spacer)

        seqs = filter_out_seqs_with_secondary_struct(seqs, my_model, spacer, filter_struct, constraint_seqs)
        
        dGs = pd.DataFrame(np.zeros((len(seqs),len(seqs)),dtype=float),index=seqs,columns=seqs) # a large dict to store all the mfes of the seqs in the library 
        dGd_inter = {}
        interim_array = np.zeros((len(seqs),len(seqs)),dtype=float)
        complements = {} # a dict with all the complementary seq pairs 
        for i, seqA in enumerate(dGs.columns):
            for j, seqB in enumerate(dGs.index):
                if i >= j : #avoid doing it for twice, j -> i == i -> j
                    dGd_inter[tuple( [seqA,seqB]  )] = (i, j)

        strands = {seq:Strand(seq, name=seq) for seq in seqs}
        # my_complexes = [Complex([strands[seqA], strands[seqB]]) for seqA, seqB in dGd_inter.keys()]
        dGd_keys = list(dGd_inter.keys())
        if len(dGd_keys) > 5e5:
            n_sections_of_5_million = int(len(dGd_keys)//5e5)
        else:
            n_sections_of_5_million = 1
        
        dGd_chunks = np.array_split(dGd_keys, n_sections_of_5_million)
        # my_complexes_chunks = np.array_split(my_complexes, n_sections_of_5_million)
        
        total_elements = len(dGd_keys)
        pbar = tqdm(total=total_elements, desc='Computing dGs for all complexes')
        
        for dGd_chunk in dGd_chunks:
            my_complexes = [Complex([strands[seqA], strands[seqB]]) for seqA, seqB in dGd_chunk]
            complex_results = complex_analysis(complexes=my_complexes, model=my_model, compute=['pfunc'])

            for idx, (seqA, seqB) in enumerate(dGd_chunk):
                i, j =  dGd_inter[(seqA, seqB)]
                interim_array[i, j] = complex_results[my_complexes[idx]].free_energy
                interim_array[j, i] = complex_results[my_complexes[idx]].free_energy
            
            pbar.update(len(my_complexes))
            
        if spacer is not None:
            seqs = [seq.replace(spacer, "") for seq in seqs]

        dGs_filled = pd.DataFrame(interim_array, index=seqs,columns=seqs) # a large dict to store all the mfes of the seqs in the library 
        dGs_filled.to_pickle(dG_file)
        
    else:
        print(f'Loading dGs from {dG_file}')
        dGs_filled = pd.read_pickle(dG_file)
        
    return dGs_filled

  
def ddG_optimization_algorithm(dGs_filled, num_seq_output, optimization_steps, score_no_opt, constraint_seqs=[]):
    seqs_to_drop = list(dGs_filled.index[np.isinf(dGs_filled).any(axis=1)])
    seqs_to_drop_comp = [get_compseq(seq) for seq in seqs_to_drop]
    drop_list = seqs_to_drop + seqs_to_drop_comp
    
    if any(c_seq in drop_list for c_seq in constraint_seqs):
        
        print(f"For some reason, the constraint sequences dG calculation has failed. Reach out to Matthew to debug.")
        raise Exception("Reach out to Matthew to debug.")
    
    dGs_filled = dGs_filled.drop(index=drop_list,columns=drop_list)
    
    ddG_df = dGs_filled.copy()
    ddG_df = set_column_complement_to_inf(ddG_df)
    
    if score_no_opt:
        optimal_sequence_set = ddG_df.columns
        maximum_min_interaction_distance = ddG_df.min().min()
        
        return optimal_sequence_set, maximum_min_interaction_distance
    
    global ddG_array
    ddG_array = ddG_df.to_numpy()
    
    index_map = {seq: i for i, seq in enumerate(ddG_df.index)}
    complement_index_map = {index_map[seq]: index_map[get_compseq(seq)] for seq in ddG_df.index}
    
    indexes = list(index_map.keys())

    
    maximum_min_interaction_distances = []
    optimal_sequence_sets = []
    optimal_sequence_dfs = []
    
    n_cpus = min(len(os.sched_getaffinity(0)), len(ddG_df.columns))
    columns = list(ddG_df.columns)
    
    if len(constraint_seqs) > 0:
        columns = [col for col in columns if col not in constraint_seqs]
        constraint_indices = [index_map[seq] for seq in constraint_seqs if seq in index_map]
    else:
        constraint_indices = []

    # args = np.array([(complement_index_map, index_map, indexes, column, num_seq_output, optimization_steps) for column in columns], dtype='object')
    args = np.array([
        (complement_index_map, index_map, indexes, column, num_seq_output, optimization_steps, constraint_indices)
        for column in columns
    ], dtype='object')
    arg_chunks = np.array_split(args, n_cpus)
    
    del ddG_df
    # del args
    del dGs_filled
    
    with ProcessPoolExecutor(max_workers=n_cpus) as executor:
        results = list(tqdm(executor.map(optimization_iteration_parallel, arg_chunks), total=len(arg_chunks), desc='Optimizing sequence selection'))
    
    results = [inner for outer in results for inner in outer]
    optimal_sequence_sets = [result[0] for result in results]
    maximum_min_interaction_distances = [result[1] for result in results]
    
    # for arg in tqdm(args, total=len(args), desc='Optimizing sequence selection'):
    #     results = optimization_iteration(ddG_array, *arg)
    #     optimal_sequence_sets.append(results[0])
    #     maximum_min_interaction_distances.append(results[1])
    
    maximum_min_interaction_distance = np.max(maximum_min_interaction_distances)
    optimal_sequence_set = optimal_sequence_sets[np.argmax(maximum_min_interaction_distances)]
            
    return optimal_sequence_set, maximum_min_interaction_distance


def set_column_complement_to_inf(ddG_df):
    # Iterate through each column
    complements = [get_compseq(column) for column in ddG_df.columns]
    comp_vals = np.array([ddG_df.loc[complement, column] for complement, column in zip(complements, ddG_df.columns)])
    array = ddG_df.to_numpy()
    result_array = array - comp_vals[:, np.newaxis]
    result_array[result_array == 0] = np.inf
    result_df = pd.DataFrame(result_array, columns=ddG_df.columns, index=ddG_df.index)

    return result_df
    

def optimization_iteration_parallel(arg_chunk):
    all_results = []
    for args in arg_chunk:
        complement_index_map, index_map, indexes, column, num_seq_output, optimization_steps, constraint_indices = args
        all_results.append(optimization_iteration(ddG_array, complement_index_map, index_map, indexes, column, num_seq_output, optimization_steps, constraint_indices))
    
    return all_results
    

def optimization_iteration(ddG_array, complement_index_map, index_map, indexes, column, num_seq_output, optimization_steps, constraint_indices):
    column_idx = index_map[column]
    maximum_min_interaction_distance = -np.inf
    optimal_sequence_set = None
    
    sorted_df = ddG_array[index_map[column], :].argsort()[::-1]
    sorted_indexes = []
    for i in sorted_df:
        if (i == column_idx) or (i == complement_index_map[column_idx]):
            continue
        if i in constraint_indices:
            continue
        if get_compseq(indexes[i]) in sorted_indexes:
            continue
        else:
            sorted_indexes.append(indexes[i])
            
    sorted_df = np.array([index_map[seq] for seq in sorted_indexes])
        
    # sorted_indexes = [index for index in sorted_indexes if index != get_compseq(column)]
    old_best_indexes = sorted_indexes[:num_seq_output-1]
    
    
    top_n_seqs = [column] + old_best_indexes[:num_seq_output-1]
        
    top_n_comp = [get_compseq(seq) for seq in top_n_seqs]
    
    
    old_seqs_plus_comp = top_n_seqs + top_n_comp
    
    seq_indices = np.array([index_map[seq] for seq in old_seqs_plus_comp])
    if len(constraint_indices) > 0:
        seq_indices = np.concatenate((constraint_indices, seq_indices))
        
    mesh = np.ix_(seq_indices, seq_indices)
    
    data_array = ddG_array[mesh]
    # data_array[data_array == 0] = np.inf
    
    reached_max_counter = 0
    for step in range(optimization_steps):

        
        data_array, worst_seq, worst_seq_comp = remove_rows_cols_based_on_condition(data_array, column_idx, seq_indices, complement_index_map, constraint_indices)
        removed_index = np.array([seq_indices[worst_seq], seq_indices[worst_seq_comp]])
        removed_array_index = np.array([worst_seq, worst_seq_comp])
        
        remaining_indices, remaining_indices_complement = get_remaining_indices_and_complements(sorted_df, seq_indices, complement_index_map)
        
        data_array, seq_indices, best_min_value = maximize_min_value(ddG_array, data_array, remaining_indices, remaining_indices_complement, removed_array_index, removed_index, seq_indices)
                
        if best_min_value > maximum_min_interaction_distance:
            maximum_min_interaction_distance = best_min_value
            optimal_sequence_set = [indexes[i] for i in seq_indices]
            reached_max_counter = 0
        
        elif best_min_value == maximum_min_interaction_distance:
            reached_max_counter +=1
        
        if reached_max_counter == optimization_steps // 4:
            break
    if len(constraint_indices) > 0:
        assert len(list(set(optimal_sequence_set))) == num_seq_output * 2 + len(constraint_indices)
    else:    
        assert len(list(set(optimal_sequence_set))) == num_seq_output * 2
        
    return (optimal_sequence_set, maximum_min_interaction_distance)


def remove_rows_cols_based_on_condition(data_array, column_idx, seq_indices, complement_index_map, constraint_indices):
    increment_worst_seq_comp = False
        
    while True:
        min_index_flat = data_array.argmin()
        
        # # Convert flat index to row and column indices
        min_row, min_col = np.unravel_index(min_index_flat, data_array.shape)
        
        if min_row < min_col:
            worst_seq = min_row
        else:
            worst_seq = min_col
        
        worst_seq_comp = np.where(seq_indices == complement_index_map[seq_indices[worst_seq]])[0][0]
        
        worst_not_column = seq_indices[worst_seq] != column_idx
        worst_not_column_comp = seq_indices[worst_seq_comp] != column_idx
        
        if len(constraint_indices) > 0:
            worst_not_constraint = seq_indices[worst_seq] not in constraint_indices
            worst_not_constraint_comp = seq_indices[worst_seq_comp] not in constraint_indices
            
        if worst_not_column and worst_not_column_comp:
            if len(constraint_indices) > 0:
                if worst_not_constraint and worst_not_constraint_comp:
                    break
            else:
                break
        
        
        if increment_worst_seq_comp is False:
            data_array_bak = data_array.copy()
            increment_worst_seq_comp = True
            
        data_array[worst_seq, :] = 100
        data_array[:, worst_seq] = 100
        data_array[worst_seq_comp, :] = 100
        data_array[:, worst_seq_comp] = 100
     
    if increment_worst_seq_comp:
        data_array = data_array_bak
    
    # removed_index = [seq_indices[worst_seq], seq_indices[worst_seq_comp]]
    # assert removed_index[0] == complement_index_map[removed_index[1]]

    return data_array, worst_seq, worst_seq_comp    
   
 
def get_remaining_indices_and_complements(sorted_df, seq_indices, complement_index_map, sample_size=100):
    remaining_indices = sorted_df[~np.isin(sorted_df, seq_indices)]
    sample_size = min([sample_size, len(remaining_indices)])
    # p = compute_exponential_decay(remaining_indices)
    remaining_indices = np.random.choice(remaining_indices, size=sample_size, replace=False)#, p=p)
    remaining_indices_complement = np.array([complement_index_map[index] for index in remaining_indices])
    return remaining_indices, remaining_indices_complement


# @njit(fastmath=True)
# def compute_exponential_decay(remaining_indices):
#     i = np.arange(remaining_indices.size)
#     p = np.exp(-i / remaining_indices.size )
#     p = p / p.sum()
#     return p
 
 
@njit(fastmath=True)
def maximize_min_value(ddG_array, data_array, remaining_indices, remaining_indices_complement, removed_array_index, removed_index, seq_indices):
    best_row_index = -1
    best_min_value = -np.inf
 
    replace_site = np.where(seq_indices == removed_index[0])[0][0]
    replace_comp_site = np.where(seq_indices == removed_index[1])[0][0]
    
    rows = np.zeros((remaining_indices.size, 4, seq_indices.size))
    for idx in np.arange(remaining_indices.size):
        index = remaining_indices[idx]
        index_complement = remaining_indices_complement[idx]
        
        seq_indices[replace_site] = index
        seq_indices[replace_comp_site] = index_complement

        rows[idx, 0, :] = ddG_array[index, seq_indices]
        rows[idx, 1, :] = ddG_array[seq_indices, index]

        rows[idx, 2, :] = ddG_array[index_complement, seq_indices]
        rows[idx, 3, :] = ddG_array[seq_indices, index_complement]

    for idx in np.arange(remaining_indices.size):
        temp_min = rows[idx].min()
        if temp_min > best_min_value:
            best_min_value = temp_min
            best_row = idx
    # all_mins = rows.min(axis=(1, 2))
    # best_row = np.where(all_mins == all_mins.max())[0][0]
    best_row_index = remaining_indices[best_row]
    
    best_comp_index = remaining_indices_complement[np.where(remaining_indices == best_row_index)[0][0]]
    
    seq_indices[replace_site] = best_row_index
    seq_indices[replace_comp_site] = best_comp_index
    
    data_array[:, replace_site] = ddG_array[seq_indices, best_row_index]
    data_array[replace_site, :] = ddG_array[best_row_index, seq_indices]
        
    data_array[:, replace_comp_site] = ddG_array[seq_indices, best_comp_index]
    data_array[replace_comp_site, :] = ddG_array[best_comp_index, seq_indices]
        
    best_min_value = data_array.min()
        
    return data_array, seq_indices, best_min_value

  
def evaluate_threshold(pseqs,my_model,output_seq, dGs_filled):
    # pseqs = pseqs + [get_compseq(seq) for seq in pseqs]
    seqs = pseqs.copy()

    mins = {}
    for s in seqs:
        mins[s] = 1e9
    
    
    
    # technically Min should be just the complementary pair..
    for i, s1 in enumerate(seqs):
        for j, s2 in enumerate(seqs):
            if dGs_filled[s1].loc[s2] < mins[s1]:
                mins[s1] = dGs_filled[s1].loc[s2]

    threshold = 0.20
    dt = 1
    while threshold <= 10:

        xs = {}
        for s in seqs:
            xs[s] = 0
        for i, s1 in enumerate(seqs):
            for j, s2 in enumerate(seqs):
                delta = dGs_filled[s1].loc[s2] - mins[s1]
                if (delta < threshold) & (delta != 0):
                    xs[s1] += 1


        smax = seqs[0]
        for s in xs.keys():
            if xs[s] > xs[smax]:
                smax = s
        toremove = smax
    
        if xs[smax] == 0: # when dGs_filled[s1].loc[s2] - mins[s1] = 0 (the same pair, complementary) and < threshold
            threshold += dt
            continue
        
        if len(seqs) == output_seq*2: # to be changed (to 12) if you want 6 seqs (exclude comps) output 
            break

    threshold = np.round(threshold, 10)

    return threshold

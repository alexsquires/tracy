from icet.tools.structure_mapping import map_structure_to_reference
from icet import StructureContainer, CrossValidationEstimator, ClusterSpace
import numpy as np
import pandas as pd

def get_structure_container(cluster_space, structures, energies, verbose = True):
    structure_container = StructureContainer(cluster_space)
    could_not_be_mapped = 0
    warning_raised = 0
    for structure, energy in zip(structures, energies):
        try:
            mapped_structure, info = map_structure_to_reference(structure, cluster_space.primitive_structure)
            if info['warnings'] == []:
                structure_container.add_structure(mapped_structure, properties = {'energy': energy})
            else:
                warning_raised += 1
        except:
            could_not_be_mapped += 1
    if verbose == True:
        print(f'{len(structure_container)} structures succesfully mapped on to primitive cell')
        if warning_raised > 0:
            print(f'{warning_raised} structures raised a warning while mapping, these structures have not been added to the structure container.')
            print('to override this behaviour, perform mapping manually')
        if could_not_be_mapped > 0:
            print(f'{could_not_be_mapped} structures failed to map on to the primitive cell, these structures cannot be added to the structure container')
    return structure_container

def get_fitting_summary(structure_container, fit_method):
    cross_validation_estimator = CrossValidationEstimator(structure_container.get_fit_data(), fit_method = fit_method)
    cross_validation_estimator.train()
    cross_validation_estimator.validate()
    return cross_validation_estimator.summary


def cutoff_convergence(primitive_cell, chemical_symbols, structures, energies, cutoff_template = [], order=2, range=np.arange(0,7), fit_methods = ['rfe']):
    fitting_data = []
    for cutoff in range:
        cutoffs = cutoff_template
        cutoffs[order-2] = cutoff
        cluster_space = ClusterSpace(primitive_cell, cutoffs, chemical_symbols)
        structure_container = get_structure_container(cluster_space, structures, energies)
        for fit_method in fit_methods:
            summary = get_fitting_summary(structure_container, fit_method)
            summary[f'cutoff_{order}'] = cutoff
            fitting_data.append(summary)
    cutoff_df = pd.DataFrame(fitting_data)
    return cutoff_df 
    

def filter_for_unique_structures(structures, cluster_space):
    filtered_structures = []
    cluster_vectors = []
    for structure in structures:
        cluster_vector = cluster_space.get_cluster_vector(structure)
        if list(cluster_vector) not in cluster_vectors:
            filtered_structures.append(structure)
            cluster_vectors.append(list(cluster_vector))
    return filtered_structures


def filter_entries_for_unique_structures(entries, cluster_space):
    filtered_structures = []
    cluster_vectors = []
    for entry in entries:
        structure = AseAtomsAdaptor.get_atoms(entry.structure)
        structure = map_structure_to_reference(structure, cluster_space.primitive_structure)[0]
        cluster_vector = cluster_space.get_cluster_vector(structure)
        if list(cluster_vector) not in cluster_vectors:
        #if np.any(np.all(cluster_vector == cluster_vectors)) == False:
            filtered_structures.append(entry)
            cluster_vectors.append(list(cluster_vector))
    return filtered_structures

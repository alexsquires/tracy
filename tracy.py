from icet.tools.structure_mapping import map_structure_to_reference
from icet import StructureContainer, ClusterSpace
from trainstation import CrossValidationEstimator
import numpy as np
import pandas as pd

from pymatgen.io.ase import AseAtomsAdaptor
from mchammer.calculators import ClusterExpansionCalculator
from mchammer.ensembles import CanonicalEnsemble as CEnsemble
from mchammer.ensembles import CanonicalAnnealing as CAnneal
from icet import ClusterExpansion
import time


def get_structure_container(
    cluster_space: ClusterSpace, structures: list, properties: list
) -> StructureContainer:
    """
    Create a StructureContainer object given ClusterSpace,
    structures and energies, filtering out structures that
    don't map on to the primitive cell defined by the ClusterSpace
    args:
        cluster_space: icet.ClusterSpace
        structures: list of ase.atoms object which matches the order of the list of energies
        properties: list of floats of target properties (assumed to be energies) to be fit
    returns:
        structure_container: icet.ClusterSpace
    """
    structure_container = StructureContainer(cluster_space)
    could_not_be_mapped = 0
    warning_raised = 0
    for structure, target_property in zip(structures, properties):
        mapped_structure, info = map_structure_to_reference(
            structure, cluster_space.primitive_structure
        )
        if info["warnings"] == []:
            structure_container.add_structure(
                mapped_structure, properties={"properties": target_property}
            )
    return structure_container


def get_fitting_summary(structure_container: StructureContainer, fit_method: str):
    """
    get the stats from a trainstation given an icet
    structure container object
    args:
        structure_container: icet.StructureContainer
        fit_method: fit method for the CrossValidationEstimator
    returns:
        cve_summary: (dict) CrossValidationEstimator stats
    """
    if 'rfe+' in fit_method:
        cross_validation_estimator = CrossValidationEstimator(
            structure_container.get_fit_data(),
            fit_method="rfe",
            final_estimator=fit_method.split('+')[1],
        )
    else:
        cross_validation_estimator = CrossValidationEstimator(
            structure_container.get_fit_data(key = 'properties'), fit_method=fit_method
        )
    cross_validation_estimator.train()
    cross_validation_estimator.validate()
    return cross_validation_estimator


def cutoff_convergence(
    cluster_space: ClusterSpace,
    structures: list,
    energies: list,
    cutoff_template: list = [],
    order: float = 2,
    cutoff_range: list = np.arange(0, 7),
    fit_methods: list = ["rfe"],
    directory: str = './'
):
    """
    get convergence with respect to cutoff distance
    args:
        cluster_space: (icet.ClusterSpace)
        structures: list of atoms objects that match the order of the list of energies
        energies: list of target properties as floats (assumed to be energies)
        cutoff_template: template for cutoffs, list of length of that you want to consider orders up to
        order: order for cutoff convergence
        cutoff_range: list of cluster cutoffs to try for the selected order
        fit_methods: list of fit methods to try
    returns:
        cutoff_df: dataframe of fitting information
    """
    fitting_data = []
    for cutoff in cutoff_range:
        cutoffs = cutoff_template
        cutoffs[order - 2] = cutoff
        cluster_space = ClusterSpace(
            cluster_space.primitive_structure, cutoffs, cluster_space.chemical_symbols
        )
        structure_container = get_structure_container(
            cluster_space, structures, energies
        )
        for fit_method in fit_methods:
            summary = get_fitting_summary(structure_container, fit_method)
            summary.write_summary(f"{directory}{fit_method}_{cutoff}_{order}.fs")
    # cutoff_df = pd.DataFrame(fitting_data)


#    return cutoff_df


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
        structure = map_structure_to_reference(
            structure, cluster_space.primitive_structure
        )[0]
        cluster_vector = cluster_space.get_cluster_vector(structure)
        if list(cluster_vector) not in cluster_vectors:
            # if np.any(np.all(cluster_vector == cluster_vectors)) == False:
            filtered_structures.append(entry)
            cluster_vectors.append(list(cluster_vector))
    return filtered_structures


def run_canonical_simulation(args):
    cluster_expansion = ClusterExpansion.read(args["path_to_ce"])
    cluster_expansion = ClusterExpansion.prune()
    start = time.time()
    temperature = args["temperature"]
    calculator = ClusterExpansionCalculator(args["supercell"], cluster_expansion)
    mc = CEnsemble(
        calculator=calculator,
        structure=args["supercell"],  
        temperature=temperature,
        data_container=f"mc_data/{args['temperature']}_{args['n_supercell']}_{args['label']}.dc",
    )
    mc.run(number_of_trial_steps=args["n_steps"])


def run_canonical_annealing_simulation(args):
    cluster_expansion = ClusterExpansion.read(args["path_to_ce"])
    cluster_expansion = ClusterExpansion.prune()
    structure = args["supercell"]
    calculator = ClusterExpansionCalculator(structure, cluster_expansion)
    mc = CAnneal(
        calculator=calculator,
        structure=structure,
        T_start=args["starting_temperature"],
        T_stop=args["finish_temperature"],
        n_steps=args["n_steps"],
        ensemble_data_write_interval=200,
        trajectory_write_interval=200,
        dc_filename=f"mc_data/anneal_{args['label']}.dc",
    )
    mc.run()                                                                                                                      1,1           Top

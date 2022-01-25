from icet.tools.structure_mapping import map_structure_to_reference
from icet import StructureContainer, ClusterSpace
from trainstation import CrossValidationEstimator
import numpy as np
import pandas as pd

from pymatgen.io.ase import AseAtomsAdaptor
from mchammer.calculators import ClusterExpansionCalculator
from mchammer.ensembles import CanonicalEnsemble as CEnsemble
from icet import ClusterExpansion
import time


def get_structure_container(
    cluster_space: ClusterSpace, structures: list, energies: list
) -> StructureContainer:
    """
    Create a StructureContainer object given ClusterSpace,
    structures and energies, filtering out structures that
    don't map on to the primitive cell defined by the ClusterSpace

    args:
        cluster_space: icet.ClusterSpace
        structures: list of ase.atoms object which matches the order of the list of energies
        energies: list of floats of target properties (assumed to be energies) to be fit
    returns:
        structure_container: icet.ClusterSpace

    """
    structure_container = StructureContainer(cluster_space)
    could_not_be_mapped = 0
    warning_raised = 0
    for structure, energy in zip(structures, energies):
        try:
            mapped_structure, info = map_structure_to_reference(
                structure, cluster_space.primitive_structure
            )
            if info["warnings"] == []:
                structure_container.add_structure(
                    mapped_structure, properties={"energy": energy}
                )
        except:
            None
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
    cross_validation_estimator = CrossValidationEstimator(
        structure_container.get_fit_data(), fit_method=fit_method
    )
    cross_validation_estimator.train()
    cross_validation_estimator.validate()
    cve_summary = cross_validation_estimator.summary
    return cve_summary


def cutoff_convergence(
    cluster_space: ClusterSpace,
    structures: list,
    energies: list,
    cutoff_template: list = [],
    order: float = 2,
    cutoff_range: list = np.arange(0, 7),
    fit_methods: list = ["rfe"],
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
        cluster_space.cutoffs = cutoffs
        structure_container = get_structure_container(
            cluster_space, structures, energies
        )
        for fit_method in fit_methods:
            summary = get_fitting_summary(structure_container, fit_method)
            summary[f"cutoff_{order}"] = cutoff
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
    start = time.time()
    temperature = args["temperature"]
    calculator = ClusterExpansionCalculator(args["supercell"], cluster_expansion)
    mc = CEnsemble(
        calculator=calculator,
        structure=args["supercell"],
        ensemble_data_write_interval=200,
        trajectory_write_interval=200,
        temperature=temperature,
        data_container=f"mc_data/{args['temperature']}_{args['n_supercell']}_{args['label']}.dc",
    )
    mc.run(number_of_trial_steps=args["n_steps"])

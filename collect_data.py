"""
This script reads information from the equivalence class file created by parse_nrlist.py and specified by the Snakefile.
It works with PyMOL to save an mmCIF file of the original structure to a location specified in the Snakemake
configuration file. Alternative conformations are removed, hydrogens are added, and the modified structure is saved to
an mmCIF file specified by the Snakefile. The script collects data from the modified structure on potential hydrogen
bonds involving atoms of interest and nucleobase atom b-factors in the equivalence class member RNA chains. It writes
this data to a csv file specified by the Snakefile. If commit_hash is set to true in the Snakemake configuration file,
the commit hash of the repo will also be written within a commented line to the csv file if no uncommitted changes
have been made to the repo.
"""

import sys
import os
import csv
import subprocess
import time
from pymol import cmd
from pymol import stored
import numpy as np
import pandas as pd
from datetime import datetime
import residue_library
import eval_H_bonding
import remove_alt_conf

# Redirect stdout and stderr to log files.
stdout = sys.stdout
stderr = sys.stderr
stdout_file = open(snakemake.log.stdout, mode='w')
stderr_file = open(snakemake.log.stderr, mode='w')
sys.stdout = stdout_file
sys.stderr = stderr_file


# This function is to be run if an error is encountered.
def error(msg):
    # Create the output files expected by Snakemake.
    subprocess.run(["touch", snakemake.output.data])
    # If the structure has already been loaded, save the modified structure if specified.
    if cmd.count_atoms('all') > 0 and snakemake.config["save_modified_mmcif"]:
        cmd.save(modified_mmcif_dir + eq_class_mem_id + ".cif")
    # If the structure has already been loaded, remove the original mmCIF if specified.
    if cmd.count_atoms('all') > 0 and snakemake.config["remove_original_mmcif"]:
        subprocess.run(["rm", original_mmcif_dir + f"{pdb_id}.cif".lower()])
    # Print the error message.
    print(msg)
    # Close files, reset stdout and stderr, and exit.
    stdout_file.close()
    stderr_file.close()
    sys.stdout = stdout
    sys.stderr = stderr
    sys.exit(0)


# Check whether the folder exists for PyMOL to save mmCIF files into from the modified structure.
# If it does not exist, create the folder.
modified_mmcif_dir = snakemake.config["modified_mmcif_dir"]
if not os.path.isdir(modified_mmcif_dir):
    try:
        os.mkdir(modified_mmcif_dir)
    except FileExistsError:
        time.sleep(5.0)
        if not os.path.isdir(modified_mmcif_dir):
            os.mkdir(modified_mmcif_dir)

# If commit_hash is set to true in the Snakemake configuration file, check if any changes have been made to the repo and
# get the hash of the current git commit. If uncommitted changes have been made to anything other than files listed in
# the acceptable_changes variable defined below, print an error message and exit.
repo_changes = []
commit_hash = ""
if snakemake.config["commit_hash"]:
    repo_changes = list(subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"],
                                                cwd=os.path.dirname(os.path.realpath(__file__)))
                        .decode('ascii').strip().split("\n"))
    acceptable_changes = ['config/config.yaml', snakemake.config["rep_set_file"], snakemake.config["add_res_file"]]
    for file in repo_changes:
        if file.split(' ')[-1] in acceptable_changes:
            repo_changes.remove(file)
    if len(repo_changes) == 0:
        commit_hash = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                              cwd=os.path.dirname(os.path.realpath(__file__))).decode('ascii').strip()
    # Print the error message, close files, reset stdout and stderr, and exit.
    else:
        print("Error: Uncommitted changes have been made to the repo.")
        stdout_file.close()
        stderr_file.close()
        sys.stdout = stdout
        sys.stderr = stderr
        sys.exit(1)

# Collect the equivalence class name, PDB ID, model info, and chain info from the string describing the equivalence
# class member.
eq_class_build = []
pdb_id_build = []
model_build = []
chain_build = []
index = 0
track = 0
for char in snakemake.wildcards.eq_class_member:
    if index < 3:
        eq_class_build.append(char)
    elif char != '_' and index == 3:
        pdb_id_build.append(char)
    elif char != '_' and index == 4:
        model_build.append(char)
    elif char != '_' and index > 4:
        if index != track:
            chain_build.append([char])
            track = index
        else:
            chain_build[-1].append(char)
    if char == '_':
        index += 1
eq_class = "".join(eq_class_build[:-1])
pdb_id = "".join(pdb_id_build)
model = "".join(model_build)
chain_list = []
for chain in chain_build:
    chain_list.append("".join(chain))

# Create an identifier for the equivalence class member.
eq_class_mem_id = f'{eq_class}_{pdb_id}_{model}{"".join(["_" + chain for chain in chain_list])}'

# Check whether the folder exists for PyMOL to save mmCIF files into that are fetched from the PDB.
# If it does not exist, create the folder.
original_mmcif_dir = snakemake.config["original_mmcif_dir"]
if not os.path.isdir(original_mmcif_dir):
    try:
        os.mkdir(original_mmcif_dir)
    except FileExistsError:
        time.sleep(5.0)
        if not os.path.isdir(original_mmcif_dir):
            os.mkdir(original_mmcif_dir)

# Change fetch_path to original_mmcif_dir so that PyMOL drops the mmCIF files into this folder when running fetch.
cmd.set('fetch_path', cmd.exp_path(original_mmcif_dir))

# Construct two strings that can be used with PyMOL to select donor and rotatable donor atoms.
donor_string = ''
rotatable_donor_string = ''
for residue in residue_library.RESIDUE_LIBRARY:
    if residue["res"] in snakemake.config["included_residues"]:
        for donor in residue['don']:
            donor_string += f'(resn {residue["res"]} and name {donor[0]}) '
            if donor[2]:
                rotatable_donor_string += f'(resn {residue["res"]} and name {donor[0]}) '
donor_string = donor_string[:-1]
rotatable_donor_string = rotatable_donor_string[:-1]

# Construct two lists of strings containing residue and atom names that describe either donor or acceptor atoms of
# particular interest. Additionally, construct two strings that can be used with PyMOL to select all possible donor or
# all possible acceptor atoms of particular interest.
donors_of_interest = []
donors_of_interest_str = ''
for atom in snakemake.config["donors_of_interest"]:
    donors_of_interest.append(atom)
    donors_of_interest_str += f'(resn {atom.split(".")[0]} and name {atom.split(".")[1]}) '
donors_of_interest_str = donors_of_interest_str[:-1]
acceptors_of_interest = []
acceptors_of_interest_str = ''
for atom in snakemake.config["acceptors_of_interest"]:
    acceptors_of_interest.append(atom)
    acceptors_of_interest_str += f'(resn {atom.split(".")[0]} and name {atom.split(".")[1]}) '
acceptors_of_interest_str = acceptors_of_interest_str[:-1]

# Prepare a string that can be used with PyMOL to identify the chains of all the equivalence class member RNAs.
mem_rna_chains = " ".join(["chain " + chain for chain in chain_list])

# Retrieve the structure that contains the equivalence class member RNA chains.
cmd.fetch(pdb_id)

# If the structure contains more than one model (or "state" in PyMOL), create a new PyMOL object of just that model. If
# only a single model exists, change the name of the object to identify the model number.
if cmd.count_states(pdb_id) > 1:
    cmd.create(f'{pdb_id}_state_{model}', selection=pdb_id,
               source_state=model,
               target_state=1)
    cmd.delete(pdb_id)
else:
    cmd.set_name(pdb_id, f'{pdb_id}_state_{model}')

# Remove atoms representing alternative conformations.
remove_status = remove_alt_conf.remove(eq_class_mem_id)
successful_completion = remove_status[0]
if not successful_completion:
    error(remove_status[1])

# Remove any hydrogens that loaded with the structure.
cmd.remove('elem H')

# Add hydrogens to non-rotatable donors that are a part of or near the equivalence class member RNA.
cmd.h_add(f'(({donor_string}) and not ({rotatable_donor_string})) within {snakemake.config["search_dist"]} of '
          f'({mem_rna_chains})')

# Randomly sample 100 atom indices in the structure and record the info of the associated atoms for later comparison.
num_atoms = cmd.count_atoms('all')
indices = np.random.default_rng().integers(low=1, high=num_atoms+1, size=100)
stored.check_one = []
for index in indices:
    cmd.iterate(f'index {index}', 'stored.check_one.append([index, name, resn, resi, chain])')

# Store a list of donors of interest from the equivalence class member RNA chains.
stored.donor_list = []
cmd.iterate(f'({mem_rna_chains}) and ({donors_of_interest_str})',
            'stored.donor_list.append([index, name, resn, resi, chain, segi])')

# Store the number of heavy atoms belonging to organic molecules or polymers near the donor of interest within a
# dictionary. Additionally, store the average b-factor of the heavy atoms that make up the nucleobase containing the
# donors of interest. Values of one for the DOI and don_can_NA keys indicate that the donor atom is a donor of interest
# and belongs to a canonical DNA or RNA residue, respectively.
don_info_dict = {"don_index": [], "don_name": [], "don_resn": [], "don_resi": [], "don_chain": [], "don_segi": [],
                 "count_1": [], "count_2": [], "b_factor": [], "DOI": [], "don_can_NA": []}
for donor in stored.donor_list:
    # Count the number of heavy atoms belonging to organic molecules or polymers near the donor. Exclude the nucleobase
    # of the donor.
    count_1 = cmd.count_atoms(f'((organic or polymer) within {snakemake.config["count_dist_1"]} of index {donor[0]}) '
                              f'and not ((sidechain and byres index {donor[0]}) or elem H)')
    count_2 = cmd.count_atoms(f'((organic or polymer) within {snakemake.config["count_dist_2"]} of index {donor[0]}) '
                              f'and not ((sidechain and byres index {donor[0]}) or elem H)')
    # Collect the b-factors for side chain atoms.
    stored.b_factors = []
    cmd.iterate(f"resn {donor[2]} and resi \\{donor[3]} and chain {donor[4]} and sidechain and not elem H",
                "stored.b_factors.append(b)")
    # If no b-factors were stored, exit with an error message. This would likely be due to the PyMOL sidechain selection
    # operator not working with a particular residue.
    if len(stored.b_factors) == 0:
        error(f"Error: No b-factors were obtained from residue {donor[4]}.{donor[2]}.{donor[3]} in {eq_class_mem_id}.")
    # Ensure that the correct number of b-factors were collected.
    if ((donor[2] in ["A", "DA"] and len(stored.b_factors) != 10) or
            (donor[2] in ["C", "DC"] and len(stored.b_factors) != 8) or
            (donor[2] in ["G", "DG"] and len(stored.b_factors) != 11)):
        error(f"Error: The correct number of b-factors were not obtained from residue {donor[4]}.{donor[2]}.{donor[3]} "
              f"in {eq_class_mem_id}.")
    # Calculate the average of the b-factors.
    b_factor_avg = sum(stored.b_factors)/len(stored.b_factors)
    # Add the donor, heavy atom counts, and average b-factor to the count dictionary.
    row = donor + [count_1, count_2, b_factor_avg, 1, 1]
    for key, value in zip(don_info_dict, row):
        don_info_dict[key].append(value)

# Create a dataframe based on the info dictionary.
don_info_df = pd.DataFrame(don_info_dict).astype("object")

# Store a list of acceptors near the donors of interest within an atom pair dictionary. Also include two keys which have
# values containing both residue and atom names for donor and acceptor atoms.
don_atom_pair_dict = {"don_index": [], "don_name": [], "don_resn": [], "don_resi": [], "don_chain": [], "don_segi": [],
                      "acc_index": [], "acc_name": [], "acc_resn": [], "acc_resi": [], "acc_chain": [], "acc_segi": [],
                      "don_resn_name": [], "acc_resn_name": []}
for donor in stored.donor_list:
    # Find the acceptors near the donor. Exclude the nucleobase of the donor.
    stored.acceptors = []
    cmd.iterate(f'(acceptors within {snakemake.config["search_dist"]} of index {donor[0]}) and (organic or polymer) '
                f'and not (sidechain and byres index {donor[0]})',
                'stored.acceptors.append([index, name, resn, resi, chain, segi])')
    # Add the acceptors to the atom pair dictionary.
    for acceptor in stored.acceptors:
        row = donor + acceptor + [f'{donor[2]}.{donor[1]}', f'{acceptor[2]}.{acceptor[1]}']
        for key, value in zip(don_atom_pair_dict, row):
            don_atom_pair_dict[key].append(value)

# Create a dataframe based on the donor atom pair dictionary.
don_atom_pair_df = pd.DataFrame(don_atom_pair_dict)

# Acquire the H-bonding geometry measurements for all acceptors near each donor of interest.
don_h_bonds_dict = {"don_index": [], "don_name": [], "don_resn": [], "don_resi": [], "don_chain": [], "don_segi": [],
                    "acc_index": [], "acc_name": [], "acc_resn": [], "acc_resi": [], "acc_chain": [], "acc_segi": [],
                    "don_resn_name": [], "acc_resn_name": [], "don_acc_distance": [], "h_acc_distance": [],
                    "don_angle": [], "h_angle": [], "h_dihedral": [], "h_name": []}
for atom_pair in don_atom_pair_df.itertuples():
    # Store the donor and acceptor atom values for the dataframe row.
    don_list = [atom_pair.don_index, atom_pair.don_name, atom_pair.don_resn, atom_pair.don_resi, atom_pair.don_chain]
    acc_list = [atom_pair.acc_index, atom_pair.acc_name, atom_pair.acc_resn, atom_pair.acc_resi, atom_pair.acc_chain]
    # Retrieve the H-bond measurements for the atom pair.
    h_bond_list = eval_H_bonding.evaluate(don_list, acc_list, eq_class_mem_id, residue_library.RESIDUE_LIBRARY,
                                          snakemake.config["single_h_donors"], snakemake.config["dual_h_donors"])
    # If the H-bond evaluation is not successful, print the error message(s) and exit.
    successful_completion = h_bond_list[0]
    if not successful_completion:
        error(h_bond_list[1])
    # Add the H-bond measurements to the dictionary.
    else:
        for h_bond in h_bond_list[1]:
            row = (don_list + [atom_pair.don_segi] + acc_list + [atom_pair.acc_segi] +
                   [atom_pair.don_resn_name, atom_pair.acc_resn_name] + h_bond)
            for key, value in zip(don_h_bonds_dict, row):
                don_h_bonds_dict[key].append(value)

# Create a dataframe based on the dictionary.
don_h_bonds_df = pd.DataFrame(don_h_bonds_dict).astype("object")

# Store a list of acceptors of interest from nucleobases containing the donors of interest.
stored.acceptor_list = []
for donor in stored.donor_list:
    cmd.iterate(f'resn {donor[2]} and resi \\{donor[3]} and chain {donor[4]} and ({acceptors_of_interest_str})',
                'stored.acceptor_list.append([index, name, resn, resi, chain, segi])')

# Store a list of donors near the acceptors of interest within an atom pair dictionary. Also include two keys which have
# values containing both residue and atom names for acceptor and donor atoms.
acc_atom_pair_dict = {"don_index": [], "don_name": [], "don_resn": [], "don_resi": [], "don_chain": [], "don_segi": [],
                      "acc_index": [], "acc_name": [], "acc_resn": [], "acc_resi": [], "acc_chain": [], "acc_segi": [],
                      "don_resn_name": [], "acc_resn_name": []}
for acceptor in stored.acceptor_list:
    # Find the donors near the acceptor. Exclude the nucleobase of the acceptor.
    stored.donors = []
    cmd.iterate(f'(donors within {snakemake.config["search_dist"]} of index {acceptor[0]}) and (organic or polymer) '
                f'and not ((sidechain and byres index {acceptor[0]}) or name OP2)',
                'stored.donors.append([index, name, resn, resi, chain, segi])')
    # Add the donors to the atom pair dictionary.
    for donor in stored.donors:
        row = donor + acceptor + [f'{donor[2]}.{donor[1]}', f'{acceptor[2]}.{acceptor[1]}']
        for key, value in zip(acc_atom_pair_dict, row):
            acc_atom_pair_dict[key].append(value)

# Create a dataframe based on the acceptor atom pair dictionary.
acc_atom_pair_df = pd.DataFrame(acc_atom_pair_dict)

# Acquire the H-bonding geometry measurements for all donors near each acceptor of interest. Values of one for the AOI
# and don_can_NA keys indicate that the acceptor atom is an acceptor of interest and belongs to a canonical DNA or RNA
# residue, respectively.
acc_h_bonds_dict = {"don_index": [], "don_name": [], "don_resn": [], "don_resi": [], "don_chain": [], "don_segi": [],
                    "acc_index": [], "acc_name": [], "acc_resn": [], "acc_resi": [], "acc_chain": [], "acc_segi": [],
                    "don_resn_name": [], "acc_resn_name": [], "don_acc_distance": [], "h_acc_distance": [],
                    "don_angle": [], "h_angle": [], "h_dihedral": [], "h_name": [], "AOI": [], "don_can_NA": []}
for atom_pair in acc_atom_pair_df.itertuples():
    # Only calculate H-bonding measurements if the donor belongs to a canonical DNA or RNA residue.
    if atom_pair.don_resn in ["A", "C", "G", "U", "DA", "DC", "DG", "DT"]:
        # Store the donor and acceptor atom values for the dataframe row.
        don_list = [atom_pair.don_index, atom_pair.don_name, atom_pair.don_resn, atom_pair.don_resi,
                    atom_pair.don_chain]
        acc_list = [atom_pair.acc_index, atom_pair.acc_name, atom_pair.acc_resn, atom_pair.acc_resi,
                    atom_pair.acc_chain]
        # Retrieve the H-bond measurements for the atom pair.
        h_bond_list = eval_H_bonding.evaluate(don_list, acc_list, eq_class_mem_id, residue_library.RESIDUE_LIBRARY,
                                              snakemake.config["single_h_donors"], snakemake.config["dual_h_donors"])
        # If the H-bond evaluation is not successful, print the error message(s) and exit.
        successful_completion = h_bond_list[0]
        if not successful_completion:
            error(h_bond_list[1])
        # Add the H-bond measurements to the dictionary.
        else:
            for h_bond in h_bond_list[1]:
                row = (don_list + [atom_pair.don_segi] + acc_list + [atom_pair.acc_segi] +
                       [atom_pair.don_resn_name, atom_pair.acc_resn_name] + h_bond + [1, 1])
                for key, value in zip(acc_h_bonds_dict, row):
                    acc_h_bonds_dict[key].append(value)
    # Set H-bonding information to NaN and don_can_NA to 0 if the donor does not belong to a canonical RNA or DNA
    # residue.
    else:
        row = list(atom_pair)[1:] + [pd.NA] * 6 + [1, 0]
        for key, value in zip(acc_h_bonds_dict, row):
            acc_h_bonds_dict[key].append(value)

# Create a dataframe based on the dictionary.
acc_h_bonds_df = pd.DataFrame(acc_h_bonds_dict).astype("object")

# Record info on the atoms associated with the previously determined 100 atom indices and check whether there are any
# changes.
stored.check_two = []
for index in indices:
    cmd.iterate(f'index {index}', 'stored.check_two.append([index, name, resn, resi, chain])')
if not stored.check_one == stored.check_two:
    error(f"Error: The indices in equivalence class member {snakemake.wildcards.eq_class_member} changed.")

# Prepare a master dataframe containing heavy atom count, b-factor, H-bonding data, and other relevant information.
master_df = don_info_df.merge(don_h_bonds_df, how='outer').merge(acc_h_bonds_df, how='outer')
master_df.loc[:, ['model', 'PDB', 'eq_class_member']] = [model, pdb_id, eq_class_mem_id]

# Write a csv containing the data stored in the master dataframe.
with open(snakemake.output.data, "w") as csv_file:
    writer = csv.writer(csv_file)
    if commit_hash:
        writer.writerow([f"# dual-H-bonding-nucleobases repo git commit hash: {commit_hash}"])
    writer.writerow([f"# representative set file: {snakemake.config['rep_set_file']}"])
    writer.writerow([f"# file created on: {datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')}"])
master_df.to_csv(snakemake.output.data, index=False, mode='a', na_rep='NaN')

# Save the modified structure if specified.
if snakemake.config["save_modified_mmcif"]:
    cmd.save(modified_mmcif_dir + eq_class_mem_id + ".cif")

# Remove the original mmCIF if specified.
if snakemake.config["remove_original_mmcif"]:
    subprocess.run(["rm", original_mmcif_dir + f"{pdb_id}.cif".lower()])

# Close files and reset stdout and stderr.
stdout_file.close()
stderr_file.close()
sys.stdout = stdout
sys.stderr = stderr

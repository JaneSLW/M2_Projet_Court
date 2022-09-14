"""This program allows to determine the accessible solvant area (ASA), the relative 
accessible solvant area (rASA) and solvant accessibility percentage from the coordinates 
derived from a PDB file.

Usage:
======
python code_src.py -n protein_name
OR
python code_src.py -f pdb_file.pdb 

protein_name: protein name as written in pdb data base
pdb_file.pdb: name of the pdb file that you have
"""

__authors__ = ("Jane Schadtler-Law")
__contact__ = ("jane.schadtler-law@etu.u-paris.fr")
__date__ = "2022-09-14"

import argparse
import csv
import math
import os.path
import sys
import warnings

import numpy as np
import Bio.PDB
from scipy.spatial.distance import pdist, squareform

H20_RADIUS = 1.7 # water radius in angstrom
POINTS_SPHERE = 92


def sphere_92(vdw_r, atm_coor_vec):
	"""Create a sphere of 92 points around an atom for accessible surface calculations.

	This function creates a sphere of 92 points uniformely distributed around the point
	defined in input. It requires a radius (Van der Waals) and its coordinates.
	The radius of a water molecule will be added in the sphere calculation. 

	Parameters 
	---------- 
	vdw_r : float
		Van der Waals radius of the atom.
	atm_coor_vec : list of floats
		Cartesian coordinates of the atom.

	Returns 
	-------
	array
		Array[92,3] of the cartesians coordinates (x,y,z) of each point of the sphere. 
	"""

	points_cartesian = []	# temporary stocking of cartesian coordinates under list format
	phi = np.pi * (3 - np.sqrt(5))	# golden angle in radian
	sphere_radius = float(vdw_r) + H20_RADIUS	

	for i in range(POINTS_SPHERE):
		lat = math.asin(-1 + 2 * float(i / (92 + 1))) #	spheric coordinates, latitude
		lon = phi * i	# longitutde

		coor_x = np.cos(lon) * np.cos(lat) * sphere_radius + atm_coor_vec[0]	# cartesian coordinates
		coor_y = np.sin(lon) * np.cos(lat) * sphere_radius + atm_coor_vec[1]
		coor_z = np.sin(lat) * sphere_radius + atm_coor_vec[2]

		points_cartesian.append([coor_x, coor_y, coor_z])

	return np.array(points_cartesian) # Creating a numpy array of lists of the cartesian coordiantes points


# Dictionnaries
vdw_rad = {"C": 1.7, "N": 1.55, "O": 1.52, "P": 1.8, "S":1.8}	# Dictionnary for Van der Waals radius in Angstrom
max_asa = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLU": 223,
 "GLN": 225, "GLY": 104, "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236,"MET": 224, 
 "PHE": 240, "PRO": 159, "SER": 155, "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174}	# MaxASA maximum possible solvent accessible surface area in Angstrom2 ; Theorical Values. (Tien et al. 2013)

# Retreiving arguments
if len(sys.argv) < 2:
	sys.exit("Usage exemple: python main.py -n 1MH1")

parser_arg = argparse.ArgumentParser(description="Calculate absolute accessible solvent area and relative area.")
parser_arg.add_argument("-f","--file", type=str, help="If you have a pdb file, file name of the protein in .pdb format.")
parser_arg.add_argument("-n", "--name", type=str, help="name of the protein as found in pdb server.")
parser_arg.add_argument("-o", "--output_file", action=argparse.BooleanOptionalAction, help="Adding an output file with results per residue", default=False)

args = parser_arg.parse_args()

if args.name:
	name_protein = args.name.strip().upper()
	repository = Bio.PDB.PDBList(pdb=name_protein) #repository = PDBList(pdb='1MH1')
	repository.retrieve_pdb_file(name_protein, file_format='pdb', pdir='.') # Getting the pdb file 
	structure_id = name_protein
	filename = "pdb" + name_protein.lower() + ".ent"
elif args.file:
	try:
		os.path.isfile(args.file)
		structure_id = args.file.replace(".pdb", "").strip().upper()
		filename = args.file
	except:
		sys.exit("Impossible de trouver le fichier pdb")
	
# Parsing PDB file
with warnings.catch_warnings():
    warnings.simplefilter("ignore")	# don't show warnings
    parser_pdb = Bio.PDB.PDBParser(PERMISSIVE=1)	# create a PDBParser object (and ignore error)
    structure = parser_pdb.get_structure(structure_id, filename)

# Extracting from structure object : atoms, their coordinates, their parent residue
atom_id = []
residues = []
residue_to_atm = []
atom_coords = []

for chain in structure[0]:
	for residue in chain:
		if Bio.PDB.is_aa(residue): # NO HETERO ATM
			residues.append(str(residue)[9:12])
			for atom in residue:
				atom_id.append(str(atom)[6]) # Saving atom ids
				atom_coords.append(atom.get_coord()) # Extracting atom coordinates
				#residue_to_atm.append([str(atom.get_parent()).split(" ")[1], (str(atom.get_parent()).split(" "))[4].replace("resseq=", "")])
				residue_to_atm.append(str(atom.get_parent()).split(" ")[-3].replace("resseq=", ""))

# Checking disordered atom and deleting them
nv_atom_id = []
nv_atom_coords = []
nv_residue_to_atm = []
for i in range(len(atom_id)):
	if atom_id[i] in list(vdw_rad.keys()):
		nv_atom_id.append(atom_id[i])
		nv_atom_coords.append(atom_coords[i])
		nv_residue_to_atm.append(residue_to_atm[i])

# Distance matrix
atom_coords_array=np.array(nv_atom_coords) 
distances_array = squareform(pdist(atom_coords_array))	# pairwise distance calculation

# Total accessible surface area
ASA_per_atm = []
res_on_surface = []
tot_acc_area = 0
ASA = 0

for atm_i in range(len(nv_atom_id)):
	neighbour_atoms = np.where((distances_array[:][atm_i] < 10) & (distances_array[:][atm_i] > 0))[0]	# row index of neighboor atoms in an np.array format 
	sphere_points_array = sphere_92(vdw_rad[nv_atom_id[atm_i]], atom_coords_array[atm_i])	# calculating sphere points coordinates
	surface_points = 0
	
	for point in range(POINTS_SPHERE):
		free_point = 0
		for neighbour_atom in neighbour_atoms:
			dist_point_center = np.sqrt((sphere_points_array[point][0] - atom_coords_array[neighbour_atom][0])**2 
				+ (sphere_points_array[point][1] - atom_coords_array[neighbour_atom][1])**2 
				+ (sphere_points_array[point][2] - atom_coords_array[neighbour_atom][2])**2)	# distance between 

			if(dist_point_center < (vdw_rad[nv_atom_id[neighbour_atom]] + H20_RADIUS)):	# threshold : neighboor atom's radius + H20(solvant) radius 
				free_point = 0
				break	# if just one neighboor is too close, than this point is considered buried regardless of the other neighboors

			free_point +=1 

		if(free_point > 0): 
			surface_points += 1	# nb of free sphere's points (max = 92)
	tot_acc_area = tot_acc_area + (4 * np.pi * (vdw_rad[nv_atom_id[atm_i]] + H20_RADIUS)**2)	# total area of solvated sphere (for accessibility to solvant percentage calculation)
	ASA_per_atm.append(surface_points * (4 * np.pi * (vdw_rad[nv_atom_id[atm_i]] + H20_RADIUS)**2 ) / POINTS_SPHERE)	# total accessible surface area (ASA) in Angstrom2 
	if surface_points > 0:
		res_on_surface.append(residue_to_atm[atm_i])	# to only keep residues on surface (and accelerate rASA calculation)

ASA = sum(ASA_per_atm)	# total ASA of protein

ASA_per_res = {}	# ASA per residue to calculate rASA 
for key, value in zip(res_on_surface, ASA_per_atm):
	ASA_per_res[key] = ASA_per_res.get(key, 0) + value

rASA = 0	# relative solvent accessible surface area (rASA)
for residue_i in range(len(residues)):
	if str(residue_i) in ASA_per_res.keys():
		if residues[residue_i] in max_asa.keys(): # because there are modified residues "HYP"
			rASA = rASA + (ASA_per_res[str(residue_i)] / max_asa[residues[residue_i]])	# rASA = ASA / max_asa_tot

accessibility_100 = 100 * ASA / tot_acc_area	# accessibility to solvant percentage = 100 x area of solvated sphere exposed to solvent / total area of solvated sphere.

# Output in terminal
print(f"""
For the protein {structure_id}
	- The total accessible surface area is {ASA:.2f} \u212B\u00B2.
	- The relative accessible surface area is {rASA:.2f} \u212B\u00B2.
	- The accessibility to solvant percentage is {accessibility_100:.2f} %.
	""")

# Output in txt file
if args.output_file:
	with open(f"results_{structure_id}.txt", "w") as res_file:
		writer = csv.writer(res_file, delimiter='\t')
		res_file.write(f"""
For the protein {structure_id}
	- The total accessible surface area is {ASA:.2f} \u212B\u00B2.
	- The relative accessible surface area is {rASA:.4f} \u212B\u00B2.
	- The accessibility to solvant percentage is {accessibility_100:.2f} %. \n

Residues on surface : 

Residues	ASA
""")
		writer.writerows(zip(ASA_per_res.keys(),residues,ASA_per_res.values()))
		#res_file.write('\n'.join(zip(map(str, nv_atom_id))))

# M2_projet_court
### 2022 M2 Biologie Informatique, Projet court, Surface accessible d'une protéine
Ce programme permet de calculer la surface accessible au solvant (accessible solvent area, ASA) à partir d'un fichier sous format pdb ou à partir de l'identifiant pdb. Il calcul ausso la surface accessible au solvant relative (rASA) et le pourcentage d'accessibilité. Il est possible d'avoir le détail par résidu de l'ASA.  

## Download this repository
```
git clone https://github.com/JaneSLW/M2_projet_court.git
cd M2_projet_court
```

## Install dependencies
### Create conda environnement
```
conda env create -f projet_court.yml
```
### Load conda environment:
```
conda activate projet_court
```
## Lauch program
### With pdb file 
If you have your own pdb file, copy paste it in working directory and replace "1bkv.pdb" by your file.
```
python src/main.py -f data/2022_09_14/peptide_LK6/6a5j.pdb
```
### With protein name
The name of the protein must correspond to its pdb ENTRY id. Replace "1BKV" by any pdb ENTRY id.
The pdb file is not required.
```
python src/main.py -n 6A5J
```
### Generate ouput file with residue's details
By default, the program won't generate a TXT file. By adding "-o" or "--output_file" you will be able to keep your results in a txt file with residues details.
```
python src/main.py -n 1BKV -o 
```
OR
```
python src/main.py -f data/2022_09_14/peptide_LK6/6a5j.pdb -o
```
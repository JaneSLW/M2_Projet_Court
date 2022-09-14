# M2_projet_court
### 2022 M2 Biologie Informatique, Projet court, Surface accessible d'une protéine
Ce programme permet de calculer la surface accessible au solvant (accessible solvent area, ASA) à partir d'un fichier sous format pdb ou à partir de l'identifiant pdb. Il calcul ausso la surface accessible au solvant relative (rASA) et le pourcentage d'accessibilité. Il est possible d'avoir le détail par résidu de l'ASA.  

## Download this repository
```
git clone https://github.com/JaneSLW/M2_projet_court/Surface_accessible_JSL.git
cd Surface_accessible_JSL
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
python main.py -f 1bkv.pdb
```
### With protein name
The name of the protein must correspond to its pdb ENTRY id. Replace "1BKV" by any pdb ENTRY id.
The pdb file is not required.
```
python main.py -n 1BKV
```
### Generate ouput file with residue's details
By default, the program won't generate a TXT file. By adding "-o" or "--output_file" you will be able to keep your results in a txt file with residues details.
```
python main.py -n 1BKV -o 
```
OR
```
python main.py -f 1bkv.pdb -o
```
from Bio.PDB import *
from rdkit import Chem
from rdkit.Chem import AllChem
import tempfile
import os
from ..default_config.chemistry import radii, polarHydrogens

"""
xyzrn.py: Read a pdb file and output it is in xyzrn for use in MSMS
Pablo Gainza - LPDI STI EPFL 2019
This file is part of MaSIF.
Released under an Apache License 2.0
"""


def output_pdb_as_xyzrn(pdbfilename, xyzrnfilename, keep_hetatms=None):
    """
    pdbfilename: input pdb filename
    xyzrnfilename: output in xyzrn format.
    keep_hetatms: list of hetatms to keep
    """
    if keep_hetatms is None:
        keep_hetatms = []
    parser = PDBParser()
    struct = parser.get_structure(pdbfilename, pdbfilename)
    outfile = open(xyzrnfilename, "w")
    for atom in struct.get_atoms():
        name = atom.get_name()
        residue = atom.get_parent()
        # Ignore hetatms.
        # if residue.get_id()[0] != " " and residue.get_id()[0][-3:] != 'RC8':
        if residue.get_id()[0] != " " and residue.get_id()[0][-3:] not in keep_hetatms:
            continue
        resname = residue.get_resname()
        reskey = residue.get_id()[1]
        chain = residue.get_parent().get_id()
        atomtype = name[0]

        color = "Green"
        coords = None
        if atomtype in radii and (resname in polarHydrogens or resname in keep_hetatms):
            if atomtype == "O":
                color = "Red"
            if atomtype == "N":
                color = "Blue"
            if atomtype == "H":
                if resname in keep_hetatms:
                    pass
                elif name in polarHydrogens[resname]:
                    color = "Blue"  # Polar hydrogens
            coords = "{:.06f} {:.06f} {:.06f}".format(
                atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]
            )
            insertion = "x"
            if residue.get_id()[2] != " ":
                insertion = residue.get_id()[2]
            full_id = "{}_{:d}_{}_{}_{}_{}".format(
                chain, residue.get_id()[1], insertion, resname, name, color
            )
        if coords is not None:
            outfile.write(coords + " " + radii[atomtype] + " 1 " + full_id + "\n")


"""
sdf_to_xyzrn.py: Read an SDF file, convert to PDB, and output in XYZRN format for use in MSMS
Adapted from Pablo Gainza's xyzrn.py - LPDI STI EPFL 2019
"""


def sdf_to_pdb(sdffilename, pdbfilename=None):
    """
    Convert SDF file to PDB format using RDKit

    Args:
        sdffilename: input SDF filename
        pdbfilename: output PDB filename (if None, returns temp file path)

    Returns:
        Path to the PDB file
    """
    # Read SDF file
    supplier = Chem.SDMolSupplier(sdffilename, removeHs=False)

    if pdbfilename is None:
        # Create a temporary PDB file
        temp_file = tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False)
        pdbfilename = temp_file.name
        temp_file.close()

    # Process molecules from SDF
    with open(pdbfilename, "w") as pdb_file:
        for i, mol in enumerate(supplier):
            if mol is None:
                continue

            # Add hydrogens if not present
            if mol.GetNumAtoms() > 0:
                mol = Chem.AddHs(mol, addCoords=True)

                # Generate 3D coordinates if not present
                if mol.GetNumConformers() == 0:
                    AllChem.EmbedMolecule(mol, randomSeed=42)
                    AllChem.UFFOptimizeMolecule(mol)

                # Write to PDB
                pdb_block = Chem.MolToPDBBlock(mol)
                pdb_file.write(pdb_block)

                # Add TER and END records between molecules if multiple
                if i < len(supplier) - 1:
                    pdb_file.write("TER\n")

        pdb_file.write("END\n")

    return pdbfilename


def output_pdb_as_xyzrn_internal(pdbfilename, xyzrnfilename, keep_hetatms=None):
    """
    Internal function to convert PDB to XYZRN format
    (Same as original output_pdb_as_xyzrn but renamed for clarity)

    Args:
        pdbfilename: input pdb filename
        xyzrnfilename: output in xyzrn format
        keep_hetatms: list of hetatms to keep
    """
    if keep_hetatms is None:
        keep_hetatms = []

    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("structure", pdbfilename)
    outfile = open(xyzrnfilename, "w")

    for atom in struct.get_atoms():
        name = atom.get_name()
        residue = atom.get_parent()

        # Handle both standard residues and hetatms from SDF
        # For SDF-derived molecules, treat them as hetatms but process them
        resname = residue.get_resname()

        # For SDF-derived molecules, we might want to process all atoms
        # Check if this is a HETATM that should be kept
        if residue.get_id()[0] != " ":
            # This is a HETATM
            if resname not in keep_hetatms and resname not in ["UNL", "LIG", "MOL"]:
                # Skip hetatms not in keep list and not common ligand names
                continue

        reskey = residue.get_id()[1]
        chain = residue.get_parent().get_id()
        atomtype = name[0]

        color = "Green"
        coords = None

        # Check if atom type is in radii dictionary
        if atomtype in radii:
            # Assign colors based on atom type
            if atomtype == "O":
                color = "Red"
            elif atomtype == "N":
                color = "Blue"
            elif atomtype == "H":
                color = "Blue"  # Default polar hydrogen color
                # Check if this hydrogen is polar (if residue is in polarHydrogens)
                if resname in polarHydrogens:
                    if name not in polarHydrogens[resname]:
                        color = "Green"  # Non-polar hydrogen
            elif atomtype == "C":
                color = "Green"
            elif atomtype == "S":
                color = "Yellow"
            elif atomtype == "P":
                color = "Orange"

            # Get coordinates
            coords = "{:.06f} {:.06f} {:.06f}".format(
                atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]
            )

            # Handle insertion code
            insertion = "x"
            if residue.get_id()[2] != " ":
                insertion = residue.get_id()[2]

            # Create full atom ID
            full_id = "{}_{:d}_{}_{}_{}_{}".format(
                chain, residue.get_id()[1], insertion, resname, name, color
            )

            # Write to output file
            if coords is not None:
                outfile.write(coords + " " + radii[atomtype] + " 1 " + full_id + "\n")

    outfile.close()


def convert_sdf_to_xyzrn(
    sdffilename, xyzrnfilename, keep_hetatms=None, keep_temp_pdb=False
):
    """
    Main function to convert SDF file to XYZRN format

    Args:
        sdffilename: input SDF filename
        xyzrnfilename: output XYZRN filename
        keep_hetatms: list of residue names to keep from hetatms (default: ['UNL', 'LIG', 'MOL'])
        keep_temp_pdb: if True, keeps the temporary PDB file (default: False)

    Returns:
        Path to the output XYZRN file
    """
    if keep_hetatms is None:
        # Default: keep common ligand residue names from SDF files
        keep_hetatms = ["UNL", "LIG", "MOL"]

    temp_pdb_path = None

    try:
        # Step 1: Convert SDF to PDB
        print(f"Converting SDF file {sdffilename} to PDB format...")
        temp_pdb_path = sdf_to_pdb(sdffilename)
        print(f"Temporary PDB file created: {temp_pdb_path}")

        # Step 2: Convert PDB to XYZRN
        print(f"Converting PDB to XYZRN format...")
        output_pdb_as_xyzrn_internal(temp_pdb_path, xyzrnfilename, keep_hetatms)
        print(f"XYZRN file created: {xyzrnfilename}")

        return xyzrnfilename

    finally:
        # Clean up temporary PDB file if requested
        if temp_pdb_path and os.path.exists(temp_pdb_path) and not keep_temp_pdb:
            os.remove(temp_pdb_path)
            print(f"Temporary PDB file removed: {temp_pdb_path}")


def convert_pdb_to_xyzrn(pdbfilename, xyzrnfilename, keep_hetatms=None):
    """
    Convert a PDB file to XYZRN format for MSMS

    Args:
        pdbfilename: input PDB filename
        xyzrnfilename: output XYZRN filename
        keep_hetatms: list of residue names to keep from hetatms (default: [])

    Returns:
        Path to the output XYZRN file
    """
    if keep_hetatms is None:
        keep_hetatms = []

    print(f"Converting PDB file {pdbfilename} to XYZRN format...")
    output_pdb_as_xyzrn_internal(pdbfilename, xyzrnfilename, keep_hetatms)
    print(f"XYZRN file created: {xyzrnfilename}")

    return xyzrnfilename

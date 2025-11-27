import os
import numpy
import logging
from subprocess import Popen, PIPE
from pathlib import Path
from precise.preprocess.default_config.global_vars import apbs_bin, pdb2pqr_bin, multivalue_bin

# Initialize logger
logger = logging.getLogger(__name__)

"""
computeAPBS.py: Wrapper function to compute the Poisson Boltzmann electrostatics for a surface using APBS.
Pablo Gainza - LPDI STI EPFL 2019
This file is part of MaSIF.
Released under an Apache License 2.0
"""

def _parse_metal_ion_from_mol2(mol2_file):
    """
    Attempts to parse a MOL2 file to see if it's a single metal ion.
    Returns a dict with ion info if it is, otherwise None.
    
    This is a simplified parser assuming a single-atom MOL2.
    """
    # Common metal atom types (atom_type field in mol2)
    COMMON_METALS = ['MN', 'MG', 'ZN', 'CA', 'FE', 'K', 'NA', 'CU', 'NI', 'CO', 'CD']
    
    # Radii from PARSE force field (approximate).
    METAL_RADII_PARSE = {
        'MN': 0.90,
        'MG': 0.72,
        'ZN': 0.88,
        'CA': 1.12,
        'K':  1.52,
        'NA': 1.16,
        'CU': 0.87,
        'FE': 0.75, # Assume Fe2+
        'NI': 0.83,
        'CO': 0.82,
        'CD': 1.09,
    }

    # Formal Charges needed for PQR
    METAL_FORMAL_CHARGE = {
        'MN': 2.0, 
        'MG': 2.0,
        'ZN': 2.0,
        'CA': 2.0,
        'FE': 2.0, # Assuming Fe2+
        'K':  1.0,
        'NA': 1.0,
        'CU': 2.0, # Assuming Cu2+
        'NI': 2.0,
        'CO': 2.0, # Assuming Co2+
        'CD': 2.0,
    }

    try:
        with open(mol2_file, 'r') as f:
            lines = f.readlines()
        
        in_atom_section = False
        for line in lines:
            if line.startswith('@<TRIPOS>ATOM'):
                in_atom_section = True
                continue
            if line.startswith('@<TRIPOS>'):
                in_atom_section = False
                continue
            
            if in_atom_section:
                parts = line.strip().split()
                if len(parts) < 9:
                    continue
                
                # MOL2 format (approx):
                # atom_id(0) atom_name(1) x(2) y(3) z(4) atom_type(5) sub_id(6) sub_name(7) charge(8)
                
                atom_name = parts[1] # e.g., MN1
                x = float(parts[2])
                y = float(parts[3])
                z = float(parts[4])
                atom_type = parts[5].strip().upper() 
                
                # Check if the atom_type (e.g., 'MN') is in our list
                if atom_type in COMMON_METALS:
                    # Found a metal!
                    res_name = atom_type # e.g., 'MN', will be the residue name in PQR
                    
                    # Get the formal charge and radius from our dictionaries
                    charge = METAL_FORMAL_CHARGE.get(res_name, float(parts[8])) 
                    radius = METAL_RADII_PARSE.get(res_name, 0.90) # Default 0.90
                    
                    return {
                        'atom_name': atom_name, 
                        'res_name': res_name,  
                        'x': x,
                        'y': y,
                        'z': z,
                        'charge': charge, 
                        'radius': radius 
                    }
                else:
                    # This is a non-metal ligand, let PDB2PQR handle it
                    return None
        
    except Exception as e:
        print(f"Could not parse MOL2 file {mol2_file}: {e}")
        logger.warning(f"Could not parse MOL2 file {mol2_file}: {e}")
        return None
        
    return None # No atom found or not a metal

def computeAPBS(vertices, pdb_file, tmp_file_base, mol2_file=None):
    """
        Calls APBS, pdb2pqr, and multivalue and returns the charges per vertex.
        
        MODIFIED: Automatically handles metal ions in mol2_file by
        running PDB2PQR without --ligand and then manually appending
        the HETATM entry to the PQR file.
    """
    print('Starting APBS computation...')
    fields = tmp_file_base.split("/")[0:-1]
    directory = "/".join(fields) + "/"
    filename_base = tmp_file_base.split("/")[-1]
    pdbname = pdb_file.split("/")[-1]
    
    metal_info = None
    if mol2_file is not None:
        # Check if the mol2 file is a metal ion
        metal_info = _parse_metal_ion_from_mol2(mol2_file)

    pdb2pqr_bin = os.environ["PDB2PQR_BIN"]
    args = [
        pdb2pqr_bin, pdbname, filename_base,
        "--ff=PARSE",
        "--whitespace",
        "--noopt",
        "--apbs-input={}.in".format(filename_base),
    ]
    
    # If it's NOT a metal (metal_info is None) 
    # and mol2_file exists, add --ligand
    if mol2_file is not None and metal_info is None:
        print(f"Detected regular ligand (non-metal): {mol2_file}")
        args.append("--ligand=" + mol2_file)
    elif metal_info is not None:
         print(f"Detected metal ion: {metal_info['res_name']}. Running PDB2PQR without --ligand.")
    
    # print("Changing directory to:", directory)
    # print('Running PDB2PQR command:', ' '.join(args))
    
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=directory)
    stdout, stderr = p2.communicate()
    
    logger.info("### PDB2PQR ###\n%s", stderr.decode("utf-8").strip())

    # If it WAS a metal ion, we ran *without* --ligand.
    # Append the ion info manually to the PQR file.
    if metal_info is not None:
        print(f"Appending {metal_info['res_name']} to PQR file...")
        try:
            pqr_output_file = os.path.join(directory, filename_base)
            last_atom_id = 0
            last_res_id = 0
            
            # Read the generated PQR file to find the last atom/residue IDs
            with open(pqr_output_file, 'r') as f:
                lines = f.readlines()
            
            for line in reversed(lines):
                if line.startswith(('ATOM', 'HETATM')):
                    try:
                        parts = line.split()
                        last_atom_id = int(parts[1])
                        last_res_id = int(parts[4])
                        break # Found the last atom, stop searching
                    except (ValueError, IndexError, TypeError):
                        # Skip any malformed lines
                        continue
            
            if last_atom_id == 0:
                logger.warning("Could not find last atom ID in PQR file. Defaulting to 0.")

            # Increment IDs for the new HETATM entry
            new_atom_id = last_atom_id + 1
            new_res_id = last_res_id + 1 # Assuming ion is a new residue
            
            mi = metal_info
            
            pqr_line = (
                f"HETATM {new_atom_id} {mi['atom_name']} {mi['res_name']} {new_res_id} "
                f"{mi['x']:.3f} {mi['y']:.3f} {mi['z']:.3f} {mi['charge']:.4f} {mi['radius']:.4f}\n"
            )
            
            # Append the new line to the PQR file
            with open(pqr_output_file, 'a') as f:
                f.write(pqr_line)
            
            print(f"Appended to PQR: {pqr_line.strip()}")
        
        except Exception as e:
            print(f"FATAL: Error appending metal ion to PQR file: {e}")
            logger.error(f"Error appending metal ion to PQR file: {e}")
            raise # Stop execution if we failed to add the ion

    
    apbs_bin = os.environ["APBS_BIN"]
    args = [apbs_bin, filename_base + ".in"]
    # print("Running APBS command:", " ".join(args))
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=directory)
    stdout, stderr = p2.communicate()
    vertfile = open(directory + "/" + filename_base + ".csv", "w")
    for vert in vertices:
        vertfile.write("{},{},{}\n".format(vert[0], vert[1], vert[2]))
    vertfile.close()

    logger.info("### APBS ###\n%s",   stderr.decode("utf-8").strip())

    multivalue_bin = os.environ["MULTIVALUE_BIN"]
    args = [
        multivalue_bin,
        filename_base + ".csv", 
        filename_base + ".dx",
        filename_base + "_out.csv",
    ]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=directory)
    stdout, stderr = p2.communicate()

    # print("Running command:", " ".join(args))

    logger.info("### MULTIVALUE ###\n%s", stderr.decode("utf-8").strip())
    
    chargefile = open(tmp_file_base + "_out.csv")
    charges = numpy.array([0.0] * len(vertices))
    for ix, line in enumerate(chargefile.readlines()):
        charges[ix] = float(line.split(",")[3])

    remove_fn = os.path.join(directory, filename_base)
    os.remove(remove_fn)
    os.remove(remove_fn+'.csv')
    os.remove(remove_fn+'.dx')
    os.remove(remove_fn+'.in')
    os.remove(remove_fn+'_out.csv')

    return charges
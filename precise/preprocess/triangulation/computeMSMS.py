import os
from subprocess import Popen, PIPE

from precise.preprocess.input_output.read_msms import read_msms
from precise.preprocess.triangulation.xyzrn import output_pdb_as_xyzrn, convert_pdb_to_xyzrn
from precise.preprocess.default_config.global_vars import msms_bin 
import random
import sys
import tempfile

# Pablo Gainza LPDI EPFL 2017-2019
# Calls MSMS and returns the vertices.
# Special atoms are atoms with a reduced radius.
def computeMSMS(pdb_file,  protonate=True, ligand_code=None):
    randnum = random.randint(1,10000000)
    file_base = tempfile.gettempdir()+"/msms_"+str(randnum)
    out_xyzrn = file_base+".xyzrn"

    if protonate:
        ligand_list = []
        if ligand_code is not None:
            ligand_list.append(ligand_code)
        # convert_pdb_to_xyzrn(pdb_file, out_xyzrn, keep_hetatms=ligand_list)
        output_pdb_as_xyzrn(pdb_file, out_xyzrn, keep_hetatms=ligand_list)
        # check if the xyzrn file was created and is not empty
        # if not os.path.exists(out_xyzrn):
        #     print("Error - could not create xyzrn file from pdb file.")
        #     sys.exit(1)
        # if os.path.getsize(out_xyzrn) == 0:
        #     print("Error - there is no valid atom in the pdb file.")
        #     sys.exit(1)
    else:
        print("Error - pdb2xyzrn is deprecated.")
        sys.exit(1)
    # Now run MSMS on xyzrn file
    FNULL = open(os.devnull, 'w')
    msms_bin = os.environ["MSMS_BIN"]
    args = [msms_bin, "-density", "3.0", "-hdensity", "3.0", "-probe",\
                    "1.5", "-if",out_xyzrn,"-of",file_base, "-af", file_base]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p2.communicate()

    vertices, faces, normals, names = read_msms(file_base)
    areas = {}
    ses_file = open(file_base+".area")
    next(ses_file) # ignore header line
    for line in ses_file:
        fields = line.split()
        areas[fields[3]] = fields[1]

    # if not os.path.exists(file_base+'.vert'):
    #     print("Error - MSMS failed to run. Check the path to MSMS in global_vars.py")
    #     print("Running command:", " ".join(args))
    #     sys.exit(1)


    # Remove temporary files. 
    os.remove(file_base+'.area')
    os.remove(file_base+'.xyzrn')
    os.remove(file_base+'.vert')
    os.remove(file_base+'.face')
    return vertices, faces, normals, names, areas


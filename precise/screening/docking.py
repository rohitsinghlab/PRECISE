import subprocess
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
from precise.preprocess.input_output.extractPDB import extractPDB


# Default tool paths (can be overridden via environment variables)
UNIDOCK = "unidock"  # Uni-Dock executable
MK_PREP_LIG = "mk_prepare_ligand.py"
MK_PREP_REC = "mk_prepare_receptor.py"


def dock_many_mols(
    smiles_list: list[str],
    pdb_path: Path,
    out_prefix: Path,
    box_center: tuple[float, float, float],
    box_size: tuple[float, float, float] = (22, 22, 22),
) -> list[float]:
    """
    Dock multiple molecules using Uni-Dock.

    Args:
        smiles_list: List of SMILES strings to dock
        pdb_path: Path to receptor PDB file
        out_prefix: Output directory prefix
        box_center: (x, y, z) coordinates of docking box center
        box_size: (x, y, z) dimensions of docking box

    Returns:
        List of docking scores (kcal/mol), NaN for failed dockings
    """
    receptor_pdbqt = out_prefix / "receptor.pdbqt"
    ligand_pdbqt = out_prefix / "ligand.pdbqt"

    if not receptor_pdbqt.exists():
        prepare_receptor(pdb_path, receptor_pdbqt)

    scores: list[float] = []

    for smiles in smiles_list:
        try:
            smiles_to_pdbqt(smiles, ligand_pdbqt)

            output_pdbqt = out_prefix / "ligand_out.pdbqt"
            run_unidock(
                receptor_pdbqt=receptor_pdbqt.name,
                ligand_pdbqt=ligand_pdbqt.name,
                outdir=out_prefix,
                box_center=box_center,
                box_size=box_size,
                output_file=output_pdbqt.name,
            )

            score = get_score(output_pdbqt)
            scores.append(score)

        except Exception as e:
            print(f"\nError docking {smiles}: {e}")
            scores.append(np.nan)

    return scores


def dock_and_score(
    smiles: str,
    pdb_path: Path,
    out_prefix: Path,
    box_center_sdf: Path,
    box_size: tuple[float, float, float] = (22, 22, 22),
) -> float:
    x, y, z = sdf_center(box_center_sdf)
    scores = dock_many_mols([smiles], pdb_path, out_prefix, (x, y, z), box_size)
    return scores[0]


def run_cmd(cmd: str, cwd: Path | None = None) -> int:
    proc = subprocess.Popen(
        cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    proc.communicate()
    return proc.returncode


def sdf_center(sdf_path: Path, mol_index: int = 0) -> tuple[float, float, float]:
    """
    Calculate the geometric center of a molecule in an SDF file.

    Args:
        sdf_path: Path to SDF file
        mol_index: Index of molecule to use (default: 0)

    Returns:
        (x, y, z) coordinates of molecular center
    """
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mols = [m for m in suppl if m is not None]

    if not mols:
        raise RuntimeError(f"Failed to read SDF: {sdf_path}")

    mol = mols[mol_index]
    conf = mol.GetConformer()
    pos = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())], dtype=float
    )
    ctr = pos.mean(axis=0)
    return float(ctr[0]), float(ctr[1]), float(ctr[2])


def smiles_to_pdbqt(smiles: str, out_pdbqt: Path) -> None:
    """
    Convert SMILES to PDBQT format for docking.

    Args:
        smiles: SMILES string
        out_pdbqt: Output PDBQT file path
    """
    out_pdbqt = out_pdbqt.with_suffix(".pdbqt")
    out_sdf = out_pdbqt.with_suffix(".sdf")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError(f"RDKit failed to parse SMILES: {smiles}")

    # Generate 3D structure
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        raise RuntimeError(f"RDKit 3D embedding failed for: {smiles}")

    # Optimize geometry
    AllChem.MMFFOptimizeMolecule(mol)

    # Write SDF
    w = Chem.SDWriter(str(out_sdf))
    w.write(mol)
    w.close()

    # Convert to PDBQT
    code = run_cmd(f"{MK_PREP_LIG} -i {out_sdf} -o {out_pdbqt} --rigid_macrocycles")
    if code != 0 or not out_pdbqt.exists():
        raise RuntimeError("mk_prepare_ligand.py failed")


def prepare_receptor(pdb_path: Path, out_pdbqt: Path) -> None:
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    chainA_pdb = out_pdbqt.parent / f"{out_pdbqt.stem}_chainA.pdb"

    # Extract chain A
    extractPDB(str(pdb_path), str(chainA_pdb), "A", None, None)
    # get out_pdbqt without suffix
    out_pdbqt_no_suffix = out_pdbqt.with_suffix("")
    # Prepare receptor PDBQT
    cmd = f"{MK_PREP_REC} -i {chainA_pdb} -o {out_pdbqt_no_suffix} -p --default_altloc A --allow_bad_res"
    code = run_cmd(cmd)
    if code != 0:
        raise RuntimeError("mk_prepare_receptor.py failed")

    if not out_pdbqt.exists():
        raise RuntimeError(f"Expected PDBQT not written: {out_pdbqt}")


def run_unidock(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    outdir: Path,
    box_center: tuple[float, float, float],
    box_size: tuple[float, float, float] = (22, 22, 22),
    output_file: str = None,
    num_modes: int = 9,
    search_mode: str = "detail",
    scoring: str = "vina",
    extra_args: str = "",
) -> None:
    """
    Run Uni-Dock docking.

    Args:
        receptor_pdbqt: Receptor PDBQT filename (not full path)
        ligand_pdbqt: Ligand PDBQT filename (not full path)
        outdir: Output directory
        box_center: (x, y, z) coordinates of docking box center
        box_size: (x, y, z) dimensions of docking box
        output_file: Output PDBQT filename (optional, defaults based on ligand name)
        num_modes: Maximum number of binding modes to generate
        search_mode: Search mode (fast, balance, or detail)
        scoring: Scoring function (vina, ad4, or vinardo)
        extra_args: Additional command-line arguments
    """
    cx, cy, cz = box_center
    sx, sy, sz = box_size

    out_arg = f"--out {output_file}" if output_file else ""
    cmd = (
        f"cd {outdir}; {UNIDOCK} "
        f"--receptor {receptor_pdbqt} "
        f"--gpu_batch {ligand_pdbqt} "
        f"{out_arg} "
        f"--center_x {cx:.2f} --center_y {cy:.2f} --center_z {cz:.2f} "
        f"--size_x {sx:.2f} --size_y {sy:.2f} --size_z {sz:.2f} "
        f"--search_mode {search_mode} "
        f"--num_modes {num_modes} "
        f"--dir ./ "
        f"--scoring {scoring} "
        f"{extra_args}".strip()
    )

    code = run_cmd(cmd)
    if code != 0:
        raise RuntimeError("Uni-Dock failed")


def get_score(file_path: Path) -> float:
    if not file_path.exists():
        raise FileNotFoundError(f"Output file not found: {file_path}")

    with open(file_path, "r") as f:
        for line in f:
            # Uni-Dock outputs scores in REMARK VINA RESULT lines
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        score = float(parts[3])
                        return score
                    except (ValueError, IndexError):
                        continue

    raise RuntimeError(f"Could not find valid score in {file_path}")

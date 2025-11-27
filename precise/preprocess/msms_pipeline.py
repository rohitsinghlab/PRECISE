import numpy as np
import os
from Bio.PDB import *
import trimesh
from pathlib import Path
import re
import pyvista as pv
import shutil


def find_tool(name: str, env_var: str, search_paths: list[Path]) -> str | None:
    """
    Find a tool by checking:
    1. Environment variable (specific binary path)
    2. Search paths (directories containing the binary)
    """
    # 1. Environment variable (specific tool path)
    if env_var in os.environ:
        path = Path(os.environ[env_var])
        if path.exists() and os.access(path, os.X_OK):
            return str(path.resolve())

    # 2. Search paths
    for p in search_paths:
        target = p / name
        if target.exists() and os.access(target, os.X_OK):
            return str(target.resolve())

    return None


def verify_executables():
    """
    Verify that all required external tools are available.
    Returns True if all tools are found, False otherwise.
    """
    search_paths = []

    # 1. Custom Tools Directory from Env Var
    if "PRECISE_TOOLS_DIR" in os.environ:
        custom_tools_dir = Path(os.environ["PRECISE_TOOLS_DIR"])
        # Check if the user pointed to the root or the bin dir
        # We'll check both the dir itself and a 'bin' subdir
        search_paths.append(custom_tools_dir / "bin")
        search_paths.append(custom_tools_dir)

    # 2. Default User Directory: ~/.precise/tools/bin
    user_tools_dir = Path.home() / ".precise" / "tools" / "bin"
    search_paths.append(user_tools_dir)

    tools_to_find = [
        ("msms", "MSMS_BIN"),
        ("apbs", "APBS_BIN"),
        ("multivalue", "MULTIVALUE_BIN"),
    ]

    # For PDB2PQR
    pdb2pqr_path = shutil.which("pdb2pqr30")

    missing_tools = []

    for tool_name, env_var in tools_to_find:
        # Special handling for pdb2pqr which might have different names
        found_path = find_tool(tool_name, env_var, search_paths)

        if pdb2pqr_path:
            os.environ["PDB2PQR_BIN"] = pdb2pqr_path
        else:
            missing_tools.append("PDB2PQR (Env: PDB2PQR_BIN)")

        if found_path:
            os.environ[env_var] = found_path
        else:
            missing_tools.append(f"{tool_name} (Env: {env_var})")

    if missing_tools:
        print("Error: The following tools are missing:")
        for tool in missing_tools:
            print(f"  - {tool}")
        print(f"Searched in: {[str(p) for p in search_paths]}")
        return False

    return True


from precise.preprocess.triangulation import (
    computeMSMS,
    fix_mesh_by_edge_length,
    computeHydrophobicity,
    computeCharges,
    assignChargesToNewMesh,
    computeAPBS,
    compute_normal,
)
from precise.preprocess.input_output import (
    extractPDB,
    save_ply_trimesh_si,
    protonate,
)
from sklearn.neighbors import KDTree
from collections import namedtuple


def calculate_shape_index_pyvista(vertices, faces):
    """
    Calculate shape index using the PyVista library from vertex and face arrays.

    Args:
        vertices (np.ndarray): NumPy array of shape (N, 3) for vertex positions.
        faces (np.ndarray): NumPy array of shape (M, 3) for face connectivity.

    Returns:
        np.ndarray: An array of shape indices for each vertex.
    """
    faces_pv = np.hstack([3 * np.ones((len(faces), 1), dtype=np.int64), faces])
    pv_mesh = pv.PolyData(vertices, faces_pv)

    K = pv_mesh.curvature(curv_type="Gaussian")
    H = pv_mesh.curvature(curv_type="Mean")

    elem = np.square(H) - K
    elem[elem < 0] = 1e-8
    k1 = H + np.sqrt(elem)
    k2 = H - np.sqrt(elem)

    si = (k1 + k2) / (k1 - k2 + 1e-8)
    si = np.arctan(si) * (2.0 / np.pi)

    return si


def process_protein(
    pdb_filename: Path,
    pdb_id: str,
    chain_id: str,
    out_path: Path,
    protonate_src_loc=None,
):
    """
    Process a protein structure to compute its molecular surface and physicochemical features.

    This pipeline extracts and protonates a specified protein chain, generates its molecular
    surface using MSMS, computes vertex-level properties (electrostatic charges via APBS,
    hydrophobicity, hydrogen bonding propensity), and regularizes the mesh to ~5000 vertices
    suitable for graph neural networks. Interface regions are identified using a buried
    surface area (BSA) approach: vertices from the isolated chain surface that are ≥2.0Å
    away from the full complex surface are marked as interface, where 2.0Å represents the
    standard water probe radius used in protein-protein interaction analysis.

    Args:
            pdb_filename (Path): Path to the input PDB file
            pdb_id (str): PDB identifier for naming files
            chain_id (str): Chain identifier to process
            out_path (Path): Output path for the PLY file with surface and features

    Returns:
            None: Saves PLY file with vertices, faces, normals, charges, hbond, hphob, iface
    """

    if not verify_executables():
        raise RuntimeError(
            "Missing required external tools. Please install them first."
        )

    ligand_code, ligand_chain = None, None

    mol2_file = None
    rdmol = None

    # pdb_dir = pdb_filename.parent
    tmp_dir = out_path.parent / "temp"
    tmp_dir.mkdir(exist_ok=True)

    protonated_file = tmp_dir / f"{pdb_id}_protonated.pdb"
    protonate(pdb_filename, protonated_file, protonate_src_loc)
    pdb_filename = protonated_file

    out_filename1 = tmp_dir / f"{pdb_id}_{chain_id}"
    extractPDB(
        pdb_filename, str(out_filename1) + ".pdb", chain_id, ligand_code, ligand_chain
    )

    # Compute MSMS of surface w/hydrogens,
    vertices1, faces1, normals1, names1, areas1 = computeMSMS(
        str(out_filename1) + ".pdb", protonate=True, ligand_code=ligand_code
    )

    # Compute "charged" vertices
    vertex_hbond = computeCharges(out_filename1, vertices1, names1, ligand_code, rdmol)

    # For each surface residue, assign the hydrophobicity of its amino acid.
    vertex_hphobicity = computeHydrophobicity(names1, ligand_code, rdmol)

    # If protonate = false, recompute MSMS of surface, but without hydrogens (set radius of hydrogens to 0).
    vertices2 = vertices1
    faces2 = faces1

    mesh = trimesh.Trimesh(vertices=vertices2, faces=faces2, process=True)
    regular_mesh = fix_mesh_by_edge_length(mesh, target_edge_length=1)

    # Compute the normals
    vertex_normal = compute_normal(regular_mesh.vertices, regular_mesh.faces)

    vertex_hbond = assignChargesToNewMesh(
        regular_mesh.vertices, vertices1, vertex_hbond
    )

    vertex_hphobicity = assignChargesToNewMesh(
        regular_mesh.vertices, vertices1, vertex_hphobicity
    )
    vertex_charges = computeAPBS(
        regular_mesh.vertices,
        str(out_filename1) + ".pdb",
        str(out_filename1),
        mol2_file,
    )
    iface = np.zeros(len(regular_mesh.vertices))
    # Compute the surface of the entire complex and from that compute the interface.
    v3, f3, _, _, _ = computeMSMS(pdb_filename, protonate=True, ligand_code=ligand_code)
    # Regularize the mesh
    mesh = trimesh.Trimesh(vertices=v3, faces=f3)
    # I believe It is not necessary to regularize the full mesh. This can speed up things by a lot.
    full_regular_mesh = mesh
    # Find the vertices that are in the iface.
    v3 = full_regular_mesh.vertices
    # Find the distance between every vertex in regular_mesh.vertices and those in the full complex.
    kdt = KDTree(v3)
    d, r = kdt.query(regular_mesh.vertices)
    d = np.square(d)  # Square d, because this is how it was in the pyflann version.
    assert len(d) == len(regular_mesh.vertices)
    # 2 here is the radius of a water probe, a BSA-based definition of interface.
    iface_v = np.where(d >= 2.0)[0]
    iface[iface_v] = 1.0
    save_ply_trimesh_si(
        out_path,
        regular_mesh.vertices,
        regular_mesh.faces,
        normals=vertex_normal,
        charges=vertex_charges,
        normalize_charges=True,
        hbond=vertex_hbond,
        hphob=vertex_hphobicity,
        iface=iface,
    )


def process_multimer(
    pdb_filename: Path,
    pdb_id: str,
    chain_ids: list,
    out_path: Path,
    protonate_src_loc=None,
):
    """
    Process a protein multimer (multiple chains) to compute its molecular surface and physicochemical features.

    Args:
            pdb_filename (Path): Path to the input PDB file
            pdb_id (str): PDB identifier for naming files
            chain_ids (list): List of chain identifiers to process
            out_path (Path): Output path for the PLY file with surface and features
            protonate_src_loc (str, optional): Source location for protonation. Defaults to None.
    """

    if not verify_executables():
        raise RuntimeError(
            "Missing required external tools. Please install them first."
        )

    ligand_code, ligand_chain = None, None
    mol2_file = None
    rdmol = None

    tmp_dir = out_path.parent / "temp"
    tmp_dir.mkdir(exist_ok=True)

    protonated_file = tmp_dir / f"{pdb_id}_protonated.pdb"
    protonate(pdb_filename, protonated_file, protonate_src_loc)
    pdb_filename = protonated_file

    # Create a combined chain ID string for filenames
    combined_chain_id = "_".join(chain_ids)
    out_filename1 = tmp_dir / f"{pdb_id}_{combined_chain_id}"

    # Extract multiple chains
    extractPDB(
        pdb_filename, str(out_filename1) + ".pdb", chain_ids, ligand_code, ligand_chain
    )

    # Compute MSMS of surface w/hydrogens
    vertices1, faces1, normals1, names1, areas1 = computeMSMS(
        str(out_filename1) + ".pdb", protonate=True, ligand_code=ligand_code
    )

    # Compute "charged" vertices
    vertex_hbond = computeCharges(out_filename1, vertices1, names1, ligand_code, rdmol)

    # For each surface residue, assign the hydrophobicity of its amino acid.
    vertex_hphobicity = computeHydrophobicity(names1, ligand_code, rdmol)

    # If protonate = false, recompute MSMS of surface, but without hydrogens (set radius of hydrogens to 0).
    vertices2 = vertices1
    faces2 = faces1

    mesh = trimesh.Trimesh(vertices=vertices2, faces=faces2, process=True)
    regular_mesh = fix_mesh_by_edge_length(mesh, target_edge_length=1)

    # Compute the normals
    vertex_normal = compute_normal(regular_mesh.vertices, regular_mesh.faces)

    vertex_hbond = assignChargesToNewMesh(
        regular_mesh.vertices, vertices1, vertex_hbond
    )

    vertex_hphobicity = assignChargesToNewMesh(
        regular_mesh.vertices, vertices1, vertex_hphobicity
    )
    vertex_charges = computeAPBS(
        regular_mesh.vertices,
        str(out_filename1) + ".pdb",
        str(out_filename1),
        mol2_file,
    )

    # Interface calculation for multimer
    # For a multimer, the "interface" might be defined against other chains NOT in the multimer
    # OR if the input PDB contains ONLY the multimer, then there is no interface to compute against "full complex"
    # if the "full complex" is the same as the multimer.
    # However, usually process_protein assumes pdb_filename is the FULL complex.
    # So we compute interface against the rest of the complex.

    iface = np.zeros(len(regular_mesh.vertices))
    # Compute the surface of the entire complex
    v3, f3, _, _, _ = computeMSMS(pdb_filename, protonate=True, ligand_code=ligand_code)

    # Regularize the mesh (optional, but kept for consistency with process_protein)
    mesh = trimesh.Trimesh(vertices=v3, faces=f3)
    full_regular_mesh = mesh

    # Find the vertices that are in the iface.
    v3 = full_regular_mesh.vertices
    # Find the distance between every vertex in regular_mesh.vertices and those in the full complex.
    kdt = KDTree(v3)
    d, r = kdt.query(regular_mesh.vertices)
    d = np.square(d)
    assert len(d) == len(regular_mesh.vertices)
    # 2 here is the radius of a water probe, a BSA-based definition of interface.
    # Note: If the extracted multimer IS the full complex, d will be 0 everywhere, so iface will be 0.
    # This is expected behavior if there's no "other" chain to interact with.
    iface_v = np.where(d >= 2.0)[0]
    iface[iface_v] = 1.0

    save_ply_trimesh_si(
        out_path,
        regular_mesh.vertices,
        regular_mesh.faces,
        normals=vertex_normal,
        charges=vertex_charges,
        normalize_charges=True,
        hbond=vertex_hbond,
        hphob=vertex_hphobicity,
        iface=iface,
    )


def parse_system_id(system_id: str):
    """
    Parse a PLINDER system ID of the form:
    <PDB ID>__<biological assembly>__<receptor chain ID>__<ligand chain ID>
    Example: '7eek__1__1.A__1.I'
    Returns a dict with components.
    """
    pattern = r"^([^_]+)__([^_]+)__([^_]+)__([^_]+)$"
    match = re.match(pattern, system_id)

    if not match:
        raise ValueError(f"Invalid system_id format: {system_id}")

    pdb_id, biological_assembly, receptor_chain_id, ligand_chain_id = match.groups()

    SystemInfo = namedtuple(
        "SystemInfo",
        ["pdb_id", "biological_assembly", "receptor_chain_id", "ligand_chain_id"],
    )

    return SystemInfo(
        pdb_id=pdb_id,
        biological_assembly=biological_assembly,
        receptor_chain_id=receptor_chain_id,
        ligand_chain_id=ligand_chain_id,
    )

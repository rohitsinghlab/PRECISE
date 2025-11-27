# import pymesh
import trimesh
import numpy as np
import logging
import pyvista as pv

# Initialize logger for this module
logger = logging.getLogger(__name__)

"""
read_ply.py: Save a ply file to disk using pymesh and load the attributes used by MaSIF. 
Pablo Gainza - LPDI STI EPFL 2019
Released under an Apache License 2.0
"""

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

def save_ply(
    filename,
    vertices,
    faces=[],
    normals=None,
    charges=None,
    vertex_cb=None,
    hbond=None,
    hphob=None,
    iface=None,
    normalize_charges=False,
):
    """ Save vertices, mesh in ply format.
        vertices: coordinates of vertices
        faces: mesh
    """
    # mesh = pymesh.form_mesh(vertices, faces)
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
    if normals is not None:
        n1 = normals[:, 0]
        n2 = normals[:, 1]
        n3 = normals[:, 2]
        mesh.add_attribute("vertex_nx")
        mesh.set_attribute("vertex_nx", n1)
        mesh.add_attribute("vertex_ny")
        mesh.set_attribute("vertex_ny", n2)
        mesh.add_attribute("vertex_nz")
        mesh.set_attribute("vertex_nz", n3)
    if charges is not None:
        mesh.add_attribute("charge")
        if normalize_charges:
            charges = charges / 10
        mesh.set_attribute("charge", charges)
    if hbond is not None:
        mesh.add_attribute("hbond")
        mesh.set_attribute("hbond", hbond)
    if vertex_cb is not None:
        mesh.add_attribute("vertex_cb")
        mesh.set_attribute("vertex_cb", vertex_cb)
    if hphob is not None:
        mesh.add_attribute("vertex_hphob")
        mesh.set_attribute("vertex_hphob", hphob)
    if iface is not None:
        mesh.add_attribute("vertex_iface")
        mesh.set_attribute("vertex_iface", iface)

    # pymesh.save_mesh(
    #     filename, mesh, *mesh.get_attribute_names(), use_float=True, ascii=True
    # )
    mesh.export(filename, encoding='ascii')
    logger.info(f"Mesh saved to {filename}")

def save_ply_trimesh(
    filename,
    vertices,
    faces,
    normals=None,
    charges=None,
    vertex_cb=None,
    hbond=None,
    hphob=None,
    iface=None,
    normalize_charges=False,
):
    """
    Saves a mesh and its custom vertex attributes to a PLY file using the TriMesh library.
    This version adopts the most robust method: preparing all attributes first and
    passing them to the Trimesh constructor directly.

    Args:
        filename (str): The output file path.
        vertices (np.ndarray): (N, 3) float array of vertex positions.
        faces (np.ndarray): (M, 3) int array of face indices.
        # ... (other args are the same)
    """
    # **THE NEW STRATEGY IS HERE**
    # 1. First, prepare a dictionary for all custom vertex attributes.
    custom_attributes = {}

    if charges is not None and len(charges) == len(vertices):
        if normalize_charges:
            charges = charges / 10.0
        custom_attributes['charge'] = charges.astype(np.float32)

    if vertex_cb is not None and len(vertex_cb) == len(vertices):
        custom_attributes['vertex_cb'] = vertex_cb.astype(np.float32)

    if hbond is not None and len(hbond) == len(vertices):
        custom_attributes['hbond'] = hbond.astype(np.float32)

    if hphob is not None and len(hphob) == len(vertices):
        custom_attributes['hphob'] = hphob.astype(np.float32)

    if iface is not None and len(iface) == len(vertices):
        custom_attributes['iface'] = iface.astype(np.int32)

    # Also handle normals, preparing them to be passed to the constructor.
    # Ensure normals are None if they are not provided or have the wrong shape.
    vertex_normals = normals if (normals is not None and len(normals) == len(vertices)) else None

    # --- Debug logging ---
    logger.debug("--- Verifying attributes before creating mesh object ---")
    logger.debug(f"File: {filename}")
    if vertex_normals is not None:
         logger.debug(f"Prepared vertex normals with shape: {vertex_normals.shape}")
    else:
        logger.debug("Vertex normals are NOT being provided to the constructor.")
    
    if custom_attributes:
        logger.debug(f"Prepared custom attributes: {list(custom_attributes.keys())}")
        for key, value in custom_attributes.items():
            logger.debug(f"  - {key}: dtype={value.dtype}, shape={value.shape}")
    else:
        logger.debug("No custom attributes are being provided to the constructor.")
    logger.debug("-------------------------------------------------------")

    # 2. Create the mesh object in one go, passing all data to the constructor.
    # This is the most reliable way to ensure all attributes are recognized.
    mesh = trimesh.Trimesh(
        vertices=vertices,
        faces=faces,
        vertex_normals=vertex_normals,
        vertex_attributes=custom_attributes,
        process=False  # Set process=False to avoid trimesh modifying the data
    )

    # 3. Export the mesh.
    mesh.export(filename, file_type='ply', encoding='ascii')
    
    logger.info(f"Mesh saved to {filename}")

def save_ply_trimesh_si(
    filename,
    vertices,
    faces,
    normals=None,
    charges=None,
    vertex_cb=None,
    hbond=None,
    hphob=None,
    iface=None,
    shape_index=True,
    normalize_charges=False,
):
    """
    Saves a mesh and its custom vertex attributes to a PLY file using the TriMesh library.
    This version adopts the most robust method: preparing all attributes first and
    passing them to the Trimesh constructor directly.

    Args:
        filename (str): The output file path.
        vertices (np.ndarray): (N, 3) float array of vertex positions.
        faces (np.ndarray): (M, 3) int array of face indices.
        # ... (other args are the same)
    """
    # **THE NEW STRATEGY IS HERE**
    # 1. First, prepare a dictionary for all custom vertex attributes.
    custom_attributes = {}

    if charges is not None and len(charges) == len(vertices):
        if normalize_charges:
            charges = charges / 10.0
        custom_attributes['charge'] = charges.astype(np.float32)

    if vertex_cb is not None and len(vertex_cb) == len(vertices):
        custom_attributes['vertex_cb'] = vertex_cb.astype(np.float32)

    if hbond is not None and len(hbond) == len(vertices):
        custom_attributes['hbond'] = hbond.astype(np.float32)

    if hphob is not None and len(hphob) == len(vertices):
        custom_attributes['hphob'] = hphob.astype(np.float32)

    if iface is not None and len(iface) == len(vertices):
        custom_attributes['iface'] = iface.astype(np.int32)

    if shape_index:
        si = calculate_shape_index_pyvista(vertices, faces)
        custom_attributes['shape_index'] = si.astype(np.float32)

    # Also handle normals, preparing them to be passed to the constructor.
    # Ensure normals are None if they are not provided or have the wrong shape.
    vertex_normals = normals if (normals is not None and len(normals) == len(vertices)) else None

    # --- Debug logging ---
    logger.debug("--- Verifying attributes before creating mesh object ---")
    logger.debug(f"File: {filename}")
    if vertex_normals is not None:
         logger.debug(f"Prepared vertex normals with shape: {vertex_normals.shape}")
    else:
        logger.debug("Vertex normals are NOT being provided to the constructor.")
    
    if custom_attributes:
        logger.debug(f"Prepared custom attributes: {list(custom_attributes.keys())}")
        for key, value in custom_attributes.items():
            logger.debug(f"  - {key}: dtype={value.dtype}, shape={value.shape}")
    else:
        logger.debug("No custom attributes are being provided to the constructor.")
    logger.debug("-------------------------------------------------------")

    # 2. Create the mesh object in one go, passing all data to the constructor.
    # This is the most reliable way to ensure all attributes are recognized.
    mesh = trimesh.Trimesh(
        vertices=vertices,
        faces=faces,
        vertex_normals=vertex_normals,
        vertex_attributes=custom_attributes,
        process=False  # Set process=False to avoid trimesh modifying the data
    )

    # 3. Export the mesh.
    mesh.export(filename, file_type='ply', encoding='ascii')
    
    logger.info(f"Mesh saved to {filename}")
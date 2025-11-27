import trimesh
import numpy as np
from numpy.linalg import norm
import logging
import pymeshlab

# Initialize logger for this module
logger = logging.getLogger(__name__)

# try:
#     import pymeshfix
#     PYMESHFIX_AVAILABLE = True
# except ImportError:
#     PYMESHFIX_AVAILABLE = False
#     logger.warning("PyMeshFix not available, using trimesh-only methods")


def fix_protein_mesh(mesh, target_vertices=3000, preserve_features=True, min_vertices=2000, max_vertices=4000):
    """
    Fix protein surface mesh while preserving biological features and reducing vertices.
    
    Parameters:
    -----------
    mesh : trimesh.Trimesh
        Input protein surface mesh
    target_vertices : int
        Target number of vertices (default 3000)
    preserve_features : bool
        If True, preserves small features important for pocket detection
    min_vertices : int
        Minimum acceptable vertices (default 2000)
    max_vertices : int
        Maximum acceptable vertices (default 4000)
    
    Returns:
    --------
    trimesh.Trimesh
        Fixed mesh with controlled vertex count
    """
    
    logger.info(f"Initial mesh: {len(mesh.vertices)} vertices, {len(mesh.faces)} faces")
    logger.info(f"Target: {target_vertices} vertices (range: {min_vertices}-{max_vertices})")
    
    # Step 1: Initial cleanup
    # Remove duplicate vertices with reasonable tolerance
    mesh.merge_vertices(digits_vertex=3)
    
    # Remove degenerate faces
    mesh.remove_degenerate_faces()
    
    # Remove duplicate faces
    mesh.remove_duplicate_faces()
    
    logger.info(f"After cleanup: {len(mesh.vertices)} vertices, {len(mesh.faces)} faces")
    
    # Step 2: Fix mesh topology issues if needed
    if PYMESHFIX_AVAILABLE and not mesh.is_watertight:
        try:
            meshfix = pymeshfix.MeshFix(mesh.vertices, mesh.faces)
            meshfix.repair(
                verbose=False,
                joincomp=True,  # Join components
                remove_smallest_components=True  # Remove tiny artifacts
            )
            vertices, faces = meshfix.return_arrays()
            mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
            logger.info("Applied PyMeshFix repair")
        except Exception as e:
            logger.warning(f"PyMeshFix repair failed: {e}, continuing with trimesh")
    
    # Step 3: Fill small holes only
    if not mesh.is_watertight:
        # Check boundary loops
        boundary_loops = mesh.facets_boundary
        if len(boundary_loops) > 0 and len(boundary_loops) < 10:
            # Only fill if we have a few small holes
            max_hole_edges = max(len(loop) for loop in boundary_loops)
            if max_hole_edges < 20:  # Small holes only
                mesh.fill_holes()
                logger.info("Filled small holes")
    
    # Step 4: Aggressive vertex reduction to target range
    current_vertices = len(mesh.vertices)
    
    if current_vertices > max_vertices:
        # Calculate target faces (roughly 2 faces per vertex)
        target_faces = target_vertices * 2
        
        # Multi-step reduction for better quality
        reduction_steps = []
        
        if current_vertices > target_vertices * 10:
            # Very high resolution - do aggressive first pass
            reduction_steps.append(target_vertices * 5)
        
        if current_vertices > target_vertices * 3:
            # Medium resolution - intermediate step
            reduction_steps.append(target_vertices * 2)
        
        # Final target
        reduction_steps.append(target_faces)
        
        for step_target in reduction_steps:
            if len(mesh.faces) <= step_target:
                break
                
            logger.info(f"Reducing to ~{step_target//2} vertices...")
            
            # Use quadric decimation for feature-preserving simplification
            try:
                mesh = mesh.simplify_quadric_decimation(
                    face_count=int(step_target),
                    aggression=7 if preserve_features else 8  # 1-10, higher = more aggressive
                )
            except:
                # Fallback to simple decimation if quadric fails
                target_fraction = step_target / len(mesh.faces)
                mesh = mesh.decimate(face_count=int(step_target))
            
            # Clean up after decimation
            mesh.remove_degenerate_faces()
            mesh.remove_unreferenced_vertices()
            
            logger.info(f"  -> {len(mesh.vertices)} vertices, {len(mesh.faces)} faces")
    
    elif current_vertices < min_vertices:
        # Need to subdivide to increase vertex count
        subdivisions_needed = int(np.ceil(np.log2(min_vertices / current_vertices)))
        subdivisions_needed = min(subdivisions_needed, 2)  # Max 2 subdivisions
        
        logger.info(f"Subdividing {subdivisions_needed} times to increase vertex count...")
        mesh = mesh.subdivide(subdivisions_needed)
        
        # Then reduce to target if needed
        if len(mesh.vertices) > max_vertices:
            target_faces = target_vertices * 2
            mesh = mesh.simplify_quadric_decimation(
                face_count=target_faces,
                aggression=7 if preserve_features else 8
            )
    
    # Step 5: Final quality improvements
    
    # Remove bad triangles (but be conservative for protein surfaces)
    face_angles = mesh.face_angles
    
    # Only remove extremely bad triangles
    max_angle_threshold = np.radians(175) if preserve_features else np.radians(170)
    min_angle_threshold = np.radians(5)  # Remove slivers
    
    bad_faces_mask = np.any(face_angles > max_angle_threshold, axis=1) | \
                     np.any(face_angles < min_angle_threshold, axis=1)
    
    if np.any(bad_faces_mask) and np.sum(~bad_faces_mask) > 100:  # Keep at least 100 faces
        num_bad = np.sum(bad_faces_mask)
        if num_bad < len(mesh.faces) * 0.1:  # Only remove if <10% are bad
            logger.info(f"Removing {num_bad} bad triangles")
            good_faces = mesh.faces[~bad_faces_mask]
            mesh = trimesh.Trimesh(vertices=mesh.vertices, faces=good_faces)
            mesh.remove_unreferenced_vertices()
    
    # Light smoothing to improve mesh quality
    if preserve_features:
        # Very light smoothing to preserve pockets
        mesh = trimesh.smoothing.filter_laplacian(
            mesh, 
            iterations=1,
            lamb=0.3  # Lower = less aggressive
        )
    else:
        mesh = trimesh.smoothing.filter_laplacian(
            mesh,
            iterations=2,
            lamb=0.5
        )
    
    # Final cleanup
    mesh.fix_normals()
    mesh.remove_duplicate_faces()
    mesh.remove_unreferenced_vertices()
    
    # Step 6: Final adjustment if still outside target range
    final_vertices = len(mesh.vertices)
    
    if final_vertices > max_vertices:
        # One more reduction
        target_faces = int(target_vertices * 2)
        mesh = mesh.simplify_quadric_decimation(face_count=target_faces, aggression=7)
        mesh.remove_unreferenced_vertices()
    elif final_vertices < min_vertices and final_vertices > min_vertices * 0.8:
        # Close enough, leave as is
        pass
    
    # Print final statistics
    logger.info(f"Final mesh statistics:")
    logger.info(f"  Vertices: {len(mesh.vertices)}")
    logger.info(f"  Faces: {len(mesh.faces)}")
    logger.info(f"  Watertight: {mesh.is_watertight}")
    logger.info(f"  Volume: {mesh.volume:.3f} cubic units")
    logger.info(f"  Surface area: {mesh.area:.3f} square units")
    
    if len(mesh.vertices) < min_vertices * 0.8 or len(mesh.vertices) > max_vertices * 1.2:
        logger.warning(f"Could not achieve target vertex range!")
    
    return mesh

def fix_mesh_by_edge_length(mesh, target_edge_length):
    """
    Repairs a mesh and remeshes it to a target edge length using PyMeshLab.
    
    Args:
        mesh (trimesh.Trimesh): The input mesh.
        target_edge_length (float): The target length for the mesh edges (e.g., 1.0 for 1 Angstrom).
        
    Returns:
        trimesh.Trimesh: The processed mesh.
    """
    
    # Create a MeshSet and add the mesh from trimesh
    print('Starting pymeshlab processing to fix mesh...')
    ms = pymeshlab.MeshSet()
    ms.add_mesh(pymeshlab.Mesh(mesh.vertices, mesh.faces))
    
    ms.meshing_remove_unreferenced_vertices()
    ms.meshing_remove_duplicate_faces()
    ms.apply_filter('meshing_remove_null_faces')
    ms.meshing_repair_non_manifold_edges()
    ms.meshing_close_holes(maxholesize=100)
    ms.meshing_remove_connected_component_by_diameter(mincomponentdiag=pymeshlab.PureValue(1.0))
    
    ms.meshing_isotropic_explicit_remeshing(
        targetlen=pymeshlab.PureValue(target_edge_length),
        iterations=10
    )
    
    ms.apply_filter('apply_coord_laplacian_smoothing', stepsmoothnum=2, cotangentweight=True)

    processed_mesh = ms.current_mesh()
    vertices = processed_mesh.vertex_matrix()
    faces = processed_mesh.face_matrix()
    
    return trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

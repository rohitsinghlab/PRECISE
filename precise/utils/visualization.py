from pathlib import Path
from typing import Optional

import numpy as np
from plyfile import PlyData, PlyElement


def values_to_colormap(values: np.ndarray, cmap_type: str = "rdbu") -> np.ndarray:
    """
    Map normalized values in [0,1] to RGB colors.

    Args:
        values: Array of values in [0, 1]
        cmap_type: Color map type - 'rdbu' (red-white-blue), 'coolwarm', 'ylorbr', 'grayscale', 'purple_white'

    Returns:
        RGB array of shape [N, 3] with uint8 values
    """
    v = np.clip(values, 0.0, 1.0).astype(np.float32)

    if cmap_type == "rdbu":
        # Red (negative) -> White (neutral) -> Blue (positive)
        # For values: 0=red (239, 67, 71), 0.5=white, 1=blue (95, 201, 219)
        red_rgb = np.array([239.0 / 255.0, 67.0 / 255.0, 71.0 / 255.0])
        white_rgb = np.array([1.0, 1.0, 1.0])
        blue_rgb = np.array([95.0 / 255.0, 201.0 / 255.0, 219.0 / 255.0])

        # Vectorized interpolation
        mask_low = v < 0.5
        mask_high = v >= 0.5

        r = np.zeros_like(v)
        g = np.zeros_like(v)
        b = np.zeros_like(v)

        # Interpolate from red to white (v: 0 -> 0.5)
        t_low = v[mask_low] * 2.0  # Scale to [0, 1]
        r[mask_low] = red_rgb[0] + t_low * (white_rgb[0] - red_rgb[0])
        g[mask_low] = red_rgb[1] + t_low * (white_rgb[1] - red_rgb[1])
        b[mask_low] = red_rgb[2] + t_low * (white_rgb[2] - red_rgb[2])

        # Interpolate from white to blue (v: 0.5 -> 1)
        t_high = (v[mask_high] - 0.5) * 2.0  # Scale to [0, 1]
        r[mask_high] = white_rgb[0] + t_high * (blue_rgb[0] - white_rgb[0])
        g[mask_high] = white_rgb[1] + t_high * (blue_rgb[1] - white_rgb[1])
        b[mask_high] = white_rgb[2] + t_high * (blue_rgb[2] - white_rgb[2])

    elif cmap_type == "coolwarm":
        # Cool (blue) -> Warm (red)
        # For values: 0=blue (95, 201, 219), 0.5=white, 1=red (239, 67, 71)
        blue_rgb = np.array([95.0 / 255.0, 201.0 / 255.0, 219.0 / 255.0])
        white_rgb = np.array([1.0, 1.0, 1.0])
        red_rgb = np.array([239.0 / 255.0, 67.0 / 255.0, 71.0 / 255.0])

        # Vectorized interpolation
        mask_low = v < 0.5
        mask_high = v >= 0.5

        r = np.zeros_like(v)
        g = np.zeros_like(v)
        b = np.zeros_like(v)

        # Interpolate from blue to white (v: 0 -> 0.5)
        t_low = v[mask_low] * 2.0  # Scale to [0, 1]
        r[mask_low] = blue_rgb[0] + t_low * (white_rgb[0] - blue_rgb[0])
        g[mask_low] = blue_rgb[1] + t_low * (white_rgb[1] - blue_rgb[1])
        b[mask_low] = blue_rgb[2] + t_low * (white_rgb[2] - blue_rgb[2])

        # Interpolate from white to red (v: 0.5 -> 1)
        t_high = (v[mask_high] - 0.5) * 2.0  # Scale to [0, 1]
        r[mask_high] = white_rgb[0] + t_high * (red_rgb[0] - white_rgb[0])
        g[mask_high] = white_rgb[1] + t_high * (red_rgb[1] - white_rgb[1])
        b[mask_high] = white_rgb[2] + t_high * (red_rgb[2] - white_rgb[2])

    elif cmap_type == "ylorbr":
        # Yellow -> Orange -> Brown
        # 0=light yellow, 1=dark brown
        r = 1.0 - 0.3 * v  # Start bright, darken
        g = 1.0 - 0.6 * v  # Fade faster
        b = 1.0 - 0.85 * v  # Fade most

    elif cmap_type == "grayscale":
        # Black to White gradient
        # 0=black, 1=white
        r = v
        g = v
        b = v

    elif cmap_type == "purple_white":
        # Purple (#6062a7) to White gradient
        # 0=purple (96, 98, 167), 1=white (255, 255, 255)
        purple_rgb = np.array([96.0 / 255.0, 98.0 / 255.0, 167.0 / 255.0])
        white_rgb = np.array([1.0, 1.0, 1.0])

        r = purple_rgb[0] + v * (white_rgb[0] - purple_rgb[0])
        g = purple_rgb[1] + v * (white_rgb[1] - purple_rgb[1])
        b = purple_rgb[2] + v * (white_rgb[2] - purple_rgb[2])

    else:
        # Default: blue -> green -> red gradient
        r = np.clip(2.0 * v - 0.0, 0.0, 1.0)
        g = 1.0 - np.abs(2.0 * v - 1.0)
        b = np.clip(2.0 * (1.0 - v) - 0.0, 0.0, 1.0)

    rgb = np.stack([r, g, b], axis=1) * 255.0
    return rgb.astype(np.uint8)


def write_ply_with_colors(
    out_path: Path,
    positions: np.ndarray,
    faces: Optional[np.ndarray],
    colors: np.ndarray,
) -> None:
    """Write PLY file with vertex colors."""
    num_vertices = positions.shape[0]
    if colors.shape != (num_vertices, 3):
        raise ValueError("colors must have shape [N, 3] to match vertex positions")

    vertex_dtype = [
        ("x", "f4"),
        ("y", "f4"),
        ("z", "f4"),
        ("red", "u1"),
        ("green", "u1"),
        ("blue", "u1"),
    ]

    vertex_data = np.empty(num_vertices, dtype=vertex_dtype)
    vertex_data["x"] = positions[:, 0]
    vertex_data["y"] = positions[:, 1]
    vertex_data["z"] = positions[:, 2]
    vertex_data["red"] = colors[:, 0]
    vertex_data["green"] = colors[:, 1]
    vertex_data["blue"] = colors[:, 2]

    elements = [PlyElement.describe(vertex_data, "vertex")]

    if faces is not None and faces.size > 0:
        face_array = np.array(
            [(list(face),) for face in faces], dtype=[("vertex_indices", "i4", (3,))]
        )
        elements.append(PlyElement.describe(face_array, "face"))

    PlyData(elements, text=True).write(str(out_path))


def minmax_normalize(values: np.ndarray) -> np.ndarray:
    """Normalize values to [0, 1] range."""
    if values.size == 0:
        return values
    vmin = float(values.min())
    vmax = float(values.max())
    if vmax - vmin < 1e-8:
        return np.ones_like(values) * 0.5  # Return middle value if constant
    return (values - vmin) / (vmax - vmin)


def create_colored_surfaces(
    ply_path: Path,
    output_dir: Optional[Path] = None,
    properties: Optional[list] = None,
) -> dict:
    """
    Create colored surface visualizations for specified properties.

    Args:
        ply_path: Path to input PLY file
        output_dir: Output directory (defaults to same directory as input)
        properties: List of properties to visualize. Options:
            - 'charge': Poisson-Boltzmann (electrostatics)
            - 'hbond': Hydrophilicity (H-bond donor/acceptor)
            - 'hphob': Hydrophobicity
            - 'shape_index': Shape index (curvature)
            - 'all': All properties
            If None, creates all properties.

    Returns:
        Dictionary mapping property names to output file paths
    """
    # Read PLY file
    plydata = PlyData.read(str(ply_path))

    # Extract vertex positions
    vertices = np.column_stack(
        [plydata["vertex"]["x"], plydata["vertex"]["y"], plydata["vertex"]["z"]]
    )

    # Extract feature fields
    pb_raw = np.array(plydata["vertex"]["charge"])  # Poisson-Boltzmann (electrostatics)
    hb = np.array(plydata["vertex"]["hbond"])  # Hydrophilicity (H-bond donor/acceptor)
    hp = np.array(plydata["vertex"]["hphob"])  # Hydrophobicity
    si = np.array(plydata["vertex"]["shape_index"])  # Shape index (curvature)

    # Extract faces
    faces = np.array([face for face in plydata["face"]["vertex_indices"]])

    # Normalize Poisson-Boltzmann to [0, 1] for coloring
    pb_lower, pb_upper = -3.0, 3.0
    pb = np.clip(pb_raw, pb_lower, pb_upper)
    pb_normalized = (pb - pb_lower) / (pb_upper - pb_lower)  # Now in [0, 1]

    # Normalize other properties to [0, 1]
    hp_normalized = minmax_normalize(hp)
    hb_normalized = minmax_normalize(hb)
    si_normalized = minmax_normalize(si)

    # Define all available features
    all_features = {
        "charge": {
            "name": "Poisson-Boltzmann (Electrostatics)",
            "values": pb_normalized,
            "cmap": "rdbu",
        },
        "hbond": {
            "name": "Hydrophilicity (H-bond)",
            "values": hb_normalized,
            "cmap": "coolwarm",
        },
        "hphob": {
            "name": "Hydrophobicity",
            "values": hp_normalized,
            "cmap": "ylorbr",
        },
        "shape_index": {
            "name": "Shape Index (Curvature)",
            "values": si_normalized,
            "cmap": "purple_white",
        },
    }

    # Determine which properties to visualize
    if properties is None or "all" in properties:
        features_to_use = all_features
    else:
        features_to_use = {k: v for k, v in all_features.items() if k in properties}

    if output_dir is None:
        output_dir = ply_path.parent
    else:
        output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    output_files = {}

    # Create colored surface for each feature
    for prop_key, feature_info in features_to_use.items():
        feature_values = feature_info["values"]
        cmap = feature_info["cmap"]
        feature_name = feature_info["name"]

        # Convert values to RGB colors
        colors = values_to_colormap(feature_values, cmap_type=cmap)

        # Save as PLY file with vertex colors
        filename = (
            feature_name.replace(" ", "_").replace("(", "").replace(")", "").lower()
        )
        output_path = output_dir / f"{ply_path.stem}_{filename}.ply"
        write_ply_with_colors(output_path, vertices, faces, colors)

        output_files[prop_key] = output_path

    return output_files


def create_binding_score_surface(
    ply_path: Path,
    vertex_scores: np.ndarray,
    output_path: Optional[Path] = None,
    cmap_type: str = "coolwarm",
) -> Path:
    """
    Create a colored surface based on binding scores.

    Args:
        ply_path: Path to input PLY file
        vertex_scores: Array of binding scores (one per vertex), in [0, 1]
        output_path: Output path for colored PLY (defaults to ply_path with _binding_scores suffix)
        cmap_type: Color map type (default: 'coolwarm')

    Returns:
        Path to output PLY file
    """
    # Read PLY file
    plydata = PlyData.read(str(ply_path))

    # Extract vertex positions
    vertices = np.column_stack(
        [plydata["vertex"]["x"], plydata["vertex"]["y"], plydata["vertex"]["z"]]
    )

    # Extract faces
    faces = np.array([face for face in plydata["face"]["vertex_indices"]])

    # Validate vertex_scores
    if len(vertex_scores) != len(vertices):
        raise ValueError(
            f"vertex_scores length ({len(vertex_scores)}) must match number of vertices ({len(vertices)})"
        )

    # Normalize scores to [0, 1]
    scores_normalized = minmax_normalize(vertex_scores)

    # Convert to colors
    colors = values_to_colormap(scores_normalized, cmap_type=cmap_type)

    # Determine output path
    if output_path is None:
        output_path = ply_path.parent / f"{ply_path.stem}_binding_scores.ply"
    else:
        output_path = Path(output_path)

    # Write colored PLY
    write_ply_with_colors(output_path, vertices, faces, colors)

    return output_path

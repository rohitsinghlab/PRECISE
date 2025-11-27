import numpy as np
import pandas as pd
import random
import time
from pathlib import Path
from typing import Callable, List, Tuple, Dict
from rich.console import Console, Group
from rich.panel import Panel
from rich.live import Live
from rich.columns import Columns
from rich.spinner import Spinner
from rich.text import Text


def get_children(Z: np.ndarray, cluster_id: int) -> list[int]:
    n = Z.shape[0] + 1
    if cluster_id < n:
        return []  # Leaf node has no children

    left, right = Z[int(cluster_id - n), 0], Z[int(cluster_id - n), 1]
    return [int(left), int(right)]


def get_all_descendants(Z: np.ndarray, cluster_id: int) -> list[int]:
    n = Z.shape[0] + 1
    if cluster_id < n:
        return [int(cluster_id)]  # Leaf node

    left, right = get_children(Z, cluster_id)
    return get_all_descendants(Z, left) + get_all_descendants(Z, right)


def get_all_children(Z: np.ndarray, ids: List[int]) -> Dict[int, Dict[str, List[int]]]:
    results = {}
    for cid in ids:
        results[cid] = {
            "immediate": get_children(Z, cid),
            "descendants": get_all_descendants(Z, cid),
        }
    return results


def get_path_to_ancestor(Z: np.ndarray, leaf_id: int, ancestor: int) -> list[int]:
    n = Z.shape[0] + 1
    current = leaf_id
    path = [current]

    # Move upward until we reach the root (2n-2) or the ancestor
    while current < 2 * n - 2:
        # Find which row merges this cluster
        rows = np.where((Z[:, 0] == current) | (Z[:, 1] == current))[0]
        if len(rows) == 0:
            break  # Reached the root (no parent)

        r = rows[0]
        parent = n + r

        path.append(int(parent))
        if parent == ancestor:
            return path
        current = parent

    return path


def cluster_and_score_per_level(
    Z: np.ndarray,
    scoring_func: Callable[[int], float],
    filtering_func: Callable[[List[Tuple[int, float]], int], List[int]],
    nodelist: List[int],
    score_assignment_dict: Dict[int, Tuple[int, float]],
    depth: int,
    csize_thres: int = 5,
    B_MAX: float = 10000,
) -> List[int]:
    n = Z.shape[0] + 1

    def is_leaf(d: int) -> bool:
        return d < n

    childmap = get_all_children(Z, nodelist)
    child_nodelist = []

    for k, dcmap in childmap.items():
        leaf_descendents = dcmap["descendants"]

        if len(leaf_descendents) <= csize_thres:
            for leaf in leaf_descendents:
                if leaf not in score_assignment_dict:
                    score_assignment_dict[leaf] = (leaf, scoring_func(leaf))
        else:
            d0, d1 = dcmap["immediate"]

            for d in [d0, d1]:
                if is_leaf(d) and d not in score_assignment_dict:
                    score_assignment_dict[d] = (d, scoring_func(d))
                else:
                    desc = get_all_children(Z, [d])
                    assert len(desc) > 0, "Non-leaf node must have children"

                    leaf = random.sample(desc[d]["descendants"], 1)[0]
                    ancestors_leaf_to_d = get_path_to_ancestor(Z, leaf, d)

                    if leaf not in score_assignment_dict:
                        score = scoring_func(leaf)
                        score_assignment_dict[leaf] = (leaf, score)
                    else:
                        score = score_assignment_dict[leaf][1]

                    if np.isnan(score):
                        score = B_MAX

                    assert d in ancestors_leaf_to_d, "Path must contain ancestor"
                    for a in ancestors_leaf_to_d:
                        if a not in score_assignment_dict:
                            score_assignment_dict[a] = (leaf, score)
                        else:
                            p_leaf, p_score = score_assignment_dict[a]
                            if p_score > score:  # Update if we found better score
                                score_assignment_dict[a] = (leaf, score)

                    child_nodelist.append((d, score_assignment_dict[d][1]))

    # Filter nodes for next level
    return filtering_func(child_nodelist, depth)


def cluster_and_score(
    Z: np.ndarray,
    scoring_func: Callable[[int], float],
    filtering_func: Callable[[List[Tuple[int, float]], int], List[int]],
    csize_thres: int = 5,
    B_MAX: float = 10000,
    verbose: bool = True,
    outprefix: Path | None = None,
    use_rich: bool = True,
    tree_depth: int = 3,
    smiles_list: List[str] = None,
) -> Dict[int, float]:
    score_assignment_dict = {}
    n = Z.shape[0] + 1
    nlist = [2 * n - 2]
    depth = 0

    console = Console() if use_rich and verbose else None

    live_context = (
        Live(console=console, refresh_per_second=4) if use_rich and console else None
    )

    level_times = []
    level_start_time = None

    try:
        if live_context:
            live_context.start()

        while len(nlist) != 0:
            if verbose:
                depth += 1

                avg_time = sum(level_times) / len(level_times) if level_times else None
                current_time = (
                    time.time() - level_start_time if level_start_time else None
                )

                if use_rich and console:
                    layout = build_visualization_layout(
                        depth,
                        nlist,
                        score_assignment_dict,
                        n,
                        smiles_list,
                        avg_time_per_level=avg_time,
                        current_level_time=current_time,
                    )

                    if live_context:
                        live_context.update(layout)
                    else:
                        console.print()
                        console.print(layout)
                        console.print()

            level_start_time = time.time()

            nlist = cluster_and_score_per_level(
                Z,
                scoring_func,
                filtering_func,
                nlist,
                score_assignment_dict,
                depth,
                csize_thres,
                B_MAX,
            )

            if level_start_time:
                level_times.append(time.time() - level_start_time)

            if outprefix is not None:
                outfile = outprefix / f"depth-{depth}.scores.csv"
                dfx = pd.DataFrame(
                    score_assignment_dict.items(), columns=["id", "vina_info"]
                )
                dfx.to_csv(outfile, index=False)
    finally:
        if live_context:
            live_context.stop()

    score_assignment_dict = {
        k: v[1] for k, v in score_assignment_dict.items() if k == v[0]
    }

    return score_assignment_dict


def build_visualization_layout(
    depth: int,
    nlist: List[int],
    score_assignment_dict: Dict[int, Tuple[int, float]],
    n: int,
    smiles_list: List[str] = None,
    avg_time_per_level: float = None,
    current_level_time: float = None,
) -> Columns:
    nlistsorted = sorted(
        [(nl, score_assignment_dict.get(nl, (np.nan, np.nan))) for nl in nlist],
        key=lambda x: x[1][1] if isinstance(x[1], tuple) else x[1],
    )

    stats_lines = []
    for nl, (sc0, sc1) in nlistsorted:
        if isinstance(sc0, (int, np.integer)) and not np.isnan(sc1):
            stats_lines.append(f"  Node {nl}: Leaf {sc0} → {float(sc1):.4f}")
        else:
            stats_lines.append(f"  Node {nl}: Not yet scored")

    stats_text = "\n".join(stats_lines) if stats_lines else "  No nodes scored yet"

    total_scored = sum(1 for k, v in score_assignment_dict.items() if k == v[0])

    scored_mols = [(k, v[1]) for k, v in score_assignment_dict.items() if k == v[0]]
    top_mols = sorted(scored_mols, key=lambda x: x[1])[:10]

    leaderboard_lines = []
    for rank, (leaf_id, score) in enumerate(top_mols, 1):
        smiles = (
            smiles_list[leaf_id]
            if smiles_list and leaf_id < len(smiles_list)
            else "N/A"
        )
        if len(smiles) > 45:
            smiles_display = smiles[:45] + "..."
        else:
            smiles_display = smiles

        if score < -10.0:
            color = "bold green"
        elif score < -8.0:
            color = "green"
        elif score < -6.0:
            color = "yellow"
        else:
            color = "white"

        leaderboard_lines.append(
            f"[{color}]{rank:2d}. [{score:7.2f}] Mol {leaf_id:4d} | {smiles_display}[/{color}]"
        )

    leaderboard_text = (
        "\n".join(leaderboard_lines)
        if leaderboard_lines
        else "[dim]No molecules scored yet[/dim]"
    )

    if scored_mols:
        scores = [s for _, s in scored_mols]
        best_score = min(scores)
        avg_score = sum(scores) / len(scores)
        score_stats = (
            f"[bold]Best:[/bold] {best_score:.2f}  [bold]Avg:[/bold] {avg_score:.2f}"
        )
    else:
        score_stats = "[dim]No scores yet[/dim]"

    # Create content with spinner
    spinner_text = "Running molecular docking..."
    if avg_time_per_level is not None and current_level_time is not None:
        spinner_text = f"Running molecular docking... [Current: {current_level_time:.1f}s | Avg: {avg_time_per_level:.1f}s/level]"
    elif avg_time_per_level is not None:
        spinner_text = (
            f"Running molecular docking... [Avg: {avg_time_per_level:.1f}s/level]"
        )

    spinner = Spinner("dots", text=Text(spinner_text, style="dim"))

    progress_text = Text()
    progress_text.append("Depth: ", style="bold cyan")
    progress_text.append(f"{depth}\n")
    progress_text.append("Exploring: ", style="bold cyan")
    progress_text.append(f"{len(nlist)} node(s)\n")
    progress_text.append("Total Scored: ", style="bold cyan")
    progress_text.append(f"{total_scored} / {n} molecules\n")
    progress_text.append_text(Text.from_markup(score_stats))
    progress_text.append("\n\n")
    progress_text.append("Current Level Nodes:\n", style="bold")
    progress_text.append_text(Text.from_markup(stats_text))

    info_content = Group(spinner, progress_text)

    info_panel = Panel(
        info_content,
        title="Search Progress",
        border_style="magenta",
        padding=(1, 2),
    )

    leaderboard_panel = Panel(
        leaderboard_text,
        title="Top Scoring Molecules",
        border_style="green",
        padding=(1, 2),
    )

    return Columns([info_panel, leaderboard_panel], equal=False, expand=True)

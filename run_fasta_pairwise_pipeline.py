from __future__ import annotations

import argparse
import csv
from itertools import combinations
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd
import matplotlib.pyplot as plt


# -------------------------
# Core functions (as used before)
# -------------------------
def read_fasta_simple(fasta_path: str | Path) -> Tuple[list[str], list[str]]:
    """
    Minimal FASTA parser.
    Returns (ids, sequences) where each sequence is concatenated across lines.
    """
    fasta_path = Path(fasta_path)
    ids: list[str] = []
    seqs: list[str] = []

    cur_id: Optional[str] = None
    cur_seq_parts: list[str] = []

    with fasta_path.open("r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    seqs.append("".join(cur_seq_parts))
                    cur_seq_parts = []
                cur_id = line[1:].strip()
                ids.append(cur_id)
            else:
                cur_seq_parts.append(line)

        if cur_id is not None:
            seqs.append("".join(cur_seq_parts))

    if len(ids) != len(seqs):
        raise ValueError(f"FASTA parse error: ids={len(ids)} seqs={len(seqs)} in {fasta_path}")

    return ids, seqs


def simple_distance(seq1: str, seq2: str) -> int:
    """
    Basic distance used in our analysis:
    Count mismatches position-wise up to the shorter length.
    """
    m = min(len(seq1), len(seq2))
    return sum(a != b for a, b in zip(seq1[:m], seq2[:m]))


def fasta_to_pairwise_distance_csv(
    fasta_path: str | Path,
    out_csv_path: str | Path,
    *,
    overwrite: bool = True,
) -> Path:
    """
    FASTA -> CSV with columns: ID1, ID2, Distance.

    Writes streaming rows to CSV so it won't allocate a huge DataFrame in memory.
    Overwrites by default.
    """
    fasta_path = Path(fasta_path)
    out_csv_path = Path(out_csv_path)
    out_csv_path.parent.mkdir(parents=True, exist_ok=True)

    if out_csv_path.exists() and not overwrite:
        return out_csv_path

    ids, seqs = read_fasta_simple(fasta_path)

    with out_csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["ID1", "ID2", "Distance"])

        if len(seqs) < 2:
            return out_csv_path

        for (i, s1), (j, s2) in combinations(enumerate(seqs), 2):
            w.writerow([ids[i], ids[j], simple_distance(s1, s2)])

    return out_csv_path


def pairwise_distance_csv_to_histograms(
    pairwise_csv_path: str | Path,
    out_full_png: str | Path,
    out_zoom_png: str | Path,
    *,
    zoom_max_distance: int = 50,
    threshold: Optional[int] = None,     # default: no threshold line
    tick_step_full: int = 5,
    tick_step_zoom: int = 1,
    show_grid: bool = False,             # default: no grid lines
    font_size: int = 18,
    axis_title_font_size: int = 22,
    main_title_font_size: int = 28,
    fig_width: float = 14.4,             # 20% narrower than 18
    fig_height: float = 6.0,
    overwrite: bool = True,
) -> Tuple[Path, Path]:
    """
    CSV -> two histogram PNGs.
    Overwrites by default.
    """
    pairwise_csv_path = Path(pairwise_csv_path)
    out_full_png = Path(out_full_png)
    out_zoom_png = Path(out_zoom_png)
    out_full_png.parent.mkdir(parents=True, exist_ok=True)
    out_zoom_png.parent.mkdir(parents=True, exist_ok=True)

    if (out_full_png.exists() or out_zoom_png.exists()) and not overwrite:
        return out_full_png, out_zoom_png

    df = pd.read_csv(pairwise_csv_path)
    if "Distance" not in df.columns:
        raise ValueError(f"CSV must include a 'Distance' column: {pairwise_csv_path}")

    distances = df["Distance"].dropna().astype(int).tolist()

    x_label = "Distance"
    y_label = "Count"

    # Empty safety
    if not distances:
        for out_path, title in [
            (out_full_png, "Pairwise Distance Histogram Full Range"),
            (out_zoom_png, f"Pairwise Distance Histogram Focused ≤{zoom_max_distance}"),
        ]:
            plt.figure(figsize=(fig_width, fig_height))
            plt.title(title, fontsize=main_title_font_size)
            plt.xlabel(x_label, fontsize=axis_title_font_size)
            plt.ylabel(y_label, fontsize=axis_title_font_size)
            plt.xticks(fontsize=font_size)
            plt.yticks(fontsize=font_size)
            if show_grid:
                plt.grid(axis="y")
            plt.tight_layout()
            plt.savefig(out_path)
            plt.close()
        return out_full_png, out_zoom_png

    # -------- Full --------
    max_full = max(distances)
    bins_full = range(0, max_full + 2)

    plt.figure(figsize=(fig_width, fig_height))
    plt.hist(distances, bins=bins_full, edgecolor="black")

    if threshold is not None:
        plt.axvline(x=threshold, color="red", linestyle="--", label=f"Threshold = {threshold}")
        plt.legend(fontsize=font_size)

    xticks_full = list(range(0, max_full + 1, max(1, tick_step_full)))
    plt.xticks(xticks_full, rotation=90, fontsize=font_size)
    plt.yticks(fontsize=font_size)

    plt.title("Pairwise Distance Histogram", fontsize=main_title_font_size)
    plt.xlabel(x_label, fontsize=axis_title_font_size)
    plt.ylabel(y_label, fontsize=axis_title_font_size)

    if show_grid:
        plt.grid(axis="y")

    plt.tight_layout()
    plt.savefig(out_full_png)
    plt.close()

    # -------- Zoom --------
    distances_zoom = [d for d in distances if d <= zoom_max_distance]
    bins_zoom = range(0, zoom_max_distance + 2)

    plt.figure(figsize=(fig_width, fig_height))
    plt.hist(distances_zoom, bins=bins_zoom, edgecolor="black")

    if threshold is not None:
        plt.axvline(x=threshold, color="red", linestyle="--", label=f"Threshold = {threshold}")
        plt.legend(fontsize=font_size)

    xticks_zoom = list(range(0, zoom_max_distance + 1, max(1, tick_step_zoom)))
    plt.xticks(xticks_zoom, rotation=90, fontsize=font_size)
    plt.yticks(fontsize=font_size)

    plt.title(f"Pairwise Distance Histogram Focused ≤{zoom_max_distance}", fontsize=main_title_font_size)
    plt.xlabel(x_label, fontsize=axis_title_font_size)
    plt.ylabel(y_label, fontsize=axis_title_font_size)

    if show_grid:
        plt.grid(axis="y")

    plt.tight_layout()
    plt.savefig(out_zoom_png)
    plt.close()

    return out_full_png, out_zoom_png


# -------------------------
# Folder conventions inferred from z.zip
# -------------------------
FASTA_SUBFOLDERS = {
    "light": "fasta/light_fasta",
    "heavy": "fasta/heavy_fasta",
    "combined": "fasta/combined_light+heavy_fasta",
    "cdr3": "fasta/CDR3_fasta",
}

PAIRWISE_SUBFOLDERS = {
    "light": "pairwise_distance/light_pairwise",
    "heavy": "pairwise_distance/heavy_pairwise",
    "combined": "pairwise_distance/combined_light+heavy_pairwise",
    "cdr3": "pairwise_distance/CDR3_pairwise",
}

PLOT_SUBFOLDERS = {
    "light": "plot/light_plots",
    "heavy": "plot/heavy_plots",
    "combined": "plot/combined_light+heavy_plots",
    "cdr3": "plot/CDR3_plots",
}


def infer_dataset_and_type(fasta_file: Path) -> Optional[tuple[str, str, str]]:
    """
    Returns (dataset_key, seq_type, variant) where:
      - dataset_key: 'dataset6'
      - seq_type: 'light' | 'heavy' | 'combined' | 'cdr3'
      - variant: 'IMGT' | 'Regular' | 'CDR3'
    """
    name = fasta_file.name

    m = __import__("re").search(r"(dataset\d+)_", name, flags=__import__("re").IGNORECASE)
    if not m:
        return None
    dataset_key = m.group(1).lower()

    if "heavy_cdr3" in name.lower():
        return dataset_key, "cdr3", "CDR3"

    if "light_chain_sequences" in name.lower():
        seq_type = "light"
    elif "heavy_chain_sequences" in name.lower():
        seq_type = "heavy"
    elif "combined_light_heavy_sequences" in name.lower():
        seq_type = "combined"
    else:
        return None

    variant = "IMGT" if "imgt" in name.lower() else "Regular"
    return dataset_key, seq_type, variant


def build_output_paths(root: Path, dataset_key: str, seq_type: str, variant: str) -> tuple[Path, Path, Path]:
    """
    Returns (csv_path, full_png_path, zoom_png_path) for a given FASTA.
    Naming matches the style in your z.zip.
    """
    # CSV
    pairwise_dir = root / PAIRWISE_SUBFOLDERS[seq_type]
    pairwise_dir.mkdir(parents=True, exist_ok=True)

    if seq_type == "cdr3":
        csv_name = f"{dataset_key}_heavy_cdr3_pairwise_distances.csv"
    else:
        if variant == "IMGT":
            csv_name = f"{dataset_key}_{seq_type}_pairwise_distances.csv"
        else:
            csv_name = f"{dataset_key}_regular_{seq_type}_pairwise_distances.csv"

    csv_path = pairwise_dir / csv_name

    # PNGs
    plot_dir = root / PLOT_SUBFOLDERS[seq_type]
    plot_dir.mkdir(parents=True, exist_ok=True)

    ds_num = dataset_key.replace("dataset", "")
    ds_prefix = f"Dataset{ds_num}"

    if seq_type == "light":
        base = f"{ds_prefix}_Light_Chain"
    elif seq_type == "heavy":
        base = f"{ds_prefix}_Heavy_Chain"
    elif seq_type == "combined":
        base = f"{ds_prefix}_Combined_Chain"
    else:
        base = f"{ds_prefix}_Heavy_CDR3"

    if seq_type != "cdr3":
        base = f"{base}_{variant}"

    full_png = plot_dir / f"{base}_full.png"
    zoom_png = plot_dir / f"{base}_zoom.png"

    return csv_path, full_png, zoom_png


def main():
    ap = argparse.ArgumentParser(
        description="Run FASTA -> pairwise CSV -> histogram PNGs using the folder layout from z.zip"
    )
    ap.add_argument(
        "--root",
        type=str,
        default=".",
        help="Project root that contains fasta/, pairwise_distance/, plot/ (default: current directory).",
    )
    ap.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing CSV/PNG files. (Default behavior is overwrite anyway.)",
    )
    ap.add_argument(
        "--no-overwrite",
        action="store_true",
        help="Do NOT overwrite existing files (skip if outputs exist).",
    )
    ap.add_argument(
        "--threshold",
        type=int,
        default=None,
        help="Optional vertical threshold line (e.g., 13). Default: no line.",
    )
    ap.add_argument("--zoom-max", type=int, default=50, help="Zoom histogram max distance (default: 50).")
    ap.add_argument("--tick-full", type=int, default=5, help="X tick step in full plot (default: 5).")
    ap.add_argument("--tick-zoom", type=int, default=1, help="X tick step in zoom plot (default: 1).")

    args = ap.parse_args()

    root = Path(args.root).resolve()
    overwrite = True
    if args.no_overwrite:
        overwrite = False
    if args.overwrite:
        overwrite = True

    fasta_root = root / "fasta"
    if not fasta_root.exists():
        raise FileNotFoundError(f"Could not find 'fasta' folder under: {root}")

    fasta_files = sorted(fasta_root.rglob("*.fasta"))
    if not fasta_files:
        print(f"No FASTA files found under: {fasta_root}")
        return

    processed = 0
    skipped = 0

    for fasta_file in fasta_files:
        info = infer_dataset_and_type(fasta_file)
        if info is None:
            skipped += 1
            continue

        dataset_key, seq_type, variant = info
        csv_path, full_png, zoom_png = build_output_paths(root, dataset_key, seq_type, variant)

        # If not overwriting, skip when outputs already exist
        if not overwrite and csv_path.exists() and full_png.exists() and zoom_png.exists():
            skipped += 1
            continue

        fasta_to_pairwise_distance_csv(fasta_file, csv_path, overwrite=True if overwrite else False)

        pairwise_distance_csv_to_histograms(
            csv_path,
            full_png,
            zoom_png,
            zoom_max_distance=args.zoom_max,
            threshold=args.threshold,
            tick_step_full=args.tick_full,
            tick_step_zoom=args.tick_zoom,
            show_grid=False,
            overwrite=True if overwrite else False,
        )

        processed += 1

    print(f"Done. Processed: {processed}, Skipped: {skipped}")
    print(f"Root: {root}")


if __name__ == "__main__":
    main()

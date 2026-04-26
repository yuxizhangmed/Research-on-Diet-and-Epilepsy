#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_cortical_surface.py

GitHub-ready script for visualizing cortical brain-region association results
on the fsaverage surface using the Desikan-Killiany atlas.

What this script does
---------------------
1. Reads a CSV file containing cortical association results.
2. Keeps significant regions according to a P-value threshold.
3. Maps region-level beta values onto fsaverage cortical vertices.
4. Generates 1 x 4 surface plots:
   Left lateral, left medial, right medial, right lateral.
5. Saves both the surface figure and a matching-check CSV.

Expected input columns
----------------------
Required:
    Block
    Beta
    P_Value

One of the following region-label columns is required:
    Brain_Region
        Example: "Area of lateraloccipital (left hemisphere)"
                 "Thickness of superiorfrontal (right hemisphere)"

    Atlas_Label
        Example: "ctx-lh-lateraloccipital"
                 "ctx-rh-superiorfrontal"

Example command
---------------
python scripts/plot_cortical_surface.py \
    --input data/brain_region_results.csv \
    --output results/cortical_surface \
    --blocks Cortical_Area Cortical_Thickness \
    --p-threshold 0.05

Notes
-----
- The script will download fsaverage through MNE if it is not already present.
- Do not upload real UK Biobank data to GitHub.
- This script is intended for reproducible figure generation from de-identified
  summary-level association results.
"""

import argparse
import gc
import os
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
import mne
from nilearn import plotting


DEFAULT_BG_COLOR = "#DEDEDE"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot cortical region-level beta values on fsaverage Desikan-Killiany surfaces."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to input CSV containing Block, Beta, P_Value, and Brain_Region or Atlas_Label.",
    )
    parser.add_argument(
        "--output",
        default="results/cortical_surface",
        help="Output directory for figures and matching-check CSV files.",
    )
    parser.add_argument(
        "--subjects-dir",
        default="data/subjects",
        help="Directory used by MNE to store fsaverage. This directory should not be committed to GitHub.",
    )
    parser.add_argument(
        "--blocks",
        nargs="*",
        default=None,
        help="Blocks to plot. If omitted, all unique Block values in the input file are plotted.",
    )
    parser.add_argument(
        "--p-threshold",
        type=float,
        default=0.05,
        help="P-value threshold for selecting regions to display. Default: 0.05.",
    )
    parser.add_argument(
        "--bg-color",
        default=DEFAULT_BG_COLOR,
        help="Figure background color. Default: #DEDEDE.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output image resolution. Default: 300.",
    )
    parser.add_argument(
        "--file-prefix",
        default="",
        help="Optional prefix added to output filenames.",
    )
    return parser.parse_args()


def check_required_columns(df: pd.DataFrame) -> str:
    required = ["Block", "Beta", "P_Value"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Input file is missing required columns: {missing}")

    if "Brain_Region" in df.columns:
        return "Brain_Region"
    if "Atlas_Label" in df.columns:
        return "Atlas_Label"

    raise ValueError(
        "Input file must contain either 'Brain_Region' or 'Atlas_Label' as the region-label column."
    )


def prepare_dataframe(input_file: str) -> Tuple[pd.DataFrame, str]:
    df = pd.read_csv(input_file)
    label_col = check_required_columns(df)

    df["Block"] = df["Block"].astype(str).str.strip()
    df[label_col] = df[label_col].astype(str).str.strip()
    df["Beta"] = pd.to_numeric(df["Beta"], errors="coerce")
    df["P_Value"] = pd.to_numeric(df["P_Value"], errors="coerce")

    df = df.dropna(subset=["Block", label_col, "Beta", "P_Value"]).copy()
    if df.empty:
        raise ValueError("No valid rows remain after removing missing Block, region, Beta, or P_Value.")

    return df, label_col


def fetch_fsaverage(subjects_dir: str) -> Dict[str, str]:
    subjects_dir_path = Path(subjects_dir)
    subjects_dir_path.mkdir(parents=True, exist_ok=True)

    mne.set_config("SUBJECTS_DIR", str(subjects_dir_path), set_env=True)
    print("Preparing fsaverage surface files...")
    mne.datasets.fetch_fsaverage(subjects_dir=str(subjects_dir_path), verbose=False)

    fsaverage_dir = subjects_dir_path / "fsaverage"
    paths = {
        "lh_annot": fsaverage_dir / "label" / "lh.aparc.annot",
        "rh_annot": fsaverage_dir / "label" / "rh.aparc.annot",
        "lh_mesh": fsaverage_dir / "surf" / "lh.pial",
        "rh_mesh": fsaverage_dir / "surf" / "rh.pial",
        "lh_sulc": fsaverage_dir / "surf" / "lh.sulc",
        "rh_sulc": fsaverage_dir / "surf" / "rh.sulc",
    }

    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing fsaverage files: {missing}")

    return {key: str(value) for key, value in paths.items()}


def load_annot(annot_path: str) -> Tuple[np.ndarray, List[str]]:
    labels, _, names = nib.freesurfer.read_annot(annot_path)
    names = [name.decode("utf-8") if isinstance(name, bytes) else str(name) for name in names]
    return labels, names


def parse_atlas_label(label: str) -> Tuple[Optional[str], str]:
    """
    Parse FreeSurfer-style labels:
        ctx-lh-lateraloccipital
        ctx-rh-superiorfrontal
    """
    x = str(label).strip().lower()

    if x.startswith("ctx-lh-"):
        return "lh", x.replace("ctx-lh-", "", 1)
    if x.startswith("ctx-rh-"):
        return "rh", x.replace("ctx-rh-", "", 1)

    return None, x


def parse_brain_region(label: str) -> Tuple[Optional[str], str]:
    """
    Parse natural-language labels:
        Area of lateraloccipital (left hemisphere)
        Thickness of superiorfrontal (right hemisphere)
        Volume of accumbens (left)
    """
    x = str(label).strip().lower()

    if "left" in x:
        hemi = "lh"
    elif "right" in x:
        hemi = "rh"
    else:
        hemi = None

    # Remove common measurement prefixes.
    region = re.sub(r"^(area of|volume of|thickness of|mean fa in|mean md in|mean icvf in|mean isovf in)\s+", "", x)
    # Remove final parenthetical hemisphere information.
    region = re.sub(r"\s*\(.*?\)\s*$", "", region)
    # Match Desikan names, which do not use spaces.
    region = region.replace(" ", "").replace("_", "").replace("-", "")

    return hemi, region


def parse_region(label: str, label_col: str) -> Tuple[Optional[str], str]:
    if label_col == "Atlas_Label":
        return parse_atlas_label(label)
    return parse_brain_region(label)


def pretty_block_name(block_name: str) -> str:
    return str(block_name).replace("_", " ")


def safe_filename(text: str) -> str:
    text = pretty_block_name(text)
    text = re.sub(r"[^\w\-. ]+", "", text)
    text = text.strip().replace(" ", "_")
    return text or "block"


def set_axes_background(ax, color: str) -> None:
    ax.set_facecolor(color)
    try:
        ax.xaxis.set_pane_color((0.87, 0.87, 0.87, 1.0))
        ax.yaxis.set_pane_color((0.87, 0.87, 0.87, 1.0))
        ax.zaxis.set_pane_color((0.87, 0.87, 0.87, 1.0))
    except Exception:
        pass


def build_hemi_map(
    sig_df: pd.DataFrame,
    label_col: str,
    hemi: str,
    vertex_labels: np.ndarray,
    name_to_idx: Dict[str, int],
) -> Tuple[np.ndarray, List[str], List[str]]:
    stat_map = np.zeros(len(vertex_labels), dtype=float)
    matched: List[str] = []
    unmatched: List[str] = []

    for _, row in sig_df.iterrows():
        raw_label = str(row[label_col])
        beta = float(row["Beta"])
        row_hemi, region = parse_region(raw_label, label_col)

        if row_hemi != hemi:
            continue

        if region in name_to_idx:
            idx = name_to_idx[region]
            stat_map[vertex_labels == idx] = beta
            matched.append(raw_label)
        else:
            unmatched.append(raw_label)

    return stat_map, sorted(set(matched)), sorted(set(unmatched))


def choose_colormap(beta_values: Iterable[float]) -> Tuple[str, float, float]:
    beta_values = np.asarray(list(beta_values), dtype=float)
    vmin = float(np.min(beta_values))
    vmax = float(np.max(beta_values))

    if vmin >= 0:
        return "Reds", 0.0, vmax
    if vmax <= 0:
        return "Blues_r", vmin, 0.0

    absmax = max(abs(vmin), abs(vmax))
    return "RdBu_r", -absmax, absmax


def plot_block(
    df: pd.DataFrame,
    label_col: str,
    block_name: str,
    p_threshold: float,
    paths: Dict[str, str],
    output_dir: str,
    bg_color: str,
    dpi: int,
    file_prefix: str = "",
) -> None:
    sub = df[(df["Block"] == block_name) & (df["P_Value"] < p_threshold)].copy()

    print("\n" + "=" * 90)
    print(f"Block: {block_name}")
    print(f"Rows with P < {p_threshold}: {sub.shape[0]}")

    if sub.empty:
        print("No significant cortical regions; skipped.")
        return

    sub = sub.sort_values([label_col, "P_Value"]).drop_duplicates(label_col, keep="first")
    print(f"Regions after removing duplicated labels: {sub.shape[0]}")

    lh_vertex_labels, lh_names = load_annot(paths["lh_annot"])
    rh_vertex_labels, rh_names = load_annot(paths["rh_annot"])
    lh_name_to_idx = {name.lower(): idx for idx, name in enumerate(lh_names)}
    rh_name_to_idx = {name.lower(): idx for idx, name in enumerate(rh_names)}

    lh_map, lh_matched, lh_unmatched = build_hemi_map(
        sub, label_col, "lh", lh_vertex_labels, lh_name_to_idx
    )
    rh_map, rh_matched, rh_unmatched = build_hemi_map(
        sub, label_col, "rh", rh_vertex_labels, rh_name_to_idx
    )

    matched_all = sorted(set(lh_matched + rh_matched))
    unmatched_all = sorted(set(lh_unmatched + rh_unmatched))

    print(f"Matched cortical regions: {len(matched_all)}")
    if unmatched_all:
        print("Unmatched labels. These may be subcortical regions or labels not present in aparc:")
        for item in unmatched_all:
            print(f"  - {item}")

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    match_df = sub[[label_col, "Beta", "P_Value"]].copy()
    match_df["Matched"] = match_df[label_col].isin(matched_all)
    match_file = output_path / f"{file_prefix}{safe_filename(block_name)}_match_check.csv"
    match_df.to_csv(match_file, index=False, encoding="utf-8-sig")
    print(f"Matching-check table saved: {match_file}")

    if len(matched_all) == 0:
        print("No renderable cortical regions for this block; skipped plotting.")
        return

    beta_values = sub.loc[sub[label_col].isin(matched_all), "Beta"].values
    cmap, vmin_plot, vmax_plot = choose_colormap(beta_values)

    fig = plt.figure(figsize=(22, 5.8), facecolor=bg_color)
    axes = [
        fig.add_subplot(1, 4, 1, projection="3d"),
        fig.add_subplot(1, 4, 2, projection="3d"),
        fig.add_subplot(1, 4, 3, projection="3d"),
        fig.add_subplot(1, 4, 4, projection="3d"),
    ]

    for ax in axes:
        set_axes_background(ax, bg_color)

    panels = [
        ("left", "lateral", axes[0], "Left Lateral"),
        ("left", "medial", axes[1], "Left Medial"),
        ("right", "medial", axes[2], "Right Medial"),
        ("right", "lateral", axes[3], "Right Lateral"),
    ]

    for hemi, view, ax, title in panels:
        stat_map = lh_map if hemi == "left" else rh_map
        mesh = paths["lh_mesh"] if hemi == "left" else paths["rh_mesh"]
        sulc = paths["lh_sulc"] if hemi == "left" else paths["rh_sulc"]

        plotting.plot_surf_stat_map(
            surf_mesh=mesh,
            stat_map=stat_map,
            hemi=hemi,
            view=view,
            bg_map=sulc,
            bg_on_data=True,
            darkness=0.6,
            axes=ax,
            cmap=cmap,
            colorbar=False,
            threshold=1e-12,
            vmin=vmin_plot,
            vmax=vmax_plot,
        )

        set_axes_background(ax, bg_color)
        ax.set_title(title, fontsize=18, fontweight="bold", pad=18)

    fig.suptitle(pretty_block_name(block_name), fontsize=24, fontweight="bold", y=1.03)
    plt.subplots_adjust(left=0.03, right=0.90, top=0.82, bottom=0.08, wspace=0.20)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin_plot, vmax=vmax_plot))
    sm.set_array([])

    cb_ax = fig.add_axes([0.915, 0.20, 0.012, 0.52])
    cb_ax.set_facecolor(bg_color)
    cbar = fig.colorbar(sm, cax=cb_ax)
    cbar.set_label("Beta value", fontsize=14)
    cbar.ax.tick_params(labelsize=11)

    figure_file = output_path / f"{file_prefix}{safe_filename(block_name)}_desikan_surface.png"
    plt.savefig(figure_file, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())

    fig.clf()
    plt.close("all")
    gc.collect()

    print(f"Figure saved: {figure_file}")


def main() -> None:
    args = parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    df, label_col = prepare_dataframe(str(input_path))
    paths = fetch_fsaverage(args.subjects_dir)

    blocks = args.blocks if args.blocks else sorted(df["Block"].dropna().unique())
    print(f"Region label column: {label_col}")
    print(f"Blocks to plot: {blocks}")

    for block in blocks:
        plot_block(
            df=df,
            label_col=label_col,
            block_name=block,
            p_threshold=args.p_threshold,
            paths=paths,
            output_dir=args.output,
            bg_color=args.bg_color,
            dpi=args.dpi,
            file_prefix=args.file_prefix,
        )

    print("\nAll cortical surface plots finished.")


if __name__ == "__main__":
    main()

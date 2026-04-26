#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_subcortical_volume.py

GitHub-ready script for visualizing subcortical volume association results
on MNI152 slices using the Harvard-Oxford subcortical atlas.

What this script does
---------------------
1. Reads a CSV file containing subcortical structure-level association results.
2. Keeps significant structures according to a P-value threshold.
3. Maps beta values onto the Harvard-Oxford subcortical atlas.
4. Draws a 1 x 4 slice figure:
   Left sagittal, coronal, right sagittal, axial.
5. Saves the figure and a matching-check CSV.

Expected input columns
----------------------
Required:
    Atlas_Label
    Beta
    P_Value

Optional:
    qFDR or FDR

Example Atlas_Label values
--------------------------
Left-Hippocampus
Right-Hippocampus
Left-Putamen
Right-Putamen
Left-Thalamus
Right-Thalamus
Left-Amygdala
Right-Amygdala
Left-Caudate
Right-Caudate
Left-Pallidum
Right-Pallidum
Left-Accumbens-area
Right-Accumbens-area

Example command
---------------
python scripts/plot_subcortical_volume.py \
    --input data/subcortical_volume_results.csv \
    --output results/subcortical_volume.png \
    --p-threshold 0.05

Notes
-----
- The script downloads/uses Nilearn datasets as needed.
- Do not upload real UK Biobank data to GitHub.
- This script is intended for reproducible figure generation from de-identified
  summary-level association results.
"""

import argparse
import re
from pathlib import Path
from typing import Dict, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from nilearn import datasets, image, plotting


DEFAULT_BG_COLOR = "#DEDEDE"


HARVARD_OXFORD_NAME_MAP: Dict[str, str] = {
    "Left-Hippocampus": "Left Hippocampus",
    "Right-Hippocampus": "Right Hippocampus",
    "Left-Putamen": "Left Putamen",
    "Right-Putamen": "Right Putamen",
    "Left-Thalamus": "Left Thalamus",
    "Right-Thalamus": "Right Thalamus",
    "Left-Amygdala": "Left Amygdala",
    "Right-Amygdala": "Right Amygdala",
    "Left-Caudate": "Left Caudate",
    "Right-Caudate": "Right Caudate",
    "Left-Pallidum": "Left Pallidum",
    "Right-Pallidum": "Right Pallidum",
    "Left-Accumbens-area": "Left Accumbens",
    "Right-Accumbens-area": "Right Accumbens",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot subcortical volume beta values on MNI152 slices."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to input CSV containing Atlas_Label, Beta, and P_Value.",
    )
    parser.add_argument(
        "--output",
        default="results/subcortical_volume.png",
        help="Output figure path. Default: results/subcortical_volume.png.",
    )
    parser.add_argument(
        "--p-threshold",
        type=float,
        default=0.05,
        help="P-value threshold for selecting subcortical structures. Default: 0.05.",
    )
    parser.add_argument(
        "--bg-color",
        default=DEFAULT_BG_COLOR,
        help="Figure background color. Default: #DEDEDE.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=600,
        help="Output image resolution. Default: 600.",
    )
    parser.add_argument(
        "--smooth-fwhm",
        type=float,
        default=1.6,
        help="FWHM used to lightly smooth the effect image. Default: 1.6.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.0015,
        help="Display threshold used to suppress minor smoothing artifacts. Default: 0.0015.",
    )
    parser.add_argument(
        "--white-bg",
        action="store_true",
        help="Use a white background instead of the default gray background.",
    )
    return parser.parse_args()


def normalize_label(label: str) -> str:
    """
    Normalize common subcortical labels to the FreeSurfer-like format used here.
    """
    x = str(label).strip()

    replacements = {
        "Left Hippocampus": "Left-Hippocampus",
        "Right Hippocampus": "Right-Hippocampus",
        "Left Putamen": "Left-Putamen",
        "Right Putamen": "Right-Putamen",
        "Left Thalamus": "Left-Thalamus",
        "Right Thalamus": "Right-Thalamus",
        "Left Amygdala": "Left-Amygdala",
        "Right Amygdala": "Right-Amygdala",
        "Left Caudate": "Left-Caudate",
        "Right Caudate": "Right-Caudate",
        "Left Pallidum": "Left-Pallidum",
        "Right Pallidum": "Right-Pallidum",
        "Left Accumbens": "Left-Accumbens-area",
        "Right Accumbens": "Right-Accumbens-area",
    }

    if x in replacements:
        return replacements[x]

    # Convert strings such as "Volume of hippocampus (left)".
    x_lower = x.lower()
    if "left" in x_lower:
        side = "Left"
    elif "right" in x_lower:
        side = "Right"
    else:
        return x

    structure = re.sub(r"^(volume of|area of|thickness of)\s+", "", x_lower)
    structure = re.sub(r"\s*\(.*?\)\s*$", "", structure).strip()
    structure = structure.replace("_", "-").replace(" ", "-")

    structure_map = {
        "hippocampus": "Hippocampus",
        "putamen": "Putamen",
        "thalamus": "Thalamus",
        "amygdala": "Amygdala",
        "caudate": "Caudate",
        "pallidum": "Pallidum",
        "accumbens": "Accumbens-area",
        "accumbens-area": "Accumbens-area",
    }

    return f"{side}-{structure_map.get(structure, structure)}"


def prepare_dataframe(input_file: str, p_threshold: float) -> pd.DataFrame:
    df = pd.read_csv(input_file)

    required = ["Atlas_Label", "Beta", "P_Value"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Input file is missing required columns: {missing}")

    df["Atlas_Label"] = df["Atlas_Label"].astype(str).str.strip()
    df["Normalized_Label"] = df["Atlas_Label"].apply(normalize_label)
    df["Beta"] = pd.to_numeric(df["Beta"], errors="coerce")
    df["P_Value"] = pd.to_numeric(df["P_Value"], errors="coerce")

    df = df.dropna(subset=["Atlas_Label", "Normalized_Label", "Beta", "P_Value"]).copy()
    if df.empty:
        raise ValueError("No valid rows remain after removing missing Atlas_Label, Beta, or P_Value.")

    df_sig = df[df["P_Value"] < p_threshold].copy()
    if df_sig.empty:
        raise ValueError(f"No subcortical structures pass P < {p_threshold}.")

    # Keep the smallest P-value if a structure is duplicated.
    df_sig = df_sig.sort_values(["Normalized_Label", "P_Value"]).drop_duplicates(
        "Normalized_Label", keep="first"
    )

    return df_sig


def choose_colormap(beta_values: np.ndarray) -> Tuple[str, float, float, bool]:
    vmin = float(np.min(beta_values))
    vmax = float(np.max(beta_values))

    if vmin >= 0:
        return "Reds", 0.0, vmax, False
    if vmax <= 0:
        return "Blues_r", vmin, 0.0, False

    absmax = max(abs(vmin), abs(vmax))
    return "RdBu_r", -absmax, absmax, True


def build_effect_image(df_sig: pd.DataFrame):
    atlas = datasets.fetch_atlas_harvard_oxford("sub-maxprob-thr25-2mm")
    atlas_img = image.load_img(atlas.maps)
    atlas_data = atlas_img.get_fdata()
    labels = list(atlas.labels)
    label_to_index = {label: idx for idx, label in enumerate(labels)}

    effect_data = np.zeros(atlas_data.shape, dtype=float)
    matched = []
    unmatched = []

    for _, row in df_sig.iterrows():
        source_label = row["Normalized_Label"]
        beta = float(row["Beta"])
        ho_label = HARVARD_OXFORD_NAME_MAP.get(source_label)

        if ho_label is None or ho_label not in label_to_index:
            unmatched.append(row["Atlas_Label"])
            continue

        idx = label_to_index[ho_label]
        effect_data[atlas_data == idx] = beta
        matched.append(row["Atlas_Label"])

    effect_img = image.new_img_like(atlas_img, effect_data)
    return effect_img, matched, unmatched


def save_match_check(df_sig: pd.DataFrame, matched: list, output_file: str) -> None:
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    match_df = df_sig[["Atlas_Label", "Normalized_Label", "Beta", "P_Value"]].copy()
    match_df["Matched"] = match_df["Atlas_Label"].isin(matched)

    if "qFDR" in df_sig.columns:
        match_df["qFDR"] = df_sig["qFDR"]
    elif "FDR" in df_sig.columns:
        match_df["FDR"] = df_sig["FDR"]

    match_file = output_path.with_suffix(".match_check.csv")
    match_df.to_csv(match_file, index=False, encoding="utf-8-sig")
    print(f"Matching-check table saved: {match_file}")


def main() -> None:
    args = parse_args()
    bg_color = "white" if args.white_bg else args.bg_color

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df_sig = prepare_dataframe(str(input_path), args.p_threshold)
    print(f"Significant subcortical structures: {len(df_sig)}")

    effect_img, matched, unmatched = build_effect_image(df_sig)
    print(f"Matched structures: {len(matched)}")
    if unmatched:
        print("Unmatched structures:")
        for item in unmatched:
            print(f"  - {item}")

    save_match_check(df_sig, matched, str(output_path))

    if len(matched) == 0:
        raise ValueError("No structures could be matched to the Harvard-Oxford subcortical atlas.")

    bg_img = datasets.load_mni152_template(resolution=1)

    effect_1mm = image.resample_to_img(
        effect_img,
        bg_img,
        interpolation="continuous",
        force_resample=True,
        copy_header=True,
    )
    effect_smooth = image.smooth_img(effect_1mm, fwhm=args.smooth_fwhm)
    bg_smooth = image.smooth_img(bg_img, fwhm=0.6)

    beta_values = df_sig.loc[df_sig["Atlas_Label"].isin(matched), "Beta"].values
    cmap, vmin_plot, vmax_plot, symmetric_cbar = choose_colormap(beta_values)

    fig = plt.figure(figsize=(22, 5.8), facecolor=bg_color)

    panel_pos = [
        [0.03, 0.12, 0.20, 0.62],
        [0.27, 0.12, 0.20, 0.62],
        [0.51, 0.12, 0.20, 0.62],
        [0.75, 0.12, 0.20, 0.62],
    ]

    axes = [fig.add_axes(pos) for pos in panel_pos]
    for ax in axes:
        ax.set_facecolor(bg_color)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    panels = [
        ("x", -18, axes[0], "Left Sagittal"),
        ("y", -8, axes[1], "Coronal"),
        ("x", 18, axes[2], "Right Sagittal"),
        ("z", 2, axes[3], "Axial"),
    ]

    for display_mode, cut_coord, ax, title in panels:
        plotting.plot_stat_map(
            effect_smooth,
            bg_img=bg_smooth,
            display_mode=display_mode,
            cut_coords=[cut_coord],
            annotate=False,
            draw_cross=False,
            black_bg=False,
            cmap=cmap,
            symmetric_cbar=symmetric_cbar,
            colorbar=False,
            threshold=args.threshold,
            vmin=vmin_plot,
            vmax=vmax_plot,
            axes=ax,
            dim=0.15,
        )
        ax.set_facecolor(bg_color)
        ax.set_title(title, fontsize=18, fontweight="bold", pad=22)

    fig.suptitle("Subcortical Volume", fontsize=24, fontweight="bold", y=0.98)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin_plot, vmax=vmax_plot))
    sm.set_array([])

    cb_ax = fig.add_axes([0.94, 0.20, 0.012, 0.52])
    cb_ax.set_facecolor(bg_color)
    cbar = fig.colorbar(sm, cax=cb_ax)
    cbar.set_label("Beta value", fontsize=14)
    cbar.ax.tick_params(labelsize=12)

    plt.savefig(
        output_path,
        dpi=args.dpi,
        bbox_inches="tight",
        facecolor=fig.get_facecolor(),
    )
    plt.close(fig)

    print(f"Figure saved: {output_path}")


if __name__ == "__main__":
    main()

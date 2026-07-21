import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from bioio import BioImage
from scipy import ndimage
from scipy.ndimage import label, binary_fill_holes
from skimage import filters, morphology
from skimage.color import label2rgb

# Configuration
DATA_DIR = "./data"
OUTPUT_DIR = "./output"
DAPI_CHANNEL = 0
HA_CHANNEL = 1
CPSF6_CHANNEL = 2
CAPSID_CHANNEL = 3
NUCLEI_DIAMETER_PX = 140
SIZE_TOLERANCE = 0.3
LABEL_IMAGE_DIR = "./output/label_images"
CHANNELS = [
    ("DAPI", DAPI_CHANNEL),
    ("HA", HA_CHANNEL),
    ("CPSF6", CPSF6_CHANNEL),
    ("Capsid", CAPSID_CHANNEL),
]

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Define experimental conditions based on file indices
CONDITION_MAPPING = {
    7: "D37_RR-VLPs", 10: "D37_RR-VLPs", 11: "D37_RR-VLPs", 12: "D37_RR-VLPs", 13: "D37_RR-VLPs",
    14: "D102-VLPs", 15: "D102-VLPs", 16: "D102-VLPs", 22: "D102-VLPs",
    17: "Uninfected", 18: "Uninfected", 19: "Uninfected", 20: "Uninfected", 21: "Uninfected"
}


def get_condition_from_filename(filename):
    """Extract image index from filename and return condition."""
    # Filenames are formatted like "10_Multichannel Z-Stack_20260622_67.vsi",
    # where the leading number is the file index used in CONDITION_MAPPING.
    numbers = re.findall(r"\d+", Path(filename).stem)
    if not numbers:
        return "Unknown"
    return CONDITION_MAPPING.get(int(numbers[0]), "Unknown")


def segment_nuclei_3d(dapi_stack):
    """
    Segment nuclei from 3D DAPI z-stack.
    
    Parameters:
    dapi_stack: 3D numpy array (z, y, x)
    
    Returns:
    labeled_3d: 3D labeled image with nucleus IDs
    """

    # Apply Gaussian smoothing to reduce noise
    dapi_smooth = ndimage.gaussian_filter(dapi_stack, sigma=[1.0, 2.0, 2.0])

    # Create initial binary mask using the triangle threshold
    threshold = filters.threshold_triangle(dapi_smooth)
    binary_mask = dapi_smooth > threshold

    # Fill holes in the binary mask
    binary_mask = binary_fill_holes(binary_mask)

    # Apply morphological operations to clean up
    binary_mask = morphology.binary_erosion(binary_mask, morphology.ball(2))
    binary_mask = morphology.binary_dilation(binary_mask, morphology.ball(2))

    # Label connected components in 3D
    labeled_3d, num_features = label(binary_mask)

    # Filter by size: keep nuclei within acceptable size range
    min_size = int(np.pi * (NUCLEI_DIAMETER_PX * (1 - SIZE_TOLERANCE) / 2) ** 2 / 10)
    max_size = int(np.pi * (NUCLEI_DIAMETER_PX * (1 + SIZE_TOLERANCE) / 2) ** 2 * 10)

    filtered_labeled = np.zeros_like(labeled_3d)
    new_label = 0

    for nucleus_id in range(1, num_features + 1):
        nucleus_voxels = np.sum(labeled_3d == nucleus_id)
        if min_size <= nucleus_voxels <= max_size:
            new_label += 1
            filtered_labeled[labeled_3d == nucleus_id] = new_label

    return filtered_labeled


def save_label_images(dapi_stack, labeled_nuclei, output_subdir):
    """
    Save one PNG per z-slice showing the DAPI signal with segmented nuclei overlaid.

    Parameters:
    dapi_stack: 3D numpy array (z, y, x) of raw DAPI intensities
    labeled_nuclei: 3D labeled image with nucleus IDs, same shape as dapi_stack
    output_subdir: directory to save the per-slice PNGs into (created if needed)
    """
    os.makedirs(output_subdir, exist_ok=True)

    dapi_min = np.percentile(dapi_stack, 1)
    dapi_max = np.percentile(dapi_stack, 99)
    dapi_normalized = np.clip((dapi_stack - dapi_min) / (dapi_max - dapi_min), 0, 1)

    for z in range(dapi_stack.shape[0]):
        overlay = label2rgb(labeled_nuclei[z], image=dapi_normalized[z], bg_label=0, alpha=0.4)
        plt.imsave(os.path.join(output_subdir, f"z{z:03d}.png"), overlay)


def extract_intensity_metrics(image_data, labeled_nuclei, nucleus_id):
    """
    Extract intensity metrics for a single nucleus across all channels.
    
    Parameters:
    image_data: 4D array (channels, z, y, x)
    labeled_nuclei: 3D labeled image
    nucleus_id: ID of nucleus to measure
    
    Returns:
    Dictionary with intensity metrics for all channels
    """
    # Get voxels belonging to this nucleus
    nucleus_mask = labeled_nuclei == nucleus_id

    metrics = {"nucleus_id": nucleus_id}

    # Extract metrics for each channel. nucleus_mask is guaranteed non-empty:
    # segment_nuclei_3d only assigns labels to regions that passed its size filter.
    for ch_name, ch_idx in CHANNELS:
        if ch_idx >= image_data.shape[0]:
            continue
        nucleus_intensities = image_data[ch_idx][nucleus_mask]
        metrics[f"{ch_name}_mean"] = np.mean(nucleus_intensities)
        metrics[f"{ch_name}_median"] = np.median(nucleus_intensities)
        metrics[f"{ch_name}_min"] = np.min(nucleus_intensities)
        metrics[f"{ch_name}_max"] = np.max(nucleus_intensities)
        metrics[f"{ch_name}_std"] = np.std(nucleus_intensities)
        metrics[f"{ch_name}_total"] = np.sum(nucleus_intensities)

    return metrics


def process_vsi_file(filepath):
    """
    Process a single VSI file: segment nuclei and extract intensity metrics.
    
    Parameters:
    filepath: Path to VSI file
    
    Returns:
    pandas DataFrame with per-nucleus measurements
    """
    print(f"Processing: {filepath}")

    try:
        # Read image using bioio, collapsing the timepoint axis to get (C, Z, Y, X)
        bio_image = BioImage(filepath)
        image_data = bio_image.get_image_data("CZYX", T=0)
        dapi_stack = image_data[DAPI_CHANNEL]

        # Segment nuclei in 3D
        labeled_nuclei = segment_nuclei_3d(dapi_stack)
        num_nuclei = np.max(labeled_nuclei)

        print(f"  Found {num_nuclei} nuclei")

        # Save per-slice label images for visual inspection
        filename = Path(filepath).name
        label_image_dir = os.path.join(LABEL_IMAGE_DIR, Path(filepath).stem)
        save_label_images(dapi_stack, labeled_nuclei, label_image_dir)

        # Extract metrics for each nucleus
        measurements = []
        condition = get_condition_from_filename(filename)

        for nucleus_id in range(1, int(num_nuclei) + 1):
            metrics = extract_intensity_metrics(image_data, labeled_nuclei, nucleus_id)
            metrics["filename"] = filename
            metrics["condition"] = condition
            measurements.append(metrics)

        return pd.DataFrame(measurements)

    except Exception as e:
        print(f"  Error processing {filepath}: {e}")
        return None


def main():
    """Main analysis pipeline."""

    # Find all VSI files in data directory
    vsi_files = list(Path(DATA_DIR).glob("*.vsi"))

    if not vsi_files:
        print(f"No VSI files found in {DATA_DIR}")
        return

    print(f"Found {len(vsi_files)} VSI files")

    # Process each file and collect results
    all_measurements = []

    for vsi_file in sorted(vsi_files):
        df = process_vsi_file(str(vsi_file))
        if df is not None and len(df) > 0:
            all_measurements.append(df)

    # Combine all measurements
    if all_measurements:
        results_df = pd.concat(all_measurements, ignore_index=True)

        # Save full results
        output_file = os.path.join(OUTPUT_DIR, "nuclei_measurements.csv")
        results_df.to_csv(output_file, index=False)
        print(f"\nSaved full measurements to: {output_file}")
        print(f"Total nuclei measured: {len(results_df)}")

        # Generate summary statistics by condition
        summary_stats = []

        for condition in results_df["condition"].unique():
            condition_data = results_df[results_df["condition"] == condition]

            summary = {
                "condition": condition,
                "num_images": condition_data["filename"].nunique(),
                "num_nuclei": len(condition_data),
            }

            # Add mean values for each channel metric
            for col in condition_data.columns:
                if col not in ["nucleus_id", "filename", "condition"]:
                    summary[f"{col}_mean"] = condition_data[col].mean()
                    summary[f"{col}_std"] = condition_data[col].std()

            summary_stats.append(summary)

        summary_df = pd.DataFrame(summary_stats)
        summary_file = os.path.join(OUTPUT_DIR, "summary_statistics.csv")
        summary_df.to_csv(summary_file, index=False)
        print(f"Saved summary statistics to: {summary_file}")

        # Display summary
        print("\nSummary by condition:")
        print(summary_df[["condition", "num_images", "num_nuclei"]])

        # Create a simple visualization: per-nucleus mean intensity by condition,
        # shown as a boxplot with individual nuclei overlaid as a swarm plot
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle("Nuclei Intensity Analysis Summary", fontsize=16)

        condition_order = sorted(results_df["condition"].unique())

        for idx, (channel, _) in enumerate(CHANNELS):
            ax = axes[idx // 2, idx % 2]
            col = f"{channel}_mean"

            sns.boxplot(
                data=results_df, x="condition", y=col, order=condition_order,
                ax=ax, showfliers=False, color="lightgray",
            )
            sns.swarmplot(
                data=results_df, x="condition", y=col, order=condition_order,
                ax=ax, size=2, color="black", alpha=0.6,
            )
            ax.set_title(channel)
            ax.set_xlabel("")
            ax.set_ylabel("Mean intensity per nucleus")
            ax.tick_params(axis="x", rotation=30)

        fig.tight_layout()
        plot_file = os.path.join(OUTPUT_DIR, "intensity_summary.png")
        fig.savefig(plot_file, dpi=150)
        plt.close(fig)
        print(f"Saved summary plot to: {plot_file}")


if __name__ == "__main__":
    main()

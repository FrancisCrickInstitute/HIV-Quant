import os
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import ndimage
from scipy.ndimage import label, binary_fill_holes
from bioio import BioImage
import matplotlib.pyplot as plt
from skimage import filters, morphology
from skimage.segmentation import watershed

# Configuration
DATA_DIR = "./data"
OUTPUT_DIR = "./output"
DAPI_CHANNEL = 0
HA_CHANNEL = 1
CPSF6_CHANNEL = 2
CAPSID_CHANNEL = 3
NUCLEI_DIAMETER_PX = 140
SIZE_TOLERANCE = 0.3

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
    try:
        # Extract number from filename (assuming format like "image_7.vsi" or similar)
        base = Path(filename).stem
        # Try to find a number in the filename
        import re
        numbers = re.findall(r'\d+', base)
        if numbers:
            idx = int(numbers[-1])
            return CONDITION_MAPPING.get(idx, "Unknown")
    except:
        pass
    return "Unknown"

def segment_nuclei_3d(dapi_stack):
    """
    Segment nuclei from 3D DAPI z-stack.
    
    Parameters:
    dapi_stack: 3D numpy array (z, y, x)
    
    Returns:
    labeled_3d: 3D labeled image with nucleus IDs
    """
    # Normalize intensity to 0-1 range
    dapi_min = np.percentile(dapi_stack, 1)
    dapi_max = np.percentile(dapi_stack, 99)
    dapi_normalized = np.clip((dapi_stack - dapi_min) / (dapi_max - dapi_min), 0, 1)
    
    # Apply Gaussian smoothing to reduce noise
    dapi_smooth = ndimage.gaussian_filter(dapi_normalized, sigma=1.5)
    
    # Create initial binary mask using Otsu's threshold
    threshold = filters.threshold_otsu(dapi_smooth)
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
    
    channel_names = ["DAPI", "HA", "CPSF6", "Capsid"]
    
    # Extract metrics for each channel
    for ch_idx, ch_name in enumerate(channel_names):
        if ch_idx < image_data.shape[0]:
            channel_data = image_data[ch_idx]
            nucleus_intensities = channel_data[nucleus_mask]
            
            if len(nucleus_intensities) > 0:
                metrics[f"{ch_name}_mean"] = np.mean(nucleus_intensities)
                metrics[f"{ch_name}_median"] = np.median(nucleus_intensities)
                metrics[f"{ch_name}_min"] = np.min(nucleus_intensities)
                metrics[f"{ch_name}_max"] = np.max(nucleus_intensities)
                metrics[f"{ch_name}_std"] = np.std(nucleus_intensities)
                metrics[f"{ch_name}_total"] = np.sum(nucleus_intensities)
            else:
                metrics[f"{ch_name}_mean"] = 0
                metrics[f"{ch_name}_median"] = 0
                metrics[f"{ch_name}_min"] = 0
                metrics[f"{ch_name}_max"] = 0
                metrics[f"{ch_name}_std"] = 0
                metrics[f"{ch_name}_total"] = 0
    
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
        # Read image using bioio
        bio_image = BioImage(filepath)
        image_data = bio_image.data
        
        # Handle different possible data shapes
        # Expected: (channels, z, y, x) or similar
        if image_data.ndim == 4:
            # Assuming format is (C, Z, Y, X)
            dapi_stack = image_data[DAPI_CHANNEL]
        elif image_data.ndim == 3:
            # If no channel dimension, assume single channel
            dapi_stack = image_data
        else:
            print(f"Unexpected image shape: {image_data.shape}")
            return None
        
        # Ensure we have proper 3D data
        if dapi_stack.ndim == 3:
            # 3D z-stack
            z_projection = np.max(dapi_stack, axis=0)
        else:
            z_projection = dapi_stack
        
        # Segment nuclei in 3D
        labeled_nuclei = segment_nuclei_3d(dapi_stack)
        num_nuclei = np.max(labeled_nuclei)
        
        print(f"  Found {num_nuclei} nuclei")
        
        # Extract metrics for each nucleus
        measurements = []
        filename = Path(filepath).name
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
        
        # Create a simple visualization
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle("Nuclei Intensity Analysis Summary", fontsize=16)
        
        # Plot mean intensities by condition for each channel
        channels = ["DAPI", "HA", "CPSF6", "Capsid"]
        conditions = results_df["condition"].unique()
        
        for idx, channel in enumerate(channels):
            ax = axes[idx // 2, idx % 2]
            
            means = []
            stds = []
            labels = []
            
            for condition in sorted(conditions):
                cond_data = results_df[results_df["condition"] == condition]
                mean_val = cond_data[f"{channel}_mean"].mean()
                std_val = cond_data[f"{channel}_mean"].std()
                means.append(mean_val)
                stds.append(std_val)
                labels.append(condition)
            
            ax.bar(labels, means, yerr=stds, capsize=5, alpha=0.7)
            ax.set_title(channel)
            ax.set_ylabel("Mean intensity")
            ax.tick_params(axis="x", rotation=30)

        fig.tight_layout()
        plot_file = os.path.join(OUTPUT_DIR, "intensity_summary.png")
        fig.savefig(plot_file, dpi=150)
        plt.close(fig)
        print(f"Saved summary plot to: {plot_file}")


if __name__ == "__main__":
    main()
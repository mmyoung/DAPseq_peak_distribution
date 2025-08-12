import pandas as pd
import numpy as np
from pybedtools import BedTool
import matplotlib.pyplot as plt

# 1. Input Files ------------------------------------------------------------

# Replace these with your actual file paths
peak_file = "/clusterfs/jgi/groups/gentech/homes/romalley/full_DAPseq_annotation/raw_data/N4_filtered-annotated-peaks_minfoldch5_minus-2000bpTSS-to-plus-2000bpTTS_111624.tsv"  # Input peak file with 'tf' column
#peak_file = "/clusterfs/jgi/groups/gentech/homes/romalley/full_DAPseq_annotation/test_peak.bed"  # Input peak file with 'tf' column
gene_file = "/clusterfs/jgi/groups/gentech/seqtech/plant_multidap_data/genomes/annotations/Arabidopsis_thaliana_Col-0_cds_primary.gff"  # File with gene coordinates

# Load the gene coordinates and peaks as Pandas DataFrames
genes = pd.read_csv(gene_file, sep="\t", header=None, names=["chrom", "source", "region", "start", "end", "score", "strand", "idx", "gene_name"])
peaks = pd.read_csv(peak_file, sep="\t", header=0)
peaks = peaks[peaks['species'] == "Arabidopsis_thaliana_Col-0"]

peaks = peaks[peaks['n_cons_species_minfrac0'] == 4]

# 2. Generate Common Bins for All Genes ------------------------------------

# Function to generate common bins (as defined earlier)
# def generate_bins(start, end, num_bins, bin_prefix):
#     bin_size = (end - start) // num_bins
#     return [(start + i * bin_size, start + (i + 1) * bin_size, f"{bin_prefix}_{i+1}") for i in range(num_bins)]

# # List to store all bins
# all_bins = []

# # Loop through each gene and align to common bins
# for _, row in genes.iterrows():
#     chrom = row["chrom"]
#     gene_start, gene_end = row["start"], row["end"]
#     strand = row["strand"]

#     # Define regions for upstream, genebody, and downstream
#     if strand == "+":
#         if gene_start > 2000:
#             upstream_start, upstream_end = gene_start - 2000, gene_start
#         else:
#             upstream_start, upstream_end = 0, gene_start
#         downstream_start, downstream_end = gene_end, gene_end + 2000
#     else:  # Reverse for "-" strand
#         if gene_start > 2000:
#             downstream_start, downstream_end = gene_start - 2000, gene_start
#         else:
#             downstream_start, downstream_end = 0, gene_start
#         upstream_start, upstream_end = gene_end, gene_end + 2000
        
#         gene_start, gene_end = gene_end, gene_start  # Flip start/end for consistency

#     # Generate bins for upstream, genebody, and downstream
#     upstream_bins = generate_bins(upstream_start, upstream_end, 20, "up_bin")
#     genebody_bins = generate_bins(gene_start, gene_end, 10, "gb_bin")
#     downstream_bins = generate_bins(downstream_start, downstream_end, 20, "down_bin")

#     # Add chromosome and bin coordinates to the list
#     # Ensure start ≤ end for all bins before appending
#     for start, end, bin_name in upstream_bins + genebody_bins + downstream_bins:
#         if start > end:  # Swap if necessary
#             start, end = end, start
#         all_bins.append((chrom, start, end, bin_name))


# # Create a DataFrame of all common bins
# bins_df = pd.DataFrame(all_bins, columns=["chrom", "start", "end", "bin_name"])

# Save common bins to a BED file (only need to generate this once)
common_bins_file = "common_bins.bed"
#bins_df.to_csv(common_bins_file, sep="\t", index=False, header=False)

# 3. Split the Peak File by 'tf' Column -------------------------------------

# Split the peaks DataFrame based on unique 'tf' values
unique_tfs = peaks['tf'].unique()

# Create a dictionary to store results for each TF
all_results = []

# 4. Process Each Subset -----------------------------------------------------

# Loop through each TF and perform the peak count
for tf in unique_tfs:
    # Subset the peaks DataFrame for the current TF
    tf_peaks = peaks[peaks['tf'] == tf]
    tf_peak_tmp = pd.DataFrame()
    tf_peak_tmp['chr'] = tf_peaks['peak_chr']
    tf_peak_tmp['start'] = tf_peaks['peak_start'] + tf_peaks['peak_summit'] - 10
    tf_peak_tmp['end'] = tf_peaks['peak_start'] + tf_peaks['peak_summit'] + 10

    # Convert TF-specific peaks to a BEDTool object
    tf_peak_tmp.to_csv("tf_peak_tmp.bed",index=False, sep="\t", header=False)
    tf_peaks_bed = BedTool("tf_peak_tmp.bed")

    # Intersect to count peaks in each bin (using pre-generated bins)
    common_bins_bed = BedTool(common_bins_file)
    overlap = common_bins_bed.intersect(tf_peaks_bed, c=True)

    # Convert the result to a DataFrame
    overlap_df = pd.read_csv(overlap.fn, sep="\t", header=None, names=["chrom", "start", "end", "bin_name", "peak_count"])

    # Summarize peak counts by bin name
    summary_df = overlap_df.groupby("bin_name")["peak_count"].sum().reset_index()

    # Sort bins numerically
    summary_df["bin_number"] = summary_df["bin_name"].str.extract(r'(\d+)').astype(int)
    summary_df = summary_df.sort_values("bin_number")

    # Add the TF label for the current subset
    summary_df["tf"] = tf

    # Calculate the total number of peaks for this TF
    total_peaks = tf_peaks.shape[0]

    # Add a column for peak frequency (normalized by total peaks for this TF)
    summary_df["peak_frequency"] = summary_df["peak_count"] / total_peaks

    # Append the result for this TF to the list of results
    all_results.append(summary_df)

# 5. Combine Results from All TFs ------------------------------------------

# Concatenate the results for all TFs into one DataFrame
final_df = pd.concat(all_results, ignore_index=True)

# Save the final concatenated results
#output_file = "peak_distribution_across_bins_all_tfs.csv"
output_file = "C4peak_distribution_across_bins_all_tfs.csv"
final_df[["tf", "bin_name", "peak_count","peak_frequency"]].to_csv(output_file, index=False)

# Print a preview of the final results
print("Peak distribution across bins for all TFs:")
print(final_df.head())
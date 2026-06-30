#Import the libraries

from pathlib import Path
import subprocess
import shutil
import os
import pandas as pd
import re
#Define the paths

metadata_file='../data/input_files/manifest_curated_nugent.tsv'
folder_results='../results'
os.makedirs(folder_results, exist_ok=True)

def create_metadata_tabulate(metadata_file, folder_results):
    cmd = ['qiime', 'metadata', 'tabulate',
           '--m-input-file', f'{metadata_file}',
           '--o-visualization', f'{folder_results}/sample-metadata-viz.qzv']
    subprocess.run(cmd, check=True)

#create_metadata_tabulate(metadata_file, folder_results)


def import_data(metadata_file, folder_results):
    cmd = ['qiime', 'tools', 'import',
           '--type', 'SampleData[PairedEndSequencesWithQuality]',
           '--input-path', f'{metadata_file}',
           '--output-path', f'{folder_results}/demux-paired-end.qza',
           '--input-format', 'PairedEndFastqManifestPhred33V2']
    subprocess.run(cmd, check=True)

#import_data(metadata_file, folder_results)


def demux_summarize(folder_results):
    cmd = ['qiime', 'demux', 'summarize',
           '--i-data', f'{folder_results}/demux-paired-end.qza',
           '--o-visualization', f'{folder_results}/demux-paired-end.qzv']
    subprocess.run(cmd, check=True)

#demux_summarize(folder_results)

def denoise_paired(folder_results):
       os.makedirs(f'{folder_results}/denoised-paired_results', exist_ok=True)
       cmd=['qiime', 'dada2', 'denoise-paired',
            '--i-demultiplexed-seqs', f'{folder_results}/demux-paired-end.qza',
            '--p-trunc-len-f', '240',
            '--p-trunc-len-r', '200',
            '--p-trim-left-f' ,'0', 
            '--p-trim-left-r', '0',
            '--p-n-threads', '32',
            '--o-base-transition-stats', f'{folder_results}/denoised-paired_results/base-transition-stats.qza',
            '--o-table', f'{folder_results}/denoised-paired_results/table.qza',
            '--o-representative-sequences', f'{folder_results}/denoised-paired_results/rep-seqs.qza',
            '--o-denoising-stats', f'{folder_results}/denoised-paired_results/denoising-stats.qza']
       subprocess.run(cmd, check=True)
       
#denoise_paired(folder_results)       


def make_denoising_visualization(folder_results):
       cmd=['qiime', 'metadata', 'tabulate',
            '--m-input-file', f'{folder_results}/denoised-paired_results/denoising-stats.qza',
            '--o-visualization',f'{folder_results}/denoised-paired_results/denoising-stats.qzv']
       subprocess.run(cmd, check=True)

#make_denoising_visualization(folder_results)

denoised_folder='../results/denoised-paired_results'

def featuretable_summaries(denoised_folder, folder_results, metadata_file):
       output_folder=f'{folder_results}/feature-table_summarise'
       os.makedirs(output_folder, exist_ok=True)
       cmd_1=['qiime', 'feature-table', 'summarize',
            '--i-table', f'{denoised_folder}/table.qza',
            '--m-metadata-file', f'{metadata_file}',
            '--o-summary', f'{output_folder}/table.qzv',
            '--o-feature-frequencies', f'{output_folder}/feature-frequencies.qza',
            '--o-sample-frequencies', f'{output_folder}/sample-frequencies.qza'
            ]
       cmd_2=['qiime', 'feature-table', 'tabulate-seqs',
              '--i-data', f'{denoised_folder}/rep-seqs.qza',
              '--o-visualization', f'{output_folder}/rep-seqs.qzv']
       subprocess.run(cmd_1, check=True)
       subprocess.run(cmd_2,check=True)
       
#featuretable_summaries(denoised_folder, folder_results, metadata_file)       


#Phylogenetic analysis


def phylogenetic_analysis(denoised_folder, folder_results):
       output_folder=f'{folder_results}/phylogenetic_analysis'
       os.makedirs(output_folder, exist_ok=True)
       cmd=['qiime', 'phylogeny', 'align-to-tree-mafft-fasttree',
            '--i-sequences', f'{denoised_folder}/rep-seqs.qza',
            '--o-alignment', f'{output_folder}/aligned-rep-seqs.qza',
            '--o-masked-alignment', f'{output_folder}/masked-aligned-rep-seqs.qza',
            '--o-tree', f'{output_folder}/unrooted-tree.qza',
            '--o-rooted-tree', f'{output_folder}/rooted-tree.qza']
       subprocess.run(cmd, check=True)
       
#phylogenetic_analysis(denoised_folder, folder_results)      


#Alpha and beta diversity analysis 

phylogenetic_folder=f'{folder_results}/phylogenetic_analysis'

def diversity_metrics(
    denoised_folder,
    folder_results,
    phylogenetic_folder,
    metadata_file,
    overwrite=True
):
    sampling_depths = [500, 3200, 6000, 10000]

    denoised_folder = Path(denoised_folder)
    folder_results = Path(folder_results)
    phylogenetic_folder = Path(phylogenetic_folder)
    metadata_file = Path(metadata_file)

    table_qza = denoised_folder / "table.qza"
    rooted_tree_qza = phylogenetic_folder / "rooted-tree.qza"

    base_output_dir = folder_results / "diversity_metrics"
    base_output_dir.mkdir(parents=True, exist_ok=True)

    for depth in sampling_depths:
        output_dir = base_output_dir / f"sampling_depth_{depth}"

        if output_dir.exists():
            if overwrite:
                shutil.rmtree(output_dir)
            else:
                print(f"Skipping depth {depth}: {output_dir} already exists")
                continue

        cmd = [
            "qiime", "diversity", "core-metrics-phylogenetic",
            "--i-phylogeny", str(rooted_tree_qza),
            "--i-table", str(table_qza),
            "--p-sampling-depth", str(depth),
            "--m-metadata-file", str(metadata_file),
            "--output-dir", str(output_dir)
        ]

        print(f"\nRunning diversity analysis at sampling depth: {depth}")
        print(" ".join(cmd))

        subprocess.run(cmd, check=True)

    print("\nAll diversity analyses completed.")

#diversity_metrics(denoised_folder, folder_results, phylogenetic_folder, metadata_file, overwrite=True)    

def alpha_group_significance(
    folder_results,
    metadata_file,
    sampling_depths=None,
    overwrite=True
):
    """
    Run QIIME 2 alpha-group-significance for all alpha diversity vectors
    generated by core-metrics-phylogenetic at multiple sampling depths.
    """

    if sampling_depths is None:
        sampling_depths = [500, 3200, 6000, 10000]

    folder_results = Path(folder_results)
    metadata_file = Path(metadata_file)

    diversity_base = folder_results / "diversity_metrics"

    alpha_metrics = {
        "faith_pd": "faith_pd_vector.qza",
        "observed_features": "observed_features_vector.qza",
        "shannon": "shannon_vector.qza",
        "evenness": "evenness_vector.qza",
    }

    if not metadata_file.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_file}")

    for depth in sampling_depths:
        depth_dir = diversity_base / f"sampling_depth_{depth}"

        if not depth_dir.exists():
            print(f"Skipping depth {depth}: folder does not exist: {depth_dir}")
            continue

        output_dir = depth_dir / "alpha_group_significance"
        output_dir.mkdir(parents=True, exist_ok=True)

        for metric_name, vector_file in alpha_metrics.items():
            alpha_vector = depth_dir / vector_file
            output_qzv = output_dir / f"{metric_name}_group_significance.qzv"

            if not alpha_vector.exists():
                print(f"Skipping {metric_name} at depth {depth}: missing {alpha_vector}")
                continue

            if output_qzv.exists():
                if overwrite:
                    output_qzv.unlink()
                else:
                    print(f"Skipping existing file: {output_qzv}")
                    continue

            cmd = [
                "qiime", "diversity", "alpha-group-significance",
                "--i-alpha-diversity", str(alpha_vector),
                "--m-metadata-file", str(metadata_file),
                "--o-visualization", str(output_qzv),
            ]

            print(f"\nRunning alpha-group-significance")
            print(f"Depth: {depth}")
            print(f"Metric: {metric_name}")
            print(" ".join(cmd))

            subprocess.run(cmd, check=True)

    print("\nAlpha-group-significance analyses completed.")
    

'''   
alpha_group_significance(
    folder_results,
    metadata_file,
    sampling_depths=[500, 3200, 6000, 10000],
    overwrite=True
)
'''

def safe_name(value):
    value = str(value)
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value.strip("_")


def get_valid_grouping_columns(metadata_file, max_unique_groups=12, min_group_size=2):
    """
    Select metadata columns suitable for QIIME2 beta-group-significance.

    A valid column must:
    - not be a sample ID, barcode, or filepath column
    - have at least 2 groups
    - not have all values unique
    - have at least one group with >= min_group_size samples
    - not have too many unique groups
    """

    metadata_file = Path(metadata_file)
    metadata = pd.read_csv(metadata_file, sep="\t", dtype=str)

    # Remove QIIME2 metadata type row if present
    first_col = metadata.columns[0]
    metadata = metadata[metadata[first_col] != "#q2:types"]

    excluded_exact = {
        "sample-id",
        "sampleid",
        "sample_id",
        "#sampleid",
        "#sample-id",
        "id",
        "barcode-sequence",
        "forward-absolute-filepath",
        "reverse-absolute-filepath",
        "sample_base_id",
    }

    excluded_patterns = [
        "filepath",
        "barcode",
        "sequence",
        "absolute",
        "path",
    ]

    valid_columns = []
    skipped_columns = {}

    for col in metadata.columns:
        col_norm = col.strip().lower()

        if col_norm in excluded_exact:
            skipped_columns[col] = "identifier, barcode, or filepath column"
            continue

        if any(pattern in col_norm for pattern in excluded_patterns):
            skipped_columns[col] = "identifier, barcode, or filepath-like column"
            continue

        values = metadata[col].dropna()

        if values.empty:
            skipped_columns[col] = "empty column"
            continue

        counts = values.value_counts()
        n_groups = len(counts)

        if n_groups < 2:
            skipped_columns[col] = "only one group"
            continue

        if n_groups == len(values):
            skipped_columns[col] = "all values are unique"
            continue

        if counts.max() < min_group_size:
            skipped_columns[col] = f"no group has at least {min_group_size} samples"
            continue

        if n_groups > max_unique_groups:
            skipped_columns[col] = f"too many groups: {n_groups}"
            continue

        valid_columns.append(col)

    print("\nValid metadata columns for beta-group-significance:")
    for col in valid_columns:
        print(f"- {col}")

    print("\nSkipped metadata columns:")
    for col, reason in skipped_columns.items():
        print(f"- {col}: {reason}")

    return valid_columns


def beta_group_significance(
    folder_results,
    metadata_file,
    sampling_depths=None,
    beta_metrics=None,
    pairwise=True,
    overwrite=True,
    max_unique_groups=12
):
    if sampling_depths is None:
        sampling_depths = [500, 3200, 6000, 10000]

    if beta_metrics is None:
        beta_metrics = {
            "unweighted_unifrac": "unweighted_unifrac_distance_matrix.qza",
            "weighted_unifrac": "weighted_unifrac_distance_matrix.qza",
            "jaccard": "jaccard_distance_matrix.qza",
            "bray_curtis": "bray_curtis_distance_matrix.qza",
        }

    folder_results = Path(folder_results)
    metadata_file = Path(metadata_file)

    metadata_columns = get_valid_grouping_columns(
        metadata_file=metadata_file,
        max_unique_groups=max_unique_groups,
        min_group_size=2
    )

    diversity_base = folder_results / "diversity_metrics"

    for depth in sampling_depths:
        depth_dir = diversity_base / f"sampling_depth_{depth}"

        if not depth_dir.exists():
            print(f"\nSkipping sampling depth {depth}: folder not found: {depth_dir}")
            continue

        beta_output_base = depth_dir / "beta_group_significance"
        beta_output_base.mkdir(parents=True, exist_ok=True)

        for beta_name, distance_file in beta_metrics.items():
            distance_matrix = depth_dir / distance_file

            if not distance_matrix.exists():
                print(f"\nSkipping {beta_name} at depth {depth}: missing {distance_matrix}")
                continue

            metric_output_dir = beta_output_base / beta_name
            metric_output_dir.mkdir(parents=True, exist_ok=True)

            for column in metadata_columns:
                column_safe = safe_name(column)

                column_output_dir = metric_output_dir / column_safe
                column_output_dir.mkdir(parents=True, exist_ok=True)

                output_qzv = column_output_dir / f"{beta_name}_{column_safe}_group_significance.qzv"

                if output_qzv.exists():
                    if overwrite:
                        output_qzv.unlink()
                    else:
                        print(f"Skipping existing file: {output_qzv}")
                        continue

                cmd = [
                    "qiime", "diversity", "beta-group-significance",
                    "--i-distance-matrix", str(distance_matrix),
                    "--m-metadata-file", str(metadata_file),
                    "--m-metadata-column", column,
                    "--o-visualization", str(output_qzv),
                ]

                if pairwise:
                    cmd.append("--p-pairwise")

                print("\nRunning beta-group-significance")
                print(f"Sampling depth: {depth}")
                print(f"Beta metric: {beta_name}")
                print(f"Metadata column: {column}")
                print(" ".join(cmd))

                try:
                    subprocess.run(cmd, check=True)
                except subprocess.CalledProcessError:
                    print(f"WARNING: QIIME failed for column '{column}' at depth {depth} using {beta_name}. Skipping.")
                    continue

    print("\nBeta-group-significance analyses completed.")

'''
beta_group_significance(
    folder_results=folder_results,
    metadata_file=metadata_file,
    sampling_depths=[500, 3200, 6000, 10000],
    pairwise=True,
    overwrite=True
)
'''




# Alpha rarefaction plotting 



def alpha_rarefaction(
    denoised_folder,
    folder_results,
    phylogenetic_folder,
    metadata_file,
    max_depth=10000,
    overwrite=True
):
    """
    Run QIIME 2 alpha-rarefaction.

    This generates an interactive .qzv visualization to evaluate whether
    sequencing depth is sufficient for alpha-diversity estimation.
    """

    denoised_folder = Path(denoised_folder)
    folder_results = Path(folder_results)
    phylogenetic_folder = Path(phylogenetic_folder)
    metadata_file = Path(metadata_file)

    table_qza = denoised_folder / "table.qza"
    rooted_tree_qza = phylogenetic_folder / "rooted-tree.qza"

    output_dir = folder_results / "alpha_rarefaction"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_qzv = output_dir / f"alpha_rarefaction_max_depth_{max_depth}.qzv"

    if output_qzv.exists():
        if overwrite:
            output_qzv.unlink()
        else:
            print(f"Skipping alpha rarefaction because output already exists: {output_qzv}")
            return

    cmd = [
        "qiime", "diversity", "alpha-rarefaction",
        "--i-table", str(table_qza),
        "--i-phylogeny", str(rooted_tree_qza),
        "--p-max-depth", str(max_depth),
        "--m-metadata-file", str(metadata_file),
        "--o-visualization", str(output_qzv),
    ]

    print("\nRunning alpha rarefaction")
    print(" ".join(cmd))

    subprocess.run(cmd, check=True)

    print(f"\nAlpha rarefaction completed: {output_qzv}")


alpha_rarefaction(
    denoised_folder=denoised_folder,
    folder_results=folder_results,
    phylogenetic_folder=phylogenetic_folder,
    metadata_file=metadata_file,
    max_depth=10000,
    overwrite=True
)


# Taxonomic analysis



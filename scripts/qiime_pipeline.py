#Import the libraries

from pathlib import Path
import subprocess
import shutil
import os
import pandas as pd
import re

# ============================================================
# BASE PATHS
# ============================================================

metadata_file = Path("../data/input_files/manifest_curated_nugent.tsv")
folder_results = Path("../results")
pretrained_model = Path('../data/classifiers/silva-138-99-nb-classifier.qza')


# ============================================================
# PATH DICTIONARY
# ============================================================

paths = {
    "metadata_dir": folder_results / "00_metadata",
    "import_dir": folder_results / "01_import",
    "denoise_dir": folder_results / "02_dada2_denoise",
    "feature_summary_dir": folder_results / "03_feature_table_summary",
    "phylogeny_dir": folder_results / "04_phylogeny",
    "diversity_dir": folder_results / "05_diversity",
    "taxonomy_dir": folder_results / "06_taxonomy",
    "barplot_dir": folder_results / "07_taxa_barplots",
    "export_dir": folder_results / "08_exported_files"
}


# ============================================================
# CREATE ONLY DIRECTORY PATHS
# ============================================================

def create_defined_folders(paths):
    for key, path in paths.items():
        if key.endswith("_dir"):
            path.mkdir(parents=True, exist_ok=True)


create_defined_folders(paths)

def create_metadata_tabulate(metadata_file):
    cmd = ['qiime', 'metadata', 'tabulate',
           '--m-input-file', f'{metadata_file}',
           '--o-visualization', paths['metadata_dir']/'sample-metadata-viz.qzv']
    subprocess.run(cmd, check=True)
    print("sample-metadata-viz.qzv created")

create_metadata_tabulate(metadata_file)


def import_data(metadata_file):
    cmd = ['qiime', 'tools', 'import',
           '--type', 'SampleData[PairedEndSequencesWithQuality]',
           '--input-path', f'{metadata_file}',
           '--output-path', paths['import_dir']/'demux-paired-end.qza',
           '--input-format', 'PairedEndFastqManifestPhred33V2']
    subprocess.run(cmd, check=True)
    print("demux-paired-end.qza created")

import_data(metadata_file)


def demux_summarize():
    cmd = ['qiime', 'demux', 'summarize',
           '--i-data', paths['import_dir']/'demux-paired-end.qza',
           '--o-visualization', paths['import_dir']/'demux-paired-end.qzv']
    subprocess.run(cmd, check=True)
    print("demux-paired-end.qzv created")

#demux_summarize()

def denoise():
    cmd=['qiime', 'dada2', 'denoise-paired',
            '--i-demultiplexed-seqs', paths['import_dir']/'demux-paired-end.qza',
            '--p-trim-left-f' ,'0',
            '--p-trim-left-r', '0',
            '--p-trunc-len-f', '0',
            '--p-trunc-len-r', '0',
            '--p-n-threads', '32',
            '--o-table', paths['denoise_dir']/'table.qza',
            '--o-representative-sequences', paths['denoise_dir']/'rep-seqs.qza',
            '--o-denoising-stats', paths['denoise_dir']/'denoising-stats.qza']
    subprocess.run(cmd, check=True)
    print("denoising step finished")

denoise()


def make_denoising_visualization():
    cmd=['qiime', 'metadata', 'tabulate',
         '--m-input-file', paths['denoise_dir']/'denoising-stats.qza',
         '--o-visualization',paths['denoise_dir']/'denoising-stats.qzv']
    subprocess.run(cmd, check=True)
    print("denoising step finished")

make_denoising_visualization()


def featuretable_summaries(metadata_file):
       cmd_1=['qiime', 'feature-table', 'summarize',
            '--i-table', paths['denoise_dir']/'table.qza',
            '--m-sample-metadata-file', f'{metadata_file}',
            '--o-visualization', paths['denoise_dir']/'table.qzv']
       cmd_2=['qiime', 'feature-table', 'tabulate-seqs',
              '--i-data', paths['denoise_dir']/'rep-seqs.qza',
              '--o-visualization', paths['denoise_dir']/'rep-seqs.qzv']
       subprocess.run(cmd_1, check=True)
       subprocess.run(cmd_2,check=True)

featuretable_summaries(metadata_file)


#Phylogenetic analysis


def phylogenetic_analysis():
       cmd=['qiime', 'phylogeny', 'align-to-tree-mafft-fasttree',
            '--i-sequences', paths['denoise_dir']/'rep-seqs.qza',
            '--o-alignment', paths['denoise_dir']/'aligned-rep-seqs.qza',
            '--o-masked-alignment', paths['phylogeny_dir']/'masked-aligned-rep-seqs.qza',
            '--o-tree', paths['phylogeny_dir']/'unrooted-tree.qza',
            '--o-rooted-tree', paths['phylogeny_dir']/'rooted-tree.qza']
       subprocess.run(cmd, check=True)

#phylogenetic_analysis()


#Alpha and beta diversity analysis


def diversity_metrics():
    sample_depth=[500,3200,6000,12000]
    for x in sample_depth:

        cmd=["qiime", "diversity", "core-metrics-phylogenetic",
             "--i-phylogeny", paths['phylogeny_dir']/'rooted-tree.qza',
             "--i-table", paths['denoise_dir']/'table.qza',
             "--p-sampling-depth", f'{x}',
            "--m-metadata-file", f'{metadata_file}',
           '--output-dir', paths['diversity_dir']/f'{x}']
        subprocess.run(cmd, check=True)

#diversity_metrics()


#Beta-group-significance analyses complete

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


def alpha_rarefaction(metadata_file):
    output_folder=paths["diversity_dir"]/'alpha-rarefaction'
    os.makedirs(output_folder, exist_ok=True)
    cmd=['qiime', 'diversity', 'alpha-rarefaction',
         '--i-table', paths['denoise_dir']/'table.qza',
        '--i-phylogeny', paths['phylogeny_dir']/'rooted-tree.qza',
        '--p-max-depth', '1000',
        '--m-metadata-file', f'{metadata_file}',
        '--o-visualization', f'{output_folder}/alpha-rarefaction.qzv']
    subprocess.run(cmd, check=True)

#alpha_rarefaction(metadata_file)

# Taxonomic analysis

def taxonomy_assigment(pretrained_model, metadata_file):
    cmd_1=['qiime', 'feature-classifier', 'classify-sklearn',
           '--i-classifier', f'{pretrained_model}',
           '--i-reads', paths['denoise_dir']/'rep-seqs.qza',
           '--o-classification', paths['taxonomy_dir']/'taxonomy.qza']
    cmd_2=['qiime', 'metadata', 'tabulate',
           '--m-input-file', paths['taxonomy_dir']/'taxonomy.qza',
           '--o-visualization', paths['taxonomy_dir']/'taxonomy.qzv']
    cmd_3=['qiime', 'taxa', 'barplot',
           '--i-table', paths['denoise_dir']/'table.qza',
           '--i-taxonomy', paths['taxonomy_dir']/'taxonomy.qza',
           '--m-metadata-file', f'{metadata_file}',
           '--o-visualization', paths['barplot_dir']/'taxa-bar-plots.qzv' ]
    subprocess.run(cmd_1, check=True)
    subprocess.run(cmd_2, check=True)
    subprocess.run(cmd_3, check=True)

#taxonomy_assigment(pretrained_model, metadata_file)






# Extract qza information


def extract():
    for x in os.listdir(paths['denoise_dir']):
        if x.endswith('rep-seqs.qza') or x.endswith('table.qza'):
            file=os.path.join(paths['denoise_dir'],x)
            cmd=['qiime', 'tools', 'extract',
                 '--input-path', f'{file}',
                 '--output-path',paths['export_dir']/f'{x}']
            subprocess.run(cmd, check=True)

#extract()


def extract_taxonomy():
    for x in os.listdir(paths['taxonomy_dir']):
        if x.endswith('.qza'):
            file=os.path.join(paths['taxonomy_dir'],x)
            cmd=['qiime', 'tools', 'extract',
                 '--input-path',f'{file}',
                 '--output-path',paths['export_dir']/'silva_taxonomy']
            subprocess.run(cmd, check=True)

#extract_taxonomy()

## Run speciate IT to assign taxonomic classification

fasta_sequences='../results/08_exported_files/rep-seqs.qza/31c77f38-a305-49c1-b877-a39b046d0c67/data/dna-sequences.fasta'
program='/home/rare/programs/speciateIT/vSpeciateDB_models/vSpeciateIT_V3V4'


def classification(fasta_sequences, program):
    output_folder='../results/09_speciate_it_classification'
    os.makedirs(output_folder, exist_ok=True)
    cmd=['classify',
         '-d', f'{program}',
         '-i', f'{fasta_sequences}',
         '-o', f'{output_folder}']
    subprocess.run(cmd, check=True)

#classification(fasta_sequences,program)

from pathlib import Path
from datetime import datetime
from uuid import uuid4
import subprocess
import shutil
import re
import pandas as pd


# ============================================================
# CONFIG
# ============================================================

PROJECT_ROOT = Path.cwd().resolve()

MANIFEST_FILE = PROJECT_ROOT / "data" / "input_files" / "manifest_curated_nugent.tsv"

# If you have a separate metadata file, change this.
# If not, the manifest will also be used as metadata.
METADATA_FILE = MANIFEST_FILE

RESULTS_ROOT = PROJECT_ROOT / "results" / "qiime2_unique_runs"

CLASSIFIERS = {
    "backbone_2024_09_full_length": PROJECT_ROOT / "data" / "classifiers" / "2024.09.backbone.full-length.nb.sklearn-1.4.2.qza",
    "silva_138_99": PROJECT_ROOT / "data" / "classifiers" / "silva-138-99-nb-classifier.qza",
}

THREADS = 32

DADA2_TRIM_LEFT_F = 0
DADA2_TRIM_LEFT_R = 0
DADA2_TRUNC_LEN_F = 0
DADA2_TRUNC_LEN_R = 0

SAMPLING_DEPTHS = [500, 3200, 6000, 10000]
ALPHA_RAREFACTION_MAX_DEPTH = 10000

RUN_ID = datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + uuid4().hex[:8]
RUN_DIR = RESULTS_ROOT / f"run_{RUN_ID}"


# ============================================================
# BASIC UTILITIES
# ============================================================

def run_cmd(cmd):
    cmd = [str(x) for x in cmd]
    print("\n$ " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def require_file(path, label):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")
    return path


def make_dir(path):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def safe_name(value):
    value = str(value)
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value.strip("_")


def initialize_unique_run():
    if RUN_DIR.exists():
        raise FileExistsError(f"RUN_DIR already exists. This should not happen: {RUN_DIR}")

    make_dir(RUN_DIR)

    print("\nUnique QIIME2 run created")
    print(f"RUN_ID: {RUN_ID}")
    print(f"RUN_DIR: {RUN_DIR}")

    return RUN_DIR


# ============================================================
# PATHS FOR THIS UNIQUE RUN
# ============================================================

def get_run_paths(run_dir):
    run_dir = Path(run_dir)

    paths = {
        "run_dir": run_dir,

        "metadata_viz_dir": run_dir / "00_metadata",
        "import_dir": run_dir / "01_import",
        "denoise_dir": run_dir / "02_dada2_denoise",
        "feature_summary_dir": run_dir / "03_feature_table_summary",
        "phylogeny_dir": run_dir / "04_phylogeny",
        "diversity_dir": run_dir / "05_diversity",
        "taxonomy_dir": run_dir / "06_taxonomy",
        "barplot_dir": run_dir / "07_taxa_barplots",

        "demux_qza": run_dir / "01_import" / "demux-paired-end.qza",
        "demux_qzv": run_dir / "01_import" / "demux-paired-end.qzv",

        "table_qza": run_dir / "02_dada2_denoise" / "table.qza",
        "rep_seqs_qza": run_dir / "02_dada2_denoise" / "rep-seqs.qza",
        "denoising_stats_qza": run_dir / "02_dada2_denoise" / "denoising-stats.qza",
        "denoising_stats_qzv": run_dir / "02_dada2_denoise" / "denoising-stats.qzv",

        "table_qzv": run_dir / "03_feature_table_summary" / "table.qzv",
        "rep_seqs_qzv": run_dir / "03_feature_table_summary" / "rep-seqs.qzv",

        "aligned_rep_seqs_qza": run_dir / "04_phylogeny" / "aligned-rep-seqs.qza",
        "masked_aligned_rep_seqs_qza": run_dir / "04_phylogeny" / "masked-aligned-rep-seqs.qza",
        "unrooted_tree_qza": run_dir / "04_phylogeny" / "unrooted-tree.qza",
        "rooted_tree_qza": run_dir / "04_phylogeny" / "rooted-tree.qza",
    }

    for key, path in paths.items():
        if key.endswith("_dir"):
            make_dir(path)

    return paths


# ============================================================
# STEP 00: METADATA VISUALIZATION
# ============================================================

def create_metadata_tabulate(metadata_file, paths):
    metadata_file = require_file(metadata_file, "Metadata file")

    output_qzv = paths["metadata_viz_dir"] / "sample-metadata-viz.qzv"

    run_cmd([
        "qiime", "metadata", "tabulate",
        "--m-input-file", metadata_file,
        "--o-visualization", output_qzv,
    ])

    print(f"\nMetadata visualization created: {output_qzv}")


# ============================================================
# STEP 01: IMPORT PAIRED-END DATA
# ============================================================

def import_paired_end_data(manifest_file, paths):
    manifest_file = require_file(manifest_file, "Manifest file")

    run_cmd([
        "qiime", "tools", "import",
        "--type", "SampleData[PairedEndSequencesWithQuality]",
        "--input-path", manifest_file,
        "--output-path", paths["demux_qza"],
        "--input-format", "PairedEndFastqManifestPhred33V2",
    ])

    print(f"\nImported paired-end data: {paths['demux_qza']}")


def demux_summarize(paths):
    run_cmd([
        "qiime", "demux", "summarize",
        "--i-data", paths["demux_qza"],
        "--o-visualization", paths["demux_qzv"],
    ])

    print(f"\nDemux summary created: {paths['demux_qzv']}")


# ============================================================
# STEP 02: DADA2 DENOISING
# ============================================================

def denoise_paired(paths):
    run_cmd([
        "qiime", "dada2", "denoise-paired",
        "--i-demultiplexed-seqs", paths["demux_qza"],

        "--p-trim-left-f", DADA2_TRIM_LEFT_F,
        "--p-trim-left-r", DADA2_TRIM_LEFT_R,
        "--p-trunc-len-f", DADA2_TRUNC_LEN_F,
        "--p-trunc-len-r", DADA2_TRUNC_LEN_R,

        "--p-n-threads", THREADS,

        "--o-table", paths["table_qza"],
        "--o-representative-sequences", paths["rep_seqs_qza"],
        "--o-denoising-stats", paths["denoising_stats_qza"],
    ])

    print("\nDADA2 denoising completed")


def make_denoising_visualization(paths):
    run_cmd([
        "qiime", "metadata", "tabulate",
        "--m-input-file", paths["denoising_stats_qza"],
        "--o-visualization", paths["denoising_stats_qzv"],
    ])

    print(f"\nDenoising stats visualization created: {paths['denoising_stats_qzv']}")


# ============================================================
# STEP 03: FEATURE TABLE SUMMARIES
# ============================================================

def feature_table_summaries(metadata_file, paths):
    metadata_file = require_file(metadata_file, "Metadata file")

    run_cmd([
        "qiime", "feature-table", "summarize",
        "--i-table", paths["table_qza"],
        "--m-sample-metadata-file", metadata_file,
        "--o-visualization", paths["table_qzv"],
    ])

    run_cmd([
        "qiime", "feature-table", "tabulate-seqs",
        "--i-data", paths["rep_seqs_qza"],
        "--o-visualization", paths["rep_seqs_qzv"],
    ])

    print("\nFeature table summaries created")


# ============================================================
# STEP 04: PHYLOGENETIC ANALYSIS
# ============================================================

def phylogenetic_analysis(paths):
    run_cmd([
        "qiime", "phylogeny", "align-to-tree-mafft-fasttree",
        "--i-sequences", paths["rep_seqs_qza"],
        "--o-alignment", paths["aligned_rep_seqs_qza"],
        "--o-masked-alignment", paths["masked_aligned_rep_seqs_qza"],
        "--o-tree", paths["unrooted_tree_qza"],
        "--o-rooted-tree", paths["rooted_tree_qza"],
    ])

    print("\nPhylogenetic analysis completed")


# ============================================================
# STEP 05: CORE DIVERSITY METRICS
# ============================================================

def diversity_metrics(metadata_file, paths, sampling_depths):
    metadata_file = require_file(metadata_file, "Metadata file")

    for depth in sampling_depths:
        output_dir = paths["diversity_dir"] / f"sampling_depth_{depth}"
        make_dir(output_dir)

        run_cmd([
            "qiime", "diversity", "core-metrics-phylogenetic",
            "--i-phylogeny", paths["rooted_tree_qza"],
            "--i-table", paths["table_qza"],
            "--p-sampling-depth", depth,
            "--m-metadata-file", metadata_file,
            "--output-dir", output_dir,
        ])

        print(f"\nDiversity metrics completed for sampling depth: {depth}")

    print("\nAll diversity metrics completed")


# ============================================================
# STEP 06: ALPHA GROUP SIGNIFICANCE
# ============================================================

def alpha_group_significance(metadata_file, paths, sampling_depths):
    metadata_file = require_file(metadata_file, "Metadata file")

    alpha_metrics = {
        "faith_pd": "faith_pd_vector.qza",
        "observed_features": "observed_features_vector.qza",
        "shannon": "shannon_vector.qza",
        "evenness": "evenness_vector.qza",
    }

    for depth in sampling_depths:
        depth_dir = paths["diversity_dir"] / f"sampling_depth_{depth}"

        if not depth_dir.exists():
            print(f"Skipping depth {depth}. Folder not found: {depth_dir}")
            continue

        output_dir = depth_dir / "alpha_group_significance"
        make_dir(output_dir)

        for metric_name, vector_file in alpha_metrics.items():
            alpha_vector = depth_dir / vector_file
            output_qzv = output_dir / f"{metric_name}_group_significance.qzv"

            if not alpha_vector.exists():
                print(f"Skipping {metric_name} at depth {depth}. Missing: {alpha_vector}")
                continue

            run_cmd([
                "qiime", "diversity", "alpha-group-significance",
                "--i-alpha-diversity", alpha_vector,
                "--m-metadata-file", metadata_file,
                "--o-visualization", output_qzv,
            ])

    print("\nAlpha group significance completed")


# ============================================================
# STEP 07: BETA GROUP SIGNIFICANCE
# ============================================================

def get_valid_grouping_columns(metadata_file, max_unique_groups=12, min_group_size=2):
    metadata_file = require_file(metadata_file, "Metadata file")
    metadata = pd.read_csv(metadata_file, sep="\t", dtype=str)

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


def beta_group_significance(metadata_file, paths, sampling_depths, pairwise=True):
    metadata_file = require_file(metadata_file, "Metadata file")

    beta_metrics = {
        "unweighted_unifrac": "unweighted_unifrac_distance_matrix.qza",
        "weighted_unifrac": "weighted_unifrac_distance_matrix.qza",
        "jaccard": "jaccard_distance_matrix.qza",
        "bray_curtis": "bray_curtis_distance_matrix.qza",
    }

    metadata_columns = get_valid_grouping_columns(
        metadata_file=metadata_file,
        max_unique_groups=12,
        min_group_size=2,
    )

    for depth in sampling_depths:
        depth_dir = paths["diversity_dir"] / f"sampling_depth_{depth}"

        if not depth_dir.exists():
            print(f"Skipping depth {depth}. Folder not found: {depth_dir}")
            continue

        beta_output_base = depth_dir / "beta_group_significance"
        make_dir(beta_output_base)

        for beta_name, distance_file in beta_metrics.items():
            distance_matrix = depth_dir / distance_file

            if not distance_matrix.exists():
                print(f"Skipping {beta_name} at depth {depth}. Missing: {distance_matrix}")
                continue

            metric_output_dir = beta_output_base / beta_name
            make_dir(metric_output_dir)

            for column in metadata_columns:
                column_safe = safe_name(column)
                column_output_dir = metric_output_dir / column_safe
                make_dir(column_output_dir)

                output_qzv = column_output_dir / f"{beta_name}_{column_safe}_group_significance.qzv"

                cmd = [
                    "qiime", "diversity", "beta-group-significance",
                    "--i-distance-matrix", distance_matrix,
                    "--m-metadata-file", metadata_file,
                    "--m-metadata-column", column,
                    "--o-visualization", output_qzv,
                ]

                if pairwise:
                    cmd.append("--p-pairwise")

                try:
                    run_cmd(cmd)
                except subprocess.CalledProcessError:
                    print(
                        f"WARNING: QIIME failed for column '{column}', "
                        f"depth {depth}, beta metric {beta_name}. Skipping."
                    )

    print("\nBeta group significance completed")


# ============================================================
# STEP 08: ALPHA RAREFACTION
# ============================================================

def alpha_rarefaction(metadata_file, paths, max_depth):
    metadata_file = require_file(metadata_file, "Metadata file")

    output_dir = paths["diversity_dir"] / "alpha_rarefaction"
    make_dir(output_dir)

    output_qzv = output_dir / f"alpha_rarefaction_max_depth_{max_depth}.qzv"

    run_cmd([
        "qiime", "diversity", "alpha-rarefaction",
        "--i-table", paths["table_qza"],
        "--i-phylogeny", paths["rooted_tree_qza"],
        "--p-max-depth", max_depth,
        "--m-metadata-file", metadata_file,
        "--o-visualization", output_qzv,
    ])

    print(f"\nAlpha rarefaction completed: {output_qzv}")


# ============================================================
# STEP 09: TAXONOMIC ASSIGNMENT
# ============================================================

def taxonomic_assignment(paths, classifiers):
    for classifier_name, classifier_path in classifiers.items():
        classifier_path = require_file(classifier_path, f"Classifier {classifier_name}")

        output_dir = paths["taxonomy_dir"] / classifier_name
        make_dir(output_dir)

        taxonomy_qza = output_dir / "taxonomy.qza"
        taxonomy_qzv = output_dir / "taxonomy.qzv"

        run_cmd([
            "qiime", "feature-classifier", "classify-sklearn",
            "--i-classifier", classifier_path,
            "--i-reads", paths["rep_seqs_qza"],
            "--o-classification", taxonomy_qza,
        ])

        run_cmd([
            "qiime", "metadata", "tabulate",
            "--m-input-file", taxonomy_qza,
            "--o-visualization", taxonomy_qzv,
        ])

        print(f"\nTaxonomy completed for classifier: {classifier_name}")

    print("\nTaxonomic assignment completed for all classifiers")


# ============================================================
# STEP 10: TAXA BARPLOTS
# ============================================================

def taxa_barplots(metadata_file, paths):
    metadata_file = require_file(metadata_file, "Metadata file")

    taxonomy_files = sorted(paths["taxonomy_dir"].glob("*/taxonomy.qza"))

    if not taxonomy_files:
        raise FileNotFoundError(f"No taxonomy.qza files found inside: {paths['taxonomy_dir']}")

    for taxonomy_qza in taxonomy_files:
        classifier_name = taxonomy_qza.parent.name

        output_dir = paths["barplot_dir"] / classifier_name
        make_dir(output_dir)

        output_qzv = output_dir / "taxa-bar-plots.qzv"

        run_cmd([
            "qiime", "taxa", "barplot",
            "--i-table", paths["table_qza"],
            "--i-taxonomy", taxonomy_qza,
            "--m-metadata-file", metadata_file,
            "--o-visualization", output_qzv,
        ])

        print(f"\nTaxa barplot completed for classifier: {classifier_name}")

    print("\nTaxa barplots completed for all classifiers")


# ============================================================
# MAIN PIPELINE
# ============================================================

def main():
    require_file(MANIFEST_FILE, "Manifest file")
    require_file(METADATA_FILE, "Metadata file")

    for classifier_name, classifier_path in CLASSIFIERS.items():
        require_file(classifier_path, f"Classifier {classifier_name}")

    run_dir = initialize_unique_run()
    paths = get_run_paths(run_dir)

    create_metadata_tabulate(METADATA_FILE, paths)

    import_paired_end_data(MANIFEST_FILE, paths)
    demux_summarize(paths)

    denoise_paired(paths)
    make_denoising_visualization(paths)

    feature_table_summaries(METADATA_FILE, paths)

    phylogenetic_analysis(paths)

    diversity_metrics(
        metadata_file=METADATA_FILE,
        paths=paths,
        sampling_depths=SAMPLING_DEPTHS,
    )

    alpha_group_significance(
        metadata_file=METADATA_FILE,
        paths=paths,
        sampling_depths=SAMPLING_DEPTHS,
    )

    beta_group_significance(
        metadata_file=METADATA_FILE,
        paths=paths,
        sampling_depths=SAMPLING_DEPTHS,
        pairwise=True,
    )

    alpha_rarefaction(
        metadata_file=METADATA_FILE,
        paths=paths,
        max_depth=ALPHA_RAREFACTION_MAX_DEPTH,
    )

    taxonomic_assignment(
        paths=paths,
        classifiers=CLASSIFIERS,
    )

    taxa_barplots(
        metadata_file=METADATA_FILE,
        paths=paths,
    )

    print("\nPipeline completed successfully")
    print(f"All outputs are inside this unique folder:\n{run_dir}")


if __name__ == "__main__":
    main()
#!/usr/bin/env python3

from pathlib import Path
import subprocess
import shutil
import numpy as np
import pandas as pd


# ============================================================
# CONFIGURATION
# Paths based on your current QIIME pipeline script
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent

metadata_file = REPO_ROOT / "data/input_files/manifest_curated_nugent.tsv"
folder_results = REPO_ROOT / "results"
denoised_folder = folder_results / "denoised-paired_results"

rep_seqs_qza = denoised_folder / "rep-seqs.qza"
table_qza = denoised_folder / "table.qza"

speciateit_dir = Path("/home/rare/programs/speciateIT")
speciateit_binary = speciateit_dir / "bin/linux/classify"

output_folder = folder_results / "speciateit_refinement"

overwrite = True
force_species = False

# If force_species = True, speciateIT will use --skip-err-thld.
# For the main analysis, keep this False.


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def run_command(cmd):
    cmd = [str(x) for x in cmd]
    print("\n" + " ".join(cmd))
    subprocess.run(cmd, check=True)


def find_speciateit_model(speciateit_dir):
    """
    Find the V3V4 speciateIT model inside the external speciateIT folder.
    """

    candidate_models = [
        speciateit_dir / "vSpeciateDB_models/vSpeciateIT_V3V4",
        speciateit_dir / "vSpeciateDB_models/V3V4",
        speciateit_dir / "vSpeciateDB_models/vSpeciateIT_v3v4",
    ]

    for model in candidate_models:
        if model.exists():
            return model

    models_root = speciateit_dir / "vSpeciateDB_models"

    if models_root.exists():
        for model in models_root.iterdir():
            if model.is_dir() and "v3v4" in model.name.lower():
                return model

    raise FileNotFoundError(
        f"Could not find the speciateIT V3V4 model inside: {models_root}\n"
        "Check the installed model folder name."
    )


def check_inputs():
    if not metadata_file.exists():
        print(f"\nWARNING: Metadata file not found: {metadata_file}")
        print("The script will skip taxa barplot.")

    if not rep_seqs_qza.exists():
        raise FileNotFoundError(f"Missing representative sequences file: {rep_seqs_qza}")

    if not table_qza.exists():
        raise FileNotFoundError(f"Missing feature table file: {table_qza}")

    if not speciateit_dir.exists():
        raise FileNotFoundError(f"speciateIT folder not found: {speciateit_dir}")

    if not speciateit_binary.exists():
        raise FileNotFoundError(
            f"speciateIT binary not found: {speciateit_binary}\n"
            "Check that speciateIT was installed correctly."
        )


def convert_speciateit_name_to_taxon(classification):
    """
    Convert speciateIT species names into QIIME-style taxonomy strings.

    Example:
    Gardnerella_vaginalis
    -> d__Bacteria; g__Gardnerella; s__Gardnerella vaginalis
    """

    if pd.isna(classification):
        return (
            "d__Bacteria; g__Unclassified; s__Unclassified",
            "Unclassified",
            "Unclassified",
        )

    raw = str(classification).strip()

    if raw == "" or raw.lower() in {"na", "none", "unknown", "unclassified"}:
        return (
            "d__Bacteria; g__Unclassified; s__Unclassified",
            "Unclassified",
            "Unclassified",
        )

    parts = raw.split("_")

    if len(parts) >= 3 and parts[0] in {"Ca", "Candidatus"}:
        genus = f"Candidatus {parts[1]}"
        species = "Candidatus " + " ".join(parts[1:])
    elif len(parts) >= 2:
        genus = parts[0]
        species = " ".join(parts)
    else:
        genus = parts[0]
        species = parts[0]

    taxon = f"d__Bacteria; g__{genus}; s__{species}"

    return taxon, genus, species


def parse_speciateit_results(speciateit_results):
    rows = []

    with open(speciateit_results, "r") as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            parts = line.split()

            if len(parts) < 2:
                continue

            feature_id = parts[0]

            if feature_id.lower() in {
                "sequence",
                "sequenceid",
                "sequence_id",
                "feature",
                "featureid",
                "feature_id",
                "asv",
                "asvid",
                "asv_id",
            }:
                continue

            classification = parts[1]
            confidence = parts[2] if len(parts) >= 3 else "0"
            decisions = parts[3] if len(parts) >= 4 else "0"

            taxon, genus, species = convert_speciateit_name_to_taxon(classification)

            rows.append(
                {
                    "Feature ID": feature_id,
                    "Taxon": taxon,
                    "Combined_Taxon": taxon,
                    "Confidence": confidence,
                    "speciateIT_classification": classification,
                    "speciateIT_confidence": confidence,
                    "speciateIT_decisions": decisions,
                    "Genus": genus,
                    "Species": species,
                }
            )

    if not rows:
        raise ValueError(f"No valid speciateIT records were parsed from: {speciateit_results}")

    taxonomy_df = pd.DataFrame(rows)

    taxonomy_df["Confidence"] = pd.to_numeric(
        taxonomy_df["Confidence"],
        errors="coerce",
    ).fillna(0)

    taxonomy_df["speciateIT_confidence"] = pd.to_numeric(
        taxonomy_df["speciateIT_confidence"],
        errors="coerce",
    ).fillna(0)

    taxonomy_df["speciateIT_decisions"] = pd.to_numeric(
        taxonomy_df["speciateIT_decisions"],
        errors="coerce",
    ).fillna(0).astype(int)

    return taxonomy_df


def read_feature_table_tsv(feature_table_tsv):
    feature_table = pd.read_csv(feature_table_tsv, sep="\t", skiprows=1)

    first_col = feature_table.columns[0]
    feature_table = feature_table.rename(columns={first_col: "Feature ID"})

    return feature_table


def create_abundance_tables(feature_table, taxonomy_df, output_folder):
    sample_columns = [col for col in feature_table.columns if col != "Feature ID"]

    feature_abundance = feature_table.copy()

    for col in sample_columns:
        feature_abundance[col] = pd.to_numeric(
            feature_abundance[col],
            errors="coerce",
        ).fillna(0)

    feature_abundance = feature_abundance.set_index("Feature ID")
    taxonomy_indexed = taxonomy_df.set_index("Feature ID")

    common_features = feature_abundance.index.intersection(taxonomy_indexed.index)

    feature_abundance = feature_abundance.loc[common_features]
    taxonomy_indexed = taxonomy_indexed.loc[common_features]

    annotated_abundance = feature_abundance.join(
        taxonomy_indexed[["Genus", "Species"]],
        how="inner",
    )

    species_count = annotated_abundance.groupby("Species")[sample_columns].sum()
    genus_count = annotated_abundance.groupby("Genus")[sample_columns].sum()

    species_relative = species_count.div(
        species_count.sum(axis=0).replace(0, np.nan),
        axis=1,
    ).fillna(0) * 100

    genus_relative = genus_count.div(
        genus_count.sum(axis=0).replace(0, np.nan),
        axis=1,
    ).fillna(0) * 100

    species_count.to_csv(output_folder / "species_count_table_species_by_samples.csv")
    species_count.T.to_csv(output_folder / "species_count_table_samples_by_species.csv")

    species_relative.to_csv(output_folder / "species_relative_abundance_species_by_samples.csv")
    species_relative.T.to_csv(output_folder / "species_relative_abundance_samples_by_species.csv")

    genus_count.to_csv(output_folder / "genus_count_table_genus_by_samples.csv")
    genus_count.T.to_csv(output_folder / "genus_count_table_samples_by_genus.csv")

    genus_relative.to_csv(output_folder / "genus_relative_abundance_genus_by_samples.csv")
    genus_relative.T.to_csv(output_folder / "genus_relative_abundance_samples_by_genus.csv")


# ============================================================
# MAIN WORKFLOW
# ============================================================

def main():
    check_inputs()

    speciateit_model = find_speciateit_model(speciateit_dir)

    if output_folder.exists() and overwrite:
        shutil.rmtree(output_folder)

    output_folder.mkdir(parents=True, exist_ok=True)

    exported_repseqs_folder = output_folder / "exported_rep_seqs"
    exported_table_folder = output_folder / "exported_table"
    speciateit_output_folder = output_folder / "speciateit_output"

    exported_repseqs_folder.mkdir(parents=True, exist_ok=True)
    exported_table_folder.mkdir(parents=True, exist_ok=True)
    speciateit_output_folder.mkdir(parents=True, exist_ok=True)

    fasta_file = exported_repseqs_folder / "dna-sequences.fasta"
    biom_file = exported_table_folder / "feature-table.biom"
    feature_table_tsv = exported_table_folder / "feature-table.tsv"

    speciateit_results = speciateit_output_folder / "MC_order7_results.txt"

    taxonomy_tsv = output_folder / "taxonomy.tsv"
    taxonomy_full_tsv = output_folder / "taxonomy_speciateIT_full.tsv"

    taxonomy_qza = output_folder / "taxonomy.qza"
    taxonomy_qzv = output_folder / "taxonomy.qzv"
    taxa_barplot_qzv = output_folder / "taxa-bar-plots.qzv"

    merged_feature_table_csv = output_folder / "merged_feature_table_speciateIT.csv"

    print("\nRepository root:")
    print(REPO_ROOT)

    print("\nMetadata file:")
    print(metadata_file)

    print("\nResults folder:")
    print(folder_results)

    print("\nDenoised folder:")
    print(denoised_folder)

    print("\nUsing speciateIT binary:")
    print(speciateit_binary)

    print("\nUsing speciateIT model:")
    print(speciateit_model)

    print("\nOutput folder:")
    print(output_folder)

    # 1. Export QIIME representative sequences to FASTA
    run_command(
        [
            "qiime",
            "tools",
            "export",
            "--input-path",
            rep_seqs_qza,
            "--output-path",
            exported_repseqs_folder,
        ]
    )

    if not fasta_file.exists():
        raise FileNotFoundError(f"Expected FASTA file was not created: {fasta_file}")

    # 2. Run speciateIT V3V4
    classify_cmd = [
        speciateit_binary,
        "-d",
        speciateit_model,
        "-i",
        fasta_file,
        "-o",
        speciateit_output_folder,
    ]

    if force_species:
        classify_cmd.append("--skip-err-thld")

    run_command(classify_cmd)

    if not speciateit_results.exists():
        raise FileNotFoundError(
            f"speciateIT did not create expected output: {speciateit_results}"
        )

    # 3. Parse speciateIT output
    taxonomy_df = parse_speciateit_results(speciateit_results)

    taxonomy_df[["Feature ID", "Taxon", "Confidence"]].to_csv(
        taxonomy_tsv,
        sep="\t",
        index=False,
    )

    taxonomy_df.to_csv(
        taxonomy_full_tsv,
        sep="\t",
        index=False,
    )

    # 4. Export QIIME feature table
    run_command(
        [
            "qiime",
            "tools",
            "export",
            "--input-path",
            table_qza,
            "--output-path",
            exported_table_folder,
        ]
    )

    if not biom_file.exists():
        raise FileNotFoundError(f"Expected BIOM file was not created: {biom_file}")

    # 5. Convert BIOM to TSV
    run_command(
        [
            "biom",
            "convert",
            "-i",
            biom_file,
            "-o",
            feature_table_tsv,
            "--to-tsv",
        ]
    )

    # 6. Merge feature table with speciateIT taxonomy
    feature_table = read_feature_table_tsv(feature_table_tsv)

    merged_feature_table = feature_table.merge(
        taxonomy_df,
        on="Feature ID",
        how="left",
    )

    merged_feature_table.to_csv(merged_feature_table_csv, index=False)

    # 7. Create species/genus abundance tables for notebook analysis
    create_abundance_tables(
        feature_table=feature_table,
        taxonomy_df=taxonomy_df,
        output_folder=output_folder,
    )

    # 8. Import speciateIT taxonomy into QIIME
    run_command(
        [
            "qiime",
            "tools",
            "import",
            "--type",
            "FeatureData[Taxonomy]",
            "--input-path",
            taxonomy_tsv,
            "--output-path",
            taxonomy_qza,
        ]
    )

    # 9. Create taxonomy visualization
    run_command(
        [
            "qiime",
            "metadata",
            "tabulate",
            "--m-input-file",
            taxonomy_qza,
            "--o-visualization",
            taxonomy_qzv,
        ]
    )

    # 10. Create taxa barplot
    if metadata_file.exists():
        run_command(
            [
                "qiime",
                "taxa",
                "barplot",
                "--i-table",
                table_qza,
                "--i-taxonomy",
                taxonomy_qza,
                "--m-metadata-file",
                metadata_file,
                "--o-visualization",
                taxa_barplot_qzv,
            ]
        )
    else:
        print("\nSkipping taxa barplot because metadata file was not found.")

    print("\nSpeciateIT refinement completed.")

    print("\nMain outputs:")
    print(f"Raw speciateIT output: {speciateit_results}")
    print(f"QIIME taxonomy TSV: {taxonomy_tsv}")
    print(f"Full speciateIT taxonomy table: {taxonomy_full_tsv}")
    print(f"Merged feature table for notebook: {merged_feature_table_csv}")

    print("\nNotebook-ready abundance tables:")
    print(f"Species count table, species x samples: {output_folder / 'species_count_table_species_by_samples.csv'}")
    print(f"Species count table, samples x species: {output_folder / 'species_count_table_samples_by_species.csv'}")
    print(f"Species relative abundance, species x samples: {output_folder / 'species_relative_abundance_species_by_samples.csv'}")
    print(f"Species relative abundance, samples x species: {output_folder / 'species_relative_abundance_samples_by_species.csv'}")
    print(f"Genus relative abundance, genus x samples: {output_folder / 'genus_relative_abundance_genus_by_samples.csv'}")
    print(f"Genus relative abundance, samples x genus: {output_folder / 'genus_relative_abundance_samples_by_genus.csv'}")

    print("\nQIIME outputs:")
    print(f"QIIME taxonomy QZA: {taxonomy_qza}")
    print(f"QIIME taxonomy QZV: {taxonomy_qzv}")
    print(f"QIIME taxa barplot: {taxa_barplot_qzv}")


if __name__ == "__main__":
    main()
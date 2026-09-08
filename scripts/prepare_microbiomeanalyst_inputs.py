#!/usr/bin/env python3
"""
Prepare MicrobiomeAnalyst input files from QIIME 2 exported outputs.

Outputs:
- asv_table_microbiomeanalyst.tsv: feature count table with #NAME as first column
- taxonomy_microbiomeanalyst.tsv: taxonomy table with #TAXONOMY as first column
- asv_id_mapping.tsv: identity map confirming original QIIME ASV feature IDs
- manifest_curated_nugent.tsv: copied sample metadata manifest
- rooted-tree.nwk: rooted phylogenetic tree extracted from QIIME 2 .qza
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import zipfile
from pathlib import Path

from biom import load_table


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BIOM = (
    PROJECT_ROOT
    / "results"
    / "08_exported_files"
    / "table.qza"
    / "7baea51d-54ab-4bcf-b314-03b787993a95"
    / "data"
    / "feature-table.biom"
)
DEFAULT_TAXONOMY = (
    PROJECT_ROOT / "results" / "10_cst_heatmap" / "refined_taxonomy_qiime.tsv"
)
DEFAULT_MANIFEST = PROJECT_ROOT / "data" / "input_files" / "manifest_curated_nugent.tsv"
DEFAULT_ROOTED_TREE = PROJECT_ROOT / "results" / "04_phylogeny" / "rooted-tree.qza"
DEFAULT_OUTPUT_DIR = (
    PROJECT_ROOT / "results" / "11_input_files_microbiomeanalyst"
)

TAXONOMY_COLUMNS = [
    "#TAXONOMY",
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
]
RANKS = TAXONOMY_COLUMNS[1:]
PREFIX_TO_RANK = {
    "d": "Kingdom",
    "k": "Kingdom",
    "p": "Phylum",
    "c": "Class",
    "o": "Order",
    "f": "Family",
    "g": "Genus",
    "s": "Species",
}
D_PREFIX_TO_RANK = {
    "D_0": "Kingdom",
    "D_1": "Phylum",
    "D_2": "Class",
    "D_3": "Order",
    "D_4": "Family",
    "D_5": "Genus",
    "D_6": "Species",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create MicrobiomeAnalyst ASV abundance and taxonomy input files."
    )
    parser.add_argument(
        "--biom",
        type=Path,
        default=DEFAULT_BIOM,
        help=f"Input BIOM table. Default: {DEFAULT_BIOM}",
    )
    parser.add_argument(
        "--taxonomy",
        type=Path,
        default=DEFAULT_TAXONOMY,
        help=f"Input QIIME taxonomy TSV. Default: {DEFAULT_TAXONOMY}",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for output TSV files. Default: {DEFAULT_OUTPUT_DIR}",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help=f"Manifest/metadata file to copy into the output directory. Default: {DEFAULT_MANIFEST}",
    )
    parser.add_argument(
        "--rooted-tree",
        type=Path,
        default=DEFAULT_ROOTED_TREE,
        help=f"Input QIIME rooted tree artifact. Default: {DEFAULT_ROOTED_TREE}",
    )
    return parser.parse_args()


def format_count(value: float) -> str:
    value = float(value)
    if math.isclose(value, round(value), rel_tol=0, abs_tol=1e-9):
        return str(int(round(value)))
    return f"{value:g}"


def read_taxonomy(taxonomy_path: Path) -> dict[str, str]:
    with taxonomy_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Feature ID", "Taxon"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            missing_list = ", ".join(sorted(missing))
            raise ValueError(f"Taxonomy file is missing required column(s): {missing_list}")

        taxonomy_by_feature = {}
        for row in reader:
            feature_id = row["Feature ID"].strip()
            taxonomy_by_feature[feature_id] = row["Taxon"].strip()

    return taxonomy_by_feature


def parse_taxon_string(taxon: str) -> dict[str, str]:
    parsed = {rank: "" for rank in RANKS}

    for position, raw_part in enumerate(taxon.split(";")):
        part = raw_part.strip()
        if not part:
            continue

        rank = RANKS[position] if position < len(RANKS) else None
        value = part

        if "__" in part:
            prefix, value = part.split("__", 1)
            prefix = prefix.strip()
            value = value.strip()
            rank = PREFIX_TO_RANK.get(prefix.lower(), D_PREFIX_TO_RANK.get(prefix))
        else:
            value = value.strip()

        if rank in parsed:
            parsed[rank] = value

    return parsed


def write_abundance_table(
    table,
    feature_ids: list[str],
    sample_ids: list[str],
    output_path: Path,
) -> None:
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["#NAME", *sample_ids])

        for feature_id in feature_ids:
            counts = table.data(feature_id, axis="observation")
            writer.writerow([feature_id, *[format_count(value) for value in counts]])


def write_taxonomy_table(
    feature_ids: list[str],
    taxonomy_by_feature: dict[str, str],
    output_path: Path,
) -> None:
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(TAXONOMY_COLUMNS)

        for feature_id in feature_ids:
            taxon = taxonomy_by_feature.get(feature_id, "")
            parsed = parse_taxon_string(taxon)
            writer.writerow([feature_id, *[parsed[rank] for rank in RANKS]])


def write_id_mapping(
    feature_ids: list[str],
    output_path: Path,
) -> None:
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Feature ID", "Original Feature ID"])
        for feature_id in feature_ids:
            writer.writerow([feature_id, feature_id])


def extract_tree_nwk(tree_qza_path: Path, output_path: Path) -> None:
    with zipfile.ZipFile(tree_qza_path) as archive:
        tree_members = [
            member
            for member in archive.namelist()
            if member.endswith("/data/tree.nwk")
        ]

        if not tree_members:
            raise ValueError(f"No data/tree.nwk file found in artifact: {tree_qza_path}")
        if len(tree_members) > 1:
            raise ValueError(
                f"Expected one data/tree.nwk file in artifact, found {len(tree_members)}"
            )

        tree_text = archive.read(tree_members[0]).decode("utf-8").strip()

    if not tree_text.endswith(";"):
        raise ValueError(f"Extracted tree is not valid Newick; missing trailing ';': {tree_qza_path}")

    output_path.write_text(f"{tree_text}\n")


def main() -> None:
    args = parse_args()
    biom_path = args.biom.resolve()
    taxonomy_path = args.taxonomy.resolve()
    manifest_path = args.manifest.resolve()
    rooted_tree_path = args.rooted_tree.resolve()
    output_dir = args.output_dir.resolve()

    if not biom_path.exists():
        raise FileNotFoundError(f"BIOM file not found: {biom_path}")
    if not taxonomy_path.exists():
        raise FileNotFoundError(f"Taxonomy file not found: {taxonomy_path}")
    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest file not found: {manifest_path}")
    if not rooted_tree_path.exists():
        raise FileNotFoundError(f"Rooted tree artifact not found: {rooted_tree_path}")

    output_dir.mkdir(parents=True, exist_ok=True)

    table = load_table(str(biom_path))
    feature_ids = [str(feature_id) for feature_id in table.ids(axis="observation")]
    sample_ids = [str(sample_id) for sample_id in table.ids(axis="sample")]
    taxonomy_by_feature = read_taxonomy(taxonomy_path)

    missing_taxonomy = sorted(set(feature_ids) - set(taxonomy_by_feature))
    if missing_taxonomy:
        raise ValueError(
            f"{len(missing_taxonomy)} BIOM feature(s) are missing from taxonomy file. "
            f"First missing feature: {missing_taxonomy[0]}"
        )

    abundance_path = output_dir / "asv_table_microbiomeanalyst.tsv"
    taxonomy_output_path = output_dir / "taxonomy_microbiomeanalyst.tsv"
    mapping_path = output_dir / "asv_id_mapping.tsv"
    manifest_output_path = output_dir / manifest_path.name
    rooted_tree_output_path = output_dir / "rooted-tree.nwk"

    write_abundance_table(table, feature_ids, sample_ids, abundance_path)
    write_taxonomy_table(feature_ids, taxonomy_by_feature, taxonomy_output_path)
    write_id_mapping(feature_ids, mapping_path)
    shutil.copy2(manifest_path, manifest_output_path)
    extract_tree_nwk(rooted_tree_path, rooted_tree_output_path)

    print(f"Wrote ASV abundance table: {abundance_path}")
    print(f"Wrote taxonomy table: {taxonomy_output_path}")
    print(f"Wrote ASV ID identity map: {mapping_path}")
    print(f"Copied manifest: {manifest_output_path}")
    print(f"Extracted rooted tree: {rooted_tree_output_path}")
    print(f"Features: {len(feature_ids)}")
    print(f"Samples: {len(sample_ids)}")


if __name__ == "__main__":
    main()

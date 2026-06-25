#Import the libraries

import os, subprocess

#Define the paths

folder_results="../results"
metadata_file="../data/input_files/metadata.tsv"


def create_metadata_tabulate(metadata_file, folder_results):
    cmd=['qiime', 'metadata', 'tabulate',
         '--m-input-file', f'{metadata_file}',
         '--o-visualization', f'{folder_results}/sample-metadata-viz.qzv']
    subprocess.run(cmd, check=True)
    
#create_metadata_tabulate(metadata_file, folder_results)    


def import_data(metadata_file, folder_results):
    cmd=['qiime', 'tools', 'import',
         '--type', 'SampleData[PairedEndSequencesWithQuality]',
         '--input-path', f'{metadata_file}',
         '--output-path', f'{folder_results}/demux-paired-end.qza',
         '--input-format', 'PairedEndFastqManifestPhred33V2']
    subprocess.run(cmd, check=True)

import_data(metadata_file, folder_results)
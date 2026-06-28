#Import the libraries

import os
import subprocess

#Define the paths

metadata_file='../data/input_files/manifest_curated_nugent.tsv'
folder_results='../results'
os.makedirs(folder_results, exist_ok=True)

def create_metadata_tabulate(metadata_file, folder_results):
    cmd = ['qiime', 'metadata', 'tabulate',
           '--m-input-file', f'{metadata_file}',
           '--o-visualization', f'{folder_results}/sample-metadata-viz.qzv']
    subprocess.run(cmd, check=True)

create_metadata_tabulate(metadata_file, folder_results)


def import_data(metadata_file, folder_results):
    cmd = ['qiime', 'tools', 'import',
           '--type', 'SampleData[PairedEndSequencesWithQuality]',
           '--input-path', f'{metadata_file}',
           '--output-path', f'{folder_results}/demux-paired-end.qza',
           '--input-format', 'PairedEndFastqManifestPhred33V2']
    subprocess.run(cmd, check=True)

import_data(metadata_file, folder_results)


def demux_summarize(folder_results):
    cmd = ['qiime', 'demux', 'summarize',
           '--i-data', f'{folder_results}/demux-paired-end.qza',
           '--o-visualization', f'{folder_results}/demux-paired-end.qzv']
    subprocess.run(cmd, check=True)

demux_summarize(folder_results)

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
       
denoise_paired(folder_results)       


def make_denoising_visualization(folder_results):
       cmd=['qiime', 'metadata', 'tabulate',
            '--m-input-file', f'{folder_results}/denoised-paired_results/denoising-stats.qza',
            '--o-visualization',f'{folder_results}/denoised-paired_results/denoising-stats.qzv']
       subprocess.run(cmd, check=True)

make_denoising_visualization(folder_results)

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
       
featuretable_summaries(denoised_folder, folder_results, metadata_file)       


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
       
phylogenetic_analysis(denoised_folder, folder_results)      

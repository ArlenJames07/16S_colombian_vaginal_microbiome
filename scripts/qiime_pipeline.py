#Import the libraries

import os
import subprocess

#Define the paths

metadata_file='../data/input_files/manifest.tsv'
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


def featuretable_summaries(folder_results):
       os.makedirs(f'{folder_results}/feature-table_summarise)


qiime feature-table summarize \
  --i-table table.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-summary table.qzv \
  --o-feature-frequencies feature-frequencies.qza \
  --o-sample-frequencies sample-frequencies.qza
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv       


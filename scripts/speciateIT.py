import os, subprocess


# Extract qza information 


denoised_folder='../results/denoised-paired_results'
results_folder='../results'

def extract(denoised_folder, results_folder):
    output_folder=f'{results_folder}/qza_extracted'
    os.makedirs(output_folder, exist_ok=True)
    for x in os.listdir(denoised_folder):
        if x.endswith('rep-seqs.qza') or x.endswith('table.qza'):
            file=os.path.join(denoised_folder,x)
            cmd=['qiime', 'tools', 'extract',
                 '--input-path', f'{file}',
                 '--output-path',f'{output_folder}/{x}']
            subprocess.run(cmd, check=True)
            
#extract(denoised_folder, results_folder)

qiime_taxonomy='../results/taxonomic_assignment'

def extract_taxonomy(qiime_taxonomy,results_folder):
    output_folder=f'{results_folder}/taxonomy_qiime'
    os.makedirs(output_folder, exist_ok=True)
    for root, dirs, files in os.walk(qiime_taxonomy):
        for x in files:
            if x.endswith('taxonomy.qzv'):
                file=os.path.join(root,x)
                basename=file.split('/')[3]
                cmd=['qiime', 'tools', 'extract',
                     '--input-path', f'{file}',
                     '--output-path',f'{output_folder}/{basename}']
                subprocess.run(cmd, check=True)
#extract_taxonomy(qiime_taxonomy,results_folder)

## Run speciate IT to assign taxonomic classification 

fasta_sequences='../results/qza_extracted/rep-seqs.qza/55afecd2-dd83-4720-a290-4e888c827966/data/dna-sequences.fasta'
program='/home/rare/programs/speciateIT/vSpeciateDB_models/vSpeciateIT_V3V4'


def classification(fasta_sequences, results_folder, program):
    output_folder=f'{results_folder}/speciate_it_classification'
    os.makedirs(output_folder, exist_ok=True)
    cmd=['classify', 
         '-d', f'{program}', 
         '-i', f'{fasta_sequences}',
         '-o', f'{output_folder}']
    subprocess.run(cmd, check=True)

#classification(fasta_sequences, results_folder, program)   
#!/software/Anaconda3-4.3.0-el7-x86_64/bin/python3

import os
import shutil
import argparse

install_dir = "/project2/onibasu/dropseqRunner_NEW"
work_dir = os.getcwd()
protocol_map = {'drop': 'Snakefile_drop.smk', '10x-v2': 'Snakefile_10x.smk', '10x-v3': 'Snakefile_10x.smk'}

REF_DICT_DROP = {'hg38': os.path.join(install_dir, 'STAR_indeces/hg38_noalt_juncGencodeV27_STAR_2.7.1_61'), 
		'cer3': os.path.join(install_dir,'STAR_indeces/sacCer3_s288c'),
		'candida': os.path.join(install_dir, 'STAR_indeces/C_albicans_SC5314'), 
		'mm10': os.path.join(install_dir, 'STAR_indeces/mm10_noaltn_juncGencodemV15_61'), 
		'danRer11': os.path.join(install_dir, 'STAR_indeces/danRer11_primary_61')}
		
REF_DICT_10X = {'hg38': os.path.join(install_dir, 'STAR_indeces/hg38_noalt_juncGencodeV27_STAR_2.7.1_99')}

def check_gzip(files):
    '''
    Checks if all given files are in gzip format
    '''
    for f in files:
        if f.split('.')[-1] != 'gz':
            return False
    return True

def make_config(args, install_dir, work_dir):
  
    if args.protocol == '10x-v2':
        barcodes = os.path.join(install_dir, 'barcodes/whitelist.v2.txt.gz')
        
        try:
          ref = REF_DICT_10X[args.ref]
        except:
          assert False, f"{args.ref} is not available for {args.protocol}"
          
    elif args.protocol == '10x-v3':
      
        barcodes = os.path.join(install_dir, 'barcodes/whitelist.v3.txt.gz')
        try:
          ref = REF_DICT_10X[args.ref]
        except:
          assert False, f"{args.ref} is not available for {args.protocol}"
          
    else:
        barcodes = None
        ref = REF_DICT_DROP[args.reff]

    config=f"""proj_dir: {work_dir}/
genome_index: {ref}/
refFlat: {os.path.join(ref, "refFlat_for_picard.refFlat")} 
scripts: {os.path.join(install_dir, "Scripts/")}
cell_num: 10000
barcode: "CCCCCCCCCCCCNNNNNNNN"
dir_log: log/
Sample: {args.sample}
Protocol: {args.protocol}
10x_barcodes: {barcodes}
"""
    
    with open('config.yaml', 'w') as f:
         f.write(config)

def make_submit_snakemake(args, install_dir, work_dir):
    '''
    Creates the SLURM that submits Snakefile on the cluster
    '''
    smk = protocol_map[args.protocol]

    cmd =f"""#!/bin/bash

#SBATCH --job-name=snakemake
#SBATCH --output=snakelog.out
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=4G
#SBATCH --tasks-per-node=4

source /project2/onibasu/dropseqRunner_NEW/miniconda3/bin/activate dropRunner
snakemake \\
    -kp \\
    --ri \\
    -j 150 \\
    --latency-wait 150 \\
    --cluster-config {install_dir}/cluster_solo.json \\
    -c "sbatch \\
        --mem={{cluster.mem}} \\
        --nodes={{cluster.n}} \\
        --tasks-per-node={{cluster.tasks}} \\
        --partition={{cluster.partition}} \\
        --account={{cluster.account}} \\
        --job-name={{cluster.name}} \\
        --output={{cluster.logfile}}" \\
    -s {os.path.join(install_dir, smk)} \\
    --configfile {os.path.join(work_dir, "config.yaml")}
"""

    with open('submit_snakemake.sbatch', 'w') as f:
        f.write(cmd)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    
    # Required
    parser.add_argument('--R1', type=str, help='Absolute path to gzipped read 1 fastq file (comma-delimited list of files if multiple.) REQUIRED')
    parser.add_argument('--R2', type=str, help='Absolute path to gzipped read 2 fastq file (comma-delimited list of files if multiple.) REQUIRED')
    parser.add_argument('--ref', type=str, help='Name of reference genome. Options: hg38, cer3, candida, mm10, danRer11. REQUIRED')
    parser.add_argument('--protocol', type=str, help='Single cell protocol used to produce data. Options: drop, 10x-v2, 10x-v3')
    # Optional/conditional
    parser.add_argument('--rerun', action='store_true', help='This flag re-runs a previously failed attempt.')
    parser.add_argument('--cluster', action='store_true', help='Provide this flag if this job should be run on the cluster.')
    parser.add_argument('--sample', type=str, help='sample name. Optional.')
    
    args = parser.parse_args()
    
    if args.rerun:
         assert os.path.isfile('submit_snakemake.sbatch'), \
            'sbatch file not found. Are you sure you ran this pipeline before?'
         os.system('sbatch submit_snakemake.sbatch')   
    
    if args.R1 == None or args.R2 == None or args.ref == None or args.protocol == None:
         raise Exception('Required arguments not provided. Please provide R1, R2, an indices folder name, and a protocol')
    
    assert os.path.isfile(args.R1), "Please provide a gzipped fastq file for read 1"
    assert os.path.isfile(args.R2), "Please provide a gzipped fastq files for read 2."
    assert args.protocol in ['drop', '10x-v2', '10x-v3'], 'Please provide a valid protocol! Options are drop, 10x-v2, or 10x-v3'
    assert args.ref in ['hg38','cer3','candida','mm10','danRer11'], 'Please provide a valid reference genome name! The options are hg38 and cer3, OR provide a STAT index directory.'
    
    if args.cluster == None:
        args.cluster = False
    if args.sample == None:
        args.sample = f'{args.protocol}-experiment'

    r1, r2 = args.R1.split(','), args.R2.split(',') 
    assert len(r1) == len(r2), \
    'Number of files in read 1 and read 2 are not the same. Please provide a read 1 and read 2 file for each experiment.'
    
    os.system('mkdir fastq')

    if check_gzip(r1) and check_gzip(r2):
        for R1,R2 in zip(r1, r2):
            f1_name = R1.split('/')[-1]
            f2_name = R2.split('/')[-1]
            
            assert 'R1' in f1_name, 'R1 filename does not contain "R1". Did you give R2 twice?'
            assert 'R2' in f2_name, 'R2 filename does not contain "R2". Did you give R1 twice?'

            if os.path.isabs(R1):
              os.system(f'ln -s {R1} fastq/{f1_name}')
            else:
              os.system(f'ln -s ../{R1} fastq/{f1_name}')
            if os.path.isabs(R2):
              os.system(f'ln -s {R2} fastq/{f2_name}')
            else:
              os.system(f'ln -s ../{R2} fastq/{f2_name}')
    else:
        msg = 'File format not recognized. Please make sure you provide gzipped fastq files (files should end with .fastq.gz)'
        raise TypeError(msg)

    make_config(args, install_dir, work_dir)
    
    if args.cluster:
        make_submit_snakemake(args, install_dir, work_dir)
        print('Submitting snakemake job..')
        os.system('sbatch submit_snakemake.sbatch')
        print('Snakemake job has been submitted to the cluster.\nType "qstat -u CNETID" to see the progress of the snakefile.')
        
    else: 
        smk = protocol_map[args.protocol]
        print('Running snakemake directly on this node. This may not finish because alignment requires >30GB of RAM.')
        os.system(f'source /project2/onibasu/dropseqRunner_NEW/miniconda3/bin/activate dropRunner; snakemake -kp --ri -s {install_dir}/{smk} --configfile {work_dir}/config.yaml')




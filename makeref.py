#!/software/Anaconda3-4.3.0-el7-x86_64/bin/python3

import os
import shutil
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=str, help='Genome fasta file in gzip format.')
    parser.add_argument('--gtf', type=str, help='Genome annotations in GTF format.')
    parser.add_argument('--readLen', type=int, help='Read length used for alignment to this genome. Default: 62')
    parser.add_argument('--outDir', type=str, help='Output directory of of genome index')
    parser.add_argument('--genomeNbases', type=int, help='A STAR argument that needs to be adjusted for smaller (non-human) genomes. Default: 14')
    parser.add_argument('--contigsParam', type=int, help='For a genome with larger number of contigs, it is recommended to scale this parameter. Default: 18')
    parser.add_argument('--cluster', action='store_true', help='Add this flag if the job should run on the cluster.')
    args = parser.parse_args()
   
    if args.fasta == None or args.gtf == None or args.outDir == None:
         raise Exception('Required arguments not provided. Please provide a fasta file, annotations, and an output directory name.')    

    assert os.path.isfile(args.fasta), 'Please provide a valid fasta file.' 
    assert os.path.isfile(args.gtf), 'Please provide a valid gtf file.' 
    assert not os.path.isdir(args.outDir), \
    'Cannot create output directory because it already exists. Please provide the location and name of a non-existing directory.'
    
    if args.genomeNbases == None:
        args.genomeNbases = 14
    if args.contigsParam == None:
        args.contigsParam = 18
    if args.readLen == None:
        args.readLen = 62
    overhang = args.readLen - 1
     
    print('Setting up directory and creating auxililary files..')
    
    os.system(f'mkdir {args.outDir}')
    
    os.system(f'source /project2/onibasu/dropseqRunner_NEW/miniconda3/bin/activate dropRunner; gtfToGenePred {args.gtf} tmp -genePredExt')
    cmd = f"""awk '{{print $12"\t"$0}}' tmp | cut -f1-11 > {args.outDir}/refFlat_for_picard.refFlat; rm tmp"""
    os.system(cmd)
     
    if os.stat(f'{args.outDir}/refFlat_for_picard.refFlat').st_size == 0:
        msg = 'WARNING: GTF to refFlat may have failed; you can safely ignore this if your GTF file doesnt have intron/exon information'
        print(msg)
    
    if args.cluster:
      
        cmd =f"""#!/bin/bash

#SBATCH --job-name=genome_index
#SBATCH --output=genome_index.log
#SBATCH --error=genome_index.error
#SBATCH --time=10:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=50G
#SBATCH --tasks-per-node=4

source /project2/onibasu/dropseqRunner_NEW/miniconda3/bin/activate dropRunner
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir {args.outDir} --genomeFastaFiles {args.fasta} --sjdbGTFfile {args.gtf} --sjdbOverhang {overhang} --genomeSAindexNbases {args.genomeNbases} --genomeChrBinNbits {args.contigsParam}
"""

        with open('generate_indices.sbatch', 'w') as f:
            f.write(cmd)
        
        os.system('sbatch generate_indices.sbatch')
        print('Index generation has been submitted to the cluster. Type qstat -u CNETID to check the status.')
    
    else:
      
      print('Genome index generation will run locally on this machine. This may not complete due to STARs large memory requirement.')
      os.system(f'source /project2/onibasu/dropseqRunner_NEW/miniconda3/bin/activate dropRunner; STAR --runThreadN 1 --runMode genomeGenerate --genomeDir {args.outDir}/ --genomeFastaFiles {args.fasta} --sjdbGTFfile {args.gtf} --sjdbOverhang {overhang} --genomeSAindexNbases {args.genomeNbases} --genomeChrBinNbits {args.contigsParam}')

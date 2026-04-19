#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ISTSD-Profiler: A parallelized pipeline for high-throughput characterization
of insertion sequence target-site duplications across microbial genomes.

Dependencies:
    Python packages: mpi4py, biopython, pandas, rpy2
    External tools: samtools, bedtools, seqkit
    R packages: ggplot2, ggseqlogo

Usage:
    mpiexec -n <num_processes> python3 ISTSD-Profiler.py \\
        -f <IS_elements.fasta> -g <genome.fasta> -n <num_splits>
"""

__version__ = "1.0.0"
__author__ = "Yong Sheng"
__email__ = "johnsheng.sjtu@vip.163.com"
__license__ = "MIT"

import pandas as pd
import datetime
import argparse
from mpi4py import MPI
from Bio.Seq import Seq
import math
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import os
import shutil
import pickle
import time
import subprocess
from rpy2.robjects import r
import gc
pandas2ri.activate() 
parser = argparse.ArgumentParser(
    description="ISTSD-Profiler: High-throughput characterization of IS element "
                "target-site duplications across microbial genomes.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="Example:\n"
           "  mpiexec -n 14 python3 ISTSD-Profiler.py \\\n"
           "    -f IS_elements.fasta -g genome.fasta -n 22"
)
parser.add_argument('-f', '--IS', metavar='string', help='Specify the FASTA file containing IS elements.')
parser.add_argument('-g', '--genome', metavar='string', help='Specify the FASTA file containing genome sequences for analysis.')
parser.add_argument('-n', type=int, metavar='int', help='Split genome sequences into N parts. Note: N should be less than or equal to the total number of sequence entries in the genome sequences.')
args = parser.parse_args()
IS_file = args.IS
genome_file = args.genome
num_parts = args.n
comm = MPI.COMM_WORLD 
rank = comm.Get_rank() 
size = comm.Get_size() 
tsd_data = {} 
min_tsd_length = 2
max_tsd_length = 14
def create_directory(directory):
    """Create a directory if it does not already exist."""
    if not os.path.exists(directory):
        os.makedirs(directory)


def serialize_and_store_data(obj, name): 
    """Serialize and store a Python object to disk using pickle."""
    folder_name = "pickle_data"
    os.makedirs(folder_name, exist_ok=True) 
    file_name = f"{folder_name}/{name}.pkl"
    try:
        with open(file_name, "wb") as f: 
            pickle.dump(obj, f) 
        print(f"[INFO] {name} was stored successfully in the folder {folder_name}.")
    except OSError:
        print(f"[ERROR] Failed to store {name} in the folder {folder_name}")


def generate_fai(fasta_file):
    """Generate FASTA index (.fai) and return a dictionary of chromosome lengths."""
    fai_file = fasta_file + ".fai"
    if not os.path.exists(fai_file):
        start_time = time.time()
        cmd = f"samtools faidx {fasta_file}" 
        subprocess.run(cmd, shell=True, check=True)
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"[INFO] Building index for {fasta_file} took {execution_time} seconds.")
    seqkit_fai_file = fasta_file + ".seqkit.fai" 
    shutil.copy(fai_file, seqkit_fai_file) 
    chromosome_lengths = {} 
    with open(fai_file, "r") as f:
        for line in f: 
            chromosome, length = line.strip().split("\t")[0:2] 
            chromosome_lengths[chromosome] = int(length)
    return chromosome_lengths


def chunk_list(input_list, size):
    """Split a list into approximately equal-sized chunks for parallel processing."""
    num_items = len(input_list)
    chunk_size = math.ceil(num_items / size) 
    extra_chunks = num_items % chunk_size 
    if extra_chunks > 0:
        chunks = [input_list[i*chunk_size:(i+1)*chunk_size] for i in range(size-1)] 
        chunks += [input_list[(size-1)*chunk_size:]] 
    else:
        chunks = [input_list[i*chunk_size:(i+1)*chunk_size] for i in range(size)]
    return chunks


def extract_is_information(file_name):
    """Parse IS element alignment information from a seqkit locate GTF output file."""
    is_info_list = []
    with open(file_name) as file:
        for line in file:
            data = line.split("\t")
            chromosome = data[0]
            is_start = int(data[3])
            is_end = int(data[4])
            direction = data[6]
            is_name = data[8].rstrip().replace('gene_id "', '').replace('";', '')
            if direction == '-':
                direction = -1
            elif direction == '+':
                direction = 1
            is_info_list.append((chromosome, is_start, is_end, is_name, direction))
    return is_info_list


def extract_sequence(chromosome, start, end, fasta_file): 
    """Extract a DNA sequence from a FASTA file using bedtools getfasta (1-based coordinates)."""
    start = int(start) 
    end = int(end)
    cmd = f"echo '{chromosome}\t{start-1}\t{end}' | bedtools getfasta -fi {fasta_file} -bed -"
    result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    lines = result.stdout.strip().split("\n")[1:] 
    sequence = "".join(lines).replace("\n", "").upper()
    return sequence


def extract_tsd_sequence(chunk, fasta_file, chromosome_lengths):
    """Identify target-site duplications (TSDs) flanking IS elements.
    TSD sequences are orientation-normalized: for IS elements on the minus
    strand, the TSD is reverse-complemented to reflect the transposase
    perspective.
    """
    if not chunk: 
        return {}, [] 
    result = ({}, []) 
    for is_info in chunk: 
        chromosome, is_start, is_end, is_name, direction = is_info
        chromosome_length = chromosome_lengths[chromosome]
        found_tsd = False 
        for i in range(min_tsd_length, max_tsd_length + 1):
            L_tsd_start = is_start - i 
            L_tsd_end = is_start - 1 
            R_tsd_start = is_end + 1
            R_tsd_end = is_end + i 
            if L_tsd_start >= 1 and R_tsd_end <= chromosome_length:
                left_tsd_sequence = extract_sequence(chromosome, L_tsd_start, L_tsd_end, fasta_file)
                right_tsd_sequence = extract_sequence(chromosome, R_tsd_start, R_tsd_end, fasta_file) 
                if left_tsd_sequence == right_tsd_sequence: 
                    found_tsd = True
                    if i not in result[0]:
                        result[0][i] = []
                    if direction == -1:
                        tsd_seq_to_store = reverse_complement(left_tsd_sequence)
                    else:
                        tsd_seq_to_store = left_tsd_sequence
                    result[0][i].append((is_name, tsd_seq_to_store, is_info))
        if not found_tsd:
            result[1].append((is_name, is_info))
    return result[0], result[1] 


def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    seq = Seq(sequence) 
    reverse_complement_sequence = str(seq.reverse_complement())
    return reverse_complement_sequence


def split_fasta_file(fasta_file, num_files):
    """Split a FASTA file into multiple parts using seqkit split."""
    temp_folder = fasta_file + ".split"
    if os.path.exists(temp_folder):
        print(f"[INFO] {temp_folder} already exists. Skipping the splitting procedure.")
    else:
        os.makedirs(temp_folder, exist_ok=True)
        seqkit_cmd = f"seqkit split -p {num_files} {fasta_file} -O {temp_folder}" 
        subprocess.run(seqkit_cmd, shell=True, check=True)
    split_files = [os.path.join(temp_folder, f) for f in os.listdir(temp_folder) if os.path.isfile(os.path.join(temp_folder, f))]
    return temp_folder, split_files


def seqkit_locate(is_fa_file, fasta_file): 
    """Run seqkit locate to map IS elements onto genome sequences."""
    is_info_file_individual = f"{fasta_file}.is_info.txt"
    seqkit_cmd = f"seqkit locate --gtf -i -F -m 0 -f {is_fa_file} {fasta_file} -o {is_info_file_individual}" 
    subprocess.run(seqkit_cmd, shell=True, check=True)


def merge_files(output_file, input_files):
    """Merge multiple text files into a single output file."""
    with open(output_file, "a") as outfile: 
        for input_file in input_files: 
            with open(input_file, "r") as infile:
                outfile.write(infile.read())


def save_tsd_sequences(tsd_sequences, output_directory):
    """Save TSD sequences to FASTA files, separated by uniqueness at each locus."""
    global tsd_data 
    position_counts = {}
    for tsd_length, sequences in tsd_sequences.items():
        for sequence_info in sequences:
            is_name, sequence, is_info = sequence_info
            chromosome = is_info[0]
            is_start = is_info[1]
            is_end = is_info[2]
            position = (chromosome, is_start, is_end) 
            if position in position_counts:
                position_counts[position].append((tsd_length, is_name, sequence))
            else:
                position_counts[position] = [(tsd_length, is_name, sequence)]
    for position, sequences in position_counts.items():
        if len(sequences) > 1: 
            output_directory_duplicated = f"{output_directory}/ISs_multiple_TSD_lengths_same_position"
            os.makedirs(output_directory_duplicated, exist_ok=True) 
            output_file_duplicated = f"{output_directory_duplicated}/{sequences[0][1]}_TSD_duplicated.fasta"  
            chromosome, pos_start, pos_end = position
            with open(output_file_duplicated, "a") as file: 
                for sequence_info in sequences:
                    tsd_length, is_name, sequence = sequence_info
                    tsd_length_bp = f"{tsd_length}bp"  
                    sequence_name = f"{is_name}_chr{chromosome}_pos{pos_start}_{pos_end}_{tsd_length_bp}TSD_duplicated"  
                    file.write(f">{sequence_name}\n{sequence}\n")   
        elif len(sequences) == 1:
            output_directory_single = f"{output_directory}/ISs_single_TSD_length_same_position"
            os.makedirs(output_directory_single, exist_ok=True)
            output_file_single = f"{output_directory_single}/{sequences[0][1]}_{sequences[0][0]}bp_TSD.fasta" 
            chromosome, pos_start, pos_end = position
            with open(output_file_single, "a") as file:
                for sequence_info in sequences:
                    tsd_length, is_name, sequence = sequence_info
                    tsd_length_bp = f"{tsd_length}bp"  
                    sequence_name = f"{is_name}_chr{chromosome}_pos{pos_start}_{pos_end}_{tsd_length_bp}TSD"
                    file.write(f">{sequence_name}\n{sequence}\n")
            for sequence_info in sequences:
                tsd_length, is_name, sequence = sequence_info
                category = f"{is_name}-{tsd_length} bp" 
                if category not in tsd_data:
                    tsd_data[category] = []
                tsd_data[category].append(sequence) 
    #with open("IS_TSD_list.txt", "a") as file:
        #for category, tsd_sequences in tsd_data.items():
            #sequences_string = " ".join([f"'{sequence}'" for sequence in tsd_sequences]) 
            #file.write(f"{category} = ({sequences_string})\n")


def analyze_tsd_length_preference(tsd_data):
    """Analyze TSD length preferences for IS elements producing multiple TSD lengths."""
    is_tsd_preference = {}
    for category, sequences in tsd_data.items():
        is_name, tsd_length = category.rsplit("-", 1) 
        tsd_length = int(tsd_length.replace(" bp", "")) 
        if is_name not in is_tsd_preference:
            is_tsd_preference[is_name] = {} 
        if tsd_length not in is_tsd_preference[is_name]:
            is_tsd_preference[is_name][tsd_length] = 0 
        is_tsd_preference[is_name][tsd_length] += len(sequences) 
    is_tsd_preference_filtered = {}
    for is_name, tsd_lengths in is_tsd_preference.items(): 
        if len(tsd_lengths) > 1: 
            is_tsd_preference_filtered[is_name] = tsd_lengths
    df = pd.DataFrame([(is_name, tsd_length, count) for is_name, tsd_lengths in is_tsd_preference_filtered.items() for tsd_length, count in tsd_lengths.items()],
                      columns=['IS', 'TSD Length', 'Sequence Count'])
    subfolder_name = 'analyze_tsd_length_preference'
    if not os.path.exists(subfolder_name):
        os.makedirs(subfolder_name)
    file_path = os.path.join(subfolder_name, 'tsd_length_preference_data.csv')
    df.to_csv(file_path, index=False)
    del is_tsd_preference
    del is_tsd_preference_filtered


def plot_tsd_motif(category, folder):
    """Generate a TSD sequence logo plot using R ggseqlogo."""
    r(f'''
    library(ggplot2)
    library(ggseqlogo)
    p <- ggplot() + geom_logo(r_is_tsd$"{category}") + theme_logo()
    ggsave(file.path("{folder}", "{category}_TSD_motif.pdf"), p)
    ''')


def analyze_is_elements(IS_file, genome_file, num_files):
    """Main pipeline: align IS elements, detect TSDs, and generate motif analyses."""
    if rank == 0:
        print(f"[INFO] Start time: {datetime.datetime.now()}")
        print(f"[INFO] Please ensure that the sequences in the {IS_file} are unique before running.")
        print(f"[INFO] Current number of processes: {size}")
        chromosome_lengths = generate_fai(genome_file)
        serialize_and_store_data(chromosome_lengths, "chromosome_lengths")
    else:
        chromosome_lengths = None 
    if rank == 0:    
        temp_folder, split_files = split_fasta_file(genome_file, num_files)
        local_files = chunk_list(split_files, size)
    else:
        local_files = None
    comm.barrier()
    chromosome_lengths = comm.bcast(chromosome_lengths, root=0)
    local_files = comm.bcast(local_files, root=0)
    local_files = local_files[rank] 
    comm.barrier()
    if rank == 0:
        start_time = time.time()
        print("[INFO] Aligning all IS elements to the genome sequences...")
    for file in local_files:
        seqkit_locate(IS_file, file)
    del local_files 
    gc.collect() 
    comm.barrier()
    if rank == 0:
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"[INFO] Aligning all IS elements to the genome sequences took {execution_time} seconds.")
        is_info_files = [f"{file}.is_info.txt" for file in split_files] 
        del split_files
        gc.collect()
        merge_files("merged_IS_alignment_info.txt", is_info_files)
        for file in is_info_files:
            os.remove(file)
        del is_info_files
        gc.collect()
        if os.path.getsize("merged_IS_alignment_info.txt") == 0:
            print("[INFO] IS elements were not detected in any genome. Skipping further analysis.")
            comm.Abort() 
        print("[INFO] Extracting ISs alignment information...")
        is_info_list = extract_is_information("merged_IS_alignment_info.txt")
        start_time = time.time()
        print("[INFO] Analyzing TSD sequences...")
        local_files_info = chunk_list(is_info_list, size)
        del is_info_list
        gc.collect()
    else:
        local_files_info = None
    comm.barrier()
    local_files_info = comm.bcast(local_files_info, root=0)
    local_files_info = local_files_info[rank]
    single_tsd_sequences, single_is_without_tsd = extract_tsd_sequence(local_files_info, genome_file, chromosome_lengths)
    del local_files_info
    del chromosome_lengths
    gc.collect()
    comm.barrier()
    gathered_tsd_sequences = comm.gather(single_tsd_sequences, root=0)
    del single_tsd_sequences
    gc.collect()
    comm.barrier() 
    if rank == 0:
        tsd_sequences = {}
        for rank_data in gathered_tsd_sequences:
            for key, value in rank_data.items():
                if key not in tsd_sequences:
                    tsd_sequences[key] = value
                else:
                    tsd_sequences[key].extend(value) 
        del gathered_tsd_sequences
        gc.collect()
        serialize_and_store_data(tsd_sequences, "tsd_sequences")
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"[INFO] Analyzing TSD sequences took {execution_time} seconds.")
        create_directory("TSD_analysis_results")
        save_tsd_sequences(tsd_sequences, "TSD_analysis_results") 
        del tsd_sequences
        gc.collect()
        serialize_and_store_data(tsd_data, "tsd_data") 
        r_is_tsd = robjects.ListVector({}) 
        TSD_motifs_folder = "TSD_motifs"
        os.makedirs(TSD_motifs_folder, exist_ok=True)
        for category, tsd_sequences in tsd_data.items():
            r_is_tsd.rx2[category] = robjects.StrVector(tsd_sequences) 
        r.assign('r_is_tsd', r_is_tsd) 
        print("[INFO] Plotting TSD motif logos...")
        for category in tsd_data.keys():
            plot_tsd_motif(category, TSD_motifs_folder)
        print("[INFO] TSD motif logo plots complete.")
        print("[INFO] Performing TSD length preference analysis...")
        analyze_tsd_length_preference(tsd_data) 
        print("[INFO] TSD length preference analysis complete.")
        print("[INFO] Congratulations! The code execution has been completed successfully.")
        print(f"[INFO] End time: {datetime.datetime.now()}")
        print("[INFO] Software was developed by Yong Sheng from Army Medical University.")
if __name__ == "__main__": 
    analyze_is_elements(IS_file, genome_file, num_parts)
    
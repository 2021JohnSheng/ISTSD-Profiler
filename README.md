# ISTSD-Profiler

ISTSD-Profiler is an MPI-enabled command-line pipeline for high-throughput
analysis of insertion sequence (IS) target-site duplications (TSDs) across
microbial genomes.

It is designed for large-scale genome datasets and supports parallel execution,
TSD detection, TSD length profiling, and sequence logo generation.

## Overview

ISTSD-Profiler performs the following tasks:

- maps IS elements to genome sequences in parallel with MPI
- detects TSDs flanking each insertion event
- summarizes TSD length preference across detected events
- generates TSD sequence logo plots through R (`ggseqlogo`)
- handles large genome FASTA files by splitting them into smaller parts

## Requirements

The recommended setup is to create the conda or mamba environment provided in
`environment.yml`:

```bash
mamba env create -f environment.yml
conda activate istsd-profiler
```

`requirements.txt` is provided as a minimal Python-only reference, while
`environment.yml` is the recommended installation route.

## Tested environment and runtime

ISTSD-Profiler v1.0.0 was tested using the bundled demo dataset in the
following environment:

- Operating system: Ubuntu 22.04.5 LTS under Windows Subsystem for Linux 2
  (WSL2), x86_64
- Environment manager: mamba/miniforge
- Python: 3.11.15
- pandas: 2.2.3
- mpi4py: 4.0.0
- Biopython: 1.84
- rpy2: 3.5.11
- SAMtools: 1.21
- BEDTools: 2.31.1
- SeqKit: 2.8.2
- R: 4.3.3
- ggplot2: 3.5.1
- ggseqlogo: 0.2

No non-standard hardware is required for the bundled demo dataset. Large-scale
analyses of millions of genome assemblies require MPI-capable computing
resources appropriate for the dataset size.

Using the commands shown below, environment creation with mamba completed
successfully in 67.39 seconds in the tested WSL2 environment with an existing
mamba package cache. Installation time may vary with network speed and package
cache state.

The bundled demo run completed successfully in 43.66 seconds in the tested WSL2
environment and produced the expected output directories and files, including
`TSD_analysis_results/`, `TSD_motifs/`,
`analyze_tsd_length_preference/tsd_length_preference_data.csv`, `pickle_data/`
and `merged_IS_alignment_info.txt`.

## Input

The script requires:

- an IS FASTA file: `-f / --IS`
- a genome FASTA file: `-g / --genome`
- the number of genome splits: `-n`

## Example data

The `example_data/` directory contains the test files used for release
validation.

- `rmdup_s_transposons_from_tncomp_finder_IDsimp.fna`: a non-redundant FASTA
  file containing 5030 IS element sequences, corresponding to the non-redundant
  IS dataset described in the manuscript.
- `merged_Enterococcus_faecalis_RefSeq_Assembly_Complete_and_Chromosome_sample_50.fasta`:
  a merged test genome FASTA file generated from 50 randomly selected
  *Enterococcus faecalis* genomes downloaded from the NCBI RefSeq database,
  restricted to assemblies annotated as Complete and Chromosome level.

The merged genome FASTA is stored with Git LFS because of its size. If you only
need the source code, clone with `GIT_LFS_SKIP_SMUDGE=1` or download the release
source archive. If you want to run the bundled example, install Git LFS and run
`git lfs pull` after cloning.

## Usage

```bash
mpiexec -n 14 python3 ISTSD-Profiler.py \
  -f IS_elements.fasta \
  -g genome.fasta \
  -n 22
```

Example using the bundled test data:

```bash
mpiexec -n 2 python3 ISTSD-Profiler.py \
  -f example_data/rmdup_s_transposons_from_tncomp_finder_IDsimp.fna \
  -g example_data/merged_Enterococcus_faecalis_RefSeq_Assembly_Complete_and_Chromosome_sample_50.fasta \
  -n 2
```

In the tested WSL2 environment described above, this bundled demo completed in
43.66 seconds using `mpiexec -n 2`.

Arguments:

- `-f, --IS`: FASTA file containing IS elements
- `-g, --genome`: FASTA file containing genome sequences
- `-n`: number of genome splits used for parallel processing

Before rerunning the pipeline in the same working directory, it is recommended
to remove or archive previous output directories and files to avoid mixing new
results with older appended outputs.

During execution, the pipeline may also create auxiliary genome-side files in
the same directory as the input genome FASTA, including `.fai`,
`.seqkit.fai`, and `.split/`.

## Output

Typical runtime outputs include:

- `TSD_analysis_results/`: main FASTA output directory
- `TSD_analysis_results/ISs_single_TSD_length_same_position/`: results for
  insertion loci with a single detected TSD length
- `TSD_analysis_results/ISs_multiple_TSD_lengths_same_position/`: results for
  insertion loci where multiple TSD lengths are detected at the same position
- `TSD_motifs/`: PDF sequence logo outputs generated for each IS-specific
  TSD-length category (for example, `IS256-8 bp`)
- `analyze_tsd_length_preference/tsd_length_preference_data.csv`: summary table
  for TSD length preference analysis
- `pickle_data/`: serialized intermediate results used for large-scale analysis,
  stepwise execution, and manual downstream loading of finished intermediate
  data without rerunning the full pipeline
- `merged_IS_alignment_info.txt`: merged IS mapping information collected from
  split genome files before downstream TSD analysis

For a quick first inspection, start with:

- `TSD_analysis_results/` for sequence-level output files
- `TSD_motifs/` for motif PDF figures
- `analyze_tsd_length_preference/tsd_length_preference_data.csv` for summary statistics

These runtime-generated files are excluded in `.gitignore` and are usually not
intended for version control.

## Files in This Directory

- `ISTSD-Profiler.py` - main release script
- `example_data/` - release test data for installation and execution checks
- `environment.yml` - recommended conda or mamba environment definition
- `requirements.txt` - minimal Python dependencies
- `LICENSE` - MIT license
- `CITATION.cff` - machine-readable citation metadata
- `CONTRIBUTING.md` - issue and contribution guidelines
- `SECURITY.md` - security contact and reporting guidance

## License

This project is released under the MIT License. See `LICENSE` for details.

## Citation

If you use ISTSD-Profiler, please cite this software repository. Citation
metadata is provided in `CITATION.cff`; the manuscript citation will be added
after publication.

## Contact

Contact information will be added after peer review.

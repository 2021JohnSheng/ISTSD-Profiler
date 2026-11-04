# IS pre-insertion loci and TSD sequence datasets

## 1. Contents

This release provides inferred insertion sequence (IS) target-site duplication (TSD) sequences and reconstructed pre-insertion local contexts in two layers:

- **Full dataset**: records retained before assembly-level dereplication.
- **Dereplicated dataset**: records associated with 80,828 representative genome assemblies selected within species.

Each layer contains TSD sequences and pre-insertion local contexts. A local context consists of 30 bp of the left flank, one inferred TSD copy and 30 bp of the right flank. Sequence orientation is standardized to the mapped IS direction.

## 2. Archives and extraction

The release contains four FASTA archives:

- `full_dataset_tsd_sequences.zip`
- `full_dataset_preinsertion_flank30bp.zip`
- `dereplicated_dataset_tsd_sequences.zip`
- `dereplicated_dataset_preinsertion_flank30bp.zip`

`Archive_manifest.tsv` records their filenames, byte counts and SHA-256 values.

The ZIP files contain FASTA records only. Extract all four archives into the same destination directory. Their internal paths begin with `is-preinsertion-loci/`, so the datasets will merge into one directory tree without overwriting one another.

TSD and local-context files are not guaranteed to form one-to-one pairs because complete flanking windows are not available for every record or combination. Some dereplicated FASTA files are empty because no records were retained after representative-assembly filtering. These empty files are intentional.

## 3. Sequence conventions

TSD filenames follow:

```text
<IS_Name>_<TSD_length>bp_TSD.fasta
```

Local-context filenames follow:

```text
<IS_Name>_<TSD_length>bp_TSD_pre-insertion_flank30bp_filtered.fasta
```

Each local-context sequence has the expected length of `30 + TSD length + 30 bp`.

FASTA headers identify the IS, contig, source interval and TSD length. Local-context headers additionally identify the reconstructed locus. Coordinates follow the original 1-based convention.

The header suffix `_resolved_from_duplicated` is a preserved provenance label for records resolved from duplicated or ambiguous source entries. It is not part of the IS name and does not denote an independent insertion.

## 4. Quality control and sequence filtering

Archive integrity, SHA-256 values, FASTA sequence lengths and agreement with the source files were checked. The dereplicated records were also compared record by record with the expected subsets after representative-assembly filtering.

Records containing characters other than A, T, C or G are retained. For sequence-composition, motif or information-content analyses, select records using:

```text
^[ATCG]+$
```

## 5. Interpretation boundaries

- Records are computational inferences, not experimentally validated transposition events.
- Records are not necessarily independent, even after assembly-level dereplication.
- Short inferred TSDs can be ambiguous and should be interpreted with their source accessions and coordinates.
- These files do not provide site-occupancy, insertion-rate, retention-probability or independent-origin estimates.

## 6. Release

This resource is released as `data-v1.0.0` in the ISTSD-Profiler GitHub repository:

https://github.com/2021JohnSheng/ISTSD-Profiler/releases/tag/data-v1.0.0

No DOI has been created for this release.

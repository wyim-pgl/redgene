# Transgenic Element Reference Database

## Overview

A curated, downloadable FASTA database of transgenic elements commonly found
in genetically modified organisms (GMOs). Designed for NGS-based GMO detection
and characterization pipelines.

No such ready-made, open-access FASTA database currently exists. Existing
resources (EUginius, JRC CCSIS, JRC GMO-Amplicons) either restrict bulk
sequence download, provide only amplicon-length fragments, or are behind
confidential business information barriers.

## Contents

- **element_master_list.tsv** - Master list of all elements with GenBank
  accessions, coordinates, boundary rules, and standardized headers
- **manual_sequences.fa** - Manually curated sequences for elements without
  clean single-accession GenBank entries
- **build_gmo_element_db.py** - Python script to download sequences from NCBI
  and build the combined reference FASTA
- **gmo_elements_db.fa** - (generated) Final combined FASTA database

## Element coverage

| Category          | Count | Examples                                    |
|-------------------|-------|---------------------------------------------|
| Promoters         | 12    | P-35S, P-FMV, P-Ubi1, P-Act1, P-NOS        |
| Terminators       | 7     | T-nos, T-ocs, T-35S, T-E9, T-g7, T-PinII   |
| Selectable markers| 5     | nptII, bar, pat, hpt, aadA                  |
| Herbicide genes   | 3     | cp4-epsps, gox, 2mepsps                     |
| Bt genes          | 8     | cry1Ab, cry1Ac, cry1F, cry2Ab2, cry3Bb1,    |
|                   |       | cry34Ab1, cry35Ab1, vip3Aa                  |
| Reporter/other    | 3     | gus, barnase, barstar                       |
| Regulatory        | 5     | hsp70 intron, CTP2, CTP4, LB, RB            |
| **Total**         | **43**|                                             |

## Usage

```bash
# Install dependency
pip install biopython

# Build the database
python build_gmo_element_db.py \
  --master element_master_list.tsv \
  --manual manual_sequences.fa \
  --output gmo_elements_db.fa \
  --email your@email.com

# Index for BWA
bwa index gmo_elements_db.fa

# Or index for minimap2
minimap2 -d gmo_elements_db.mmi gmo_elements_db.fa
```

## FASTA header format

```
>element_class|header_id|organism|accession|boundary_rule|version
```

Example:
```
>promoter|P-35S-CaMV|Cauliflower_mosaic_virus|V00141.1|functional_35S_promoter_region|canonical_v1
```

## Adding new elements

1. Add a row to `element_master_list.tsv` with the GenBank accession and
   coordinates
2. If no clean GenBank accession exists, add the sequence to
   `manual_sequences.fa` and set accession to "MANUAL" in the master list
3. Re-run `build_gmo_element_db.py`

## Boundary rules

Each element has an explicitly declared boundary rule that defines exactly
what portion of the source sequence is included. This is critical because:

- **Promoters**: CaMV 35S has full-length, core, and enhanced variants.
  The Ubi1 promoter MUST include exon1+intron1 for functionality.
- **Terminators**: Must be clearly separated from their associated coding
  sequences. P-35S and T-35S must never be conflated.
- **Coding sequences**: Full CDS only. Transit peptide fusions are represented
  separately.
- **Bt genes**: Commercial events often use truncated or synthetic variants.
  We provide one canonical full-length CDS per gene.

## Sources

- Won Yim, University of Nevada Reno (curated list, 2025)
- Debode et al. 2019, Sci Rep 9:15595 (element names and enrichment strategy)
- Debode et al. 2013, Eur Food Res Technol 236:659 (primer/element definitions)
- EUginius database (https://euginius.eu) - element nomenclature
- JRC GMOMETHODS / GMO-Amplicons - method target definitions
- NCBI GenBank - source sequences

## Citation

If you use this database, please cite:
[manuscript in preparation]

## License

CC-BY 4.0

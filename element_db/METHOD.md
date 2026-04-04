# GMO Element Database Construction Methods

## Overview

A comprehensive reference database of transgenic element sequences was
constructed for use in GMO positive control characterization. The database
integrates curated element sequences from published literature with detection
method target sequences from the European Union Reference Laboratory for
Genetically Modified Food and Feed (EUginius).

## Data Sources

### 1. Curated Element Sequences (`gmo_elements_db.fa`)

Transgenic element sequences were manually curated from NCBI GenBank based on
published literature and established GMO detection references (Debode et al.
2019, *Scientific Reports*). Elements include:

- **Promoters**: CaMV 35S (V00141.1), enhanced 35S, FMV (NC_003554.1), NOS
  (V00087.1), OCS (AJ311874.1), SSuAra (X13611.1), TA29 (X52283.1), Ubi1
  (S94464.1), Act1 (S44221.1), CsVMV (U59751.1)
- **Terminators**: NOS (V00087.1), OCS, CaMV 35S (V00141.1), E9 (X00806.1),
  g7, PinII, pCAMBIA 35S (AF234296.1)
- **Selectable markers**: nptII (V00618.1), bar (X17220.1), pat (M22827.1),
  hpt (U40398.1), aadA (X02340.1)
- **Herbicide tolerance genes**: CP4-EPSPS (L29358.1), GOX, 2mEPSPS, GAT
- **Insect resistance genes**: cry1Ab, cry1Ac, cry1F, cry2Ab2, cry3Bb1,
  cry34Ab1, cry35Ab1, vip3A
- **Other trait genes**: CP4-EPSPS variants, barnase, barstar, dmo,
  HPPD, aad-1, aad-12

Coordinates were defined from full-length GenBank accessions using annotated
functional boundaries (e.g., promoter region, full CDS, terminator region).
Elements without public accessions were sourced from published literature
and stored as manual sequences.

**Script**: `gmo_db.py`
**Input**: `element.list` (element definitions with accession coordinates),
`manual_sequences.fa` (sequences without public accessions)

### 2. EUginius Detection Method Sequences

#### 2.1 Method Scraping (`EUginius_scraper.py`)

All detection method entries were scraped from the EUginius database
(https://www.euginius.eu/) using method codes listed in
`EUginius_method_code.txt`. For each method, the following data were extracted:

- Method code, name, type, and description
- Target scope (element-specific, construct-specific, event-specific,
  taxon-specific)
- Target GMO name(s) and target DNA element(s)
- Oligonucleotide sequences: forward primer, reverse primer, probe
- Amplicon size and sequence (when provided)

**Output**: `euginius_primers.tsv` (structured method data),
`euginius_primers.fasta` (primer/probe/amplicon sequences)

#### 2.2 Amplicon Sequence Recovery

EUginius provides amplicon sequences for many methods, but some contain
ambiguous bases (NNN) or are not provided. Amplicon sequences were recovered
through two strategies:

**Strategy A: Local primer matching** (`local_match_amplicons.py`)

Primer sequences were matched against the curated element database
(`gmo_elements_db.fa`) and known NCBI accessions from `element.list` using
regex-based matching that supports IUPAC degenerate bases. When both forward
and reverse primers matched within a sequence at a distance consistent with
the expected amplicon size (within 30% tolerance), the intervening region was
extracted as the amplicon.

**Strategy B: Clean amplicon sequences** (`fetch_amplicons.py`)

Methods with unambiguous amplicon sequences (no N bases, >20 bp) provided
directly by EUginius were used as-is.

#### 2.3 Full Source Sequence Recovery (`fetch_full_sequences.py`)

To provide broader sequence context beyond short amplicon regions, full source
sequences were recovered for each amplicon:

1. **Local database matching**: Each amplicon was checked for exact substring
   match (forward or reverse complement) against the curated element database.
   When found, the full curated element sequence was used (e.g., full promoter,
   full CDS).

2. **NCBI BLAST**: Amplicons without local matches were queried against NCBI nt
   using megaBLAST (E-value < 1e-10, identity >= 95%). BLAST hits were filtered
   to **exclude vector/plasmid sequences** using keyword detection ("vector",
   "plasmid", "cloning", "pBI", "pCAMBIA", "pBluescript", etc.) and a maximum
   length threshold of 8,000 bp. Among non-vector hits, the shortest accession
   (most likely to represent an individual gene or element) was preferred.

3. **Fallback**: When no suitable non-vector BLAST hit was found, the original
   amplicon sequence was retained.

**Output**: `euginius_fullseq.fa` (full source sequences),
`euginius_fullseq.tsv` (summary with accessions and source types)

### 3. Database Merging (`format_final_db.py`)

The final combined database was assembled with the following priority:

1. **Curated elements** (`gmo_elements_db.fa`) - highest priority
2. **Full source sequences** (`euginius_fullseq.fa`) - preferred for EUginius
   entries
3. **Amplicon sequences** (`euginius_amplicons.fa`) - fallback when full
   sequence unavailable

Duplicate sequences were removed by exact sequence matching. For methods with
both full-length and amplicon sequences, the full-length version was used.

**Output**: `gmo_combined_db.fa` (merged FASTA), `gmo_combined_db.tsv`
(summary table)

## FASTA Header Format

All sequences follow a 6-field pipe-delimited header:

```
>class|name|organism_or_target|accession|description|version
```

- **class**: sequence category (promoter, terminator, cds, construct-specific,
  element-specific, event-specific, taxon-specific)
- **name**: element name or method code
- **organism_or_target**: source organism or target GMO name(s)
- **accession**: NCBI accession or "EUginius"
- **description**: functional description or amplicon size
- **version**: versioning tag (canonical_v1, euginius_v1)

## Software and Versions

- Python 3.11
- Biopython (Bio.Blast, Bio.Entrez, Bio.SeqIO)
- Requests + BeautifulSoup4 (web scraping)
- NCBI BLAST+ (via NCBIWWW.qblast, megaBLAST)
- NCBI Entrez (efetch for sequence retrieval)

## References

- Debode F, Hulin J, Charlier C, et al. (2019) Detection and identification of
  transgenic events by next generation sequencing combined with enrichment
  technologies. *Scientific Reports* 9:15595.
- EUginius - European GMO Initiative for a Unified Database System.
  https://www.euginius.eu/
- Grohmann L, Broll H, Dagand E, et al. (2016) Guidelines for the
  single-laboratory validation of qualitative real-time PCR methods.
  *Trends in Food Science & Technology* 69:B:150-163.

## File Inventory

| File | Description |
|------|-------------|
| `element.list` | Master element definitions (class, name, accession, coordinates) |
| `manual_sequences.fa` | Sequences without public accessions |
| `EUginius_method_code.txt` | List of EUginius method codes to scrape |
| `gmo_db.py` | Script to build curated element FASTA from element.list |
| `EUginius_scraper.py` | Scraper for EUginius method detail pages |
| `fetch_amplicons.py` | NCBI BLAST-based amplicon recovery |
| `local_match_amplicons.py` | Local primer-to-sequence matching |
| `fetch_full_sequences.py` | Full source sequence recovery (BLAST + filter) |
| `format_final_db.py` | Database merging script |
| `gmo_elements_db.fa` | Curated element sequences |
| `euginius_primers.tsv` | Scraped EUginius method data |
| `euginius_amplicons.fa` | Recovered amplicon sequences |
| `euginius_amplicons.tsv` | Amplicon recovery summary |
| `euginius_fullseq.fa` | Full source sequences (non-vector) |
| `euginius_fullseq.tsv` | Full sequence recovery summary |
| `gmo_combined_db.fa` | Final merged database |
| `gmo_combined_db.tsv` | Final database summary |

## Results Summary

### Database Composition (as of 2026-04-03)

| Source | Sequences | Description |
|--------|-----------|-------------|
| Curated elements | 29 | Manually curated from GenBank (promoters, terminators, CDS) |
| EUginius full sequences | 47 | Full source sequences recovered via local match or NCBI BLAST |
| EUginius amplicons (fallback) | 54 | Amplicon-only when full sequence unavailable |
| **Total** | **130** | **120,259 bp total** |

### Full Sequence Recovery Results

From 106 EUginius amplicon entries:

| Strategy | Count | Description |
|----------|-------|-------------|
| Local match | 25 | Amplicon found within curated element (used full element sequence) |
| NCBI BLAST | 51 | Full source accession recovered (vector/plasmid hits filtered out) |
| No hit / vector only | 30 | No suitable non-vector hit; amplicon retained as fallback |

### NCBI Accessions Used for Full Sequences

#### Element-specific (genes, promoters, terminators)

| Method Code | Target | Accession | Length | Description |
|-------------|--------|-----------|--------|-------------|
| QL-ELE-00-016 | cry1Ab/1Ac | Y09787.1 | 1,857 bp | *B. thuringiensis* cryIA(c) gene |
| QL-ELE-00-025 | pat | JX139722.1 | 493 bp | *G. max* transgenic A2704-12 GM cassette |
| QL-ELE-AINT 2-5'/2-3' | I-actin | X63830.1 | 623 bp | *O. sativa* Act1 gene |
| QL-ELE-AgroBS2 | Ti-ocs flanking | HF950132.1 | 1,130 bp | *Boechera divaricarpa* GSS |
| QL-ELE-Cry1Ac | cry1Ac | HQ161057.1 | 1,867 bp | *O. sativa* Kefeng6 genomic sequence |
| QL-ELE-P-SSuAra | P-rbcS | XM_002888530.2 | 1,170 bp | *A. lyrata* rbcS mRNA |
| QL-ELE-P-TA29 | P-TA29 | X52283.1 | 6,254 bp | Tobacco TA-29 anther-specific gene |
| QL-ELE-T-35S | T-35S-CaMV | FJ154952.1 | 1,632 bp | *B. napus* T45 border junction |
| QL-ELE-cat-F1/R1 | cat (CAT) | NG_047580.1 | 851 bp | *S. sciuri* catA gene |
| QL-ELE-cry1A 4-5'/4-3' | cry1Ab | AF465640.1 | 204 bp | *Z. mays* transgenic Bt CryIA(b) |
| QL-ELE-epsps 1-5'/3-3' | CP4-EPSPS | AF464188.1 | 1,529 bp | *G. max* CP4EPSPS gene, complete cds |
| QL-TAX-AT-001 | nos | OP712404.1 | 1,133 bp | *A. tumefaciens* nos gene |

#### Construct-specific

| Method Code | Target Events | Accession | Length | Description |
|-------------|---------------|-----------|--------|-------------|
| QL-CON-00-008 | MON events (CP4-EPSPS) | FN550387.1 | 88 bp | *Z. mays* ctp2/cp4epsps |
| QL-CON-00-012 | Papaya 55-1/63-1 | KR076687.1 | 896 bp | CMV 3a/CP genes |
| QL-CON-00-014 | Multiple (P-nos::nptII) | MF521566.1 | 3,304 bp | *Petunia* DFR gene |
| QL-CON-00-015 | P-35S::bar junction | JX139719.1 | 985 bp | *O. sativa* GM cassette |
| QL-CON-35S-F/nptII-R | P-35S::nptII | KT184682.1 | 888 bp | *Z. mays* MON863 nptII gene |
| QL-CON-CP1/pFMV2 | P-FMV::CP4-EPSPS | BT022026.1 | 1,804 bp | *A. thaliana* At2g45300 gene |
| QL-CON-gat/T-pinII | gat/T-pinII junction | KP784699.1 | 917 bp | *Z. mays* 98140 GAT gene |

#### Event-specific

| Method Code | Event | Accession | Length | Description |
|-------------|-------|-----------|--------|-------------|
| QL-EVE-16-0-1RB | Papaya 16-0-1 | KU376441.1 | 234 bp | Transgenic papaya genomic seq |
| QL-EVE-18-2-4RB | Papaya 18-2-4 | KU376439.1 | 309 bp | Transgenic papaya genomic seq |
| QL-EVE-319FalF | Falcon GS 40/90 | AY096774.1 | 216 bp | *Populus* transgene insertion site |
| QL-EVE-558F/558R | e871 B. subtilis | LT622644.1 | 2,668 bp | *B. subtilis* recA with CmR |
| QL-EVE-816AvaF | Falcon GS 40/90 | XM_009128664.3 | 1,834 bp | *B. rapa* A70 mRNA |
| QL-EVE-Kef6 | Kefeng6 | HM124448.1 | 506 bp | *O. sativa* Kefeng6 RB junction |
| QL-EVE-MDB110 | OXY-235 | EU099579.1 | 1,109 bp | *B. napus* OXY-235 3' junction |
| QL-EVE-KK179 | KK179 | HF949564.1 | 1,158 bp | *Boechera* GSS |
| QL-EVE-OS-KM2 | KMD1 | EU980363.1 | 556 bp | *O. sativa* KMD1 transgenic seq |
| QL-EVE-Oxy-3F | OXY-235 | EU099579.1 | 1,109 bp | *B. napus* OXY-235 3' junction |
| QL-EVE-SS-AquAd | AquAdvantage Salmon | S65567.1 | 3,350 bp | *M. americanus* opAFP promoter |
| QT-EVE-HF-1F | Huafan No 1 | KC820135.1 | 753 bp | *A. thaliana* NBR1 gene |
| QT-EVE-OS-002 | LLRICE62 | JQ406881.1 | 628 bp | *O. sativa* LLRICE62 5' junction |
| QT-EVE-ZM-029 | MON87419 | HF949548.1 | 1,133 bp | *Boechera* GSS |

#### Taxon-specific (endogenous reference genes)

| Method Code | Species | Accession | Length | Description |
|-------------|---------|-----------|--------|-------------|
| QL-TAX-Acc1 | *Medicago sativa* | L25042.1 | 7,175 bp | ACCase mRNA |
| QL-TAX-CP-001 | *Carica papaya* | AY803756.1 | 719 bp | Chymopapain gene |
| QL-TAX-CaMVF | CaMV | HE978785.1 | 1,391 bp | CaMV coat protein gene |
| QL-TAX-DC-001 | *Dianthus caryophyllus* | AB727362.2 | 4,604 bp | Anthocyanidin synthase promoter |
| QL-TAX-PS-001 | *Pisum sativum* | HM346476.1 | 284 bp | Lectin 6 gene |
| QL-TAX-SS-AM5F | *Salmo salar* | AY614009.1 | 4,673 bp | Growth hormone I gene |
| QL-TAX-FMVorf7 | FMV | NC_003554.1 | 7,743 bp | FMV complete genome |
| QT-TAX-BV-001 | *Beta vulgaris* | AY026353.1 | 1,296 bp | Glutamine synthetase GS2 |
| QT-TAX-CP-001 | *Carica papaya* | AY803756.1 | 719 bp | Chymopapain gene |
| QT-TAX-FatA | *Brassica napus* | XM_010437010.2 | 1,476 bp | Oleoyl-ACP thioesterase |
| QT-TAX-OS-002 | *Oryza sativa* | XM_052303253.1 | 2,740 bp | Phospholipase D alpha 2 |
| QT-TAX-SPS | *Oryza sativa* | KT225496.1 | 721 bp | Sucrose phosphate synthase |
| QT-TAX-ZM-003 | *Zea mays* | K03285.1 | 360 bp | Adh1 gene 5' end |
| QT-TAX-ZM-004 | *Zea mays* | EU963452.1 | 928 bp | mRNA clone |
| QT-TAX-wx012 | *Triticum aestivum* | LC373576.1 | 4,431 bp | Waxy protein gene |
| QT-TAX-zSSIIb | *Zea mays* | AF019297.1 | 2,480 bp | Starch synthase IIb |

### Entries Retained as Amplicon Only (No Full Sequence Available)

30 methods had no suitable non-vector BLAST hit and retained their original
amplicon sequences. These include:

- **Synthetic genes**: cry1A.105, cry2Ab2, cry1F, cry3A, cry3Bb1 (synthetic/modified,
  not in GenBank as standalone genes)
- **Short elements**: bar (60 bp amplicon, only vector hits), barstar, barnase
- **Event-specific junctions**: Several event junctions with no deposited genomic
  sequences (J101, J163, Huanong No.1, Liberator, DAS1507, etc.)
- **Other**: chloramphenicol resistance short amplicon, DsRed2, some P-35S/T-ocs
  variants

### Vector Filtering Criteria

BLAST hits were classified as vectors and excluded when:
1. Title contained keywords: "vector", "plasmid", "cloning", "pBI", "pCAMBIA",
   "pBluescript", "pUC", "pGEM", "pET-", "pBR322", "pGreen", "pBIN", "pSOUP",
   "binary vector", "shuttle vector", "expression vector", "T-DNA",
   "transformation vector", "construct", "recombinant plasmid"
2. Subject length exceeded 8,000 bp (likely vector/genome rather than gene)

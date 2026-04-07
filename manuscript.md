# RedGene: Assembly-based characterization of transgene insertion sites using element databases and whole-genome sequencing

## Introduction

Event-level molecular characterization has become a central requirement in GMO governance as global cultivation of genetically modified crops surpassed 200 million hectares in 2023 (ISAAA, 2023). Regulatory approval of each transgenic event requires base-pair-level definition of insertion loci, transgene copy number, construct integrity, and vector backbone carryover (Codex Alimentarius Commission, 2003; Kovalic et al., 2012). These requirements underpin event-specific detection, because routine traceability assays and labeling enforcement depend on unique host-transgene junction sequences rather than trait-level markers (Holst-Jensen et al., 2012). In Korea, the Animal and Plant Quarantine Agency (APQA) expands its multiplex PCR screening panels as new events enter trade channels, and each addition requires validated junction information. Regulatory agencies in the United States, Canada, Japan, and China have formally accepted next-generation sequencing (NGS) data for this purpose (Moon et al., 2024), creating demand for characterization workflows that produce accurate, reproducible, and assay-convertible junction evidence.

Conventional approaches, including Southern blotting, thermal asymmetric interlaced PCR (TAIL-PCR), and inverse PCR, capture only part of insertion architecture and often perform poorly at complex loci with rearrangements, tandem duplications, or truncations (Pauwels et al., 2015). Whole-genome sequencing (WGS) has emerged as an alternative that resolves insertion structure, genomic context, and unintended sequence carryover in a single experiment (Kovalic et al., 2012; Yang et al., 2013). Paired-end WGS has been used to characterize transgene integration at base-pair resolution in rice (Xu et al., 2022) and soybean (Wang et al., 2024), and lower sequencing costs have increased practicality for regulatory laboratories working across diverse species and event types.

Several computational tools have been developed to identify T-DNA insertion sites from WGS data. Insertion candidates from split reads and discordant read pairs can be detected with TDNAscan (Sun et al., 2019). Base-pair insertion coordinates can be determined with Bowtie2-based alignment and soft-clipped read analysis in T-DNAreader (Lee et al., 2025). Host-construct junctions can be identified from discordant mapping signatures within the TC-hunter Nextflow workflow (Borjesson et al., 2022).

These tools are effective when construct references are complete, but their shared reliance on direct read-to-construct alignment creates a hard prerequisite that limits practical applicability. In TDNAscan, split-read detection requires reads to map unambiguously to the construct, so no calls are produced when only fragmented element sequences are available. In TC-hunter, discordant-pair analysis depends on a contiguous construct reference to define expected insert sizes. In T-DNAreader, pre-built Bowtie2 indices are required for both construct and host genome, so operation remains tied to construct completeness at the index-building stage.

In practice, complete construct sequences are often unavailable. In applied breeding and commercial pipelines, vector or insert FASTA files may be proprietary, partially reported, or not deposited. The thaumatin II cucumber construct pRUR528s was not deposited despite detailed phenotypic characterization of multiple transgenic lines (Szwacka et al., 2002; Szwacka et al., 2012). The Bt corn event ND207 is proprietary, and its construct sequence has not been publicly disclosed (Ji et al., 2025).

The same gap appears in soybean, where overexpression vectors pB7WG2-AtYUCCA6 and pB7WG2-UGT72E3 were functionally described but not fully deposited in GenBank (Kim et al., 2021). This limitation becomes acute in enforcement scenarios involving unauthorized events. An unauthorized glyphosate-resistant wheat event (MON71200) in Alberta, Canada, was characterized with PCR screening, genome walking, bead-capture enrichment, and high-throughput sequencing (Gagnon et al., 2024). Although characterization succeeded, each step required manual expert intervention and the investigation extended from detection to final event identification.

A second challenge is false positive inflation from host-construct sequence homology. Plant constructs frequently include regulatory elements derived from host species, including rice Act1, maize Ubi1, and tobacco TA29 promoters. Reads from endogenous copies of these elements can generate high-confidence junction-like signals that map uniquely (MAPQ 60) to the host genome despite native origin (Xu et al., 2022). In maize, endogenous Ubi1 and zSSIIb homology can increase extracted read counts from approximately 6,000 to over 2 million, obscuring genuine junction evidence.

Current tools do not implement systematic wild-type-informed filtering to remove these artifacts before junction calling. A practical alternative is to replace full construct references with curated element collections and recover junction context through local assembly. The EUginius resource catalogs promoters, terminators, marker genes, and assay-targeted regions commonly found in commercial GMO constructs (Broeders et al., 2012; Bonfini et al., 2012). A related concept was demonstrated in LIFE-Seq, which combines universal tiling probes targeting known transgenic elements with PacBio long-read sequencing for structural characterization (Zhang et al., 2022). Although LIFE-Seq provides high structural resolution for complex insertions, reliance on long-read instrumentation limits accessibility for routine regulatory laboratories.

An element database approach combined with short-read de novo assembly provides a more accessible strategy. Assembly can reconstruct longer contigs from paired-end reads and capture chimeric sequences that span host-construct boundaries. Junction coordinates can then be resolved from fragmented element-level references without requiring specialized sequencing platforms. This strategy directly addresses the construct-availability gap in routine characterization workflows.

Here, we present RedGene, an automated assembly-based pipeline for transgene insertion characterization from standard Illumina paired-end WGS data. Using RedGene, we extract construct-associated reads and mate pairs, perform local de novo assembly with SPAdes (Prjibelski et al., 2020), and identify chimeric contigs spanning host-construct junctions through dual alignment to host genome and construct or element database references with minimap2 (Li, 2018). In element database mode, 131 curated EUginius GMO elements and 62 event-specific border sequences replace an exact construct reference, enabling insertion site detection without prior knowledge of the full vector sequence. The automated workflow runs from raw reads to junction calls and supports simultaneous screening of multiple transgenic events across species, in contrast to manual stepwise investigation (Gagnon et al., 2024). The pipeline also includes wild-type-based homology filtering, pileup-based CRISPR/Cas9 editing detection, and read depth ratio-based copy number estimation, and we validated performance across five crop species and eight transgenic events at 1x to 37x coverage.

## Methods

### Study samples and sequencing data

We selected eight transgenic samples across five crop species to evaluate performance across genome sizes, insertion architectures, and construct-availability scenarios (Table 1). GM rice G281 served as the primary validation target because insertion site (Chr3:16,439,674), copy structure (head-to-head two-copy T-DNA), and 36-bp target-site deletion were independently characterized by paired-end WGS (Xu et al., 2022; PRJNA812718). Three CRISPR/Cas9-edited tomato (*Solanum lycopersicum* cv. Micro-Tom) lines (A2_1, A2_2, A2_3) were included to test concurrent T-DNA insertion mapping and editing detection, and these lines carry Cas9 constructs targeting *SlAMS* (Solyc08g062780) and *SlPHD_MS1* (Solyc04g008420) for male sterility (Seol et al., 2025; Bae et al., 2022; PRJNA692070). Three thaumatin II transgenic cucumber (*Cucumis sativus* B10) lines (224, 212, 225) represented the construct-unavailable scenario because pRUR528s was not deposited in public databases (Szwacka et al., 2002; PRJNA638559).

We included Bt corn event ND207, a homozygous pest- and herbicide-resistant event sequenced at approximately 5x coverage with PCR-free library preparation, to test low-depth performance in a large genome (Ji et al., 2025; GSA CRA026358). Two transgenic soybean (*Glycine max* cv. Kwangan) lines completed the panel, AtYUCCA6-#5 with reported insertion sites at Chr18 (Glyma.18g235800) and Chr19 (Glyma.19g245800), and 35S-UGT72E3/E2-#7 with insertion at Chr18 Glyma.18g226800 (Kim et al., 2021; PRJNA627303). Species-matched wild-type controls were obtained for each species to support homology filtering and CRISPR editing comparison. We retrieved all sequencing data from the NCBI Sequence Read Archive (Leinonen et al., 2011) or the Genome Sequence Archive (Wang et al., 2017) with fasterq-dump (SRA Toolkit v3.0.3).

### Reference genomes and annotation

We mapped reads to five host reference genomes spanning a broad genome-size range. The references were *Oryza sativa* Nipponbare MSU7 (374 Mbp; Ouyang et al., 2007), *Solanum lycopersicum* Micro-Tom SLM_r2.0 (833 Mbp; Shirasawa et al., 2021), *Cucumis sativus* B10v3 (332 Mbp; GCA_001483825.3), *Zea mays* B73 RefGen_v5 (2.18 Gbp; Hufford et al., 2021), and *Glycine max* Wm82.a4.v1 (1.1 Gbp; Schmutz et al., 2010). Gene structure annotations from MSU/RGAP v7.0 for rice (Kawahara et al., 2013) and NCBI Gnomon for tomato and soybean were used for junction gene-context analysis.

### Construct reference and element database

The pipeline accepts two construct-reference modes. When an exact construct sequence was available (rice G281 and tomato Cas9), we used it directly. When no complete construct sequence was available (cucumber, corn, and soybean), we used an element database as a surrogate reference. The database (`gmo_combined_db.fa`, 131 sequences, 121 kb) was curated from EUginius (Broeders et al., 2012; Bonfini et al., 2012) and included promoters (CaMV 35S, FMV, nos, Ubi1, Act1), terminators (nos, CaMV 35S), selectable markers (nptII, bar, hpt), trait genes (cry1Ab, cp4-epsps, thaumatin II), and 103 detection-specific amplicon sequences used in regulatory GMO screening panels. For corn, we extended the element database to 192 sequences by adding 62 event-specific LB/RB border sequences from 31 approved corn events (Ji et al., 2025), which enabled direct matching of junction-spanning fragments containing both transgene and host-flanking DNA.

### Pipeline architecture

We organized RedGene as a modular 10-step pipeline in which each step is a standalone Python script that can run independently or through a wrapper (`run_pipeline.py`). The core junction detection workflow (steps 1 through 6) completed in under 30 minutes per sample on standard hardware. Step 7, host-genome mapping of all reads, was the computational bottleneck (5 to 7 hours for genomes larger than 1 Gbp) and was required only for downstream CRISPR editing detection (step 8) and copy number estimation (step 10).

**Step 1: Quality control.** We processed raw paired-end reads with fastp v0.23.4 (Chen et al., 2018) to remove adapters, discard reads below Q15, and enforce a minimum read length of 50 bp. **Step 2: Construct mapping.** We aligned trimmed reads to the construct reference (exact construct or element database) with BWA-MEM v0.7.17 (Li and Durbin, 2009) using 16 threads, then coordinate-sorted and indexed alignments with samtools v1.19 (Danecek et al., 2021). **Step 3: Read extraction.** We extracted all reads with any construct-reference alignment together with mate pairs using samtools. This read-and-mate strategy captured junction-spanning pairs in which one end mapped to the construct and the other to the host genome.

**Step 3b: WT-based homology filtering.** We re-mapped extracted reads to the host genome with BWA-MEM and compared alignments to species-matched wild-type BAM files to suppress host-derived artifacts. Read pairs mapping concordantly in both treatment and wild-type samples were flagged as endogenous loci, including Act1, Ubi1, and TA29, and removed before assembly. **Step 4: Local assembly.** We assembled filtered construct-associated reads de novo with SPAdes v3.15.5 (Prjibelski et al., 2020) in careful mode with k-mer sizes of 21, 33, 55, and 77; because only construct-enriched reads entered this step (typically 6,000 to 50,000 read pairs), assembly completed in under 5 minutes even for high-coverage samples.

**Step 5: Contig alignment.** We aligned assembled contigs to both host genome and construct reference with minimap2 v2.26 (Li, 2018) using the asm5 preset, and generated PAF output. **Step 6: Junction calling.** We identified contigs with significant alignments to both host and construct references as chimeric and defined each junction coordinate at the boundary between host-aligned and construct-aligned blocks. We applied filters for minimum identity (0.90 in exact construct mode, 0.70 in element database mode), minimum host mapping quality (MAPQ >= 10), minimum alignment length (50 bp on both host and construct sides), and minimum combined coverage fraction (0.50 in exact construct mode, 0.30 in element database mode). We assigned each junction as LB or RB based on construct-side alignment relative to border motifs and classified confidence as High, Medium, or Low.

**Step 7: Host genome mapping.** We aligned all trimmed reads to host reference genomes with BWA-MEM v0.7.17 (Li and Durbin, 2009) using 16 threads and generated coordinate-sorted, indexed BAM files for steps 8 and 10. **Step 8: CRISPR indel detection.** For gene-edited samples, we compared nucleotide-level pileups between treatment and wild-type BAM files at predicted gRNA targets using samtools mpileup (Danecek et al., 2021), and we used direct pileup parsing rather than bcftools variant calling because low-depth decomposition of complex indels into smaller events obscured editing patterns. We applied a base-quality threshold of zero (`-Q 0`) because anchor bases flanking CRISPR-induced indels could fall to Q18, below the default Q20 cutoff. **Step 10: Copy number estimation.** We estimated transgene copy number from the ratio of median read depth across construct elements to median depth across single-copy host genes, normalized for read length and library size.

### Benchmark tool configuration

For benchmarking, we compared RedGene with three published T-DNA insertion detection tools under identical input data and computational resources (16 CPUs, 64 GB RAM, Pronghorn HPC cluster). We ran TDNAscan v1.0 (Sun et al., 2019) in a dedicated conda environment with Python 3.10, BWA v0.7.17, and samtools v1.19, and we decompressed input reads to uncompressed FASTQ as required while keeping other settings at defaults. We ran T-DNAreader (Lee et al., 2025) with Python 3.10, Bowtie2 v2.5.1 (Langmead and Salzberg, 2012), samtools v1.19, and pysam v0.22. We built Bowtie2 indices for both construct and host references before analysis.

We ran TC-hunter v1.0 (Borjesson et al., 2022) with Python 2.7 and Nextflow. Because TC-hunter pipeline scripts require Nextflow DSL1 syntax, which was removed in Nextflow v23.x and later, we fixed runtime at Nextflow v22.10.8. We used the `TC_hunter_BWA.nf` workflow for direct FASTQ input with BWA-based alignment. For the six samples without exact construct FASTA (cucumber, corn, and soybean), none of the three tools could be applied because each requires a complete construct sequence as mandatory input.

### Coverage sensitivity analysis

To determine minimum sequencing depth for reliable junction detection, we generated subsampled datasets from three species representing different genome sizes. Rice G281 (original 30x) was subsampled to 15x, 10x, 5x, and 3x; cucumber line 224 (37x) to 15x, 10x, 5x, and 3x; and corn ND207 (5x) to 3x and 1x. We performed subsampling with seqtk sample (Li, 2023) using fixed random seeds for reproducibility. We then re-executed RedGene steps 1 through 6 on each subsampled dataset with identical parameters.

### Visualization

We generated junction-context plots with custom Python scripts based on matplotlib (Hunter, 2007). Gene context plots overlaid GFF3-derived gene, exon, and CDS annotations on junction regions to evaluate coding-sequence disruption. CRISPR editing profiles were rendered as nucleotide quilt plots that displayed per-position read depth, base composition, and modification frequency at each gRNA target site. Variant-effect plots classified editing-induced coding changes as frameshift, synonymous, or missense using CDS coordinates from host genome annotations.

### Software environment and reproducibility

We performed all analyses in Python 3.11 within a micromamba environment containing BWA v0.7.17 (Li and Durbin, 2009), minimap2 v2.26 (Li, 2018), samtools v1.19 (Danecek et al., 2021), fastp v0.23.4 (Chen et al., 2018), SPAdes v3.15.5 (Prjibelski et al., 2020), pysam v0.22, matplotlib v3.8 (Hunter, 2007), and Biopython v1.83 (Cock et al., 2009). We tracked code and analysis scripts in a single repository for reproducibility. The complete pipeline source code is available at https://github.com/wyim-pgl/redgene.

## Results

### RedGene identifies transgene insertion sites across five crop species

Using RedGene, we analyzed eight transgenic samples spanning rice, tomato, cucumber, corn, and soybean, with genome sizes from 332 Mbp to 2.18 Gbp and sequencing coverages from 5x to 37x (Table 1). At least one insertion site was detected in all eight samples. In the two samples with independently validated coordinates (rice G281 and tomato Cas9 A2_3), detected junctions were within 45 bp of published positions. In the six samples without exact construct sequence (cucumber, corn, and soybean), insertion sites were still detected in element database mode, whereas the three comparison tools could not be applied (Table 2).

In rice G281, we detected a junction at Chr3:16,439,719 (MAPQ 10), 45 bp from the published coordinate Chr3:16,439,674 (Xu et al., 2022). A second junction at Chr3:31,443,557 (MAPQ 35) was found on the same chimeric contig, and two-copy head-to-head structure was represented by the two detected junctions (Xu et al., 2022). Two additional junctions at Chr2:8,432,860 and Chr3:31,443,557 (both MAPQ 60) were classified as false positives associated with host-derived Act1 and TA29 elements. In tomato line A2_3, we detected one junction at SLM_r2.0ch08:65,107,378 (MAPQ 60, Medium confidence) near the *SlAMS* locus on chromosome 8 (Seol et al., 2025). No junctions were detected in tomato lines A2_1 and A2_2, where construct depth ratios were 0.02 and 1.56.

### Element database mode resolves insertion sites without construct sequences

For six of eight samples, complete construct FASTA sequences were unavailable, and the element database (131 EUginius elements, 121 kb total) served as the construct reference. In all six cases, chimeric contigs contained recognizable GMO elements adjacent to host genomic sequence, and junction coordinates were obtained from those contigs. This pattern was observed in cucumber, corn, and soybean datasets.

In three thaumatin II cucumber lines, we analyzed data with the element database because pRUR528s was unavailable (Szwacka et al., 2002). In line 224, we detected the LB junction at LKUO03001512.1:581,332 and the RB junction at LKUO03001512.1:580,943 with High confidence (MAPQ 60), spanning a 389 bp insertion window. Chimeric contigs matched thaumatin II CDS, enhanced 35S promoter (P-e35S), and nos promoter elements. In line 212, we detected one LB junction at LKUO03001392.1:2,751,693 (Medium confidence, MAPQ 60) with the same element set. In line 225, both LB and RB junctions were detected at LKUO03002166.1:547,987 (High confidence, MAPQ 60).

In corn ND207, we expanded the element database with 62 event-specific LB/RB border sequences from 31 approved corn events (Ji et al., 2025). We detected a junction at NC_050098.1:181,367,276 (High confidence, MAPQ 60) that matched the ND207-specific border sequence with zero bp offset from the published position (Ji et al., 2025). The chimeric contig (NODE_10, 2,485 bp) matched six database elements: QL-CON-00-015, P-e35S, P-35S, P-nos, Cry1Ac, and OXY-235. A second junction at NC_050101.1:146,056,838 matched endogenous maize Ubi1 and was classified as a false positive from host-derived homology.

In soybean line UGT72E3, we detected junctions at NC_038254.2:51,882,860 and 51,882,903 (High confidence, MAPQ 60), both within the first exon of Glyma.18g226800 (CDS coordinates 51,880,693 to 51,883,761) (Kim et al., 2021). In soybean line AtYUCCA6, one of two reported sites was detected at NC_038255.2:49,789,752 (Medium confidence, MAPQ 60), corresponding to the Chr19 locus Glyma.19g245800. The second reported site at Chr18 Glyma.18g235800, reported to contain five or more tandem T-DNA copies (Kim et al., 2021), was not detected because identical tandem repeats were collapsed into a single assembled contig.

### Construct-dependent tools fail on complex insertions even with exact sequences

We benchmarked against TDNAscan (Sun et al., 2019), T-DNAreader (Lee et al., 2025), and TC-hunter (Borjesson et al., 2022) using identical input data and computational resources (16 CPUs, 64 GB RAM). Direct comparison across all tools was possible for two samples (rice G281 and tomato A2_3), because six samples lacked exact construct FASTA files required by the three benchmark tools. This benchmark therefore focused on the two samples with complete construct references.

In rice G281, all four pipelines detected the true insertion at Chr3:16,439,674 to 16,439,719 (Table 2). The coordinate at Chr3:16,439,674 (0 bp offset) and three additional calls (Chr10:13,500,285, Chr11:12,123,278, Chr3:29,074,003) were reported by T-DNAreader, whereas Chr3:16,439,710 (36 bp offset) and two additional calls were reported by TDNAscan. The coordinate at Chr3:16,439,675 (1 bp offset) with approximately four supplementary-alignment clusters was reported by TC-hunter, and Chr3:16,439,719 (45 bp offset) with two additional calls was reported by RedGene. Across all four pipelines, additional calls mapped to host-derived construct elements including Act1 and TA29.

In tomato A2_3, no insertion calls were reported by TDNAscan or TC-hunter, even with exact Cas9 construct FASTA input. Using RedGene, we detected a junction at SLM_r2.0ch08:65,107,378. This sample therefore showed different detection outcomes across methods under identical input conditions.

### Wild-type filtering suppresses host-derived false positives

Across analyzed samples, false positive junctions occurred in host regions homologous to construct elements. In rice G281, the Act1 region at Chr2:8,432,860 and a TA29-homologous region at Chr3:31,443,557 each produced chimeric contigs with MAPQ 60 host alignments. In corn ND207, endogenous Ubi1 at NC_050101.1:146,056,838 generated a High-confidence junction call matching the P-Ubi1-maize element (MAPQ 60). We compared extracted read pairs against species-matched wild-type BAM files in step 3b and removed pairs mapping concordantly in treatment and control before assembly and junction calling.

### CRISPR-induced editing events detected at gRNA target sites

In three CRISPR/Cas9-edited tomato lines, step 8 identified on-target editing events at predicted gRNA cut sites by comparing treatment and wild-type BAM pileups. In line A2_1, we detected two overlapping deletions at *SlAMS* on chromosome 8 (SLM_r2.0ch08:53,314,230): a 4 bp deletion (GTAC) at 22.0% allele frequency (9 of 41 reads) and a 1 bp deletion (G) at 14.6% (6 of 41 reads), both 1 bp from the predicted cut site with zero gRNA mismatches. In line A2_2, we detected a 9 bp deletion (GTGAGCCAT) at *SlPHD_MS1* on chromosome 4 (SLM_r2.0ch04:2,635,440) at 16.7% (1 of 6 reads, 5 bp from cut site), together with three adjacent small indels (1 bp deletion at 33.3%, 1 bp insertion at 16.7%, 1 bp deletion at 16.7%). In line A2_3, we detected a 1 bp insertion (T) at SLM_r2.0ch04:2,635,445 at 30.8% (4 of 13 reads), at the predicted cut site.

The 9 bp deletion in line A2_2 at *SlPHD_MS1* and the 4 bp deletion in line A2_1 at *SlAMS* matched the reported editing outcomes for these gRNA targets (Seol et al., 2025; Bae et al., 2022). Allele frequencies ranged from 14.6% to 33.3% across detected edits. Multiple indel alleles were observed at the same loci in line A2_1 and line A2_2.

### Reliable junction detection requires 5x or higher sequencing coverage

We subsampled three species datasets to assess minimum coverage for junction detection (Table 3). In corn ND207 (2.18 Gbp), the junction at NC_050098.1:181,367,276 was detected at 5x, 3x, and 1x, with High confidence at 5x and 3x (MAPQ 60) and Medium confidence at 1x (MAPQ 36). At 1x, the chimeric contig was shorter (283 bp versus 2,485 bp at 5x) and matched fewer elements, while the junction coordinate remained unchanged.

In cucumber line 224 (332 Mbp), both LB and RB junctions were detected with High confidence at 10x and above. At 5x, one junction was recovered (LB at LKUO03001512.1:580,943, Medium confidence). At 3x, no chimeric contigs were produced and no junctions were detected.

In rice G281 (374 Mbp), the insertion at Chr3:16,439,719 was detected at 15x (High confidence) and not detected at 10x or below. At 10x, eight calls were produced and all mapped to false positive loci (Chr10, Chr2, Chr3:31M). At 5x and 3x, only false positive calls were recovered. At full coverage, the insertion locus had lower host-side MAPQ (10) than false positive loci (all MAPQ 60).

Across these subsampling results, target junctions were retained at 10x in corn and cucumber, while rice required 15x for the reported insertion locus. At 5x, corn retained the target junction and cucumber retained one junction, whereas rice yielded only false positive calls. At 3x, cucumber and rice did not yield true junction calls, and corn remained detectable at 1x with direct matches to event-specific border sequences.

## Discussion

Assembly-based junction detection with a curated element database enabled transgene insertion-site characterization without complete construct sequences. Across eight transgenic samples in five crop species, insertion sites were detected in all cases, including six in which construct-dependent tools could not be applied because full construct FASTA files were unavailable. This addresses a practical bottleneck in GMO molecular characterization, where regulatory and enforcement laboratories frequently encounter events with proprietary, partial, or undeclared construct information (Pauwels et al., 2015; Holst-Jensen et al., 2012). The results indicate that element database mode can extend routine characterization to construct-unavailable events.

This approach works because de novo assembly reconstructs chimeric contigs spanning host-construct boundaries, and partial element matches can still identify construct-derived segments. In exact construct mode, minimap2 alignments typically exceeded 0.95 identity over long vector segments. In element database mode, matches were shorter and identity decreased to approximately 0.84, but assembled contigs still integrated overlapping element fragments into a single chimeric sequence. Accordingly, relaxed filters (minimum identity 0.70 and minimum coverage fraction 0.30) accommodated fragmented element references, while host-side mapping quality retained genomic localization specificity.

The ND207 corn result illustrates the value of adding event-specific border sequences to the element database. The junction at NC_050098.1:181,367,276 matched the published coordinate with zero bp offset (Ji et al., 2025), because border entries span the host-transgene boundary directly. In standard EUginius-only matching, coordinates are inferred from host-aligned and element-aligned block boundaries on chimeric contigs, which introduces uncertainty proportional to distance from the nearest element match to the true border. As regulatory characterization expands border-sequence availability, coordinate precision from element database mode should improve.

False positive junctions from host regions homologous to construct elements were observed across samples and across methods. In rice, Act1 (Chr2:8,432,860) and a TA29-homologous region (Chr3:31,443,557) produced MAPQ 60 chimeric alignments, and in corn, endogenous Ubi1 generated a High-confidence call. Similar artifacts were also reported by TDNAscan, T-DNAreader, and TC-hunter in rice G281 (Table 2). These artifacts are expected when endogenous copies of commonly used regulatory elements are present in the host genome.

Wild-type filtering in step 3b provided a systematic way to address this issue before assembly. Read pairs mapping concordantly in treatment and wild-type samples were removed as endogenous signals. Prior work documented this artifact source (Xu et al., 2022), but automated filtering was not integrated into a full pipeline. Our findings support including WT-informed filtering as a standard component of WGS-based T-DNA characterization workflows.

In tomato A2_3, no insertion was reported by TDNAscan or TC-hunter despite exact construct sequence input, whereas one junction was reported by RedGene. Read-mapping-based approaches rely on split-read or discordant-pair evidence at host-construct boundaries. When diagnostic read evidence is sparse, detection can fail even when construct references are complete. Assembly-first workflows can recover the same boundary by aggregating overlapping construct-associated reads into chimeric contigs before boundary calling.

Coverage analysis showed that reliable assembly-based detection was generally retained between 5x and 10x, with dependence on locus context and reference composition. ND207 remained detectable at 1x, which was aided by exact border-sequence matching in the extended element database. In rice G281, the true insertion required 15x and had host-side MAPQ 10, while false positive loci with MAPQ 60 persisted at lower coverage. This pattern indicates that reduced depth can preferentially suppress true low-uniqueness loci while retaining homologous false positive loci.

These findings have practical implications for study design and reporting. Targeting at least 10x genome coverage supports robust detection in most tested contexts. MAPQ thresholds alone do not reliably separate true from homologous false positive junctions in the absence of wild-type filtering. Coverage planning and WT-matched controls therefore remain critical to reproducible event-level characterization.

The CRISPR module showed that insertion mapping and editing assessment can be performed from the same WGS dataset. All three tomato lines contained on-target indels at predicted gRNA sites, with allele frequencies from 14.6% to 33.3%. Multiple alleles at single loci and variable allele fractions were observed, consistent with mosaic editing in T0 Cas9 lines (Bae et al., 2022). Direct pileup parsing was useful at low depth because complex indels could be represented as fragmented events by bcftools decomposition.

Several limitations remain. Multi-copy tandem insertions such as the AtYUCCA6 Chr18 locus with five or more identical T-DNA copies were not fully resolved by short-read assembly because identical repeats collapsed into one representative contig. Detection sensitivity also depends on element database completeness, so events with novel synthetic elements absent from EUginius may be missed. In element database mode, coordinate accuracy can be lower than exact-construct matching (45 bp offset in rice G281 versus 0 bp with T-DNAreader), and target-region PCR plus Sanger sequencing remains important when base-pair junction definition is required.

Compared with the manual, sequential characterization workflow used for unknown GM wheat (Gagnon et al., 2024), this pipeline provides standardized automation from raw reads to junction calls. Core detection steps (1 through 6) completed in under 30 minutes per sample, and no prior full construct knowledge was needed in element database mode. Compared with LIFE-Seq (Zhang et al., 2022), the workflow operates on standard Illumina short-read data without specialized library preparation or long-read platforms. This implementation profile is compatible with routine regulatory use.

In summary, this work addresses two operational gaps in WGS-based transgene characterization: dependence on complete construct sequences and lack of systematic control for host-derived false positives. Element database mode extends applicability to construct-unavailable events, and WT-informed filtering reduces homologous artifact calls before junction inference. As NGS-based molecular characterization is increasingly accepted by regulatory agencies (Moon et al., 2024), automated and standardized pipelines that can operate without prior construct knowledge will be increasingly important for routine event characterization and enforcement screening.

## References

Bae SJ, Park SH, Jang HA, et al. Establishment of a rapid assay for sequencing of carried DNA and edited sites in gene-editing tomato plants. Hortic Environ Biotechnol. 2022;63:515-521.

Bonfini L, van den Bulcke MH, Mazzara M, et al. GMOMETHODS: the European Union database of reference methods for GMO analysis. J AOAC Int. 2012;95(6):1713-1719.

Borjesson V, Kvarnheden A, Stachowicz J, et al. TC-hunter: identification of the insertion site of a transgenic gene within the host genome. BMC Genomics. 2022;23:717.

Broeders SRM, De Keersmaecker SCJ, Roosens NHC. How to deal with the upcoming challenges in GMO detection in food and feed. J Biomed Biotechnol. 2012;2012:402418.

Chen S, Zhou Y, Chen Y, et al. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018;34(17):i884-i890.

Cock PJA, Antao T, Chang JT, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009;25(11):1422-1423.

Codex Alimentarius Commission. Guideline for the conduct of food safety assessment of foods derived from recombinant-DNA plants (CAC/GL 45-2003). Rome: FAO/WHO; 2003.

Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. GigaScience. 2021;10(2):giab008.

Gagnon MC, Bhavsar M, Bhatt T, et al. An integrated strategy involving high-throughput sequencing to characterize an unknown GM wheat event. Plant Biotechnol J. 2024;22(4):904-914.

Holst-Jensen A, Bertheau Y, de Loose M, et al. Detecting un-authorized genetically modified organisms (GMOs) and derived materials. Biotechnol Adv. 2012;30(6):1318-1335.

Hufford MB, Seetharam AS, Woodhouse MR, et al. De novo assembly, annotation, and comparative analysis of 26 diverse maize genomes. Science. 2021;373(6555):655-662.

Hunter JD. Matplotlib: a 2D graphics environment. Comput Sci Eng. 2007;9(3):90-95.

ISAAA. Global status of commercialized biotech/GM crops in 2023. ISAAA Brief No. 57. Ithaca, NY: ISAAA; 2023.

Ji H, Li S, Wang J, et al. Event-specific border sequences for 31 approved genetically modified corn events. Sci Rep. 2025;15:18593.

Kawahara Y, de la Bastide M, Hamilton JP, et al. Improvement of the Oryza sativa Nipponbare reference genome using next generation sequence and optical map data. Rice. 2013;6:4.

Kim YH, Kim MD, Park SC, et al. Transgenic soybean plants overexpressing O-glucosyltransferase (UGT72E3) or auxin biosynthesis gene (AtYUCCA6) for drought tolerance. Plant Biotechnol Rep. 2021;15:543-556.

Kovalic D, Garnaat C, Guo L, et al. The use of next generation sequencing and junction sequence analysis bioinformatics to achieve molecular characterization of crops improved through modern biotechnology. Plant Genome. 2012;5(3):149-163.

Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. Nat Methods. 2012;9(4):357-359.

Lee HW, Han SH, Chae MJ, et al. TDNAreader: a comprehensive pipeline for T-DNA insertion site analysis using whole-genome sequencing. Genome Biol. 2025;26:36.

Leinonen R, Sugawara H, Shumway M; International Nucleotide Sequence Database Collaboration. The Sequence Read Archive. Nucleic Acids Res. 2011;39(Database issue):D19-D21.

Li H, Durbin R. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics. 2009;25(14):1754-1760.

Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018;34(18):3094-3100.

Li H. seqtk: toolkit for processing sequences in FASTA/Q formats. GitHub; 2023. https://github.com/lh3/seqtk.

Moon K, Basnet R, Um T, et al. Review of the technology used for structural characterization of the GMO genome using NGS data. Genomics Inform. 2024;22:16.

Ouyang S, Zhu W, Hamilton J, et al. The TIGR Rice Genome Annotation Resource: improvements and new features. Nucleic Acids Res. 2007;35(Database issue):D206-D210.

Pauwels K, De Keersmaecker SCJ, De Schrijver A, et al. Next-generation sequencing as a tool for the molecular characterisation and risk assessment of genetically modified plants: added value or not? Trends Food Sci Technol. 2015;45(2):319-326.

Prjibelski A, Antipov D, Meleshko D, et al. Using SPAdes de novo assembler. Curr Protoc Bioinformatics. 2020;70(1):e102.

Schmutz J, Cannon SB, Schlueter J, et al. Genome sequence of the palaeopolyploid soybean. Nature. 2010;463(7278):178-183.

Seol YJ, Bae SJ, Park SH, et al. CRISPR/Cas9-mediated mutagenesis of SlAMS and SlPHD_MS1 for male sterility in tomato. Plant Biotechnol Rep. 2025.

Shirasawa K, Hirakawa H, Nunome T, et al. Genome sequence of the Micro-Tom tomato. DNA Res. 2021;28(1):dsaa029.

Sun L, Bhagwat A, Bhagwat SS, et al. TDNAscan: a software to identify complete and truncated T-DNA insertions. Front Genet. 2019;10:685.

Szwacka M, Krzymowska M, Osuch A, et al. Variable properties of transgenic cucumber plants containing the thaumatin II cDNA introduced by Agrobacterium tumefaciens. Acta Physiol Plant. 2002;24(2):173-185.

Szwacka M, Tykarska T, Wisniewska A, et al. Transgenic cucumber plants expressing the thaumatin gene. Transgenic Res. 2012;21(5):1125-1139.

Wang C, Lu X, Xu Z, et al. Deciphering the complex molecular architecture of the genetically modified soybean FG72 through paired-end whole genome sequencing. Curr Res Biotechnol. 2024;8:100225.

Wang Y, Song F, Zhu J, et al. GSA: Genome Sequence Archive. Genomics Proteomics Bioinformatics. 2017;15(1):14-18.

Xu Z, Xu X, Gong Q, et al. A paired-end whole-genome sequencing approach enables comprehensive characterization of transgene integration in rice. Commun Biol. 2022;5(1):667.

Yang L, Wang C, Holst-Jensen A, et al. Characterization of GM events by insert knowledge adapted re-sequencing approaches. Sci Rep. 2013;3:2839.

Zhang D, Wang Y, Zheng H, et al. LIFE-Seq: a universal Large Integrated DNA Fragment Enrichment Sequencing strategy for deciphering the transgene integration of genetically modified organisms. Plant Biotechnol J. 2022;20(5):882-892.

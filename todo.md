# RedGene Pipeline - TODO

## Current: Corn ND207 (DONE — steps 1-6)
- [x] Download GSA CRA026358 PCR-Free data (ND207_PCRFree_R1/R2.fq.gz)
- [x] Download B73 RefGen_v5 reference (2.18 Gbp)
- [x] Create corn border sequence DB (62 LB/RB from Table S1, Sci Rep 2025)
- [x] Create combined construct DB (element_db 130 + corn borders 62 = 192 seqs)
- [x] Step 1: QC (fastp)
- [x] Steps 2-6: construct map → extract → assembly → contig map → junction detection
- [x] Fix: --min-identity 0.70 for element_db/combined_db in run_pipeline.py
- [x] Junction found: **Chr3:181,367,276** (ND207-LB + Bt construct, High confidence)
- [x] Update results.md with findings
- [ ] Step 7: Host mapping (BWA, ~5-7h for 2.18 Gbp genome)
- [ ] Step 10: Copy number estimation

---

## Next: Soybean G2-EPSPS/GAT

**Rationale**: Larger crop genome (~1.1 Gbp) for host complexity testing

**Weakness**: Accession visibility is somewhat limited — original Frontiers 2016 paper did not deposit WGS data publicly

**Paper**: Yang et al. (2016) Front Plant Sci 7:1009. DOI:10.3389/fpls.2016.01009
- Events: GE-J16 (Chr19:50,543,767–50,543,792), ZH10-6 (Chr17:7,980,527–7,980,541)
- Platform: HiSeq 2500, 125bp PE, ~21x coverage
- Reference used: Wm82.a2.v1 (Phytozome)

**Data**:
- Reference genome: GCF_000004515.6 Glycine max Wm82.a4.v1 (~1.1 Gbp, 20 chr + scaffolds)
- SRA: Not publicly deposited (CAAS, China) — need alternative source or contact authors
- Construct: G2-EPSPS (glyphosate resistance) + GAT (glufosinate acetyltransferase)

### Tasks
- [x] Download soybean reference genome (Wm82.a4.v1) — 946 MB
- [x] Index reference (BWA + samtools faidx) — BWA index 21 min
- [ ] Obtain WGS data (contact authors or find alternative accession)
- [ ] Create construct reference DB (G2-EPSPS + GAT + nptII + vector elements)
- [ ] Run pipeline steps 1-6
- [ ] Run step 7 (host mapping)
- [ ] Analyze junction results — expected: Chr19 and Chr17 insertions
- [ ] Document in results.md

---

## Next: Cucumber (PRJNA638559)

**Rationale**: Public data for species expansion validation (~350 Mbp genome, smallest test species)

**Weakness**: Insertion ground truth is relatively weak

**Paper**: Szwacka et al. (2021) PMC8139995. Thaumatin II transgenic cucumber lines.
- Construct: pRUR528 (thaumatin II under 35S CaMV + nptII under nos promoter)
- **Line 212**: Chr6, intergenic, 1 T-DNA copy, 1304bp deletion
- **Line 224**: Chr2, promoter of G6936 (disrupts G6838), 1 T-DNA copy, 361bp deletion
- **Line 225**: Chr2, intergenic, 2 T-DNA copies + vector backbone, 95bp deletion

**Data**:
- Reference genome: GCA_001483825.3 CucSat-B10v3 (332 Mbp)
- SRA runs (HiSeq 2000, 100bp PE, ~37x each):
  - SRR12081227 — thaumatin line 212 (63.4M reads)
  - SRR12081931 — thaumatin line 224 (63.4M reads)
  - SRR12082195 — thaumatin line 225 (63.3M reads)

### Tasks
- [x] Download cucumber reference genome (B10v3) — 332 MB
- [x] Index reference (BWA + samtools faidx) — BWA index 6 min
- [x] Download SRA data — 3 runs, ~3.9 GB each (23 GB total)
- [x] Add thaumatin II (J01209.1) to element_db — now 131 sequences
- [x] Run pipeline steps 1-6 for each line — all 3 completed
- [x] Fix s06_junction.py coverage calculation bug (union merge for element_db)
- [x] Junction detection:
  - Line 212: LKUO03001392.1:2,751,693 (LB, Medium)
  - Line 224: LKUO03001512.1:580,628–581,332 (LB+RB, High)
  - Line 225: LKUO03002166.1:547,987 (LB+RB, High)
- [x] Document initial results in results.md
- [x] Coverage sensitivity tests (15x, 10x, 5x, 3x): ≥10x reliable, 5x partial, 3x fails
- [ ] Run step 7 (host mapping) for copy number estimation

---

## Completed Species

### Rice G281
- Insertion: Chr3:16,439,674 (bp resolution)
- Copy number: 2 (head-to-head tandem)
- Coverage sensitivity: reliable at >= 15x (374 Mbp genome)

### Tomato Micro-Tom (PRJNA692070)
- WT (SRR13450615): No transgene (negative control)
- Cas9 A2_1 (SRR13450618): Cas9 removed, heterozygous edits (SlAMS 6bp del, SlPHD_MS1 9bp del)
- Cas9 A2_2 (SRR13450617): Cas9 present, homozygous edits, T-DNA on chr08
- Cas9 A2_3 (SRR13450616): Cas9 present, heterozygous edits, T-DNA at chr08:65,107,378
- Coverage sensitivity: reliable at >= 10x (833 Mbp genome)

# element_db/Makefile
# Builds gmo_combined_db_v2.fa (4-way src tag + cd-hit-est dedup @ 95% identity)
# and gmo_combined_db_manifest.tsv (T1 audit-header consumer).
#
# cd-hit-est parameter rationale (Issue #10 M-6)
# ----------------------------------------------
#   -c 0.95   Cluster sequences with >=95% identity.  Chosen because 90% is
#             too permissive (merges distinct payload variants such as
#             nptII vs hpt fragments that share the Tn5 transposon backbone)
#             while 99% is too strict (keeps trivially different CDS copies
#             that collapse later at BLAST time anyway).  95% is the
#             accepted community default for GMO element curation.
#   -n 10     Word size matched to -c.  cd-hit-est requires -n 10 when
#             -c is in [0.95, 1.0]; any lower -n triggers a runtime error.
#             Source: cd-hit-est user guide §2.1.2.
#   -M 8000   8 GB RAM cap — comfortably above our ~300MB working-set but
#             stays safe on the Pronghorn login nodes used for ad-hoc builds.
#   -T 4      4 threads — build is I/O bound past this point.
.PHONY: all clean

DB      := gmo_combined_db_v2.fa
RAW     := _raw_combined.fa
MANI    := gmo_combined_db_manifest.tsv
CDH_C   := 0.95    # see header comment for rationale (Issue #10 M-6)

SRCS := common_payload.fa gmo_combined_db.fa cas9_sgrna.fa euginius_missing.fa payload_cds.fa crl_amplicons.fa

all: $(DB) $(MANI)

$(RAW): $(SRCS) tag_and_concat.py
	@echo "[build] tagging source + concat..."
	python3 tag_and_concat.py $(SRCS) > $@.tmp
	mv $@.tmp $@

$(DB): $(RAW)
	@echo "[build] cd-hit-est -c $(CDH_C)..."
	cd-hit-est -i $< -o $@.tmp -c $(CDH_C) -n 10 -M 8000 -T 4
	mv $@.tmp $@
	mv $@.tmp.clstr $@.clstr  # preserve cluster listing for audit

$(MANI): $(DB)
	@printf "name\tmd5\tbuild_date\tseq_count\n" > $@.tmp
	@md5=$$(md5sum $< | awk '{print $$1}'); \
	 count=$$(grep -c '^>' $<); \
	 date=$$(date -I); \
	 printf "element_db\t%s\t%s\t%s\n" "$$md5" "$$date" "$$count" >> $@.tmp
	mv $@.tmp $@

clean:
	rm -f $(DB) $(RAW) $(MANI) $(DB).clstr

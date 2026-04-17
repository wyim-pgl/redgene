# element_db/Makefile
# Builds gmo_combined_db_v2.fa (4-way src tag + cd-hit-est dedup @ 95% identity)
# and gmo_combined_db_manifest.tsv (T1 audit-header consumer).
.PHONY: all clean

DB      := gmo_combined_db_v2.fa
RAW     := _raw_combined.fa
MANI    := gmo_combined_db_manifest.tsv
CDH_C   := 0.95

SRCS := common_payload.fa gmo_combined_db.fa cas9_sgrna.fa euginius_missing.fa payload_cds.fa

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

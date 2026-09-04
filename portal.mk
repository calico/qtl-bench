# Download every public input the benchmark is built from.
#
#   make -f portal.mk all      QTL summary stats + fine-mapping + expression + annotation
#   make -f portal.mk eqtl     one modality
#
# Everything lands in $(PORTAL), which is data (gitignored); this file is code.

PORTAL  ?= data/gtex11/portal

BASE    := https://storage.googleapis.com/adult-gtex/bulk-qtl/v11
CIS     := $(BASE)/single-tissue-cis-qtl
SUSIE   := $(BASE)/susie-qtl
EXPR    := https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq
GENCODE := https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48

GCT     := GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct
GTF     := gencode48_basic_nort_protein.gtf

.PHONY: all eqtl sqtl paqtl expr genes

all: eqtl sqtl paqtl expr genes

eqtl: $(PORTAL)/GTEx_Analysis_v11_eQTL/.unpacked \
      $(PORTAL)/GTEx_Analysis_v11_eQTL/.susie_unpacked

sqtl: $(PORTAL)/GTEx_Analysis_v11_sQTL/.unpacked \
      $(PORTAL)/GTEx_Analysis_v11_sQTL/.susie_unpacked

paqtl: $(PORTAL)/GTEx_Analysis_v11_apaQTL/.unpacked \
       $(PORTAL)/GTEx_Analysis_v11_apaQTL/.susie_unpacked

expr: $(PORTAL)/$(GCT)

genes: $(PORTAL)/$(GTF)

# ── QTL summary statistics and fine-mapping ───────────────────────────────────
# signif_pairs come from the cis-QTL tarballs, SuSiE_summary from the susie ones.

$(PORTAL)/%_SuSiE.tar:
	mkdir -p $(PORTAL)
	wget -c -O $@ $(SUSIE)/$*_SuSiE.tar

$(PORTAL)/%.tar:
	mkdir -p $(PORTAL)
	wget -c -O $@ $(CIS)/$*.tar

$(PORTAL)/%/.susie_unpacked: $(PORTAL)/%_SuSiE.tar
	tar -xf $< -C $(PORTAL)
	touch $@

$(PORTAL)/%/.unpacked: $(PORTAL)/%.tar
	tar -xf $< -C $(PORTAL)
	touch $@

# ── expression ────────────────────────────────────────────────────────────────

$(PORTAL)/$(GCT):
	mkdir -p $(PORTAL)
	curl -sL $(EXPR)/$(GCT).gz | gunzip -c > $@

# ── gene annotation ───────────────────────────────────────────────────────────
# GENCODE v48 basic, restricted to protein-coding transcripts and stripped of
# readthroughs, which would otherwise contaminate UTR3 annotation with
# spuriously extended 3' UTRs and pollute the splice-site and PAS sets.
# Piped rather than staged: the intermediates are ~1 GB each. Result is 890 MB.

$(PORTAL)/$(GTF):
	mkdir -p $(PORTAL)
	curl -sL $(GENCODE)/gencode.v48.basic.annotation.gtf.gz | gunzip -c \
	| grep -F 'transcript_type "protein_coding";' \
	| grep -v 'readthrough_transcript\|readthrough_gene' > $@

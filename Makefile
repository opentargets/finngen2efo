# Makefile --- 

# Copyright (C) 2020 Open Targets

# Author: David Ochoa <ochoa@ebi.ac.uk>

#################################
# Constants
#################################

EFORELEASETAG = v3.39.0

#################################
# Paths (DIRECTORIES)
#################################

ROOTDIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
TEMPDIR= $(ROOTDIR)/temp
SRCDIR= $(ROOTDIR)/src

#################################
# Paths (Tools)
#################################

# bins
RSCRIPT ?= $(shell which Rscript)
CURL ?= $(shell which curl)
ROBOT ?= $(shell which robot)
CUT ?= $(shell which cut)
ONTOMA ?= $(shell which ontoma)
SED ?= $(shell which sed)
JQ ?= $(shell which jq)

RMIRROR ?= http://cran.uk.r-project.org

EFOFILE ?= $(TEMPDIR)/efo_otar_slim_$(EFORELEASETAG).owl
EFOCACHE ?= $(TEMPDIR)/efo_otar_slim_$(EFORELEASETAG).owl
MAPPINGSFILE ?= $(TEMPDIR)/all_mappings.tsv
MANIFESTFILE_R2 ?= $(TEMPDIR)/finngen_manifest_r2.tsv
MANIFESTFILE_R3 ?= $(TEMPDIR)/finngen_manifest_r3.tsv
MANIFESTFILE_R4 ?= $(TEMPDIR)/finngen_manifest_r4.tsv
MANIFESTFILE_R5 ?= $(TEMPDIR)/finngen_manifest_r5.tsv
MANIFESTFILE_R6 ?= $(TEMPDIR)/finngen_manifest_r6.tsv
ONTOMARESULTS_R2 ?= $(TEMPDIR)/ontoma_results_r2.tsv
ONTOMARESULTS_R3 ?= $(TEMPDIR)/ontoma_results_r3.tsv
ONTOMARESULTS_R4 ?= $(TEMPDIR)/ontoma_results_r4.tsv
ONTOMARESULTS_R5 ?= $(TEMPDIR)/ontoma_results_r5.tsv
ONTOMARESULTS_R6 ?= $(TEMPDIR)/ontoma_results_r6.tsv
PHENOTYPES_R2 ?= $(TEMPDIR)/finngen_phenotypes_r2.tsv
PHENOTYPES_R3 ?= $(TEMPDIR)/finngen_phenotypes_r3.tsv
PHENOTYPES_R4 ?= $(TEMPDIR)/finngen_phenotypes_r4.tsv
PHENOTYPES_R5 ?= $(TEMPDIR)/finngen_phenotypes_r5.tsv
PHENOTYPES_R6 ?= $(TEMPDIR)/finngen_phenotypes_r6.tsv

#######################


#### Phony targets
.PHONY: R-deps

# ALL
all: R-deps create-temp

# DEPENDENCIES
R-deps:
#cran
	$(R) --slave -e "library('tidyverse')" >/dev/null 2>&1 || \
	$(R) -e "install.packages(c('tidyverse'), repos=c('$(RMIRROR)'))"

# CREATES TEMPORARY DIRECTORY
create-temp:
	mkdir -p $(TEMPDIR)

test: $(MAPPINGSFILE)

manifest: $(MANIFESTFILE_R6)

$(EFOFILE):
	$(CURL) -L https://github.com/EBISPOT/efo/releases/download/$(EFORELEASETAG)/efo_otar_slim.owl > $@

$(MANIFESTFILE_R2): create-temp
	$(CURL) https://storage.googleapis.com/finngen-public-data-r2/summary_stats/r2_manifest.tsv > $@

$(MANIFESTFILE_R3): create-temp
	$(CURL) https://storage.googleapis.com/finngen-public-data-r3/summary_stats/r3_manifest.tsv > $@

$(MANIFESTFILE_R4): create-temp
	$(CURL) https://storage.googleapis.com/finngen-public-data-r4/summary_stats/R4_manifest.tsv > $@

$(MANIFESTFILE_R5): create-temp
	$(CURL) https://storage.googleapis.com/finngen-public-data-r5/summary_stats/R5_manifest.tsv > $@

$(MANIFESTFILE_R6):
	$(CURL) https://storage.googleapis.com/finngen-public-data-r6/summary_stats/R6_manifest.tsv > $@

$(PHENOTYPES_R2): create-temp
	$(CURL) http://r2.finngen.fi/api/phenos | $(JQ) -r '.[]| @json' | $(JQ) -r '[.phenocode, .phenostring] | @tsv' > $@

$(PHENOTYPES_R3): create-temp
	$(CURL) http://r3.finngen.fi/api/phenos  | $(JQ) -r '.[]| @json' | $(JQ) -r '[.phenocode, .phenostring] | @tsv' > $@

$(PHENOTYPES_R4): create-temp
	$(CURL) http://r4.finngen.fi/api/phenos  | $(JQ) -r '.[]| @json' | $(JQ) -r '[.phenocode, .phenostring] | @tsv' > $@

$(PHENOTYPES_R5): create-temp
	$(CURL) https://r5.finngen.fi/api/phenos  | $(JQ) -r '.[]| @json' | $(JQ) -r '[.phenocode, .phenostring] | @tsv' > $@

$(PHENOTYPES_R6): create-temp
	$(CURL) https://r6.finngen.fi/api/phenos  | $(JQ) -r '.[]| @json' | $(JQ) -r '[.phenocode, .phenostring] | @tsv' > $@

$(ONTOMARESULTS_R2): create-temp $(MANIFESTFILE_R2)
	$(CUT) -f2 $(MANIFESTFILE_R3) | $(SED) '/^$$/d' | $(ONTOMA) - $@

$(ONTOMARESULTS_R3): create-temp $(MANIFESTFILE_R3)
	$(CUT) -f2 $(MANIFESTFILE_R3) | $(SED) '/^$$/d' | $(ONTOMA) - $@

$(ONTOMARESULTS_R4): create-temp $(MANIFESTFILE_R4)
	$(CUT) -f2 $(MANIFESTFILE_R4) | $(SED) '/^$$/d' | $(ONTOMA) - $@

$(ONTOMARESULTS_R5): create-temp $(MANIFESTFILE_R5)
	$(CUT) -f2 $(MANIFESTFILE_R5) | $(ONTOMA) -s - $@

$(ONTOMARESULTS_R6): $(MANIFESTFILE_R6)
	$(CUT) -f2 $(MANIFESTFILE_R6) | \
	$(SED) '/^$$/d' | \
	$(ONTOMA) \
	--cache-dir $(TEMPDIR) \
	--outfile $@

$(MAPPINGSFILE): $(EFOFILE)
	$(ROBOT) query \
	--input $(TEMPDIR)/efo_otar_slim.owl \
	--query $(SRCDIR)/efo_all_mappings.sparql \
	$@

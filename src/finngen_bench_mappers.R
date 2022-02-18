library(tidyverse)
library(fs)


#manifest <- read_tsv("../temp/finngen_manifest.tsv")
phenotypes <- read_tsv("../temp/finngen_phenotypes_r4.tsv",
                       col_names = c("phenocode", "phenostring"))


ontoma <- read_tsv("../temp/ontoma_results_r4.tsv")

zooma <- read_csv("../data/zooma_finngen_r4.csv") %>%
  setNames(c("phenostring", "label", "score", "term")) %>%
  filter(score != "Did not map")

manual_mappings <- read_csv("../data/finngen_df3_efo_manual_mappings.csv")

## mk_mappings_old <- read_tsv("../data/finngen_r4_mk_mapping.tsv",
##                         col_names = c("queryId", "query", "term", "label", "score"))

mk_mappings <- jsonlite::stream_in(file("~/Downloads/mapped_traits_v9.json")) %>%
  rename(queryId = traitId, query = traitName, term = efoId, label = efoName)


predictions <- bind_rows(
  phenotypes %>%
  left_join(mk_mappings %>%
            filter(score == 1) %>%
            select(phenocode = queryId, term, label) %>%
            distinct(),
            by = c("phenocode")) %>%
  mutate(method = "ot_score_1"),
  phenotypes %>%
  left_join(mk_mappings %>%
            filter(score >= 0.996) %>%
            select(phenocode = queryId, term, label) %>%
            distinct(),
            by = c("phenocode")) %>%
  mutate(method = "ot_score_0.996"),
  phenotypes %>%
  left_join(mk_mappings %>%
            filter(score >= 0.95) %>%
            select(phenocode = queryId, term, label) %>%
            distinct(),
            by = c("phenocode")) %>%
  mutate(method = "ot_score_0.95"),
  phenotypes %>%
  left_join(mk_mappings %>%
            filter(score >= 0.8) %>%
            select(phenocode = queryId, term, label) %>%
            distinct(),
            by = c("phenocode")) %>%
  mutate(method = "ot_score_0.8"),
  phenotypes %>%
  left_join(mk_mappings %>%
            filter(score >= 0.7) %>%
            select(phenocode = queryId, term, label) %>%
            distinct(),
            by = c("phenocode")) %>%
  mutate(method = "ot_score_0.7"),
  phenotypes %>%
  left_join(mk_mappings %>%
            filter(score >= 0.6) %>%
            select(phenocode = queryId, term, label) %>%
            distinct(),
            by = c("phenocode")) %>%
  mutate(method = "ot_score_0.6"),
  phenotypes %>%
  left_join(mk_mappings %>%
            filter(score >= 0.5) %>%
            select(phenocode = queryId, term, label) %>%
            distinct(),
            by = c("phenocode")) %>%
  mutate(method = "ot_score_0.5"),
  phenotypes %>%
  left_join(mk_mappings %>%
            filter(score >= 0.4) %>%
            select(phenocode = queryId, term, label) %>%
            distinct(),
            by = c("phenocode")) %>%
  mutate(method = "ot_score_0.4"),
    phenotypes %>%
  left_join(mk_mappings %>%
            filter(score >= 0.3) %>%
            select(phenocode = queryId, term, label) %>%
            distinct(),
            by = c("phenocode")) %>%
  mutate(method = "ot_score_0.3"),
  phenotypes %>%
  left_join(mk_mappings %>%
            filter(score >= 0.2) %>%
            select(phenocode = queryId, term, label) %>%
            distinct(),
            by = c("phenocode")) %>%
  mutate(method = "ot_score_0.2"),
  phenotypes %>%
  left_join(mk_mappings %>%
            filter(score >= 0.1) %>%
            select(phenocode = queryId, term, label) %>%
            distinct(),
            by = c("phenocode")) %>%
  mutate(method = "ot_score_0.1"),
  phenotypes %>%
  left_join(zooma %>%
            filter(score == "Good") %>%
            select(phenostring, term, label) %>%
            distinct(),
            by = c("phenostring")) %>%
  mutate(method = "zooma_good"),
  phenotypes %>%
  left_join(zooma %>%
            select(phenostring, term, label) %>%
            distinct(),
            by = c("phenostring")) %>%
  mutate(method = "zooma_any"),
  phenotypes %>%
  left_join(ontoma %>%
            filter(quality == "match") %>%
            mutate(term = path_file(term)) %>%
            select(phenostring = query, term, label) %>%
            distinct(),
            by = c("phenostring")) %>%
  mutate(method = "ontoma_match"),
  phenotypes %>%
  left_join(ontoma %>%
            mutate(term = path_file(term)) %>%
            select(phenostring = query, term, label) %>%
            distinct(),
            by = c("phenostring")) %>%
  mutate(method = "ontoma_any")
) %>%
  distinct()

bench <- predictions %>%
  ## filtering to universe of R2/R3
  inner_join(manual_mappings %>%
             filter(valid) %>%
             select(phenocode = NAME) %>%
             distinct(),
             by = "phenocode") %>%
  ## adding manual mappings as gold
  left_join(manual_mappings %>%
            filter(valid) %>%
            mutate(efo_cls = path_file(efo_cls)) %>%
            select(phenocode = NAME, term = efo_cls, valid) %>%
            distinct(),
            by = c("phenocode", "term")) %>%
  mutate(valid = replace_na(valid, FALSE)) %>%
  mutate(predicted = !is.na(label))

bench %>%
  group_by(method) %>%
  summarise(valid = sum(valid),
            predicted = sum(predicted),
            precision = sum(valid) / sum(predicted),
            recall = sum(valid) / n()) %>%
  rowwise() %>%
  mutate(F = mean(c(precision, recall))) %>%
  arrange(desc(precision))

## failExamples <- bench %>%
##   filter(method == "ot_score_1") %>%
##   filter(!valid) %>%
##   filter(predicted) %>%
##   left_join(
##     manual_mappings %>%
##     filter(valid) %>%
##     mutate(efo_cls = path_file(efo_cls)) %>%
##     select(phenocode = NAME, sandraterm = efo_cls, sandralabel = efo_label),
##     by = c("phenocode")) %>%
##   filter(!is.na(sandralabel)) %>%
##   select(-valid, -predicted, -method)

## failExamples <- bench %>%
##   filter(method == "mk_score_any") %>%
##   filter(!valid) %>%
##   filter(predicted) %>%
##   select(-valid, -predicted, -method) %>%
##   inner_join(
##     bench %>%
##     filter(method == "ontoma_any") %>%
##     filter(valid) %>%
##     filter(predicted) %>%
##     rename(ontomaterm = term, ontomalabel = label) %>%
##     select(-valid, -predicted, -method),
##     by = c("phenocode", "phenostring")
##   )


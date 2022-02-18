library(tidyverse)
library(fs)

############################################
## Load files
############################################

### endpoints and manifest from FinnGen
endpoints <- read_csv("../data/FINNGEN_ENDPOINTS_DF4_2019-06-27_Public.csv",
                      guess_max = 4000,
                      comment = "")
#manifest <- read_tsv("../temp/finngen_manifest.tsv")
phenotypes <- read_tsv("../temp/finngen_phenotypes_r4.tsv",
                       col_names = c("phenocode", "phenostring"))

### Mappings extracted from EFO
ICDmappings <- read_delim("../temp/all_mappings.tsv",
                          delim = "\t",
                          comment = "@en",
                          escape_backslash = FALSE,
                          escape_double = FALSE
                          ) %>%
    setNames(c("xref", "cls", "label")) %>%
    filter(str_detect(xref, "ICD")) %>%
    separate(xref, c("xref_source", "xref"), ":") %>%
    filter(xref_source != "ICDO") %>%
    mutate(xref_source = str_replace(xref_source, "ICD-10|ICD10CM|ICD10", "ICD_10")) %>%
    mutate(xref_source = str_replace(xref_source, "ICD9CM|ICD9", "ICD_9")) %>%
    mutate(xref = str_replace(xref, "\\.", "")) %>%
    mutate(cls = str_replace_all(cls, "[<>]", ""))

zooma <- read_tsv("../temp/ontoma_results_r4.tsv")

# up to r3
manual_mappings <- read_csv("../data/finngen_df3_efo_manual_mappings.csv")

mk_mappings <- jsonlite::stream_in(file("~/Downloads/mapped_traits_v9.json")) %>%
  rename(queryId = traitId, query = traitName, term = efoId, label = efoName)

############################################
## ALL ICD9/10 mapping (messy with regexp)
############################################

icdMap_temp <- endpoints %>%
    filter(row_number() != 1) %>%
    select(NAME, LONGNAME, HD_ICD_10, HD_ICD_9) %>%
    gather("xref_source", "xref", -NAME, -LONGNAME) %>%
    mutate(xref_source = str_replace(xref_source, "HD_", "")) %>%
    filter(!is.na(xref) & xref != "$!$") %>%
    ## simplifying xrefs
    mutate(simple_xref = xref) %>%
    mutate(simple_xref = str_replace_all(simple_xref, "%", "")) %>%
    mutate(simple_xref = str_replace_all(simple_xref, " ", "")) %>%
    mutate(simple_xref = str_replace_all(simple_xref,
                                         "\\[([A-Z0-9]+?)\\|([A-Z0-9]+?)\\]",
                                         "\\[\\1<or>\\2\\]")) %>%
    mutate(simple_xref = str_replace_all(simple_xref,
                                         "\\[([A-Z0-9]+?)\\-([A-Z0-9]+?)\\]",
                                         "\\[\\1<to>\\2\\]")) %>%
    mutate(simple_xref = str_split(simple_xref, "[\\|&]")) %>%
    unnest(simple_xref) %>%
    distinct() %>%
    mutate(clause = str_replace_all(simple_xref, "(?:[A-Z0-9]+(\\[[0-9A-Z]+((<to>)|(<or>))[0-9A-Z]+\\])?[A-Z0-9]*)", "\\1")) %>%
    mutate(clause = str_replace_all(clause, "[\\[\\]]", "")) %>%
    mutate(clause = str_split(clause, "((<to>)|(<or>))")) %>%
    mutate(clause = ifelse(lengths(clause) > 2, simple_xref, clause)) %>%
    mutate(clause_label = ifelse(lengths(clause) == 2, list(c("start", "end")), "")) %>%
    unnest() %>%
    spread(clause_label, clause) %>%
    mutate(end = as.numeric(str_replace(end, "\\*", "")))
partA <- icdMap_temp %>%
    filter(str_count(simple_xref, "or") == 1) %>%
    rowwise() %>%
    mutate(vector = list(c(start, end))) %>%
    ungroup()
partB <- icdMap_temp %>%
    filter(str_count(simple_xref, "to") == 1 & start != "A") %>%
    rowwise() %>%
    mutate(vector = list(as.character(seq(start, end)))) %>%
    ungroup()
partC <- icdMap_temp %>%
    filter(str_count(simple_xref, "to") != 1 & start == "A" & str_count(simple_xref, "or") != 1) %>%
    rowwise() %>%
    mutate(vector = list(c(NA))) %>%
    ungroup()
icdMap <- bind_rows(partA, partB, partC) %>%
    unnest(vector) %>%
    mutate(simple_xref = str_replace(simple_xref,
                                     "\\[[0-9A-Z]+((<to>)|(<or>))[0-9A-Z]+\\]",
                                     vector)) %>%
    distinct() %>%
    left_join(ICDmappings, by = c("simple_xref" = "xref", "xref_source")) %>%
    select(-V1, -vector, -start, -end)

###########################################
### create combined
###########################################

result <- endpoints %>%
  select(NAME, LONGNAME) %>%
  inner_join(phenotypes %>% select(phenocode),
             by = c("NAME" = "phenocode")) %>%
  left_join(icdMap %>%
            select(NAME, xref_source, xref, simple_xref, cls, label),
            by = c("NAME")) %>%
  bind_rows(zooma %>%
            filter(!is.na(term)) %>%
            mutate(source = str_replace_all(source, " ", "_")) %>%
            select(LONGNAME = query,
                   cls = term,
                   label,
                   xref_source = source) %>%
            left_join(phenotypes %>%
                      select(NAME = phenocode, LONGNAME = phenostring),
                      by = "LONGNAME")) %>%
  bind_rows(
    mk_mappings %>%
    rename(NAME = queryId,
           LONGNAME = query,
           cls = term,
           quality = score) %>%
    mutate(xref_source = "ot_mappings")) %>%
  mutate(label = str_to_lower(label)) %>%
  mutate(cls = path_file(cls)) %>%
  group_by(NAME, LONGNAME, cls, label) %>%
  summarize_all(~paste(unique(na.omit(.)), collapse = ';')) %>%
  ungroup() %>%
  group_by(NAME, LONGNAME) %>%
  mutate(size = n()) %>%
  ungroup() %>%
  filter(size == 1 | (!is.na(cls) & !is.na(label) & !is.na(xref_source))) %>%
  rename(efo_cls = cls) %>%
  rename(efo_label = label) %>%
  select(-size)

#######################################
## Combine with manual file from DF2
#######################################

manual <- manual_mappings %>%
  filter(valid) %>%
  mutate(efo_cls = path_file(efo_cls))

merged <- result %>%
    filter((!NAME %in% unique(manual$NAME))) %>%
    mutate(valid = FALSE) %>%
    mutate(release = "r4") %>%
    bind_rows(manual) %>%
    arrange(NAME)


merged %>%
    write_csv("~/Desktop/df4_test.csv")


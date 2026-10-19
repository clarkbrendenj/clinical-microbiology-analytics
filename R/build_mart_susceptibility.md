Build the ARMD Susceptibility Mart with R
================
Brenden Clark

- [1. Import and stage the source
  data](#1-import-and-stage-the-source-data)
- [2. Consolidate culture records and reconcile AST
  results](#2-consolidate-culture-records-and-reconcile-ast-results)
- [3. Apply reviewed organism
  mappings](#3-apply-reviewed-organism-mappings)
- [4. Classify resistance phenotypes from AST
  evidence](#4-classify-resistance-phenotypes-from-ast-evidence)
- [5. Build reporting cultures and reconcile their
  AST](#5-build-reporting-cultures-and-reconcile-their-ast)
- [6. Select the first qualifying culture within each reporting
  category](#6-select-the-first-qualifying-culture-within-each-reporting-category)
- [7. Apply antibiotic panels and build the susceptibility
  mart](#7-apply-antibiotic-panels-and-build-the-susceptibility-mart)
- [8. Validate against the Tableau CSV and export the R
  mart](#8-validate-against-the-tableau-csv-and-export-the-r-mart)

This report transforms the original ARMD culture data into the cleaned
susceptibility dataset used by the Tableau dashboard.

## 1. Import and stage the source data

### 1.1 Load packages and define the source

Load the required packages and specify the original cultures CSV and
ARMD release version.

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.2     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.4     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(lubridate)

raw_data_file = "../data/raw/microbiology_cultures_cohort.csv"
source_release = "2025-04-11"
```

### 1.2 Define the reporting period

Use January 1, 2023 through December 31, 2023 to match the Tableau
dataset. The reporting window will be applied after resolving
culture-level fields.

``` r
analysis_window_start = as.Date("2023-01-01")
analysis_window_end_exclusive = analysis_window_start %m+% months(12)
```

### 1.3 Import the raw data

Import all fields as text to preserve identifiers and control data types
explicitly. Convert blank CSV fields to missing values while preserving
the original susceptibility text.

``` r
raw_armd =
  read_csv(
    raw_data_file,
    col_types = cols(.default = col_character()),
    na = "",
    trim_ws = FALSE
  ) %>%
  janitor::clean_names()
```

### 1.4 Clean and standardize fields

Rename identifiers, parse timestamps, trim text, and handle the source’s
literal “Null” values. Retain all source rows, including records with
missing AST results, for subsequent reconciliation.

``` r
armd_staged =
  raw_armd %>%
  transmute(
    source_version = source_release,
    patient_key = anon_id,
    encounter_key = pat_enc_csn_id_coded,
    culture_order_key = order_proc_id_coded,
    order_time_jittered_utc_raw = order_time_jittered_utc,
    culture_time_shifted = ymd_hms(
      order_time_jittered_utc,
      tz = "UTC",
      quiet = TRUE
    ),
    ordering_mode = na_if(str_trim(ordering_mode), "Null"),
    culture_description = na_if(str_trim(culture_description), "Null"),
    was_positive = as.integer(na_if(str_trim(was_positive), "Null")),
    organism = na_if(str_trim(organism), "Null"),
    antibiotic = na_if(str_trim(antibiotic), "Null"),
    susceptibility_raw = susceptibility,
    susceptibility = str_trim(na_if(susceptibility, "Null"))
  )
```

### 1.5 Validate the import

Check the row count, timestamp parsing, missing identifiers, and
specimen categories. ARMD timestamps are shifted research dates.

``` r
import_check =
  armd_staged %>%
  summarise(
    source_rows = n(),
    timestamp_parse_failures = sum(
      !is.na(order_time_jittered_utc_raw) &
        !(order_time_jittered_utc_raw %in% c("", "Null")) &
        is.na(culture_time_shifted)
    ),
    missing_patient_keys = sum(is.na(patient_key)),
    missing_order_keys = sum(is.na(culture_order_key))
  )

print(import_check)
```

    ## # A tibble: 1 × 4
    ##   source_rows timestamp_parse_failures missing_patient_keys missing_order_keys
    ##         <int>                    <int>                <int>              <int>
    ## 1     2241050                        0                    0                  0

``` r
glimpse(armd_staged)
```

    ## Rows: 2,241,050
    ## Columns: 13
    ## $ source_version              <chr> "2025-04-11", "2025-04-11", "2025-04-11", …
    ## $ patient_key                 <chr> "JC2744063", "JC1713666", "JC1669304", "JC…
    ## $ encounter_key               <chr> "131368600230", "131300064625", "131272997…
    ## $ culture_order_key           <chr> "928257722", "697566032", "620809641", "51…
    ## $ order_time_jittered_utc_raw <chr> "2023-12-23 22:29:00+00:00", "2020-12-27 0…
    ## $ culture_time_shifted        <dttm> 2023-12-23 22:29:00, 2020-12-27 00:40:00,…
    ## $ ordering_mode               <chr> "Inpatient", "Inpatient", "Inpatient", "In…
    ## $ culture_description         <chr> "URINE", "URINE", "BLOOD", "URINE", "URINE…
    ## $ was_positive                <int> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, …
    ## $ organism                    <chr> "KLEBSIELLA PNEUMONIAE", "KLEBSIELLA PNEUM…
    ## $ antibiotic                  <chr> "Ertapenem", "Ertapenem", "Ertapenem", "Er…
    ## $ susceptibility_raw          <chr> "Susceptible", "Susceptible", "Susceptible…
    ## $ susceptibility              <chr> "Susceptible", "Susceptible", "Susceptible…

``` r
armd_staged %>%
  count(culture_description, sort = TRUE) %>%
  print(n = Inf)
```

    ## # A tibble: 3 × 2
    ##   culture_description       n
    ##   <chr>                 <int>
    ## 1 URINE               1516950
    ## 2 BLOOD                435774
    ## 3 RESPIRATORY          288326

## 2. Consolidate culture records and reconcile AST results

The source contains repeated culture details across organism and
antibiotic rows. Create one record per culture, retain the original
organism labels, and reconcile each organism-antibiotic endpoint before
applying organism mappings or resistance classifications.

### 2.1 Create one record per culture

A culture is identified by its source release, patient, encounter, and
order. For each culture, retain an order-level field only when its
source values agree. If a field contains different values, including a
mixture of a value and missing data, set that field to missing and flag
the culture for review. The flag identifies a conflict; it does not
automatically exclude the culture.

``` r
culture_id = c(
  "source_version",
  "patient_key",
  "encounter_key",
  "culture_order_key"
)

organism_id = c(culture_id, "organism")

cultures =
  armd_staged %>%
  group_by(across(all_of(culture_id))) %>%
  summarise(
    source_row_count = n(),
    has_order_level_conflict = n_distinct(
      culture_time_shifted,
      ordering_mode,
      culture_description,
      was_positive
    ) > 1,
    across(
      c(
        culture_time_shifted,
        ordering_mode,
        culture_description,
        was_positive
      ),
      ~ if (n_distinct(.x) == 1) first(.x) else .x[NA_integer_]
    ),
    .groups = "drop"
  ) %>%
  mutate(
    cohort_id = case_when(
      culture_description == "BLOOD" ~ "blood",
      culture_description == "URINE" ~ "urine",
      culture_description == "RESPIRATORY" ~ "respiratory",
      TRUE ~ NA_character_
    )
  )
```

`across()` applies the same consistency check to each of the four
culture fields. In the conflicting case, `.x[NA_integer_]` returns a
missing value with the field’s original type, preserving date-time,
character, and integer types. `n_distinct()` includes missing values in
its count by default.

### 2.2 Retain the original organism observations

Create one record per culture and original organism label, then attach
the consolidated culture details. These represent source organism
observations; ARMD does not supply verified biological isolate
identifiers. The explicit join keys prevent unrelated columns from
affecting the match.

``` r
culture_organisms =
  armd_staged %>%
  filter(!is.na(organism)) %>%
  select(all_of(organism_id)) %>%
  distinct() %>%
  inner_join(
    cultures,
    by = culture_id,
    na_matches = "never",
    relationship = "many-to-one"
  )
```

`na_matches = "never"` prevents missing identifiers from matching one
another, consistent with the existing SQL joins.
`relationship = "many-to-one"` requires each organism observation to
match at most one culture record.

### 2.3 Reconcile repeated AST results

Group antibiotic records within each culture and original organism
label. Repeated identical results resolve to one result. Different
results, including a valid result combined with a missing value, remain
conflicts. An endpoint containing only missing values retains a missing
status.

| Source values within an endpoint | Status   | Retained result |
|----------------------------------|----------|-----------------|
| Susceptible, Susceptible         | resolved | Susceptible     |
| Susceptible, Resistant           | conflict | Missing         |
| Susceptible, missing             | conflict | Missing         |
| Missing, missing                 | missing  | Missing         |

Rows without an antibiotic name do not define an AST endpoint. They
remain represented in the culture and organism data frames above.

``` r
ast_endpoints =
  armd_staged %>%
  filter(
    !is.na(organism),
    !is.na(antibiotic)
  ) %>%
  group_by(across(all_of(c(organism_id, "antibiotic")))) %>%
  summarise(
    source_result_rows = n(),
    missing_result_rows = sum(is.na(susceptibility)),
    result_state_count = n_distinct(susceptibility),
    candidate_result = first(susceptibility),
    .groups = "drop"
  ) %>%
  mutate(
    resolution_status = case_when(
      result_state_count > 1 ~ "conflict",
      is.na(candidate_result) ~ "missing",
      TRUE ~ "resolved"
    ),
    susceptibility = if_else(
      resolution_status == "resolved",
      candidate_result,
      NA_character_
    ),
    repeated_result_rows = source_result_rows - result_state_count
  ) %>%
  select(-candidate_result) %>%
  semi_join(
    culture_organisms,
    by = organism_id,
    na_matches = "never"
  )
```

The candidate result is retained only when all source values agree, so
the resolved output is independent of source row order. The
`semi_join()` keeps endpoints linked to an existing organism observation
without adding columns or multiplying rows.

### 2.4 Inspect culture and AST summaries

Check how many cultures have conflicting fields, whether encounter
identifiers are missing, and how many AST endpoints are resolved,
missing, or conflicting. The source row total in `culture_check` should
match the imported 2,241,050 rows. Culture counts are expected to be
smaller because each source culture may have many organism and
antibiotic rows.

``` r
culture_check =
  cultures %>%
  summarise(
    culture_count = n(),
    source_rows_accounted_for = sum(source_row_count),
    cultures_with_conflicting_fields = sum(has_order_level_conflict),
    missing_encounter_keys = sum(is.na(encounter_key)),
    missing_culture_times = sum(is.na(culture_time_shifted))
  )

ast_check =
  ast_endpoints %>%
  count(resolution_status, name = "endpoint_count")

print(culture_check, width = Inf)
```

    ## # A tibble: 1 × 5
    ##   culture_count source_rows_accounted_for cultures_with_conflicting_fields
    ##           <int>                     <int>                            <int>
    ## 1        751075                   2241050                                0
    ##   missing_encounter_keys missing_culture_times
    ##                    <int>                 <int>
    ## 1                      0                     0

``` r
print(ast_check)
```

    ## # A tibble: 3 × 2
    ##   resolution_status endpoint_count
    ##   <chr>                      <int>
    ## 1 conflict                    9047
    ## 2 missing                     1991
    ## 3 resolved                 1588327

``` r
culture_organisms %>%
  count(cohort_id, name = "culture_organism_count", sort = TRUE) %>%
  print(n = Inf)
```

    ## # A tibble: 3 × 2
    ##   cohort_id   culture_organism_count
    ##   <chr>                        <int>
    ## 1 urine                        97184
    ## 2 respiratory                  22411
    ## 3 blood                        13651

The resulting `cultures`, `culture_organisms`, and `ast_endpoints` data
frames will support organism mapping and phenotype classification in the
next stage.

## 3. Apply reviewed organism mappings

The source may describe the same organism using several labels,
including older names and resistance or morphotype qualifiers. Join the
existing dbt reference tables to attach standardized identities and
reporting groups. Using the same reference files keeps the R and SQL
reporting definitions aligned.

Retain the original `organism` field and one record per original source
label on each culture. Resistance will be evaluated for those
observations before labels are consolidated into reporting populations.

### 3.1 Load the organism reference tables

Use the organism dictionary and Gram classification CSVs already stored
in `dbt/seeds/`. The relative path below assumes the same working
directory as your existing `../data/raw/` path. Change `seed_dir` if
necessary.

``` r
library(tidyverse)
library(lubridate)

seed_dir = "../dbt/seeds"

organism_mapping =
  read_csv(
    file.path(seed_dir, "organism_mapping.csv"),
    col_types = cols(
      .default = col_character(),
      include_primary_summary = col_integer()
    ),
    na = "",
    trim_ws = FALSE
  ) %>%
  mutate(
    across(
      c(organism_reporting_group, source_qualifier),
      ~ na_if(str_trim(.x), "")
    ),
    has_reviewed_mapping = TRUE
  )

organism_reporting_groups =
  read_csv(
    file.path(seed_dir, "organism_reporting_groups.csv"),
    col_types = cols(.default = col_character()),
    na = "",
    trim_ws = FALSE
  )
```

`organism_standardized` identifies the reviewed organism name.
`organism_reporting_group` identifies the group used in the final
report; for example, related source identifications can share a
reporting complex. `source_qualifier` retains information such as a
source MRSA designation. Those source qualifiers remain traceable, but
AST evidence will determine the resistance classifications used by this
pipeline.

### 3.2 Attach standardized identities and Gram groups

Match each original organism label exactly to the dictionary, then
attach the reporting group’s Gram classification. The `left_join()`
operations retain unmatched observations for review. Each reference
table must supply at most one match, preventing the joins from
multiplying records.

``` r
mapped_organisms =
  culture_organisms %>%
  filter(!is.na(cohort_id)) %>%
  left_join(
    organism_mapping,
    by = c("organism" = "organism_source"),
    na_matches = "never",
    relationship = "many-to-one"
  ) %>%
  mutate(
    has_reviewed_mapping = coalesce(has_reviewed_mapping, FALSE),
    mapping_status = if_else(
      has_reviewed_mapping,
      mapping_status,
      "unmapped"
    ),
    include_primary_summary = coalesce(include_primary_summary, 0L)
  ) %>%
  left_join(
    organism_reporting_groups,
    by = "organism_reporting_group",
    na_matches = "never",
    relationship = "many-to-one"
  )
```

`has_reviewed_mapping` indicates whether the source label matched the
dictionary. The `include_primary_summary` flag records the reference
table’s reporting policy; it does not yet determine the final patient
population. No organism observations are consolidated at this stage.

### 3.3 Review mapping coverage

The previous stage contained 133,246 culture-organism observations
across the three cohorts. This mapping step should retain that count.
Review the unmatched labels and the records enabled for primary
reporting that lack a Gram classification. An unmatched label may be
outside the intended report scope; review its identity and reporting
purpose before extending the dictionary.

``` r
mapping_check =
  mapped_organisms %>%
  summarise(
    culture_organism_records = n(),
    distinct_source_labels = n_distinct(organism),
    matched_records = sum(has_reviewed_mapping),
    unmapped_records = sum(!has_reviewed_mapping),
    primary_reporting_records = sum(include_primary_summary == 1L),
    primary_records_missing_gram = sum(
      include_primary_summary == 1L & is.na(gram_stain_group)
    )
  )

unmapped_labels =
  mapped_organisms %>%
  filter(!has_reviewed_mapping) %>%
  count(organism, name = "culture_organism_count", sort = TRUE)

print(mapping_check, width = Inf)
```

    ## # A tibble: 1 × 6
    ##   culture_organism_records distinct_source_labels matched_records
    ##                      <int>                  <int>           <int>
    ## 1                   133246                    314          126405
    ##   unmapped_records primary_reporting_records primary_records_missing_gram
    ##              <int>                     <int>                        <int>
    ## 1             6841                    126405                            0

``` r
print(unmapped_labels, n = 20)
```

    ## # A tibble: 273 × 2
    ##    organism                               culture_organism_count
    ##    <chr>                                                   <int>
    ##  1 MYCOBACTERIUM AVIUM COMPLEX                               799
    ##  2 ACHROMOBACTER XYLOSOXIDANS                                489
    ##  3 STREPTOCOCCUS ANGINOSUS GROUP                             314
    ##  4 STREPTOCOCCUS MITIS GROUP                                 269
    ##  5 MYCOBACTERIUM ABSCESSUS                                   164
    ##  6 ACINETOBACTER SPECIES                                     152
    ##  7 CORYNEBACTERIUM PROPINQUUM                                141
    ##  8 PROTEUS VULGARIS GROUP                                    137
    ##  9 STAPHYLOCOCCUS LUGDUNENSIS                                130
    ## 10 CORYNEBACTERIUM PSEUDODIPHTHERITICUM                      128
    ## 11 VIRIDANS GROUP STREPTOCOCCI                               126
    ## 12 STREPTOCOCCUS PYOGENES (GROUP A)                          122
    ## 13 PROVIDENCIA RETTGERI                                      109
    ## 14 MYCOBACTERIUM CHELONAE                                    105
    ## 15 ACINETOBACTER BAUMANNII COMPLEX                            96
    ## 16 BACTEROIDES FRAGILIS                                       96
    ## 17 STREPTOCOCCUS DYSGALACTIAE                                 85
    ## 18 PSEUDOMONAS SPECIES - NOT P.AERUGINOSA                     83
    ## 19 CORYNEBACTERIUM SPECIES                                    82
    ## 20 PROVIDENCIA STUARTII                                       73
    ## # ℹ 253 more rows

``` r
mapped_organisms %>%
  count(
    mapping_status,
    include_primary_summary,
    name = "culture_organism_count",
    sort = TRUE
  ) %>%
  print(n = Inf)
```

    ## # A tibble: 2 × 3
    ##   mapping_status include_primary_summary culture_organism_count
    ##   <chr>                            <int>                  <int>
    ## 1 reviewed                             1                 126405
    ## 2 unmapped                             0                   6841

The resulting `mapped_organisms` data frame retains the culture
identifiers and original labels alongside standardized organism
information. The next stage will combine these observations with
`ast_endpoints` to evaluate resistance phenotypes.

## 4. Classify resistance phenotypes from AST evidence

Use the reconciled antibiotic results to evaluate CRE, extended-spectrum
cephalosporin resistance (ESCE-R), oxacillin resistance in *S. aureus*,
vancomycin resistance in *E. faecalis* and *E. faecium*, and multidrug
resistance in *P. aeruginosa*. The rules below reproduce the existing
project’s definitions using a small drug lookup and shared
transformations.

The classification remains at the culture and original organism-label
level. Source name qualifiers remain available for tracing records; the
resistance classifications use AST evidence.

### 4.1 Select the reviewed reporting scope and identify applicable rules

The SQL applies its reviewed reporting-scope filter after
classification. Here, apply the same filter first: classification uses
each original organism observation independently, so the excluded
observations cannot change the retained observations’ evidence. The
6,841 unmatched observations remain available in `mapped_organisms` and
`unmapped_labels` for review.

For the supplied dictionary’s 41 source labels, eligibility based on
standardized identities agrees with the SQL’s standardized-name and
source-label fallback checks. This lets the R code use the standardized
identities directly. Recheck this equivalence if the dictionary is
extended.

``` r
library(tidyverse)
library(lubridate)

organism_scope =
  mapped_organisms %>%
  filter(
    include_primary_summary == 1L,
    mapping_status == "reviewed",
    !is.na(organism_reporting_group)
  ) %>%
  mutate(
    standardized_name = coalesce(organism_standardized, ""),
    cre_standard =
      str_detect(
        standardized_name,
        "^(Citrobacter |Enterobacter |Klebsiella )"
      ) |
      standardized_name %in% c("Escherichia coli", "Serratia marcescens"),
    cre_proteus_morganella =
      str_starts(standardized_name, "Proteus ") |
      standardized_name == "Morganella morganii",
    esce =
      standardized_name %in% c(
        "Escherichia coli",
        "Klebsiella aerogenes",
        "Klebsiella oxytoca",
        "Klebsiella pneumoniae"
      ) |
      str_starts(standardized_name, "Enterobacter "),
    oxacillin = standardized_name == "Staphylococcus aureus",
    vre = standardized_name %in% c(
      "Enterococcus faecalis",
      "Enterococcus faecium"
    ),
    mdr = standardized_name == "Pseudomonas aeruginosa"
  ) %>%
  select(-standardized_name)

rule_eligibility =
  organism_scope %>%
  select(
    all_of(organism_id),
    cre_standard,
    cre_proteus_morganella,
    esce,
    oxacillin,
    vre,
    mdr
  ) %>%
  pivot_longer(
    cols = -all_of(organism_id),
    names_to = "rule",
    values_to = "eligible"
  ) %>%
  filter(eligible) %>%
  select(-eligible)
```

### 4.2 Define the antibiotics used by each rule

Each lookup row describes one rule and antimicrobial category. The pipe
character separates drug names within the compact definition;
`separate_longer_delim()` turns them into one row per antibiotic.

The Proteus/Morganella CRE panel excludes imipenem-based evidence. The
MDR *P. aeruginosa* definition evaluates six categories; its
aminoglycoside category uses tobramycin alone, matching the project’s
SQL.

``` r
phenotype_panels =
  tribble(
    ~rule, ~category, ~antibiotic,
    "cre_standard", "cre",
    "Ertapenem|Imipenem|Meropenem|Meropenem/Vaborbactam|Imipenem/Relebactam",
    "cre_proteus_morganella", "cre",
    "Ertapenem|Meropenem|Meropenem/Vaborbactam",
    "esce", "esce",
    "Cefepime|Ceftriaxone|Cefotaxime|Ceftazidime|Ceftazidime/Avibactam|Ceftolozane/Tazobactam",
    "oxacillin", "oxacillin",
    "Oxacillin|Cefoxitin",
    "vre", "vre",
    "Vancomycin",
    "mdr", "es_ceph",
    "Cefepime|Ceftazidime|Ceftazidime/Avibactam|Ceftolozane/Tazobactam",
    "mdr", "fluoroquinolone",
    "Ciprofloxacin|Levofloxacin",
    "mdr", "aminoglycoside",
    "Tobramycin",
    "mdr", "carbapenem",
    "Imipenem|Meropenem|Imipenem/Relebactam",
    "mdr", "piperacillin_tazobactam",
    "Piperacillin/Tazobactam",
    "mdr", "cefiderocol",
    "Cefiderocol"
  ) %>%
  separate_longer_delim(antibiotic, delim = "|")
```

### 4.3 Summarize the observed AST evidence

Attach the applicable rules to observed antibiotic endpoints. Some drugs
contribute to more than one rule, so the first join intentionally
permits multiple matches in both tables. The second join keeps only
rules that apply to the specific original organism observation.

A categorized endpoint must be resolved and contain Susceptible,
Intermediate, or Resistant. Resistant is positive evidence for the CRE,
ESCE-R, oxacillin, and vancomycin rules. Both Intermediate and Resistant
are positive evidence for the MDR rule.

``` r
rule_evidence =
  ast_endpoints %>%
  inner_join(
    phenotype_panels,
    by = "antibiotic",
    relationship = "many-to-many"
  ) %>%
  inner_join(
    rule_eligibility,
    by = c(organism_id, "rule"),
    na_matches = "never",
    relationship = "many-to-one"
  ) %>%
  mutate(
    is_categorized =
      resolution_status == "resolved" &
      susceptibility %in% c("Susceptible", "Intermediate", "Resistant"),
    is_positive = if_else(
      rule == "mdr",
      susceptibility %in% c("Intermediate", "Resistant"),
      susceptibility %in% "Resistant"
    )
  ) %>%
  group_by(across(all_of(c(organism_id, "rule", "category")))) %>%
  summarise(
    observed_endpoint_count = n(),
    positive_count = sum(is_categorized & is_positive),
    nonpositive_count = sum(is_categorized & !is_positive),
    invalid_count = sum(!is_categorized),
    .groups = "drop"
  )
```

This summarizes observed endpoints only. An untested drug is not
automatically an invalid result. Antibiotics outside these phenotype
panels remain in `ast_endpoints` for susceptibility reporting later.

### 4.4 Assign an evidence status to each applicable category

Create each required rule/category combination, including the six MDR
categories. Categories with no observed endpoints become `not_tested`.
Positive evidence takes precedence over invalid results; otherwise,
invalid results take precedence over observed nonpositive evidence.

``` r
required_categories =
  phenotype_panels %>%
  distinct(rule, category)

phenotype_category_status =
  rule_eligibility %>%
  inner_join(
    required_categories,
    by = "rule",
    relationship = "many-to-many"
  ) %>%
  left_join(
    rule_evidence,
    by = c(organism_id, "rule", "category"),
    na_matches = "never",
    relationship = "one-to-one"
  ) %>%
  mutate(
    across(
      c(
        observed_endpoint_count,
        positive_count,
        nonpositive_count,
        invalid_count
      ),
      ~ replace_na(.x, 0L)
    ),
    category_status = case_when(
      positive_count > 0 ~ "positive",
      invalid_count > 0 ~ "not_evaluable",
      nonpositive_count > 0 ~ "observed_non_positive",
      TRUE ~ "not_tested"
    )
  )
```

### 4.5 Derive the five phenotype classifications

For the four single-panel phenotypes, positive evidence produces a
positive classification. Observed nonpositive evidence without invalid
endpoints produces `not_detected`. Invalid evidence or no usable testing
produces `not_evaluable`.

For MDR *P. aeruginosa*, at least three positive categories establish
positive evidence. If the positive and unknown categories together
cannot reach three, classify the observation as `not_detected`.
Otherwise, retain `not_evaluable`. Unknown categories include both
invalid and untested categories.

``` r
ordinary_phenotype_status =
  phenotype_category_status %>%
  filter(rule != "mdr") %>%
  transmute(
    across(all_of(organism_id)),
    phenotype = case_when(
      category == "cre" ~ "ast_cre_status",
      category == "esce" ~ "ast_esce_enterobacterales_status",
      category == "oxacillin" ~ "ast_oxacillin_resistance_status",
      category == "vre" ~ "ast_vancomycin_resistance_status"
    ),
    status = case_when(
      category_status == "positive" ~ "positive",
      category_status == "observed_non_positive" ~ "not_detected",
      TRUE ~ "not_evaluable"
    )
  )

mdr_summary =
  phenotype_category_status %>%
  filter(rule == "mdr") %>%
  group_by(across(all_of(organism_id))) %>%
  summarise(
    mdr_positive_category_count = sum(category_status == "positive"),
    mdr_unknown_category_count = sum(
      category_status %in% c("not_evaluable", "not_tested")
    ),
    .groups = "drop"
  ) %>%
  mutate(
    ast_mdr_p_aeruginosa_status = case_when(
      mdr_positive_category_count >= 3 ~ "positive",
      mdr_positive_category_count + mdr_unknown_category_count < 3 ~
        "not_detected",
      TRUE ~ "not_evaluable"
    )
  )

phenotype_columns = c(
  "ast_cre_status",
  "ast_esce_enterobacterales_status",
  "ast_oxacillin_resistance_status",
  "ast_vancomycin_resistance_status",
  "ast_mdr_p_aeruginosa_status"
)

phenotype_status =
  bind_rows(
    ordinary_phenotype_status,
    mdr_summary %>%
      transmute(
        across(all_of(organism_id)),
        phenotype = "ast_mdr_p_aeruginosa_status",
        status = ast_mdr_p_aeruginosa_status
      )
  ) %>%
  mutate(phenotype = factor(phenotype, levels = phenotype_columns)) %>%
  pivot_wider(
    id_cols = all_of(organism_id),
    names_from = phenotype,
    values_from = status,
    names_expand = TRUE,
    values_fill = "not_applicable"
  )

phenotype_organisms =
  organism_scope %>%
  left_join(
    phenotype_status,
    by = organism_id,
    na_matches = "never",
    relationship = "one-to-one"
  ) %>%
  mutate(
    across(
      all_of(phenotype_columns),
      ~ replace_na(.x, "not_applicable")
    )
  )
```

An organism outside a rule’s scope receives `not_applicable`. Keeping
`not_detected`, `not_evaluable`, and `not_applicable` separate prevents
incomplete or irrelevant testing from being interpreted as negative
evidence.

### 4.6 Review the phenotype counts

Expect 126,405 culture-organism observations after applying the existing
reviewed reporting scope. Each should have five nonmissing status
fields. These remain culture-level source observations, before
reporting-window selection, label consolidation, or first-culture
selection per patient.

``` r
phenotype_population_check =
  phenotype_organisms %>%
  summarise(
    culture_organism_records = n(),
    distinct_source_labels = n_distinct(organism),
    missing_status_records = sum(
      if_any(all_of(phenotype_columns), is.na)
    )
  )

phenotype_check =
  phenotype_organisms %>%
  select(all_of(phenotype_columns)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "phenotype",
    values_to = "status"
  ) %>%
  count(phenotype, status, name = "culture_organism_count") %>%
  pivot_wider(
    names_from = status,
    values_from = culture_organism_count,
    values_fill = 0L
  )

print(phenotype_population_check, width = Inf)
```

    ## # A tibble: 1 × 3
    ##   culture_organism_records distinct_source_labels missing_status_records
    ##                      <int>                  <int>                  <int>
    ## 1                   126405                     41                      0

``` r
print(phenotype_check, width = Inf)
```

    ## # A tibble: 5 × 5
    ##   phenotype                        not_applicable not_detected not_evaluable
    ##   <chr>                                     <int>        <int>         <int>
    ## 1 ast_cre_status                            40725        59614         25868
    ## 2 ast_esce_enterobacterales_status          49638        54160         15451
    ## 3 ast_mdr_p_aeruginosa_status              115939         6572          2819
    ## 4 ast_oxacillin_resistance_status          116289         7541           153
    ## 5 ast_vancomycin_resistance_status         124316         1596           243
    ##   positive
    ##      <int>
    ## 1      198
    ## 2     7156
    ## 3     1075
    ## 4     2422
    ## 5      250

The resulting `phenotype_organisms` data frame is ready for assigning
reporting categories and consolidating labels within the correct scope.
CRE and MDR classifications are needed here because they determine which
observations can enter the main susceptibility-reporting populations.

## 5. Build reporting cultures and reconcile their AST

The source organism observations now have resistance classifications.
Assign each observation to its reporting category, limit the population
to the reporting period, and combine source labels that belong to the
same reporting culture and category.

The overall *S. aureus* population also needs a separate culture-level
view that combines labels across phenotypes. First-culture selection per
patient follows in the next stage.

### 5.1 Assign categories and apply the reporting period

Assign the reporting category in the same priority order as the SQL:
CRE, MDR *P. aeruginosa*, oxacillin resistance, vancomycin resistance,
ESCE-R, and then `base`. The `base` category means that no positive
phenotype above took precedence; it does not establish susceptibility or
exclude incomplete phenotype testing.

Apply the inclusive start and exclusive end dates in UTC. The source
labels on each culture share the consolidated culture timestamp, so
filtering them before combining labels preserves the reporting scope.

``` r
window_organisms =
  phenotype_organisms %>%
  mutate(
    reporting_phenotype = case_when(
      ast_cre_status == "positive" ~ "cre",
      ast_mdr_p_aeruginosa_status == "positive" ~ "mdr_p_aeruginosa",
      ast_oxacillin_resistance_status == "positive" ~ "oxacillin_resistant",
      ast_vancomycin_resistance_status == "positive" ~ "vancomycin_resistant",
      ast_esce_enterobacterales_status == "positive" ~ "esce_r",
      TRUE ~ "base"
    )
  ) %>%
  filter(
    !is.na(culture_time_shifted),
    culture_time_shifted >= as.POSIXct(analysis_window_start, tz = "UTC"),
    culture_time_shifted < as.POSIXct(
      analysis_window_end_exclusive,
      tz = "UTC"
    )
  )

reporting_culture_id = c(
  culture_id,
  "cohort_id",
  "organism_reporting_group",
  "reporting_phenotype"
)

overall_culture_id = c(
  culture_id,
  "cohort_id",
  "organism_reporting_group"
)
```

### 5.2 Combine source labels within the correct reporting scope

When source labels are combined, retain their original names as a list
for traceability. Combine phenotype evidence using the precedence
`positive`, `not_evaluable`, `not_detected`, and `not_applicable`.

`reporting_cultures` combines labels only within the same culture,
organism reporting group, and phenotype. CRE and MDR groups remain
available here; they will not form main susceptibility-reporting rows.

`overall_cultures` combines all labels for *S. aureus* on the same
culture, regardless of phenotype.

``` r
combine_phenotype_status = function(x) {
  case_when(
    any(x == "positive") ~ "positive",
    any(x == "not_evaluable") ~ "not_evaluable",
    any(x == "not_detected") ~ "not_detected",
    TRUE ~ "not_applicable"
  )
}

reporting_cultures =
  window_organisms %>%
  group_by(across(all_of(reporting_culture_id))) %>%
  summarise(
    culture_time_shifted = first(culture_time_shifted),
    gram_stain_group = first(gram_stain_group),
    source_label_count = n_distinct(organism),
    source_labels = list(sort(unique(organism))),
    across(all_of(phenotype_columns), combine_phenotype_status),
    .groups = "drop"
  )

overall_cultures =
  window_organisms %>%
  filter(organism_reporting_group == "Staphylococcus aureus") %>%
  group_by(across(all_of(overall_culture_id))) %>%
  summarise(
    culture_time_shifted = first(culture_time_shifted),
    gram_stain_group = first(gram_stain_group),
    source_label_count = n_distinct(organism),
    source_labels = list(sort(unique(organism))),
    source_reporting_phenotypes = list(sort(unique(reporting_phenotype))),
    across(all_of(phenotype_columns), combine_phenotype_status),
    .groups = "drop"
  )
```

### 5.3 Define a shared reconciliation for reporting-level AST

There are two distinct reconciliation scopes in this project. Stage 2
resolved repeated values within one original organism label. Here,
combine those already-reconciled endpoints across source labels within a
reporting group.

Any source endpoint marked as conflicting keeps the reporting result in
conflict. Different resolved results also conflict. One consistent
resolved result can be retained alongside a separate label’s missing
endpoint. An endpoint absent from a label contributes no AST record.

| Endpoints from different source labels        | Reporting result     |
|-----------------------------------------------|----------------------|
| Resolved Susceptible and resolved Susceptible | Resolved Susceptible |
| Resolved Susceptible and missing              | Resolved Susceptible |
| Resolved Susceptible and resolved Resistant   | Conflict             |
| Any source endpoint already marked conflict   | Conflict             |
| Missing endpoints only                        | Missing              |

Define one function so both the stratified and overall populations use
this reconciliation consistently.

``` r
reconcile_reporting_ast = function(data, group_columns) {
  data %>%
    group_by(across(all_of(c(group_columns, "antibiotic")))) %>%
    summarise(
      source_endpoint_count = n(),
      source_conflict_count = sum(resolution_status == "conflict"),
      source_missing_count = sum(resolution_status == "missing"),
      distinct_resolved_result_count = n_distinct(
        susceptibility[resolution_status == "resolved"],
        na.rm = TRUE
      ),
      candidate_result = first(
        susceptibility[resolution_status == "resolved"],
        default = NA_character_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      resolution_status = case_when(
        source_conflict_count > 0 | distinct_resolved_result_count > 1 ~
          "conflict",
        distinct_resolved_result_count == 1 ~ "resolved",
        TRUE ~ "missing"
      ),
      susceptibility = if_else(
        resolution_status == "resolved",
        candidate_result,
        NA_character_
      ),
      reconciliation_reason = case_when(
        source_conflict_count > 0 ~ "source_endpoint_conflict",
        distinct_resolved_result_count > 1 ~ "discordant_source_results",
        distinct_resolved_result_count == 1 & source_missing_count > 0 ~
          "resolved_with_missing_source_endpoint",
        distinct_resolved_result_count == 1 ~ "resolved_consistent_result",
        TRUE ~ "missing_only"
      )
    ) %>%
    select(-candidate_result)
}
```

### 5.4 Reconcile AST separately for stratified and overall cultures

Attach the reporting identities to the original endpoints using all
culture and original-organism keys. This retains each endpoint’s link to
the original label that supplied its phenotype evidence.

For stratified reporting, retain the main susceptibility source
categories and reconcile within phenotype. MSSA candidates will be
derived from qualifying `base` *S. aureus* cultures in the next stage.
For overall *S. aureus*, reconcile across all its labels on the culture.

``` r
window_ast =
  ast_endpoints %>%
  inner_join(
    window_organisms %>%
      select(
        all_of(organism_id),
        cohort_id,
        organism_reporting_group,
        reporting_phenotype
      ),
    by = organism_id,
    na_matches = "never",
    relationship = "many-to-one"
  )

stratified_ast =
  window_ast %>%
  filter(
    reporting_phenotype %in% c(
      "base",
      "esce_r",
      "oxacillin_resistant",
      "vancomycin_resistant"
    )
  ) %>%
  reconcile_reporting_ast(reporting_culture_id)

overall_ast =
  window_ast %>%
  filter(organism_reporting_group == "Staphylococcus aureus") %>%
  reconcile_reporting_ast(overall_culture_id)
```

### 5.5 Inspect the reporting populations and AST statuses

Review the number of source observations in the window, the resulting
phenotype culture groups, and the separate overall *S. aureus* cultures.
The overall and stratified views overlap; their counts are not additive
unique-culture or unique-patient totals.

``` r
reporting_scope_check = tibble(
  population = c(
    "Original organism observations in the reporting window",
    "Consolidated culture-organism-phenotype groups",
    "Overall Staphylococcus aureus cultures"
  ),
  record_count = c(
    nrow(window_organisms),
    nrow(reporting_cultures),
    nrow(overall_cultures)
  )
)

reporting_ast_check =
  bind_rows(
    stratified_ast %>%
      count(resolution_status, name = "endpoint_count") %>%
      mutate(scope = "stratified"),
    overall_ast %>%
      count(resolution_status, name = "endpoint_count") %>%
      mutate(scope = "overall_s_aureus")
  ) %>%
  select(scope, resolution_status, endpoint_count) %>%
  arrange(scope, resolution_status)

print(reporting_scope_check, width = Inf)
```

    ## # A tibble: 3 × 2
    ##   population                                             record_count
    ##   <chr>                                                         <int>
    ## 1 Original organism observations in the reporting window        10862
    ## 2 Consolidated culture-organism-phenotype groups                10835
    ## 3 Overall Staphylococcus aureus cultures                          531

``` r
reporting_cultures %>%
  count(cohort_id, reporting_phenotype, name = "culture_group_count") %>%
  print(n = Inf)
```

    ## # A tibble: 16 × 3
    ##    cohort_id   reporting_phenotype  culture_group_count
    ##    <chr>       <chr>                              <int>
    ##  1 blood       base                                 850
    ##  2 blood       cre                                    3
    ##  3 blood       esce_r                                91
    ##  4 blood       oxacillin_resistant                   60
    ##  5 blood       vancomycin_resistant                  12
    ##  6 respiratory base                                 714
    ##  7 respiratory cre                                    2
    ##  8 respiratory esce_r                                24
    ##  9 respiratory mdr_p_aeruginosa                      34
    ## 10 respiratory oxacillin_resistant                   42
    ## 11 urine       base                                8104
    ## 12 urine       cre                                   18
    ## 13 urine       esce_r                               832
    ## 14 urine       mdr_p_aeruginosa                      10
    ## 15 urine       oxacillin_resistant                   22
    ## 16 urine       vancomycin_resistant                  17

``` r
print(reporting_ast_check, n = Inf)
```

    ## # A tibble: 6 × 3
    ##   scope            resolution_status endpoint_count
    ##   <chr>            <chr>                      <int>
    ## 1 overall_s_aureus conflict                      35
    ## 2 overall_s_aureus missing                       47
    ## 3 overall_s_aureus resolved                    6079
    ## 4 stratified       conflict                     622
    ## 5 stratified       missing                      115
    ## 6 stratified       resolved                  163793

The next stage will select each patient’s first qualifying culture
independently within each published category, including separate MRSA,
MSSA, and overall *S. aureus* populations.

## 6. Select the first qualifying culture within each reporting category

Each patient can contribute one first qualifying culture per source
release, specimen cohort, organism reporting group, and published
category. Categories are selected independently. A patient can therefore
contribute to both an ordinary and a resistant category, or to MRSA and
MSSA.

Overall *S. aureus* has its own first-culture selection across
phenotypes. The reporting window was already applied in the previous
stage.

### 6.1 Reproduce the deterministic culture ordering

Order cultures by their shifted timestamp. When timestamps tie, the SQL
uses a hashed culture key to choose consistently. Reproduce that key
from the same ordered identifiers so that R selects the same culture.
The hash helper uses the `digest` package. If it is unavailable, install
it once with `install.packages("digest")` in the console.

This helper will also be reused when creating the final mart key.

``` r
make_dbt_key = function(data, key_columns) {
  key_text =
    data %>%
    select(all_of(key_columns)) %>%
    mutate(
      across(
        everything(),
        ~ replace_na(as.character(.x), "_dbt_utils_surrogate_key_null_")
      )
    ) %>%
    unite("key_text", everything(), sep = "-") %>%
    pull(key_text)

  map_chr(
    key_text,
    ~ digest::digest(.x, algo = "md5", serialize = FALSE)
  )
}

culture_keys =
  window_organisms %>%
  select(all_of(culture_id)) %>%
  distinct() %>%
  mutate(
    culture_key = make_dbt_key(pick(all_of(culture_id)), culture_id)
  )

stratified_sources =
  reporting_cultures %>%
  left_join(
    culture_keys,
    by = culture_id,
    na_matches = "never",
    relationship = "many-to-one"
  ) %>%
  mutate(source_reporting_phenotype = reporting_phenotype)
```

### 6.2 Identify eligible MSSA cultures

MSSA requires a `base` *S. aureus* culture with an oxacillin-resistance
status of `not_detected`, plus susceptible oxacillin or cefoxitin marker
evidence. All observed marker endpoints must be resolved Susceptible;
missing, conflicting, Intermediate, or Resistant marker endpoints block
qualification.

The grouping below starts with observed marker rows, so every resulting
group contains at least one marker. `all()` then requires every observed
marker in that group to be resolved Susceptible. A culture with no
marker endpoint cannot enter this marker table.

``` r
mssa_marker_cultures =
  stratified_ast %>%
  filter(
    organism_reporting_group == "Staphylococcus aureus",
    reporting_phenotype == "base",
    antibiotic %in% c("Oxacillin", "Cefoxitin")
  ) %>%
  group_by(across(all_of(reporting_culture_id))) %>%
  summarise(
    all_observed_markers_susceptible = all(
      resolution_status == "resolved" &
      susceptibility %in% "Susceptible"
    ),
    .groups = "drop"
  ) %>%
  filter(all_observed_markers_susceptible) %>%
  select(all_of(reporting_culture_id))

mssa_candidates =
  stratified_sources %>%
  filter(
    organism_reporting_group == "Staphylococcus aureus",
    reporting_phenotype == "base",
    ast_oxacillin_resistance_status == "not_detected"
  ) %>%
  inner_join(
    mssa_marker_cultures,
    by = reporting_culture_id,
    na_matches = "never",
    relationship = "one-to-one"
  ) %>%
  mutate(
    reporting_phenotype = "oxacillin_susceptible",
    organism_name = str_c(organism_reporting_group, " — MSSA")
  )
```

`source_reporting_phenotype` stays `base` for these MSSA cultures. This
preserves the link to the correct stratified AST records while
publishing the category as `oxacillin_susceptible`.

### 6.3 Define the remaining published categories

Include ordinary non-*S. aureus* organisms, ESCE-R observations without
positive CRE evidence, VRE in the two defined *Enterococcus* species,
and MRSA in *S. aureus*. CRE and MDR *P. aeruginosa* groups do not enter
these main susceptibility populations.

``` r
ordinary_candidates =
  stratified_sources %>%
  filter(
    (
      reporting_phenotype == "base" &
      organism_reporting_group != "Staphylococcus aureus"
    ) |
    (
      reporting_phenotype == "esce_r" &
      ast_cre_status != "positive"
    ) |
    (
      reporting_phenotype == "vancomycin_resistant" &
      organism_reporting_group %in% c(
        "Enterococcus faecalis",
        "Enterococcus faecium"
      )
    ) |
    (
      reporting_phenotype == "oxacillin_resistant" &
      organism_reporting_group == "Staphylococcus aureus"
    )
  ) %>%
  mutate(
    organism_name = case_when(
      reporting_phenotype == "esce_r" ~
        str_c(organism_reporting_group, " — ESCE-R"),
      reporting_phenotype == "vancomycin_resistant" ~
        str_c(organism_reporting_group, " — VRE"),
      reporting_phenotype == "oxacillin_resistant" ~
        str_c(organism_reporting_group, " — MRSA"),
      TRUE ~ organism_reporting_group
    )
  )

stratified_candidates =
  bind_rows(ordinary_candidates, mssa_candidates)
```

### 6.4 Select first cultures independently within each category

Sort eligible cultures by time and the deterministic culture key, then
retain the first row in each patient-category group. Each published
category has at most one candidate row per culture after consolidation,
so the culture key resolves timestamp ties within that group.

The overall *S. aureus* selection uses all its eligible cultures,
regardless of MRSA or MSSA qualification.

``` r
patient_category_id = c(
  "source_version",
  "cohort_id",
  "patient_key",
  "organism_reporting_group",
  "reporting_phenotype"
)

first_stratified =
  stratified_candidates %>%
  arrange(culture_time_shifted, culture_key) %>%
  group_by(across(all_of(patient_category_id))) %>%
  slice_head(n = 1) %>%
  ungroup()

first_overall =
  overall_cultures %>%
  left_join(
    culture_keys,
    by = culture_id,
    na_matches = "never",
    relationship = "many-to-one"
  ) %>%
  mutate(
    reporting_phenotype = "overall",
    source_reporting_phenotype = NA_character_,
    organism_name = organism_reporting_group
  ) %>%
  arrange(culture_time_shifted, culture_key) %>%
  group_by(across(all_of(patient_category_id))) %>%
  slice_head(n = 1) %>%
  ungroup()

selected_cultures =
  bind_rows(first_stratified, first_overall)
```

### 6.5 Attach AST for each selected culture and scope

For stratified selections, use the original source phenotype when
matching AST. This is particularly important for MSSA, whose published
category differs from its `base` source category. Overall *S. aureus*
uses the separate overall reconciliation.

Keep `selected_cultures` as the eligible population independently of the
AST join. A selected culture without an antibiotic endpoint must still
contribute to the eligible reporting denominator.

``` r
selected_columns = c(
  culture_id,
  "cohort_id",
  "organism_reporting_group",
  "reporting_phenotype",
  "organism_name",
  "gram_stain_group",
  "culture_time_shifted",
  "culture_key"
)

stratified_join_id = c(
  culture_id,
  "cohort_id",
  "organism_reporting_group",
  "source_reporting_phenotype"
)

first_stratified_ast =
  first_stratified %>%
  select(all_of(selected_columns), source_reporting_phenotype) %>%
  inner_join(
    stratified_ast %>%
      rename(source_reporting_phenotype = reporting_phenotype),
    by = stratified_join_id,
    na_matches = "never",
    relationship = "one-to-many"
  )

first_overall_ast =
  first_overall %>%
  select(all_of(selected_columns)) %>%
  inner_join(
    overall_ast,
    by = overall_culture_id,
    na_matches = "never",
    relationship = "one-to-many"
  )

selected_ast =
  bind_rows(first_stratified_ast, first_overall_ast)

selected_culture_id = c(
  culture_id,
  "cohort_id",
  "organism_reporting_group",
  "reporting_phenotype"
)
```

### 6.6 Check selections against the supplied Tableau populations

Review the candidate and selected counts. Then compare the nine *S.
aureus* category counts against the eligible counts in the supplied
Tableau CSV. These expected values apply to the current 2025-04-11
release and 2023 reporting window.

MRSA and MSSA counts need not sum to the overall count because each
population selects its first culture independently. For example, the
supplied urine populations contain 22 MRSA, 38 MSSA, and 58 overall
patient-category observations.

``` r
first_selection_check = tibble(
  population = c(
    "Stratified candidate cultures",
    "First stratified cultures",
    "Overall Staphylococcus aureus candidate cultures",
    "First overall Staphylococcus aureus cultures"
  ),
  record_count = c(
    nrow(stratified_candidates),
    nrow(first_stratified),
    nrow(overall_cultures),
    nrow(first_overall)
  )
)

staph_reference_counts = tribble(
  ~cohort_id, ~reporting_phenotype, ~expected_count,
  "blood", "overall", 150L,
  "blood", "oxacillin_resistant", 59L,
  "blood", "oxacillin_susceptible", 91L,
  "respiratory", "overall", 292L,
  "respiratory", "oxacillin_resistant", 34L,
  "respiratory", "oxacillin_susceptible", 258L,
  "urine", "overall", 58L,
  "urine", "oxacillin_resistant", 22L,
  "urine", "oxacillin_susceptible", 38L
)

staph_selection_check =
  selected_cultures %>%
  filter(organism_reporting_group == "Staphylococcus aureus") %>%
  count(cohort_id, reporting_phenotype, name = "selected_count") %>%
  full_join(
    staph_reference_counts,
    by = c("cohort_id", "reporting_phenotype"),
    relationship = "one-to-one"
  ) %>%
  mutate(
    selected_count = replace_na(selected_count, 0L),
    matches_reference =
      !is.na(expected_count) & selected_count == expected_count
  ) %>%
  arrange(cohort_id, reporting_phenotype)

print(first_selection_check, width = Inf)
```

    ## # A tibble: 4 × 2
    ##   population                                       record_count
    ##   <chr>                                                   <int>
    ## 1 Stratified candidate cultures                           10763
    ## 2 First stratified cultures                                9203
    ## 3 Overall Staphylococcus aureus candidate cultures          531
    ## 4 First overall Staphylococcus aureus cultures              500

``` r
print(staph_selection_check, n = Inf, width = Inf)
```

    ## # A tibble: 9 × 5
    ##   cohort_id   reporting_phenotype   selected_count expected_count
    ##   <chr>       <chr>                          <int>          <int>
    ## 1 blood       overall                          150            150
    ## 2 blood       oxacillin_resistant               59             59
    ## 3 blood       oxacillin_susceptible             91             91
    ## 4 respiratory overall                          292            292
    ## 5 respiratory oxacillin_resistant               34             34
    ## 6 respiratory oxacillin_susceptible            258            258
    ## 7 urine       overall                           58             58
    ## 8 urine       oxacillin_resistant               22             22
    ## 9 urine       oxacillin_susceptible             38             38
    ##   matches_reference
    ##   <lgl>            
    ## 1 TRUE             
    ## 2 TRUE             
    ## 3 TRUE             
    ## 4 TRUE             
    ## 5 TRUE             
    ## 6 TRUE             
    ## 7 TRUE             
    ## 8 TRUE             
    ## 9 TRUE

``` r
selected_cultures %>%
  count(cohort_id, reporting_phenotype, name = "selected_culture_count") %>%
  print(n = Inf)
```

    ## # A tibble: 17 × 3
    ##    cohort_id   reporting_phenotype   selected_culture_count
    ##    <chr>       <chr>                                  <int>
    ##  1 blood       base                                     678
    ##  2 blood       esce_r                                    83
    ##  3 blood       overall                                  150
    ##  4 blood       oxacillin_resistant                       59
    ##  5 blood       oxacillin_susceptible                     91
    ##  6 blood       vancomycin_resistant                      12
    ##  7 respiratory base                                     380
    ##  8 respiratory esce_r                                    19
    ##  9 respiratory overall                                  292
    ## 10 respiratory oxacillin_resistant                       34
    ## 11 respiratory oxacillin_susceptible                    258
    ## 12 urine       base                                    6841
    ## 13 urine       esce_r                                   672
    ## 14 urine       overall                                   58
    ## 15 urine       oxacillin_resistant                       22
    ## 16 urine       oxacillin_susceptible                     38
    ## 17 urine       vancomycin_resistant                      16

The next stage will expand the selected populations to approved
organism-antibiotic pairs, apply scoped exclusions, and calculate the
final mart’s counts, percentages, and display flags.

## 7. Apply antibiotic panels and build the susceptibility mart

Start with the selected cultures, expand each population to its approved
antibiotics, and then attach AST. This keeps untested antibiotics in the
eligible denominator. Apply the same reviewed panel and exclusion files
as the SQL pipeline.

### 7.1 Load approved antibiotic pairs and exclusion rules

Only panel rows marked for the primary dashboard and reviewed are
approved. An exclusion can apply to one cohort or phenotype, or to all
of them using `*`. These are the project’s existing reporting policies;
they do not change the AST evidence retained in earlier stages.

``` r
approved_pairs =
  read_csv(
    file.path(seed_dir, "reportable_organism_antibiotic_pairs.csv"),
    col_types = cols(
      .default = col_character(),
      pair_order = col_integer(),
      include_primary_dashboard = col_integer()
    ),
    na = "",
    trim_ws = FALSE
  ) %>%
  filter(
    include_primary_dashboard == 1L,
    review_status == "reviewed"
  ) %>%
  select(organism_reporting_group, cohort_id, antibiotic, pair_order)

exclusion_rules =
  read_csv(
    file.path(seed_dir, "organism_antibiotic_exclusions.csv"),
    col_types = cols(.default = col_character()),
    na = "",
    trim_ws = FALSE
  ) %>%
  transmute(
    organism_reporting_group,
    antibiotic,
    cohort_scope = cohort_id,
    phenotype_scope = reporting_phenotype
  ) %>%
  distinct()
```

### 7.2 Determine which pairs each reporting population can publish

The same organism can have different exclusions in its ordinary and
resistant categories. Resolve those rules after the published phenotype
is known. Make a distinct list of excluded combinations, then remove
them with `anti_join()`. Overlapping rules therefore cannot multiply the
reporting rows.

The panel joins intentionally expand one population to several
antibiotics and reuse the same organism panel across phenotypes.

``` r
population_id = c(
  "source_version",
  "cohort_id",
  "organism_reporting_group",
  "reporting_phenotype"
)

mart_id = c(population_id, "antibiotic")

population_panels =
  selected_cultures %>%
  select(all_of(population_id)) %>%
  distinct() %>%
  inner_join(
    approved_pairs,
    by = c("organism_reporting_group", "cohort_id"),
    na_matches = "never",
    relationship = "many-to-many"
  )

excluded_pairs =
  population_panels %>%
  inner_join(
    exclusion_rules,
    by = c("organism_reporting_group", "antibiotic"),
    na_matches = "never",
    relationship = "many-to-many"
  ) %>%
  filter(
    cohort_scope == "*" | cohort_id == cohort_scope,
    phenotype_scope == "*" | reporting_phenotype == phenotype_scope
  ) %>%
  select(all_of(mart_id)) %>%
  distinct()

reporting_pairs =
  population_panels %>%
  anti_join(excluded_pairs, by = mart_id, na_matches = "never")
```

### 7.3 Expand selected cultures and attach their antibiotic results

Every eligible culture receives one row for each permitted antibiotic.
The `left_join()` then attaches any observed endpoint. A missing joined
status means the antibiotic was untested. An endpoint whose status is
`missing` means it was present in the source but lacked a usable result;
these are counted separately.

``` r
eligible_grid =
  selected_cultures %>%
  select(all_of(selected_culture_id), organism_name) %>%
  inner_join(
    reporting_pairs,
    by = population_id,
    na_matches = "never",
    relationship = "many-to-many"
  )

endpoint_id = c(selected_culture_id, "antibiotic")

reporting_endpoints =
  eligible_grid %>%
  left_join(
    selected_ast %>%
      select(
        all_of(endpoint_id),
        resolution_status,
        susceptibility,
        reconciliation_reason
      ),
    by = endpoint_id,
    na_matches = "never",
    relationship = "one-to-one"
  )
```

### 7.4 Calculate result counts

Summarise to one row per reporting population and antibiotic.
`categorized_count` includes only resolved Susceptible, Intermediate,
and Resistant results. Inconclusive, Synergism, missing, conflicting,
and untested endpoints remain separate counts.

Using `%in%` below produces `FALSE` for missing values, so absent
endpoints cannot turn a count into `NA`.

``` r
mart_counts =
  reporting_endpoints %>%
  group_by(
    source_version,
    cohort_id,
    organism_reporting_group,
    organism_name,
    reporting_phenotype,
    antibiotic,
    pair_order
  ) %>%
  summarise(
    eligible_count = n(),
    tested_endpoint_count = sum(!is.na(resolution_status)),
    categorized_count = sum(
      resolution_status %in% "resolved" &
      susceptibility %in% c("Susceptible", "Intermediate", "Resistant")
    ),
    susceptible_count = sum(
      resolution_status %in% "resolved" &
      susceptibility %in% "Susceptible"
    ),
    resistant_count = sum(
      resolution_status %in% "resolved" &
      susceptibility %in% "Resistant"
    ),
    intermediate_count = sum(
      resolution_status %in% "resolved" &
      susceptibility %in% "Intermediate"
    ),
    inconclusive_count = sum(
      resolution_status %in% "resolved" &
      susceptibility %in% "Inconclusive"
    ),
    synergism_count = sum(
      resolution_status %in% "resolved" &
      susceptibility %in% "Synergism"
    ),
    missing_result_count = sum(resolution_status %in% "missing"),
    conflict_result_count = sum(resolution_status %in% "conflict"),
    untested_count = sum(is.na(resolution_status)),
    .groups = "drop"
  )
```

### 7.5 Add percentages, display flags, and final metadata

Coverage uses all eligible observations as its denominator. Percent
susceptible and percent resistant use only categorized S/I/R results. A
zero denominator produces `NA`.

The helper rounds nonnegative percentages to one decimal with ties
rounded upward, matching the SQL. Base R’s `round()` uses a different
tie rule. Scaling the counts by 1,000 before dividing also avoids
introducing an extra rounding step at percentage half-ties.

The display flag requires at least 30 categorized results. Keep all rows
and their calculated percentages, including those flagged
`insufficient_data`.

``` r
percentage_1dp = function(numerator, denominator) {
  floor(
    1000 * as.double(numerator) /
      na_if(as.double(denominator), 0) + 0.5
  ) / 10
}

mart_susceptibility_r =
  mart_counts %>%
  mutate(
    susceptibility_mart_key = make_dbt_key(pick(all_of(mart_id)), mart_id),
    .before = 1
  ) %>%
  mutate(
    categorized_result_coverage_pct = percentage_1dp(
      categorized_count,
      eligible_count
    ),
    susceptible_pct = percentage_1dp(susceptible_count, categorized_count),
    resistant_pct = percentage_1dp(resistant_count, categorized_count),
    display_status = if_else(
      categorized_count >= 30L,
      "display",
      "insufficient_data"
    ),
    analysis_window_start = .env$analysis_window_start,
    analysis_window_end_exclusive = .env$analysis_window_end_exclusive
  ) %>%
  left_join(
    organism_reporting_groups %>%
      select(organism_reporting_group, gram_stain_group),
    by = "organism_reporting_group",
    na_matches = "never",
    relationship = "many-to-one"
  ) %>%
  arrange(
    source_version,
    cohort_id,
    organism_reporting_group,
    reporting_phenotype,
    pair_order,
    antibiotic
  )
```

The resulting `mart_susceptibility_r` uses the same 26 columns and
column order as the supplied Tableau CSV. Its row order is for
readability; the final comparison will match rows by their keys.

### 7.6 Check the mart’s structure and count reconciliation

For the supplied release, window, and seed files, expect 396 rows, 26
columns, and 73 reporting populations. Every failure count below should
be zero. The display-status counts should be 187 `display` and 209
`insufficient_data`.

``` r
mart_build_check =
  mart_susceptibility_r %>%
  summarise(
    mart_rows = n(),
    duplicate_mart_keys = n() - n_distinct(susceptibility_mart_key),
    eligible_balance_failures = sum(
      eligible_count != tested_endpoint_count + untested_count
    ),
    categorized_balance_failures = sum(
      categorized_count !=
        susceptible_count + intermediate_count + resistant_count
    ),
    endpoint_balance_failures = sum(
      tested_endpoint_count !=
        categorized_count + inconclusive_count + synergism_count +
        missing_result_count + conflict_result_count
    ),
    missing_gram_groups = sum(is.na(gram_stain_group))
  ) %>%
  mutate(
    mart_columns = ncol(mart_susceptibility_r),
    reporting_populations = nrow(
      mart_susceptibility_r %>%
        select(all_of(population_id)) %>%
        distinct()
    )
  )

print(mart_build_check, width = Inf)
```

    ## # A tibble: 1 × 8
    ##   mart_rows duplicate_mart_keys eligible_balance_failures
    ##       <int>               <int>                     <int>
    ## 1       396                   0                         0
    ##   categorized_balance_failures endpoint_balance_failures missing_gram_groups
    ##                          <int>                     <int>               <int>
    ## 1                            0                         0                   0
    ##   mart_columns reporting_populations
    ##          <int>                 <int>
    ## 1           26                    73

``` r
mart_susceptibility_r %>%
  count(display_status, name = "mart_rows") %>%
  print(n = Inf)
```

    ## # A tibble: 2 × 2
    ##   display_status    mart_rows
    ##   <chr>                 <int>
    ## 1 display                 187
    ## 2 insufficient_data       209

``` r
mart_susceptibility_r %>%
  select(
    cohort_id,
    organism_name,
    antibiotic,
    eligible_count,
    categorized_count,
    susceptible_pct,
    display_status
  ) %>%
  print(n = 20, width = Inf)
```

    ## # A tibble: 396 × 7
    ##    cohort_id organism_name                                 
    ##    <chr>     <chr>                                         
    ##  1 blood     Citrobacter freundii complex                  
    ##  2 blood     Citrobacter freundii complex                  
    ##  3 blood     Citrobacter freundii complex                  
    ##  4 blood     Citrobacter freundii complex                  
    ##  5 blood     Citrobacter freundii complex                  
    ##  6 blood     Citrobacter freundii complex                  
    ##  7 blood     Citrobacter freundii complex                  
    ##  8 blood     Citrobacter koseri                            
    ##  9 blood     Citrobacter koseri                            
    ## 10 blood     Citrobacter koseri                            
    ## 11 blood     Citrobacter koseri                            
    ## 12 blood     Citrobacter koseri                            
    ## 13 blood     Citrobacter koseri                            
    ## 14 blood     Citrobacter koseri                            
    ## 15 blood     Citrobacter koseri                            
    ## 16 blood     Coagulase-negative staphylococci (unspeciated)
    ## 17 blood     Coagulase-negative staphylococci (unspeciated)
    ## 18 blood     Coagulase-negative staphylococci (unspeciated)
    ## 19 blood     Enterobacter cloacae complex                  
    ## 20 blood     Enterobacter cloacae complex                  
    ##    antibiotic                    eligible_count categorized_count
    ##    <chr>                                  <int>             <int>
    ##  1 Ciprofloxacin                              8                 8
    ##  2 Gentamicin                                 8                 8
    ##  3 Trimethoprim/Sulfamethoxazole              8                 8
    ##  4 Cefepime                                   8                 2
    ##  5 Meropenem                                  8                 8
    ##  6 Ertapenem                                  8                 8
    ##  7 Levofloxacin                               8                 7
    ##  8 Ceftriaxone                                5                 5
    ##  9 Ciprofloxacin                              5                 5
    ## 10 Gentamicin                                 5                 5
    ## 11 Trimethoprim/Sulfamethoxazole              5                 5
    ## 12 Meropenem                                  5                 5
    ## 13 Ertapenem                                  5                 5
    ## 14 Levofloxacin                               5                 5
    ## 15 Piperacillin/Tazobactam                    5                 5
    ## 16 Vancomycin                                57                57
    ## 17 Trimethoprim/Sulfamethoxazole             57                57
    ## 18 Erythromycin                              57                57
    ## 19 Gentamicin                                22                22
    ## 20 Meropenem                                 22                22
    ##    susceptible_pct display_status   
    ##              <dbl> <chr>            
    ##  1            87.5 insufficient_data
    ##  2            87.5 insufficient_data
    ##  3            87.5 insufficient_data
    ##  4            50   insufficient_data
    ##  5           100   insufficient_data
    ##  6           100   insufficient_data
    ##  7            71.4 insufficient_data
    ##  8           100   insufficient_data
    ##  9           100   insufficient_data
    ## 10           100   insufficient_data
    ## 11           100   insufficient_data
    ## 12           100   insufficient_data
    ## 13           100   insufficient_data
    ## 14           100   insufficient_data
    ## 15            80   insufficient_data
    ## 16           100   display          
    ## 17            47.4 display          
    ## 18            35.1 display          
    ## 19            95.5 insufficient_data
    ## 20           100   insufficient_data
    ## # ℹ 376 more rows

Matching these totals checks the structure, but does not yet establish
that every result matches SQL. The next stage will compare all rows and
fields with the supplied Tableau CSV, then export the verified R mart.

## 8. Validate against the Tableau CSV and export the R mart

The structure and count checks passed in Step 7. Now compare the actual
values with the existing SQL-generated Tableau CSV. Match rows using
their five identifying fields, so differences in row order do not affect
the comparison. Compare all remaining fields, including the mart hash,
counts, percentages, display flags, dates, and Gram group.

### 8.1 Import the existing Tableau export with explicit types

Set `tableau_reference_file` to your existing Tableau data-source CSV.
The path below is an example; there is no need to move or duplicate the
file. Use the same working directory as the previous sections.

Keep `source_version` as text even though it looks like a date. Parse
the two analysis-window fields as dates and the measures as numbers. The
supplied export uses the literal lowercase `null` for missing
percentages, so recognize that value as `NA` during import.

``` r
tableau_reference_file = "../data/exports/mart_susceptibility.csv"

tableau_reference =
  read_csv(
    tableau_reference_file,
    col_types = cols(
      .default = col_character(),
      pair_order = col_double(),
      eligible_count = col_double(),
      tested_endpoint_count = col_double(),
      categorized_count = col_double(),
      susceptible_count = col_double(),
      resistant_count = col_double(),
      intermediate_count = col_double(),
      inconclusive_count = col_double(),
      synergism_count = col_double(),
      missing_result_count = col_double(),
      conflict_result_count = col_double(),
      untested_count = col_double(),
      categorized_result_coverage_pct = col_double(),
      susceptible_pct = col_double(),
      resistant_pct = col_double(),
      analysis_window_start = col_date(format = "%Y-%m-%d"),
      analysis_window_end_exclusive = col_date(format = "%Y-%m-%d")
    ),
    na = c("", "null"),
    trim_ws = FALSE
  )

reference_parse_problems = problems(tableau_reference)

if (nrow(reference_parse_problems) > 0L) {
  print(reference_parse_problems, n = 20, width = Inf)
  stop("Resolve the reference CSV parsing problems before comparing marts.")
}
```

### 8.2 Check column names and row identities

Both tables must contain the same column names. The natural row identity
is `source_version`, `cohort_id`, `organism_reporting_group`,
`reporting_phenotype`, and `antibiotic`, already stored in `mart_id`.
Check these fields for missing values and duplicates before joining.

``` r
column_check = tibble(
  reference_columns = ncol(tableau_reference),
  r_columns = ncol(mart_susceptibility_r),
  same_column_names = setequal(
    names(tableau_reference),
    names(mart_susceptibility_r)
  ),
  same_column_order = identical(
    names(tableau_reference),
    names(mart_susceptibility_r)
  )
)

print(column_check, width = Inf)
```

    ## # A tibble: 1 × 4
    ##   reference_columns r_columns same_column_names same_column_order
    ##               <int>     <int> <lgl>             <lgl>            
    ## 1                26        26 TRUE              TRUE

``` r
if (!column_check$same_column_names) {
  print(list(
    only_in_reference = setdiff(
      names(tableau_reference), names(mart_susceptibility_r)
    ),
    only_in_r = setdiff(
      names(mart_susceptibility_r), names(tableau_reference)
    )
  ))
  stop("The R mart and reference CSV have different columns.")
}

r_keys =
  mart_susceptibility_r %>%
  select(all_of(mart_id))

reference_keys =
  tableau_reference %>%
  select(all_of(mart_id))

key_check = tibble(
  dataset = c("R mart", "Tableau CSV"),
  row_count = c(nrow(r_keys), nrow(reference_keys)),
  duplicate_key_rows = c(
    nrow(r_keys) - nrow(distinct(r_keys)),
    nrow(reference_keys) - nrow(distinct(reference_keys))
  ),
  missing_key_rows = c(
    sum(!complete.cases(r_keys)),
    sum(!complete.cases(reference_keys))
  )
)

print(key_check, width = Inf)
```

    ## # A tibble: 2 × 4
    ##   dataset     row_count duplicate_key_rows missing_key_rows
    ##   <chr>           <int>              <int>            <int>
    ## 1 R mart            396                  0                0
    ## 2 Tableau CSV       396                  0                0

``` r
if (any(key_check$duplicate_key_rows > 0L | key_check$missing_key_rows > 0L)) {
  stop("Resolve missing or duplicate row identities before comparing values.")
}
```

### 8.3 Find missing or extra rows

The two `anti_join()` operations identify any row identities present in
only one table. Both results should have zero rows.

``` r
missing_from_r =
  reference_keys %>%
  anti_join(r_keys, by = mart_id, na_matches = "never")

extra_in_r =
  r_keys %>%
  anti_join(reference_keys, by = mart_id, na_matches = "never")

print(missing_from_r, n = 20, width = Inf)
```

    ## # A tibble: 0 × 5
    ## # ℹ 5 variables: source_version <chr>, cohort_id <chr>, organism_reporting_group <chr>, reporting_phenotype <chr>, antibiotic <chr>

``` r
print(extra_in_r, n = 20, width = Inf)
```

    ## # A tibble: 0 × 5
    ## # ℹ 5 variables: source_version <chr>, cohort_id <chr>, organism_reporting_group <chr>, reporting_phenotype <chr>, antibiotic <chr>

### 8.4 Compare every remaining field

Convert the values to a common text representation after numbers and
dates have been parsed. This allows different field types to appear in
one long comparison table. Numeric values such as `100` and `100.0`
compare equally; dates use their ISO representation. Text values retain
their spelling and whitespace.

The join compares fields only for matching row identities. Missing and
extra rows are checked separately above. Two missing values match, while
a missing value on only one side is a difference.

``` r
value_columns = setdiff(names(tableau_reference), mart_id)

r_values =
  mart_susceptibility_r %>%
  mutate(
    across(where(is.numeric), as.double),
    across(all_of(value_columns), as.character)
  ) %>%
  pivot_longer(
    cols = all_of(value_columns),
    names_to = "field",
    values_to = "r_value"
  )

reference_values =
  tableau_reference %>%
  mutate(across(all_of(value_columns), as.character)) %>%
  pivot_longer(
    cols = all_of(value_columns),
    names_to = "field",
    values_to = "reference_value"
  )

field_comparison =
  r_values %>%
  inner_join(
    reference_values,
    by = c(mart_id, "field"),
    na_matches = "never",
    relationship = "one-to-one"
  ) %>%
  mutate(
    values_match =
      coalesce(r_value == reference_value, FALSE) |
      (is.na(r_value) & is.na(reference_value))
  )

field_differences =
  field_comparison %>%
  filter(!values_match) %>%
  select(all_of(mart_id), field, r_value, reference_value)

comparison_check = tibble(
  reference_rows = nrow(tableau_reference),
  r_rows = nrow(mart_susceptibility_r),
  missing_from_r = nrow(missing_from_r),
  extra_in_r = nrow(extra_in_r),
  nonkey_fields_compared = length(value_columns),
  field_values_compared = nrow(field_comparison),
  field_mismatches = nrow(field_differences)
) %>%
  mutate(
    matches_reference =
      missing_from_r == 0L &
      extra_in_r == 0L &
      field_mismatches == 0L &
      reference_rows > 0L
  )

print(comparison_check, width = Inf)
```

    ## # A tibble: 1 × 8
    ##   reference_rows r_rows missing_from_r extra_in_r nonkey_fields_compared
    ##            <int>  <int>          <int>      <int>                  <int>
    ## 1            396    396              0          0                     21
    ##   field_values_compared field_mismatches matches_reference
    ##                   <int>            <int> <lgl>            
    ## 1                  8316                0 TRUE

``` r
print(field_differences, n = 20, width = Inf)
```

    ## # A tibble: 0 × 8
    ## # ℹ 8 variables: source_version <chr>, cohort_id <chr>, organism_reporting_group <chr>, reporting_phenotype <chr>, antibiotic <chr>, field <chr>, r_value <chr>, reference_value <chr>

For the supplied CSV, expect 396 rows in both tables, zero missing or
extra rows, and zero field mismatches. The comparison covers 21 nonkey
fields across 396 rows: 8,316 field values, in addition to the five
identifying fields. `matches_reference` should be `TRUE`.

If differences appear, review `field_differences`, `missing_from_r`, and
`extra_in_r`. The export below stops when the comparison fails.

### 8.5 Export the verified R mart

Write a separate R-generated CSV after the comparison passes. Keep the
reference column order and use the same lowercase `null` missing-value
marker as the supplied export. Change `output_file` if you prefer a
different output folder.

``` r
if (!isTRUE(comparison_check$matches_reference)) {
  stop("The marts differ. Review the comparison results before exporting.")
}

output_file = "../data/processed/mart_susceptibility_r.csv"

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

mart_susceptibility_r %>%
  select(all_of(names(tableau_reference))) %>%
  write_csv(output_file, na = "null")

export_check = tibble(
  output_file = output_file,
  exported_rows = nrow(mart_susceptibility_r),
  exported_columns = ncol(mart_susceptibility_r),
  matches_tableau_reference = comparison_check$matches_reference
)

print(export_check, width = Inf)
```

    ## # A tibble: 1 × 4
    ##   output_file                                 exported_rows exported_columns
    ##   <chr>                                               <int>            <int>
    ## 1 ../data/processed/mart_susceptibility_r.csv           396               26
    ##   matches_tableau_reference
    ##   <lgl>                    
    ## 1 TRUE

When these checks pass, this R Markdown document demonstrates the
complete transformation from the original cultures CSV to the same
cleaned susceptibility mart used by Tableau, using the project’s shared
reference tables and tidyverse transformations.

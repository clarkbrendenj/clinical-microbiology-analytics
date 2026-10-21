# Metric Definitions

## Purpose and scope

This document defines the source population, analytical grains, organism mappings, resistance phenotypes, deduplication policies, susceptibility denominators, and reporting exclusions used in the **Clinical Microbiology Analytics Lakehouse** portfolio project.

The source is the **Stanford Antibiotic Resistance Microbiology Dataset (ARMD), release `2025-04-11`**. The pipeline uses the source's categorical antimicrobial susceptibility testing (AST) interpretations. It does not reconstruct minimum inhibitory concentrations (MICs), breakpoint versions, or unobserved isolate identities.

The output is a retrospective research summary designed to demonstrate reproducible microbiology data engineering and analysis. It is **not an official clinical antibiogram**, a validated CLSI M39 implementation, or clinical decision support. Clinical reporting references informed selected organism and antibiotic display decisions; they do not establish clinical equivalence.

## Source release and cohorts

| ARMD `culture_description` | `cohort_id` |
|---|---|
| `BLOOD` | `blood` |
| `URINE` | `urine` |
| `RESPIRATORY` | `respiratory` |

Cohorts are assigned by exact source values in `int_cohort_culture_organism`. No substring-based specimen classification is used.

Source timestamps are shifted. They support ordering within the supplied extract but do not define a shared clinical calendar year across patients. The release date and the selected reporting period are separate concepts.

### Optional 12-month reporting window

`analysis_window_start` in `dbt/dbt_project.yml` defines an inclusive shifted-date start; `dbt/macros/analysis_window.sql` computes an exclusive end exactly 12 months later. The current value is `2023-01-01`, selecting **2023-01-01 ≤ timestamp < 2024-01-01**. When the variable is absent, null, or blank, the date predicate is unrestricted; reporting models still require a non-null culture timestamp.

The filter is applied before ranking in the susceptibility branches and before patient-level aggregation in `mart_mdro`. A patient can therefore have a different selected first culture in a restricted window than in the full-release run. Filtering an already aggregated mart cannot reproduce this behavior.

Both marts publish `analysis_window_start` and `analysis_window_end_exclusive`; these are null in full-release mode. The mart keys do not include the dates, and each build replaces the configured mart tables. Retain the window fields and separate export identifiers when preserving more than one period; do not append different windows and treat a mart key as globally unique.

`int_reporting_culture_organism` remains a full-release upstream model. Notebook 14 uses it to explore alternative windows without inheriting the currently active window. Label the selected annual result a **2023 shifted-date research antibiogram**.

## Core data grains

| Model | Analytical grain | Purpose |
|---|---|---|
| `fct_culture` | Source release / patient / encounter / culture order | Prevent long-format AST rows from inflating culture counts |
| `fct_culture_organism` | Source culture–organism combination | Preserve source organism identity without inventing a biological isolate ID |
| `int_susceptibility_reconciliation` | Culture–organism–antibiotic endpoint | Classify repeated AST as resolved, missing, or conflicting |
| `fct_susceptibility` | Unambiguous source endpoint | Retain resolved source interpretations |
| `int_culture_organism_phenotype` | Source culture–organism combination | Add AST-derived resistance statuses |
| `int_reporting_culture_organism` | Release / cohort / culture / reporting organism / phenotype | Consolidate reviewed source labels after classification; full-release scope |
| `int_culture_reporting_observation` | Release / cohort / culture / reporting organism | Combine source keys across phenotypes within the window; eligible patient–organism universe |
| `int_culture_reporting_susceptibility` | AST scope / reporting source group / antibiotic | Shared reconciliation for separate overall and stratified source groups |
| `int_stratified_reporting_candidates` | Qualifying culture / organism / reporting category | Apply period and category eligibility, including MSSA marker requirements |
| `int_first_stratified_reporting_culture` | Release / cohort / patient / organism / reporting category | First qualifying culture independently within each category |
| `int_first_stratified_reporting_susceptibility` | Selected stratified culture / antibiotic | Attach reconciled AST from the matching phenotype group |
| `int_overall_first_reporting_isolate` | Release / cohort / patient / *S. aureus* | First *S. aureus* culture irrespective of phenotype |
| `int_overall_first_reporting_susceptibility` | Selected overall *S. aureus* culture / antibiotic | Attach AST reconciled across source labels on that culture |

ARMD does not include a verified isolate identifier. Model keys identify analytical records at these grains; they are not reconstructed biological isolate IDs. All downstream reporting grains refer to one configured window per build.

## Source susceptibility reconciliation

### First-stage reconciliation

`int_susceptibility_reconciliation` groups repeated source records for a culture-organism–antibiotic endpoint.

| `resolution_status` | Rule | Analytical treatment |
|---|---|---|
| `resolved` | Repeated source rows agree on one non-null susceptibility category, with no missing-category rows in that source endpoint | Retain the interpretation; eligible for `fct_susceptibility` |
| `missing` | The endpoint exists but all recorded source interpretations are `NULL` | Preserve as missing |
| `conflict` | Source rows contain different non-null categories **or** a combination of a category and a `NULL` result | Preserve the conflict; do not select a winner |

The source-level implementation treats a `NULL` as a distinct recorded result state when comparing repeated rows. For example, `Susceptible` plus `NULL` for the **same source culture-organism–drug endpoint** is classified as a conflict. This conservative rule differs from the second-stage reconciliation of **different source organism records** on a reporting culture, described below. Both stages retain their unresolved-result counts for audit.

No "latest result," "most resistant," or similar tie-break is applied because the extract does not provide enough version history to defend such a choice.

### Second-stage reconciliation within the reporting source group

Several source organism labels can map to the same reporting organism on one culture. The pipeline assigns phenotypes before consolidation and retains the source culture–organism keys. Mapping collisions alone do not remove records.

`int_culture_reporting_susceptibility` applies the same reconciliation rule to two deliberately separate scopes:

- **`overall`:** all retained source keys across phenotype labels on one *S. aureus* culture, using `int_culture_reporting_observation`.
- **`stratified`:** source keys within one culture/reporting-organism/phenotype group. AST from a different phenotype group on the same culture is not borrowed to fill a missing result.

For each source group and antibiotic:

1. Any already-conflicting source endpoint produces `conflict`.
2. Different resolved interpretations across source endpoints produce `conflict`.
3. One consistent resolved interpretation is retained, even when another source endpoint is missing; the reason is `resolved_with_missing_source_endpoint`.
4. With no resolved interpretation, the endpoint is `missing`.

The selected first-culture models attach the matching AST scope to their selected culture. A later tested culture is not substituted for a missing drug. The first-stage category-plus-null conflict rule and this cross-source-label reconciliation rule answer different questions and retain their audit reasons.

## Source AST categories

The source interpretation values retained by the pipeline are:

- `Susceptible`
- `Resistant`
- `Intermediate`
- `Inconclusive`
- `Synergism`

SQL `NULL` represents a missing interpretation. **Only `Susceptible`, `Resistant`, and `Intermediate`** enter the categorized-result denominator of the susceptibility mart. Inconclusive, synergism, missing, conflicting, and untested endpoints are tracked separately.

These are **source-defined categories**. An `Intermediate` result is not silently recoded to `Resistant`, and categorical results are not reinterpreted with a new MIC breakpoint.

## Organism terminology and reporting scope

### Version-controlled mappings

`organism_mapping.csv` maps reviewed source labels to:

- `organism_standardized`: the reviewed analytical identity;
- `organism_reporting_group`: the grouping used for reporting, which need not always be a single species;
- `source_qualifier`: source terminology retained separately from species identity;
- `mapping_status` and `include_primary_summary`: explicit review and population-selection fields.

For example, `STAPH AUREUS {MRSA}` maps to `Staphylococcus aureus`, with its source MRSA qualifier preserved. That qualifier is **not** the definitive basis for the project's AST-derived MRSA classification.

Reviewed taxonomy decisions include retaining *Citrobacter koseri* and *C. farmeri* separately from *C. freundii* complex, and mapping historical `ZZZENTEROBACTER AEROGENES` to *Klebsiella aerogenes*. Unreviewed source labels remain unmapped rather than being assigned guessed identities.

The current mapping seed has **41 reviewed source-label entries**. Reporting scope reflects the reviewed mapping, the available drug panel, cohort eligibility, and explicit exclusions; it is not determined by organism frequency alone.

### Reporting organism groups

`organism_reporting_groups.csv` contains **20 reviewed reporting groups**, including *Escherichia coli*, *Klebsiella* and *Citrobacter* species or complexes, *Staphylococcus aureus*, CoNS (unspeciated), *Enterococcus* species (unspeciated), *Morganella morganii*, and *Stenotrophomonas maltophilia*. A reporting group may represent a species, species complex, or an intentionally pooled unspeciated category.

CoNS is pooled as **coagulase-negative staphylococci (unspeciated)**; its results are not labeled as species-specific *S. epidermidis* or *S. lugdunensis* susceptibility.

Some source organisms remain outside the primary reporting mart, including organisms requiring separate AST frameworks or species-specific interpretive review. Not appearing in the mart does not mean an organism is clinically insignificant.

### Gram classification

`organism_reporting_groups.csv` defines `gram_stain_group` at the **reporting-group** level:

- `gram_positive`
- `gram_negative`

Both final marts join this reference using their base reporting-organism field. A phenotype-specific name therefore inherits the Gram classification of its underlying organism. Gram status is descriptive metadata; it does not change source grain, patient eligibility, or susceptibility denominators.

## Resistance phenotype derivation

`int_culture_organism_phenotype` evaluates resistance from **reconciled categorical AST**, not solely from phrases embedded in source organism names. The source qualifier is kept for concordance review.

Each applicable phenotype is assigned one of:

| Status | Meaning |
|---|---|
| `positive` | Sufficient qualifying resistant/non-susceptible AST evidence |
| `not_detected` | Available defining AST supports absence under the project's rule |
| `not_evaluable` | Required evidence is incomplete, conflicting, or otherwise insufficient |
| `not_applicable` | Organism is outside the phenotype's defined scope |

When collapsing source records on a culture, status precedence is `positive` → `not_evaluable` → `not_detected` → `not_applicable`. `not_evaluable` is not equivalent to susceptible or phenotype-negative.

### CRE (carbapenem-resistant Enterobacterales)

The project applies a CDC/NHSN-inspired **AST-based phenotypic CRE rule** to reviewed eligible Enterobacterales groups, including relevant *Citrobacter*, *Enterobacter*, *Escherichia*, *Klebsiella*, *Serratia*, *Proteus*, and *Morganella* mappings.

A resolved `Resistant` result to a qualifying carbapenem establishes a positive phenotype. The standard reviewed endpoints include ertapenem, imipenem, meropenem, meropenem/vaborbactam, and imipenem/relebactam where present. For *Proteus* and *Morganella*, imipenem is excluded from the project's qualifying rule because of species-specific interpretive concerns.

A missing test, unresolved endpoint, or lack of a qualifying resistant result must not automatically be interpreted as demonstrated susceptibility. The ARMD extract does not provide molecular carbapenemase results for this analysis.

### Extended-spectrum cephalosporin-resistant Enterobacterales (ESCE-R)

The project uses an **ESCE-R AST phenotype**, not a confirmed ESBL producer designation. Applicable reviewed groups include *E. coli*, *K. aerogenes*, *K. oxytoca*, *K. pneumoniae*, and *Enterobacter* groups.

Qualifying agents are cefepime, ceftriaxone, cefotaxime, ceftazidime, ceftazidime/avibactam, and ceftolozane/tazobactam. A resolved `Resistant` interpretation for any qualifying agent establishes ESCE-R under this rule.

ESCE-R and CRE statuses are derived independently. For **mutually exclusive ESCE-R display**, CRE-positive cultures take precedence, so a CRE-positive first culture is not also listed in the ESCE-R subgroup. CRE reporting records are excluded from the main susceptibility categories and remain available to the separate MDRO calculation. The same patient may contribute another eligible non-CRE culture to a main reporting category.

The source does not contain an ESBL confirmatory screen, resistance gene, or numeric MIC values suitable for proving ESBL production or rebuilding breakpoints.

### Oxacillin-resistant *Staphylococcus aureus* (MRSA phenotype)

For reviewed *S. aureus*, resistance to oxacillin or cefoxitin is positive evidence for the oxacillin-resistant phenotype (`oxacillin_resistant`). The dashboard label is **methicillin resistant (MRSA)**. This is an AST-derived project classification, not an unverified copy of the source MRSA label.

The separate **methicillin susceptible (MSSA)** subgroup (`oxacillin_susceptible`) requires an unambiguous resolved `Susceptible` oxacillin/cefoxitin marker, no conflicting or non-susceptible defining marker, and the applicable phenotype status `not_detected`. A record with unknown defining AST is not automatically assigned to MSSA.

MRSA and MSSA independently select their first qualifying culture within the window. The same patient may enter both populations on different qualifying observations, so their counts can sum to more than the overall *S. aureus* count. Unknown defining AST is not treated as MSSA.

### Vancomycin-resistant *Enterococcus* (VRE)

The VRE rule applies to reviewed *E. faecalis* and *E. faecium* only. A resolved `Resistant` vancomycin result is VRE-positive. Unspecified *Enterococcus* remains a separate reporting group and is not automatically assigned to either species-specific VRE subgroup.

### MDR *Pseudomonas aeruginosa*

The project's CDC/NHSN-derived MDR rule counts antimicrobial categories with a resolved `Intermediate` or `Resistant` result. The six reviewed categories are:

1. Extended-spectrum cephalosporins (including relevant antipseudomonal agents)
2. Fluoroquinolones
3. Tobramycin (aminoglycoside category)
4. Carbapenems
5. Piperacillin/tazobactam
6. Cefiderocol

Three or more positive categories establish MDR. Fewer than three positive categories are `not_detected` **only when the remaining unknown categories cannot raise the total to three**; otherwise the status is `not_evaluable`.

MDR *P. aeruginosa* reporting records are summarized in `mart_mdro` and excluded from the main susceptibility categories. Eligible ordinary *P. aeruginosa* records can contribute to the independently selected `base` population. Exclusion is at the reporting-record/category level; it is not a blanket exclusion of every observation from a patient who ever has MDR evidence.

## Susceptibility population: first qualifying culture within each category

### Independent category selection

`int_stratified_reporting_candidates` restricts source reporting records to the configured window and creates the ordinary (`base`), ESCE-R, VRE, MRSA, and MSSA populations under their category-specific eligibility rules. There is no published *S. aureus* `base` row: it has overall, MRSA, and MSSA rows instead. CRE and MDR *P. aeruginosa* categories are not eligible for the main susceptibility mart.

`int_first_stratified_reporting_culture` ranks separately by:

`source_version × cohort_id × patient_key × organism_reporting_group × reporting_phenotype`.

The ordering is `culture_time_shifted`, then `culture_key`, then `reporting_culture_organism_key`. Exactly one first qualifying culture is retained per category. Selection is independent of the antibiotic being summarized. Defining AST establishes phenotype eligibility (including MSSA marker requirements); the presence or result of each displayed drug does not trigger a new first-culture selection.

A patient can contribute a base culture and a different ESCE-R or VRE culture, or separately qualifying MRSA and MSSA cultures. These are overlapping patient populations. Do not sum them to estimate distinct patients, and do not require one category's denominator to be bounded by another displayed category's denominator. Each category is bounded by the eligible patient–organism universe within the same cohort and window.

`base` means the source record did not enter a higher-priority reporting category. It may include not-evaluable phenotype evidence and must not be labeled proven susceptible or wild type.

### Overall *Staphylococcus aureus*

`int_overall_first_reporting_isolate` is now restricted to *S. aureus*. It selects the earliest culture per release/cohort/patient/organism from the window-filtered `int_culture_reporting_observation`, irrespective of phenotype, using timestamp and culture key as the ordering. Its source keys span phenotype labels on that selected culture.

This overall population and the independently selected MRSA/MSSA populations answer different questions. The MRSA/MSSA rows are not filters of the overall row's selected cultures.

### Evidence of expected overlap

In Notebook 13's saved **full-release** baseline, urine *E. coli* has 34,714 base patients, 3,500 ESCE-R patients, and 1,119 patients in both categories. Blood *S. aureus* has 1,737 overall patients, 546 MRSA patients, and 1,206 MSSA patients. The sum 546 + 1,206 exceeding 1,737 is compatible with independent category selection; it is not evidence of duplicate records within a category. These figures describe the baseline, not the 2023 window.

## Reviewed antibiotic-pair policy

### Cohort-specific allow-list

**Seed:** `reportable_organism_antibiotic_pairs.csv`

The approved pair key is:

`organism_reporting_group × cohort_id × antibiotic`.

Only rows with `review_status = 'reviewed'` and `include_primary_dashboard = 1` enter the reporting grid. `review_status = 'proposed'` and `include_primary_dashboard = 0` retain candidate drugs without publishing them. Specimen-specific drugs, including urinary agents, can be approved for urine without automatically appearing in blood or respiratory reporting.

The reviewed grid is constructed **before** joining AST endpoints. An approved antibiotic stays represented when an eligible organism was not tested for that drug, so testing coverage can be measured explicitly.

### Explicit exclusions override approvals

**Seed:** `organism_antibiotic_exclusions.csv`

The exclusion key is `organism_reporting_group × cohort_id × reporting_phenotype × antibiotic`. A literal `*` in either scope column matches every cohort or category respectively. Scope fields are required; null or blank is not a wildcard. Existing species-level exclusions retain `*` in both fields.

The mart applies these rules after expanding approved pairs to the selected reporting populations, before attaching AST endpoints and aggregating counts. Any matching exclusion wins. A cohort/category-specific rule does not affect the same drug in other populations; for example, MRSA oxacillin exclusions preserve overall and MSSA oxacillin rows. Review reasons include:

| `exclusion_reason` | Meaning |
|---|---|
| `expected_resistant` | Expected/intrinsic-resistance phenotype makes routine display inappropriate |
| `ampc_reporting_suppression` | Project policy suppresses selected agents for reviewed inducible-AmpC-risk organisms |
| `inappropriate_surrogate` | Result is not appropriate for the proposed surrogate interpretation |
| `no_current_breakpoint` | No applicable current interpretive breakpoint for the reporting use |
| `phenotype_reporting_suppression` | Project-reviewed exclusion of a drug within a resistance-defined reporting category |
| `source_category_limitation` | Conventional percent-susceptible presentation is excluded because of the reviewed source-category interpretation limitation |
| `project_panel_exclusion` | Project-owner decision to omit a combination from a particular routine reporting panel |

On October 13, 2026, the project owner approved all 28 reviewed low-susceptibility rows for exclusion (27 distinct organism/cohort/drug triples because E. faecium–blood–daptomycin occurs in both base and VRE categories). Their exact cohorts and categories are enumerated in the seed with `review_basis = Project owner review 2026-10-13`. This records a project reporting decision, not a claim that every pair is intrinsically resistant or prohibited by a clinical guideline. No runtime percent-susceptible cutoff is used: the reviewed rules continue to apply if later counts change.

These are **reporting suppressions**, not edits to the original organism, AST result, or phenotype calculation. A drug may remain in the source facts and downstream reconciliation while being omitted from the reporting mart.

The AmpC policy is deliberately more restrictive than simply displaying every available categorical AST result. In the validated approved-pair audit, ceftriaxone and piperacillin/tazobactam for *E. cloacae* complex account for **six approved cohort-specific pairs** that are omitted from the mart by design.

### Inactive proposals

The current seed retains **11 inactive proposed rows** pending more detailed review of source coverage, selective testing, or organism-specific interpretation. Their inactivity must not be presented as proof that the corresponding organism–drug combinations are intrinsically resistant or never clinically relevant.

## Susceptibility reporting mart

**Model:** `mart_susceptibility`

**Grain:**

`source_version × cohort_id × organism_reporting_group × reporting_phenotype × antibiotic`.

`organism_name` is the display label. `organism_reporting_group` is the base identity used to join the approved pair and Gram-reference seeds. `reporting_phenotype` identifies overall *S. aureus* or an independently selected reporting category.

### Reporting categories

| `reporting_phenotype` | Population |
|---|---|
| `base` | First qualifying ordinary reporting culture per patient/cohort/organism; excludes *S. aureus* |
| `overall` | First *S. aureus* culture per patient/cohort irrespective of phenotype |
| `esce_r` | Independently selected first ESCE-R reporting culture; CRE takes display precedence |
| `oxacillin_resistant` | Independently selected first AST-derived MRSA reporting culture |
| `oxacillin_susceptible` | Independently selected first culture meeting the MSSA marker requirements |
| `vancomycin_resistant` | Independently selected first VRE reporting culture for reviewed *E. faecalis* or *E. faecium* |

CRE and MDR *P. aeruginosa* do not appear as main susceptibility categories. Their patient-level evidence remains in `mart_mdro`. A missing category can mean no qualifying selected population in that window; a retained row marked `insufficient_data` means a population exists but its drug-specific categorized count is below 30.

The category is part of the mart key. Patient overlap across categories is allowed, but duplicate keys within a category are not.

### Counts and identities

| Field | Definition |
|---|---|
| `susceptibility_mart_key` | Deterministic key for one row at the mart grain |
| `eligible_count` | Number of independently selected first qualifying cultures in the row's category; for `overall`, the first *S. aureus* cultures irrespective of phenotype |
| `tested_endpoint_count` | Eligible cultures with an AST endpoint for the antibiotic, regardless of whether it was resolved |
| `untested_count` | Eligible cultures with no recorded endpoint for that antibiotic |
| `categorized_count` | Resolved endpoints labeled `Susceptible`, `Intermediate`, or `Resistant` |
| `susceptible_count` | Resolved `Susceptible` endpoints |
| `intermediate_count` | Resolved `Intermediate` endpoints |
| `resistant_count` | Resolved `Resistant` endpoints |
| `inconclusive_count` | Resolved `Inconclusive` endpoints |
| `synergism_count` | Resolved `Synergism` endpoints |
| `missing_result_count` | Endpoints classified `missing` after reconciliation |
| `conflict_result_count` | Endpoints classified `conflict` after reconciliation |
| `gram_stain_group` | Reference-data classification of the base reporting group |
| `pair_order` | Review-defined ordering of antibiotic rows in the reporting panel |
| `analysis_window_start` | Inclusive shifted-date start, or null for full-release mode |
| `analysis_window_end_exclusive` | Start plus 12 months, excluded from the window; null for full-release mode |

The required row-level reconciliations are:

```text
eligible_count = tested_endpoint_count + untested_count

categorized_count = susceptible_count + intermediate_count + resistant_count

tested_endpoint_count = categorized_count + inconclusive_count
                      + synergism_count + missing_result_count
                      + conflict_result_count
```

The same `eligible_count` is used for every **approved, nonsuppressed drug** within an organism/cohort/**reporting phenotype** population. It can differ across phenotypes; no antibiotic-specific first-isolate selection is applied.

### Percentages

All percentages are calculated from the relevant mart row's counts and rounded to one decimal place.

| Field | Formula | Interpretation |
|---|---|---|
| `categorized_result_coverage_pct` | `100 × categorized_count / eligible_count` | Proportion of the eligible population with interpretable S/I/R AST for this drug; **not** a susceptibility rate |
| `susceptible_pct` | `100 × susceptible_count / categorized_count` | Proportion of interpretable source results labeled susceptible |
| `resistant_pct` | `100 × resistant_count / categorized_count` | Proportion explicitly labeled resistant; intermediate is not reclassified |

Rates with no usable denominator are null rather than fabricated zeroes. Selective antimicrobial testing may make `susceptible_pct` unrepresentative of the full eligible organism population even when the categorized count meets the display threshold.

### Sample-size flag and exploratory display

The mart retains `display_status = 'display'` when `categorized_count >= 30`; otherwise it is `insufficient_data`. This flag does not remove rows or change their calculated percentages.

The Tableau and ggplot/matplotlib antibiograms intentionally show available susceptibility percentages **without filtering on `display_status`**, so low-count organisms and phenotype populations remain visible. Interpret each percentage with its drug-specific `categorized_count`; values flagged `insufficient_data` have limited sample support. Approved-pair exclusions still apply, and a result with no categorized denominator remains blank rather than zero.

The 30-result flag is a **project reporting convention**, not certification of CLSI M39 compliance. These unfiltered figures are exploratory portfolio outputs, not clinical reporting. A genuine 0% susceptibility differs from both a missing denominator and a percentage based on few results.

## MDRO reporting mart

**Model:** `mart_mdro`

The MDRO mart reports **CRE** and **MDR *P. aeruginosa*** independently of the category-specific susceptibility table.

**Grain:** `source_version × cohort_id × organism_base_name × mdro_phenotype`.

Its patient-level denominator is constructed from the reviewed base-organism population **across eligible culture observations in the configured window**, not only from the patient's selected first overall culture and not by summing phenotype-specific categories.

For each CRE or MDR *P. aeruginosa* calculation, the model first requires a non-null timestamp within the configured window and excludes culture records where the phenotype status is `not_applicable`, groups the remaining records by source version, cohort, patient, and base reporting organism, and assigns **one patient-level status**. Precedence across eligible observations is:

1. `positive` if **any** observation is positive;
2. otherwise `not_evaluable` if **any** observation is not evaluable;
3. otherwise `not_detected` if any observation supports not detected.

A patient therefore contributes at most once to each cohort/base-organism/MDRO analysis, even with multiple cultures or phenotype-qualified source labels. This is an **any-positive-across-observations** metric, not a first-isolate measure or an incidence rate.

| Field | Definition |
|---|---|
| `denominator_count` | Eligible distinct patient/base-organism groups for the MDRO phenotype and cohort |
| `positive_count` | Patient groups classified positive for the specified phenotype |
| `not_detected_count` | Patient groups with `not_detected` evidence and neither positive nor not-evaluable evidence |
| `not_evaluable_count` | Patient groups with no positive observation and at least one not-evaluable observation, even if other observations are not detected |
| `evaluable_count` | `positive_count + not_detected_count` |
| `phenotype_pct` | `round(100 × positive_count / denominator_count, 1)`; percentage of all eligible patient groups classified phenotype-positive |
| `evaluable_pct` | `round(100 × evaluable_count / denominator_count, 1)`; percentage of eligible patient groups with a determinate positive or not-detected status |
| `gram_stain_group` | Base-organism Gram classification from the same reference seed |

Required reconciliations:

```text
denominator_count = positive_count + not_detected_count + not_evaluable_count

evaluable_count = positive_count + not_detected_count
```

Both `phenotype_pct` and `evaluable_pct` use **`denominator_count`**, which includes `not_evaluable_count`. Consequently, `evaluable_pct` reports **AST phenotype evaluability/ascertainment**, not the percentage that is phenotype-positive among evaluable patients. The SQL does **not** calculate `100 × positive_count / evaluable_count` as a separate mart field.

Division is protected by `nullif(denominator_count, 0)`, and the two percentages are rounded to one decimal place. These MDRO measures are distinct from `mart_susceptibility`'s S/I/R percentages and should be labeled separately in Tableau.

## Historical worked susceptibility example

**ARMD release:** `2025-04-11`  
**Cohort:** `blood`  
**Organism:** *Escherichia coli*  
**Antibiotic:** ceftriaxone

**Historical scope:** this saved example belongs to the earlier full-release, overall-population/subset implementation. It is retained to explain the percentage arithmetic. The current implementation has a `base` *E. coli* population and independently selected ESCE-R population; these numbers are not current 2023 results.

| Metric | Overall *E. coli* | ESCE-R subgroup |
|---|---:|---:|
| `eligible_count` | 3,022 | 467 |
| `categorized_count` | 3,016 | 466 |
| `susceptible_count` | 2,542 | 1 |
| `resistant_count` | 465 | 461 |
| `intermediate_count` | 9 | 4 |
| `categorized_result_coverage_pct` | 99.8% | 99.8% |
| `susceptible_pct` | 84.3% | 0.2% |
| `resistant_pct` | 15.4% | 98.9% |

For the overall row:

```text
2,542 susceptible + 465 resistant + 9 intermediate = 3,016 categorized

100 × 2,542 / 3,016 = 84.3% susceptible (rounded)
```

In that historical implementation, the ESCE-R row was a subset of the selected overall first cultures. The current first-within-category methodology supersedes that selection rule while retaining the percentage formulas. Category counts must not be added as unique patients in either example.

## Validation and reproducibility

The current dbt tests cover unique keys, reference relationships, accepted values, first qualifying culture selection, window boundaries, reconciled AST states, endpoint count identities, reporting-population denominators, approved-pair completeness, suppression rules, and Gram metadata.

The custom assertions in `dbt/tests/` are:

- `assert_first_stratified_qualifying_culture`: independently reconstructs the first culture in each qualifying category.
- `assert_stratified_denominators_within_population`: bounds each category by the window-eligible patient–organism universe, not by a different displayed category.
- `assert_overall_susceptibility_denominator`: compares overall *S. aureus* mart denominators with the selected overall population.
- `assert_reporting_mart_population_complete`: compares the full expected approved-pair grid with the actual mart in both directions.
- `assert_reporting_window_selection`: checks that selected overall and stratified culture timestamps lie within the configured window; the date predicate is unrestricted in full-release mode.
- `assert_no_mdro_in_main_susceptibility`: checks exclusion of CRE/MDR reporting categories and their selected-record evidence.
- `assert_no_suppressed_or_unapproved_pairs`: enforces the allow-list and exclusion overrides.
- `assert_scoped_exclusions_have_approved_pairs`: checks that cohort/category-specific exclusions identify an active reviewed antibiotic pair.

Notebook 13 supplies the phenotype inventory, overlap checks, *S. aureus* population inspection, MDRO summary, AST-scope review, and display-coverage results. Notebook 14 documents preselection availability by shifted year and month. Their saved outputs include both full-release and 2023 evidence; their Markdown labels the scope of each.

### Recorded build and current exports

The saved full `dbt build` results were generated on **October 13, 2026, at 11:30 p.m. Eastern** (`2026-10-14T03:30:57Z`): **22 models, 4 reference seeds, and 229 passing data tests** — **255 successful operations**, with no error, warning, or skipped result statuses. The saved manifest uses the 2023 window; its model SQL and reference-seed checksums match the current files. These are recorded results, not a new warehouse build performed during this documentation update.

The current exports in `data/exports/` record `analysis_window_start = 2023-01-01` and `analysis_window_end_exclusive = 2024-01-01`:

| Export | Rows | Status |
|---|---:|---|
| `mart_susceptibility.csv` | 396 | 187 rows meet the 30-result flag; 209 are flagged `insufficient_data` |
| `mart_mdro.csv` | 32 | Patient-level CRE and MDR *P. aeruginosa* summaries |

The reference seeds contain 41 organism-mapping rows, 20 organism/Gram groups, 311 active approved cohort-specific pairs, 11 inactive proposals, and **61 exclusion rules: 33 broad plus 28 cohort/category-specific rules**. These rules are applied to the category-specific reporting grid; their counts are not expected mart-row totals. The MDRO model does not depend on the susceptibility exclusion seed.

### Current sample-size coverage by cohort

| Cohort | Mart rows | `categorized_count >= 30` | `insufficient_data` | Meeting flag (%) |
|---|---:|---:|---:|---:|
| Blood | 147 | 56 | 91 | 38.1% |
| Respiratory | 114 | 31 | 83 | 27.2% |
| Urine | 135 | 100 | 35 | 74.1% |
| **Total** | **396** | **187** | **209** | **47.2%** |

These percentages measure the fraction of reporting rows meeting the 30-result flag, not susceptibility or patient-level testing coverage. The exploratory antibiograms retain available percentages from both flag states. Low counts do not establish absence of an organism or resistance phenotype.

### Saved notebook evidence and historical scope

- Notebook 07 examines source vocabulary and reference-policy checks.
- Notebook 10 checks overall *S. aureus* selection against the earliest culture in the window-filtered `int_culture_reporting_observation` population.
- Notebook 11 checks the reporting grid, counts, percentages, denominators, exclusions, sample-size flags, and Gram classification. Saved tables should be read with their recorded window and reference-policy scope.
- Notebook 12 independently reconstructs patient-level MDRO status from upstream observations within the distinct start/end pair recorded in `mart_mdro`. Its saved result reports **`pass`, 32 rows, and the bounds `2023-01-01` and `2024-01-01`**. Filtering precedes both CRE/MDR branches and patient aggregation. The reconstruction uses the mart's recorded bounds; it does not independently prove the intended configuration.
- Notebook 13 preserves historical full-release and pre-exclusion-update comparisons: **455 susceptibility rows** in the full-release baseline, **424 susceptibility rows** in the earlier 2023 export, and **34 MDRO rows** in the full-release baseline. These are distinct historical snapshots, not expected row counts for the current exports. The 28 scoped exclusion rules reduce the earlier 2023 susceptibility result from 424 to the current 396 rows.
- Notebook 14 explores shifted-year and monthly coverage before first-culture selection. Its saved monthly check has observations in all 12 shifted months for each cohort. That demonstrates monthly representation, not completeness of source capture or a clinical calendar-year population.

All snapshot counts depend on the source release, model version, reference seeds, and reporting window. Preserve window labels in exports and distinguish historical results from the current build.

## Interpretation limitations

- **Not a clinical antibiogram.** The project uses selected reporting conventions but has not established full CLSI M39/M100 compliance, clinical laboratory authorization, or clinical decision-support validity.
- **No verified isolate identifier.** An analytical culture/reporting-organism key is not proof that multiple biological isolates have been resolved.
- **Shifted dates and project period.** The release does not represent an unshifted calendar-year institutional antibiogram. First-qualifying-culture selection is a documented project policy within cohort, category, and configured window; overall selection is retained for *S. aureus*.
- **Source-provided S/I/R.** MICs, original testing methods, breakpoint editions, confirmatory AST, and resistance genes cannot be reconstructed from the available categorical fields.
- **Phenotypes are derived, not diagnoses.** ESCE-R is not equivalent to confirmed ESBL; CRE, MRSA, VRE, and MDR classifications depend on the qualifying AST present in the source.
- **Incomplete and selective testing.** Missing drugs, conflicts, and selectively tested agents can bias susceptibility percentages. `categorized_result_coverage_pct` should be interpreted alongside `%S`.
- **Mixed taxonomic detail.** Reporting groups include species, complexes, and unspeciated categories. Their pooled rates should not be interpreted as species-specific susceptibility distributions.
- **Specimen cohort stratification.** Blood, respiratory, and urine analyses are distinct; a patient can appear in more than one cohort. Cohort results should not be added without a new deduplication policy.
- **Overlapping phenotype populations.** A patient can enter multiple independently selected categories; adding their denominators does not produce a unique-patient total. `base` is not proof of wild-type susceptibility.
- **MDRO rates use another denominator.** `mart_mdro` uses patient-level phenotype status across observations; `mart_susceptibility` uses the first qualifying culture within each reporting category, plus overall *S. aureus*. Their counts and percentages must not be directly equated.
- **Reporting controls are explicit project decisions.** An excluded organism–drug pair is not necessarily untested, intrinsically resistant, or inappropriate in every clinical circumstance.

Differences across cohorts are descriptive and should not be interpreted as treatment effects or estimates of resistance rates in a general hospital population.

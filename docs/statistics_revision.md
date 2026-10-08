# Statistics revision for the manuscript

## Scope and execution status

The plot scripts already performed statistical comparisons. This revision fixes
sample grouping, count denominators, statistical evidence exports and the LEfSe
implementation, and adds PERMDISP. It does not rerun Kraken2 or change reference
databases. Earlier figures and p-values should be regenerated before adopting the
revised methods paragraph.

Validation includes 19 regression checks and an integration run on the supplied
`Taxonomy_FISH_Cumulative_Reads.xlsx` with `metadata.xlsx`. The FISH integration
uses 99 permutations for verification, not the manuscript default of 999.
It processes 26 sample columns and excludes three specimens with unrecorded
captivity time from temporal inference. Domain rows are absent from that input.
LEfSe input preparation was verified; the official R-dependent analysis was not
executed in the validation environment. Full microbiome reanalysis requires
`taxonomic_classification_clean.xlsx`, not the per-sample classification totals.

## Corrections

| Stage | Revised behavior | Evidence output |
|---|---|---|
| Metadata | Exact canonical sample matching; confirmed NCM2wild/Ncm2 alias; no prefix matching or quartile bins | `*_samples.csv` |
| Time | Recorded 0d, 7d, 30d; a wild specimen with missing Time is a documented derived 0d baseline; lab specimens with missing Time are excluded | Inclusion and derivation reason per specimen |
| Duplicates | Identical duplicate metadata rows represent one sample; conflicting metadata or colliding count columns raise errors | Explicit exception |
| Rarefaction | Multivariate hypergeometric sampling without replacement; 50 grid depths, 10 repetitions, seed 42; no read-population expansion | Curves, resampling SD, common-depth richness, execution record |
| Abundance | Full count denominator at the selected rank; low-abundance taxa retained as Other; groups average per-sample proportions | Proportions CSV and summary workbook |
| Alpha diversity | Observed richness, Shannon (natural log), bias-corrected Simpson 1-D and Chao1 | Per-sample indices, KW and pairwise tables |
| Alpha inference | Kruskal-Wallis; BH across four indices within each rank; bilateral Mann-Whitney U after significant adjusted omnibus; Bonferroni within each index | Full-precision statistics and adjusted p-values |
| Beta diversity | Bray-Curtis on per-sample proportions; PCoA with Lingoes correction when needed | Original distances, coordinates and correction diagnostics |
| Community inference | Official scikit-bio PERMANOVA and secondary ANOSIM, 999 permutations, seed 42; original Bray-Curtis distances | Global F/R, R², p-values and sample counts |
| Dispersion | Official scikit-bio PERMDISP with spatial medians on the Lingoes-corrected distance representation | Dispersion statistic and permutation p-value |
| Pairwise ANOSIM | BH correction across group pairs | Full-precision pairwise CSV |
| Taxon abundance | ANOVA and Tukey retained as exploratory; BH across displayed taxa gates Tukey; corrected compact letters | Omnibus and posthoc CSVs |
| LEfSe | Official SegataLab executables, KW alpha .05, Wilcoxon setting .05, LDA threshold 2, 30 bootstraps, .67 bootstrap fraction, strict one-against-one multiclass strategy | Official `.res`, mapped results, command log and parameters |
| Cladogram | Correlation-derived parent links removed; rank-only tables cannot define a lineage tree | Requests for unsupported cladograms fail explicitly |
| Provenance | Actual Python, OS, package versions, timestamp, input SHA256, git commit/dirty state and parameters | `*_run.json` |

No detections and missing reports are different conditions. Missing numeric
counts raise an error. Samples with zero totals at the analyzed rank have no
diversity evidence and are excluded from normalized abundance and group tests.
Unrecorded groups never enter temporal inference, even when displayed as Unknown.

The complete denominator is all assigned counts at the selected rank, not all
sequenced reads and not unclassified reads. Never combine parent and descendant
rows across ranks into a single abundance denominator. Alpha indices are
calculated from unrarefied rank counts; rarefaction is a separate diagnostic.
Unequal depth and classifier/reference coverage must be considered when
interpreting richness, particularly sparse diet databases. The automatic common
rarefaction depth is the minimum positive rank total, which can be very small.

Unrestricted permutations and independent-sample tests assume different,
independently sampled specimens. If the same individual was measured repeatedly,
or multiple specimens share experimental tanks, record subject/tank IDs and use
an appropriate blocked or hierarchical design before interpreting temporal tests.
A minimum of two samples per group is enforced; very small groups still limit
precision and asymptotic rank-test inference. A non-significant dispersion test
does not prove equal dispersion; interpret it together with sample sizes and
PERMANOVA effect size.

No ANCOM-BC or ALDEx2 execution was recovered or added. Do not describe the
exploratory ANOVA/Tukey screen as compositional differential-abundance inference.
Official LEfSe cannot retroactively validate old custom scores. Subclasses are
not supplied here, so subclass consistency is not claimed. The official CLI
exposes no seed option; exact replay of its stochastic LDA step is not promised.

## Reanalysis commands (fish-compatible)

Run from the repository root after checking out this revision. An isolated
environment keeps the previous installation available for provenance.

```fish
python -m venv .venv-statistics
source .venv-statistics/bin/activate.fish
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python scripts/run_manuscript.py \
  --data results/final_tables/taxonomic_classification_clean.xlsx \
  --metadata data/metadata.xlsx \
  --output results/manuscript_review \
  --rank genus --abundance-rank phylum --category Time
```

For Bash, activate with `source .venv-statistics/bin/activate` instead.
The core revision was verified on Python 3.12.14 with numpy 2.3.5,
pandas 2.2.3, scipy 1.17.0, statsmodels 0.15.0, scikit-bio 0.7.3,
matplotlib 3.10.8 and openpyxl 3.1.5. These are validation versions,
not inferred versions of the original analyses.

After installing official [SegataLab LEfSe](https://github.com/SegataLab/lefse)
and its R/rpy2 dependencies in an appropriate environment, repeat the same
workflow with `--run-lefse`. Check `workflow_status.csv` and the LEfSe run record.
Prepared inputs alone are not completed LEfSe results. If Domain is absent,
rebuild the taxonomy matrix from the original reports with the revised
`03_generate_table.py`, which now retains rank D and K. Do not infer Domain by
adding every phylum row without a verified parent mapping.

For diet tables, repeat with the FISH, MACROINV and PLANTS matrices in separate
output directories. Their normalized percentages are database-specific signals,
not proportions of ingested food biomass, and should not be pooled as if they
shared a validated common abundance denominator.

## Methods text after completing the revised reanalysis

Use the following text only for analyses actually rerun with this revision and
after confirming the independent-specimen design. Exclude the LEfSe sentence
until an official completed result is available. Runtime versions are cited in
the supplementary software table generated from the matching run records.

**2.12 Statistical analysis.** Sampling depth was assessed using taxonomic
rarefaction curves generated by repeated subsampling without replacement from
per-sample counts at the analyzed rank, using 50 grid depths, 10 repetitions per
depth and a random seed of 42. Curves report mean observed richness and
resampling standard deviation. Relative abundances were calculated using all
assigned counts at the selected rank; low-abundance taxa were combined as
Other, and group summaries were arithmetic means of per-sample proportions.
Samples were grouped by documented captivity time (0, 7 or 30 days). A wild
specimen without a recorded captivity time was assigned to the zero-captivity
baseline, with this derivation retained in the sample audit; lab specimens with
unrecorded time were excluded from temporal comparisons. Samples with no
detections at the analyzed rank were excluded from normalized analyses.
Observed richness, Shannon diversity (natural logarithm), bias-corrected Simpson
diversity (1-D) and bias-corrected Chao1 were calculated from unrarefied counts.
Independent time-point groups were compared using Kruskal-Wallis tests, with
Benjamini-Hochberg correction across the four indices within each rank.
Significant adjusted omnibus tests were followed by two-sided Mann-Whitney U
comparisons with Bonferroni correction across group pairs within each index.
Bray-Curtis dissimilarities calculated from per-sample relative abundances were
visualized by principal coordinates analysis, applying a Lingoes correction when
negative eigenvalues were present. Differences in community composition were
tested by PERMANOVA (Anderson 2001) on the original Bray-Curtis distances, with 999 unrestricted
permutations and a random seed of 42; pseudo-F, R² and permutation p-values were
reported. Multivariate dispersion was assessed by PERMDISP using spatial
medians in the corrected coordinate representation (Anderson 2006) and the same permutation
settings. ANOSIM was retained as a secondary comparison, with
Benjamini-Hochberg correction across pairwise group tests. Taxon-level ANOVA
and Tukey HSD screens were treated as exploratory, with Benjamini-Hochberg
correction of omnibus tests across displayed taxa. Analyses were performed in
Python using NumPy, pandas, SciPy, statsmodels, scikit-bio and Matplotlib;
openpyxl supported workbook input/output. Software versions, code revision,
sample inclusion and execution parameters are reported in the supplementary
methods.

**Optional LEfSe sentence, after official execution:** LEfSe (Segata et al. 2011)
was performed using the official implementation as an exploratory screen, with
a Kruskal-Wallis alpha of 0.05, a Wilcoxon setting of 0.05, an LDA score threshold
of 2.0, normalization to 10⁶ per sample, 30 bootstrap iterations and a 0.67
bootstrap fraction; the one-against-one multiclass strategy was used. No
subclasses were specified. Results were interpreted cautiously because
differential-abundance methods can yield discordant results (Nearing et al. 2022).

## Figure placement

| Results section | Main figure | Supporting material |
|---|---|---|
| 3.1 Sequencing and sampling depth | Rarefaction, with depth units and rank stated | QC/assembly tables; classify counts; full curves |
| 3.2 Community composition | Domain overview and phylum/genus relative-abundance bars | Alpha boxplots and per-sample counts |
| 3.3 Beta diversity | Bray-Curtis PCoA | Full global, dispersion and pairwise statistics |
| 3.2/3.3 Exploratory taxa | Official LEfSe LDA bars only after execution | All LEfSe outputs and exploratory taxon abundance screens |
| 3.4 Diet results | Separate FISH, MACROINV and PLANTS detection/abundance displays | Complete database-specific tables and sample evidence |
| 3.5 AMR/virulence results | Separate detection heatmaps for resistance and virulence | Full hits and per-database presence matrices; missing CARD report as NA |
| 3.6 Xenobiotic results, if included | Sample-by-pathway or KO heatmap | Distinct gene/KO counts; pathway associations without completeness claims |

For small groups use boxplots with all specimens visible. Select figures by
the biological question and validated data, not appearance or significance.
Boxplot whiskers follow the usual 1.5×IQR rule. Mean ± SE overlays, if retained,
must be distinguished from medians and IQR in captions. Compact letters share
only non-rejected pair comparisons; equal letters do not demonstrate equality.

## Primary references

- Anderson MJ. 2001. A new method for non-parametric multivariate analysis of variance. *Austral Ecology* 26:32–46. https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x
- Anderson MJ. 2006. Distance-based tests for homogeneity of multivariate dispersions. *Biometrics* 62:245–253. https://doi.org/10.1111/j.1541-0420.2005.00440.x
- Segata N et al. 2011. Metagenomic biomarker discovery and explanation. *Genome Biology* 12:R60. https://doi.org/10.1186/gb-2011-12-6-r60
- Nearing JT et al. 2022. Microbiome differential abundance methods produce different results across 38 datasets. *Nature Communications* 13:342. https://doi.org/10.1038/s41467-022-28034-z
- scikit-bio PERMANOVA/PERMDISP documentation: https://scikit.bio/docs/latest/generated/skbio.stats.distance.html
- Official LEfSe source: https://github.com/SegataLab/lefse

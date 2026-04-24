# Replication package: Comment on Ciacci (2024)

Re-examines results from:

> Ciacci, R. (2024). "Banning the purchase of sex increases cases of rape:
> evidence from Sweden." *Journal of Population Economics* 37 (2), Article 36.
> Ciacci, R. (2025). Follow-up working paper.

Produces the tables and figures in the 2nd-round comment manuscript
(`manuscript 2nd round.docx`, `Appendix RR2.pdf`).


## Folder layout

```
code/
    replication.do            Tables 2, 3, 4, B1, B2, C2, C3, C5 and
                              Figures 1, 3, B1, B2, B4, B5, B6.
    replication_fuzzydid.do   Table C1 (Fuzzy DiD). Split out because
                              `fuzzydid` with breps=99 takes ~30 min.
    run_all.bat               Windows helper that calls both do-files.
data/
    Sweden_tables.dta                        main monthly region panel.
    Sweden_tables3.dta                       adds airport-flights instruments
                                             (used only by Table C5).
    yearly national data rape reports.dta    used in Figure 1.
    annual data nordic countries.dta         used in Figure 3.
output/                                      generated .pdf figures and .tex tables.
```


## Data sources

- Rape reports and monthly police personnel by region: Swedish National
  Council for Crime Prevention (BRA, <https://bra.se>).
- Regional population: Statistics Sweden (SCB).
- Sex-purchase prosecutions: Swedish Prosecution Authority (Åklagarmyndigheten).
- Nordic cross-country rape counts: national statistical agencies of
  Denmark, Finland, Norway, Sweden.

The `.dta` files here were copied verbatim from Ciacci's 2024 replication
archive and from the 2025 follow-up archive. Nothing is downloaded at
run time; Block 1 of `replication.do` stops with an error if any raw file
is missing.


## How to run

1. Install Stata 18 (MP or SE). The setup block auto-installs the SSC
   packages `rdrobust`, `reghdfe`, `outreg2`, `tuples`, `fuzzydid`.
2. **Set the project root.** Open each do-file and edit the line
   marked `>>> EDIT THIS LINE <<<` near the top (Block 0):

   ```stata
   global root "C:/Users/c4041171/Dropbox/Ciacci replication/replication_220426"
   ```

   Replace the path with the absolute path to wherever you placed the
   replication package on your machine. Use forward slashes (`/`) — they
   work on Windows, macOS, and Linux. All other paths in the do-files
   (`${root}/data/...`, `${root}/output/...`) are derived from this one
   variable, so this is the only line that needs editing.
3. Run each do-file in Stata (or use `run_all.bat` on Windows).

### Runtime (per file, run alone)

| Script | `B_REPS=9` | `B_REPS=99` |
| --- | --- | --- |
| `replication.do` | ~3 min | ~5 min |
| `replication_fuzzydid.do` | ~3 min | ~30 min |

`global B_REPS = ...` at the top of each script controls replications.
Use low values for development, 99 for production.


## Generated exhibits

### Main paper
- **Figure 1** — yearly national rape reports.
- **Figure 3** — Nordic countries comparison.
- **Table 2** — RDiT estimates; cluster bootstrap on cols 3 & 6.
- **Table 3** — RDiT alternative specifications (six `rdrobust` cols).
- **Table 4** — sex-purchase prosecutions by year.

### Appendix
- **Figure B1** — RDiT residual plot, whole sample, wide window.
- **Figure B2** — RDiT residual plot, restricted sample, narrow window.
- **Figure B4** — RDiT residual plot, Ciacci (2025) season coding.
- **Figure B5** — four `rdplot` panels (uniform vs. triangular × raw vs. month-residualised).
- **Figure B6** — placebo-threshold treatment effects, two kernels.
- **Table B1** — RDiT replication with optimal bandwidth
- **Table B2** — D-RDiT counterpart of Table B1.
- **Table C1** — Fuzzy DiD (from `replication_fuzzydid.do`).
- **Table C2** — event study with out-of-window dummy.
- **Table C3** — event study with alternative clustering.
- **Table C5** — IV, three specifications, first-stage F reported.


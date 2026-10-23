# Open OnDemand manual check (step 5, about 30 minutes)

The batch jobs run the same R as the OOD app, but not the app itself. Do these steps once per kit run in a browser, as a learner would, and write the answers into `$RUN/records/listings/ood_manual_check.txt` on Negishi (collect.sh picks it up). Copy the template at the end of this page.

1. Open `https://gateway.negishi.rcac.purdue.edu`, then Interactive Apps, Bioinformatics Apps, **RStudio (bioconductor)**. Fill the form exactly as `learners/setup.md` says: partition cpu, account (`OOD_ACCOUNT`), QoS normal, 8 hours, 16 cores.
   - Record the app's exact name, every field on the form with its default, and whether there is a memory field.
   - If the form differs from `learners/fig/ood_rstudio_resources.png`, take a new screenshot at the same size.
2. Record how long it takes from Launch to the "Connect to RStudio Server" button.
3. In the RStudio console:

   ```r
   R.version.string; BiocManager::version(); .libPaths()
   packageVersion("Seurat"); packageVersion("SeuratObject")
   Sys.getenv(c("http_proxy", "https_proxy", "R_LIBS_USER", "TMPDIR"))
   ```

4. In the RStudio Terminal tab (a shell inside the app's session), run one cold episode and the egress check in the app's own R:

   ```bash
   bash ~/singlecell-rnaseq/negishi-run/run_R_chain.sh --run 2026-10-20 --pass cold --episodes 03 --r current
   bash ~/singlecell-rnaseq/negishi-run/run_R_chain.sh --run 2026-10-20 --episodes env --r current
   ```

   These overwrite the batch results for `c03` and `env-egress` in the same run, which is intended: the app is the reference. If the batch results differ from what the app gives, say so in the notes.
5. Paste the episode 3 `load-data` and `vlnplot` chunks into the console by hand. Note anything a learner would see that the batch run cannot show: pop-ups, package update prompts, the Plots pane size.

Template for `ood_manual_check.txt`:

```text
date:
app name on the menu:
form fields and defaults:
memory field present (y/n), value:
launch to connect (min):
R.version.string:
Bioconductor:
.libPaths():
Seurat / SeuratObject:
proxy variables:
run_R_chain c03 --r current: PASS/WARN/FAIL, minutes:
run_R_chain env --r current: failed calls:
notes:
```

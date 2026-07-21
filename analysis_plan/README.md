# Analysis plans — multiPPIseq DMS pipeline

Lab-notebook / knowledge-base for the multiPPIseq DMS pipeline. Each
`*_analysis_plan.md` is the present-tense **methods spec** for one analysis and
doubles as (a) a session-onboarding doc to load into a fresh session and (b) a
source for the manuscript methods. Where an analysis has accumulated decisions, a
sibling `*_decisions.md` holds the dated "why/when" log so the methods file stays
clean.

**Conventions**
- Keep the analysis name as a specific filename prefix (e.g. `dict_build_*`) so
  files for one analysis sort together in this flat folder.
- Methods file = "what we do now" (rewrite in place, prune dead-ends). Decisions
  log = "why / when it changed" (dated, append-only).
- If a fact helps a new session act correctly today → methods file. If it explains
  why we got here → decisions log.
- All analysis runs on the cluster via `qsub` — never on the login node.

## Part 1 — Dictionary building
| File | Role |
|---|---|
| [dict_build_analysis_plan.md](dict_build_analysis_plan.md) | Methods + onboarding. Start with its "Status / orientation" block. **Canonical:** Python lineage (`step1_v2.py` → bowtie2 → `step2.py` → `step3.py`), driven by `qsub/dict_build/` steps 00–12. Current batch: `multiPPIseq_DMS_dict_071526`, both libtypes, off union-merged reads. |
| [dict_build_decisions.md](dict_build_decisions.md) | Dated decision log. |

## Part 2 — cDNA edit count
| File | Role |
|---|---|
| [cdna_edit_count_analysis_plan.md](cdna_edit_count_analysis_plan.md) | Methods + onboarding. Applies a barcode→variant dictionary to cDNA reads and counts edits per variant. Driven by `qsub/cdna_edit_count/` → `script/cdna_edit_count/`. |
| [cdna_edit_count_decisions.md](cdna_edit_count_decisions.md) | Dated decision log. |

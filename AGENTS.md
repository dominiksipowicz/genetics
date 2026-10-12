# Repository Guidelines

## Project Structure & Module Organization

Core scripts in the repo root (`gene_stats.py`, `compare.py`, `compare2.py`, `simple_start.py`) parse and compare Nebula VCF exports. Person-specific folders (`Anna/`, `Babcia/`, `Dom/`, `Jadwiga/`) contain raw VCF/CRAM files alongside generated CSV or text summaries; keep names stable because scripts load them via hard-coded constants. `README.md` lists CSV headers, and `package.json` is currently unused but reserved for future tooling.

## Build, Test, and Development Commands

- `python -m venv .venv && source .venv/bin/activate`: set up a Python 3 virtualenv.
- `/opt/homebrew/bin/python3.11 family_genetics_analysis.py` quick python run
- `pip install pandas cyvcf2`: install required runtime dependencies.
- `python gene_stats.py`: parse the configured VCF and write `parsed_vcf_data.csv`.
- `python compare.py`: regenerate shared ID reports and variant counts.
- `python compare2.py`: recompute PASS SNP overlap metrics for regression checks.

## Coding Style & Naming Conventions

Use PEP 8 defaults: 4-space indentation, snake_case files and functions, and CapWords for classes if added. Group path constants near the top of each module and comment why switches are made when targeting a different sample. Favor pure helper functions returning dataframes over inline prints so new analyses can be chained or unit tested.

## Testing Guidelines

There is no automated suite yet; when adding logic, craft small VCF fixtures (strip to a few variants) and exercise them through ad-hoc scripts or `pytest` tests under a new `tests/` directory. Before opening a PR, rerun `python gene_stats.py` and `python compare2.py` on representative inputs and note record counts or Jaccard deltas in the PR description.

## Commit & Pull Request Guidelines

Follow the existing short, imperative commit style (“add Babcia”, “compare genes”) and keep generated artifacts in the same commit as the code that produced them. PRs should describe the datasets touched, summarize new outputs, and link relevant issues or discussions. Include screenshots or metrics when variant counts shift so reviewers can confirm the change quickly.

## Security & Data Handling

Treat genomic data as sensitive: avoid committing fresh Nebula downloads without approval, scrub personal identifiers from derived tables, and extend `.gitignore` for temporary exports. Use Git LFS or external storage for large binaries, and rotate any exposed credentials or signed URLs immediately.

# AGENTS.md

## Project

NGS Data Analysis Handbook — Zensical (Material for MkDocs fork) static site deployed to GitHub Pages (`https://ilypopv.github.io/NGS-Handbook/`). Cheat-sheets for QC, variant calling, pangenomics, phylogenetics, 16S. No application code; `src/ngs_handbook/__init__.py` is empty placeholder.

## Stack & Structure

- **Builder:** Zensical `>=0.0.57` (Rust-based MkDocs). Config is `zensical.toml` — not `mkdocs.yml`. `pyproject.toml:12` defines `uv_build` backend; requires Python `>=3.12` (`uv.lock`).
- **Content:** `docs/` — `index.md` + 5 chapters: `QC/`, `VarCall/`, `Pangenome/`, `Phylogenetics/` (8 pages `04_01`–`04_08`), `16S_amplicon_analysis/` (3 pages). Assets in `docs/assets/`, `docs/imgs/`, `docs/stylesheets/extra.css`.
- **Theme overrides:** `overrides/main.html` + `overrides/partials/`; `zensical.toml:40` sets `custom_dir = "overrides"`. Don't edit generated `site/`.
- **Nav:** Single source of truth is `zensical.toml:14-36` (`nav` array). Adding/renaming a page requires updating `nav` there.
- **Envs:** `envs/*.yaml` — per-chapter conda envs (`qc.yaml`, `varcall.yaml`, `panacota.yaml`, `phylo.yaml`), channels `bioconda`/`conda-forge`/`defaults` + `pip`/`ipykernel`. Data fixtures in `data/` (zips).

## Commands

```bash
uv sync --group dev          # install zensical from [dependency-groups.dev]
uv run zensical serve        # local preview with live reload
uv run zensical build --clean # production build -> site/ (CI uses `pip install zensical && zensical build --clean`)
uv run zensical build --strict # fail on warnings (broken links/nav)
```

No lint/typecheck/test suite in repo. Only verification is conda dry-run and site build.

Validate conda envs (same as CI):

```bash
conda env create -f envs/qc.yaml --dry-run
conda env create -f envs/varcall.yaml --dry-run
conda env create -f envs/panacota.yaml --dry-run
conda env create -f envs/phylo.yaml --dry-run
```

## Conventions & Gotchas

- **Never create `mkdocs.yml`** — project uses `zensical.toml`. `extra_css = ["stylesheets/extra.css"]` is relative to `docs/`.
- **`site/` is gitignored** (`.gitignore:2`) and is the Pages artifact. Don't commit it; `build --clean` wipes it.
- **Empty package:** `src/ngs_handbook` has no code; `pyproject.toml:13` script `ngs-handbook = "ngs_handbook:main"` is currently unimplemented (empty `__init__.py:1`).
- **Python version:** Enforced `>=3.12` — use `uv` not `pip` locally to respect `uv.lock`.
- **Images/links:** Use relative paths from `docs/` (e.g., `imgs/...`). Zensical `pymdownx` extensions enabled (admonition, superfences, tabbed, emoji via `zensical.extensions.emoji`).
- **Favicon same-origin collision:** `https://ilypopv.github.io` (root Jekyll site) and `/NGS-Handbook/` share origin `ilypopv.github.io`, so browsers cache favicons per-origin. Tabs honor `zensical.toml:42` (`assets/favicon.svg` -> relative `link rel=icon`), but Favorites/bookmarks often probe `https://ilypopv.github.io/favicon.ico` or `images/favicon.ico` (root site) and ignore project path. No fully reliable fix on `github.io` without a custom domain (costs money). Best-effort mitigation in place: `overrides/main.html:4` injects absolute `https://ilypopv.github.io/NGS-Handbook/assets/*` icons + `manifest.webmanifest` with `scope/start_url=/NGS-Handbook/` + copies `docs/favicon.ico` + `docs/assets/favicon-{192,512}.png` + `favicon.ico` (multi-res ICO built from `assets/favicon.svg`). Keep all those in sync; don't reintroduce nested `<html><head>` wrapper that was previously in `overrides/main.html`.
- **Branch:** CI deploys on push to `master` OR `main` (`.github/workflows/docs.yml:5`).

## CI

- `docs.yml`: `actions/setup-python@v5` (python `3.x`) -> `pip install zensical` -> `zensical build --clean` -> `upload-pages-artifact` (`site/`) -> `deploy-pages@v4`. Permissions `pages: write`, `id-token: write`.
- `test-envs.yml`: On any `push`, matrix over 4 yaml files, `conda-incubator/setup-miniconda@v2` + `conda env create -f envs/${{matrix.env}} --dry-run` (no solve/install, just validation).

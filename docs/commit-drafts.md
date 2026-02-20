# Commit Drafts — task/15-export-notebook-figures-to-svg-format

## Commit 1 — uv project setup

**Files:** `pyproject.toml`, `uv.lock`, `.python-version`, `.gitignore`

```plaintext
Initialize uv project with Python 3.12

Add pyproject.toml with all notebook dependencies
pinned to Python 3.12. Use pandas>=2.0,<2.3 for
Colab compatibility, pyjanitor<0.30 to match, and
openpyxl for Excel reads. Exclude .venv/ from
version control in .gitignore.

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>
```

---

## Commit 2 — Notebook adaptation and SVG export (combined)

**Files:** `src/notebooks/site_differentation.ipynb`

```plaintext
Adapt notebook and add dual-format figure export

Wrap google.colab import in try/except to set an
IN_COLAB flag. Guard drive.mount and os.chdir with
that flag so neither runs outside Colab. Remove
!pip install and !sudo apt install cells, now managed
by uv. Fix calibration data paths missing ../../
prefix. Remove deprecated parse_dates=True from
read_excel calls and use datetime64[ns] precision
in the basic_preprocessing dtype map. Add save_fig
helper that routes output to figures/png/ and
figures/svg/. Replace 77 savefig() calls so every
render produces both formats.

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>
```

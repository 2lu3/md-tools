# Plan: Split dit/ and pytool/ from md-tools

Closes #4.

## Goal

Move `dit/` and `pytool/` out of the `2lu3/md-tools` monorepo into
standalone repositories with a fresh snapshot (no git history transfer).

| Source | Destination |
|--------|-------------|
| `md-tools/dit/` | https://github.com/2lu3/dit |
| `md-tools/pytool/` | https://github.com/2lu3/pytool |

## Decisions

- **History**: fresh snapshot only (no `git filter-repo` / subtree split)
- **Base snapshot**: `origin/main` of md-tools at split time
- **CI**: each new repo gets `.github/workflows/increment_version.yml`
  using `astral-sh/setup-uv@v5` and `uv version --bump patch`, committing
  `pyproject.toml` and `uv.lock`
- **md-tools**: becomes a hub README pointing at the new repos; remove
  `dit/`, `pytool/`, and `.github/workflows`

## Steps

1. Bootstrap empty `2lu3/dit` and `2lu3/pytool` from `origin/main` snapshot
2. Add LICENSE, root-level CI, update READMEs; push to `main`
3. Verify `uv sync` and `pip install .` (Python >= 3.10) in each
4. On md-tools branch `4-chore-dit-pytool`: delete packages/workflows,
   update hub README, commit plan, open PR

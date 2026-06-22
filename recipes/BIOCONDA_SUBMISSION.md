# Getting Pan_TE onto bioconda (`conda install -c bioconda pan_te`)

This directory holds the conda recipes that make Pan_TE installable from bioconda. Because
bioconda only accepts **public, versioned, license-bearing** source that builds on its CI,
some steps are upstream actions you must take (creating releases, adding LICENSE files,
opening PRs). Everything in `recipes/` is authored and ready; the checklist below covers
what is left.

## What's here

| Recipe | Builds | Notes |
|--------|--------|-------|
| `look4ltrs/` | `look4ltrs` (C++/CMake) + ships `build_ltr_library.py` | repo already has a LICENSE (AGPL-1.0) |
| `mdl-repeat/` | `mdl-repeat` (C/Make, `PORTABLE=1`) | **needs a LICENSE added upstream** |
| `te-looker/` | `dtr`, `te-discover`, `te-refine`, `te-seed` (Rust/cargo) | **needs a public repo + LICENSE** |
| `pan_te/` | the pipeline scripts + `pan_te-setup-data` | depends on the three above **plus** `fastga` (already on bioconda) |
| `build_local.sh` | local smoke test from local source | verifies builds before any release exists |
| `pan_te/pan_te-setup-data` | post-install fetch of Dfam libs + ClassifyTE models | data assets can't ship in conda |

`fastga` is already on bioconda, so `pan_te` just depends on it — no recipe needed.

## Dependency / submission order

`pan_te` depends on `look4ltrs`, `mdl-repeat`, `te-looker`. Bioconda must have the three
tool packages **merged first**, then `pan_te`:

```
look4ltrs  ┐
mdl-repeat ├─►  pan_te
te-looker  ┘
```

## Pre-submission checklist (the external steps)

1. **Make te-looker a public GitHub repo.** Assumed URL in the recipe:
   `https://github.com/unavailable-2374/te-looker`. If you use a different URL, update
   `te-looker/meta.yaml` (`source.url`, `about.home`) and `pan_te/meta.yaml` description.

2. **Add a `LICENSE` file** to each repo that lacks one and set the recipe's `license:` to
   the matching SPDX id. Status today:
   - Look4LTRs — ✅ has `LICENSE` (AGPL-1.0-only).
   - mdl-repeat — ❌ add one; recipe placeholder is `GPL-3.0-or-later` (a guess).
   - te-looker — ❌ add one; recipe placeholder is `GPL-3.0-or-later` (a guess).
   - Pan_TE — ❌ add one; recipe placeholder is `GPL-3.0-or-later` (a guess).
   > A package whose AGPL dependency (look4ltrs) is statically pulled in may impose
   > license-compatibility expectations — pick the Pan_TE license deliberately.

3. **Cut a tagged release** for each repo so bioconda has a stable tarball:
   ```bash
   git tag v1.0.0 && git push origin v1.0.0      # Pan_TE, Look4LTRs, mdl-repeat
   git tag v0.3.0 && git push origin v0.3.0      # te-looker (matches te-core crate version)
   ```
   Use GitHub's auto-generated `archive/refs/tags/vX.Y.Z.tar.gz`. Pick versions you're happy
   to keep; bioconda bumps require new tags.

4. **Fill in every `sha256`.** For each `url:` in the four `meta.yaml` files:
   ```bash
   curl -sL <the-url-with-version-substituted> | sha256sum
   ```
   Replace the `0000…0000` placeholder. There are: 1 in look4ltrs, 1 in mdl-repeat,
   1 in te-looker, and **2 in pan_te** (Pan_TE tarball + Look4LTRs tarball).

5. **Smoke-test locally first** (catches build.sh / dependency bugs before the PR):
   ```bash
   conda install -n base -c conda-forge conda-build conda-verify   # one-time
   cd recipes
   ./build_local.sh lint     # render all recipes
   ./build_local.sh build    # build all 4 from LOCAL source into recipes/local_channel
   conda create -n pan_te_test -c file://$PWD/local_channel -c conda-forge -c bioconda pan_te
   ```
   `build_local.sh` rewrites each recipe's `source:` to a local `path:` so it works even
   though no release/sha256 exists yet. Override source locations with `PANTE_SRC`,
   `LOOK4LTRS_SRC`, `MDLREPEAT_SRC`, `TELOOKER_SRC` env vars.

6. **Open the bioconda PRs.** Fork `bioconda/bioconda-recipes`, then for each tool copy
   `recipes/<name>/` into the fork's `recipes/<name>/`:
   ```bash
   git clone https://github.com/<you>/bioconda-recipes
   cp -r recipes/look4ltrs  bioconda-recipes/recipes/look4ltrs
   cp -r recipes/mdl-repeat bioconda-recipes/recipes/mdl-repeat
   cp -r recipes/te-looker  bioconda-recipes/recipes/te-looker
   # open these 3 PRs; after they merge, then:
   cp -r recipes/pan_te bioconda-recipes/recipes/pan_te
   ```
   Bioconda's CI (linux-64 + osx-64) and `bioconda-utils lint` run on each PR. Common
   follow-ups: confirm the build is reproducible on osx (add `skip: True  # [osx]` to a
   recipe if a dependency is linux-only), and ensure `license_file` actually exists in the
   tarball.

## After install — the data assets

The Dfam library and ClassifyTE models are **not** in the conda package (redistribution
limits). The package ships `pan_te-setup-data`:

```bash
conda activate pan_te                       # or whatever env you installed into
pan_te-setup-data --dfam-list               # how to pick a Dfam partition
pan_te-setup-data \
    --classifyte-dir ~/Pan_TE_data/ClassifyTE \
    --dfam-url https://www.dfam.org/releases/Dfam_3.9/families/dfam39_full.0.h5.gz
# then:
Pan_TE --genome genome.fa --threads 80 --model-dir ~/Pan_TE_data/ClassifyTE
```

ClassifyTE also needs its own Python-3.8 env (`environment-classifyte.yml`) because its
pinned ML deps conflict with the main env; `run_Classifier` invokes it via
`conda run -n ClassifyTE_env`.

## Known constraints

- **osx-64**: the Rust (`te-looker`) and C/C++ tools should cross-build, but RepeatMasker /
  RepeatModeler behavior differs; if osx CI is troublesome, ship the tool packages on osx
  but mark `pan_te` linux-only.
- **`license_file`** must be present inside each release tarball — verify after tagging.
- The `0000…` sha256 placeholders are intentional and **must** be replaced (step 4); they
  are not valid hashes.

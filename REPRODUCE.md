# Reproducing the CliMA Land / Emerald model runs

The model used in this study is **CliMA Land** (the `Emerald` Julia package,
version 0.3.0). The paper's Julia environment (`src/Project.toml` +
`src/Manifest.toml`) pins every dependency exactly. All dependencies except `Emerald`
are registered Julia packages and resolve automatically from the General registry.

`Emerald` was a git dependency pinned to:

| Field | Value |
|-------|-------|
| Package | `Emerald` (uuid `56ca4449-d59c-41e3-9ae8-0f89b45aadb6`) |
| Version | `0.3.0` |
| `git-tree-sha1` | `d246c63b8f087e1af59cd3f26f3b1ffc6fab4300` |
| Original source (now unavailable) | `https://github.com/Yujie-W/Emerald.git`, rev `wyujie` |

The original development repository's history is periodically rewritten, and the pinned
commit is no longer retrievable from it. The **exact, byte-for-byte source** is therefore
archived here:

- **GitHub:** <https://github.com/braghiere/Land>, tag `shift-dangermond-clima-land`
- **Zenodo:** DOI [10.5281/zenodo.21479295](https://doi.org/10.5281/zenodo.21479295)
  (this version); concept DOI `10.5281/zenodo.21479294` (all versions)

The archived tag's pristine commit reproduces `git-tree-sha1`
`d246c63b8f087e1af59cd3f26f3b1ffc6fab4300` exactly (verifiable with
`git rev-parse 'shift-dangermond-clima-land^{tree}'` / `~1`).

## Steps

```julia
using Pkg

# 1. Get the exact Emerald v0.3.0 source (either clone the tag or download the Zenodo archive):
#    git clone --branch shift-dangermond-clima-land https://github.com/braghiere/Land.git emerald-src
#    (or unzip the Zenodo archive to ./emerald-src)

# 2. Activate the paper environment:
Pkg.activate("src")

# 3. Point Emerald at the archived source (byte-identical to git-tree d246c63b):
Pkg.develop(path="/absolute/path/to/emerald-src")

# 4. Resolve everything else from the pinned Manifest:
Pkg.instantiate()
```

Then run the Julia scripts under `julia_scripts/`. The final paper results use the
`*_prescribed_lai_ci` scripts (e.g. `target_function_free_lai_ci.jl`,
`1_gpp_reg_loop_clima_fit_prescribed_lai_ci.jl`).

> Note: `src/Manifest.toml` is intentionally left recording the original
> (now-unavailable) Emerald source, as a faithful record of provenance. Use the
> archived source above via `Pkg.develop` to reproduce the environment.

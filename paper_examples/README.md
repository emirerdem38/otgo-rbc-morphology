# Paper examples (OTGO cell morphologies)

Minimal MATLAB drivers that reproduce the **main analysis types** from:

> Erdem et al., *How Morphology Affects Optical Trapping of Red Blood Cells in Health and Disease*

These scripts are intentionally coarse (few sample points / short trajectories, `fast_demo = true` → `p = 32`) so they run as demos. Set `fast_demo = false` and increase `N` / `Nsteps` to approach the paper production settings.

A smoke test is included: `smoke_test_paper_examples.m` (verified on MATLAB R2024a).

## Requirements

- MATLAB (R2018b+ recommended)
- This repository’s `otgo1.0.0/` package (includes the extended `Cell` / `ParticleCELL` classes)
- For Brownian dynamics: `otgo1.0.0/diffusion_tensors.mat`

## Quick start

In MATLAB:

```matlab
cd paper_examples
demo_01_plot_morphology
demo_02_static_displacement
demo_03_static_rotation
demo_04_brownian_dynamics
```

Each demo calls `setup_paths` automatically.

## Morphologies

| Paper name   | `shp` string |
|-------------|--------------|
| Dumbbell    | `'dumbbell'` |
| Four-bumped | `'fourBump'` |
| Eight-bumped| `'eightBump'`|
| Neck        | `'neck'`     |
| Oblate      | `'oblate75'` |

Recommended surface density `p`: 64 (dumbbell, fourBump), 80 (eightBump, neck, oblate75).

## Optical parameters (paper)

- Power `5 mW`, `NA = 1.3`, ray set `20 × 40`, scattering cutoff `10`
- `n_medium = 1.33`, `n_RBC = 1.38`

See `paper_optical_params.m`.

## Upstream OTGO

The geometrical-optics engine is based on:

Callegari et al., *Computational toolbox for optical tweezers in geometrical optics* (2015).

Cite that work for the core OTGO framework; this repository adds parametric RBC-like morphologies and paper demos.

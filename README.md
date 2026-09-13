# Optical trapping of deformed red blood cells (OTGO)

MATLAB code accompanying:

> Erdem, Kabacaoğlu, Volpe, Callegari & Biancofiore,  
> *How Morphology Affects Optical Trapping of Red Blood Cells in Health and Disease*  
> (Biomedical Optics Express).

This repository provides an extended [OTGO](https://github.com/softmatterlab) geometrical-optics optical-tweezers toolbox with parametric RBC-like morphologies (dumbbell, four-/eight-bumped, neck, oblate), plus short demos that reproduce the main analysis types in the paper (static force/torque sweeps and Brownian dynamics).

## Contents

| Path | Description |
|------|-------------|
| `otgo1.0.0/` | OTGO 1.0.0 with extended `Cell` / `ParticleCELL` classes and morphology generators |
| `paper_examples/` | Minimal drivers matching the paper optical setup |
| `otgo_examples1.0.0/` | Additional OTGO usage examples |

## Requirements

- MATLAB (tested with R2024a)
- No extra toolboxes beyond what OTGO already uses

## Quick start

```matlab
cd paper_examples
demo_01_plot_morphology
demo_02_static_displacement
demo_03_static_rotation
demo_04_brownian_dynamics
```

See `paper_examples/README.md` for morphology name strings and parameters (`5 mW`, `NA = 1.3`, ray set `20×40`, etc.).

Run `smoke_test_paper_examples.m` for a short end-to-end check.

## Citation

If you use this code, please cite the Biomedical Optics Express article above and the original OTGO toolbox:

> Callegari, A., Mijalkov, M., Gököz, A. B. & Volpe, G.  
> Computational toolbox for optical tweezers in geometrical optics.  
> *J. Opt. Soc. Am. B* **32**, B11–B19 (2015).

## License

Research code released for reproducibility of the paper. Upstream OTGO components remain under their original terms (Callegari et al., 2015).

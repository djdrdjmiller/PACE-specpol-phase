# Discussion starting point: 3D effects on the OCI ratio

For a fixed liquid reference, cancellation in a ratio depends on equal *fractional* perturbations, not equal absolute reflectance changes. Write `rho_b_3D = rho_b_1D (1 + epsilon_b)`. Then

`R_3D / R_1D = (1 + epsilon_2260) / (1 + epsilon_1615)`

and for small perturbations `delta ln R = epsilon_2260 - epsilon_1615`. Equal fractional brightening cancels. A larger fractional loss at 1615 raises the ratio; a larger loss at 2260 lowers it. There is no universal sign for arbitrary cloud geometry. These equations isolate 3D perturbations, holding atmospheric correction and reference fixed; they are not a correction model.

## What the actual 265 K liquid LUT tells us

Inspected `SingleScatterAlbedo` values in the supplied 3 m/s ocean LUT:

| CER (µm) | SSA at 1.616 µm | SSA at 2.26 µm |
|---|---:|---:|
| 6 | 0.99500 | 0.98886 |
| 10 | 0.99186 | 0.98133 |
| 14 | 0.98903 | 0.97462 |
| 20 | 0.98482 | 0.96520 |
| 30 | 0.97802 | 0.95013 |

At these radii, 2260 absorbs more strongly per scattering interaction. That favors shorter surviving photon paths and generally different depth/horizontal sensitivity from 1615. It does not mean that a 3D reflectance bias, a ratio bias, or its magnitude follows directly from these SSA values. The number/directions of scatterings and illumination/view geometry matter too.

Cloud sides, shadows, horizontal transport and within-footprint variability can therefore perturb the bands differently. Low Sun over the Arctic and broken/multilayer clouds are relevant regimes to examine, not proof that a specific high ratio is caused by 3D effects. Existing cloud studies establish 3D bias and wavelength-dependent sampling for related MODIS bands; their numerical biases should not be copied to OCI's 2260/1615 combination. See [Platnick (2000), vertical photon transport](https://atmosphere-imager.gsfc.nasa.gov/sites/default/files/ModAtmo/Platnick%20%282000%29.pdf), [Zhang et al. (2012)](https://doi.org/10.1029/2012JD017655), and [the CAMP2Ex retrieval comparison](https://acp.copernicus.org/articles/22/8259/2022/).

## Normalization adds another dependence

Our liquid reference uses retrieved CER/COT, not known true cloud properties. To first order:

`delta ln Rnorm = delta ln rho2260 - delta ln rho1615 + delta ln T1615 - delta ln T2260 - delta ln Rliquid(CER,COT,geometry)`.

OCI2260 CER/COT is inferred partly from the same numerator band. A 3D perturbation can change the observed ratio **and** the reconstructed reference, with partial cancellation or amplification. Consequently a ratio near unity is not proof that 3D effects are absent. HARP2 CER and OCI2130 comparisons on common samples are useful sensitivity tests, but remain coupled through OCI COT, geometry, and reference assumptions.

## Diagnostics before adopting a threshold

Compare confirmed liquid cases across ratio-reference variants and stratify by solar/view zenith, relative azimuth, COT, surface regime, and independently measured cloud heterogeneity/edge distance. Inspect both individual-band reflectances, not just their ratio. Compare native OCI variability inside the L1C footprint when available; the current cache does not yet contain a footprint-level heterogeneity metric. Keep atmospheric transmission, ice/snow surface contributions, time mismatch and 3D effects as separate hypotheses. The currently global ocean LUT especially limits interpretation over Arctic surfaces when clouds are optically thin.

A useful later quantitative test is an all-liquid 3D scene experiment spanning these geometries, processed through the same 1D OCI retrieval and reference normalization. That would estimate false LTMP detections from scene structure rather than absorbing them into a case-fitted threshold. No such 3D simulation or correction has been implemented in this package.

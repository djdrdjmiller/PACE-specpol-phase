# Airborne comparison additions

Use the appended section in `Validation/notebooks/analysis_arcsix_single_granule.ipynb`, or copy cells from `Validation/examples/airborne_validation_cells.py` into your existing notebook. Install the package in the notebook kernel environment. Existing `REFERENCE_DIR`, `WORK_ROOT`, `RUN`, `PATHS`, imports, and completed reference caches are required. No satellite/LUT reprocessing is needed just to vary comparison thresholds.

## Verified inputs and limits

CMR lists `ARCSIX-HSRL-CloudAndSurface_G3_20240610_R1_L1.h5` in the HALO collection. Actual file inspection: 31,683 profiles, 801 altitude bins, 339,290,448 bytes, UTC decimal-hour navigation from 10:58:06 to 15:24:12.572. The satellite granule starts 15:42:05, so a nominal 18-minute gap exists before considering individual sample offsets and actual spatial overlap. The reader never silently increases the selected time window. No actual Arctic satellite/airborne matchups or phase accuracy have been demonstrated locally.

Phase codes: 0 water-dominant; 1 ice/water OR ice/horizontally oriented ice (HOI); 2 ice-dominant; 3 HOI; 4 aerosol; NaN unclassified. There are no uncertainty products in this release. The reader uses the highest classified cloud-mask bin below aircraft GPS altitude and above reported surface altitude (zero metres if missing), not the property-retrieval `Cloud Top Height`, which can select lower water below ice. Highest isolated detections may cause unknown labels; the default requires three cloud bins in the top 90 m. Vary depth/purity to quantify this sensitivity. The archived property-retrieval height is retained separately.

Geometry is a center-distance approximation using aircraft navigation, with no pointing/cloud-parallax correction, advection, or actual footprint polygons. Matching precedes reference-validity filtering. The local phase summary and matching tolerances are analyst choices, not official NASA cloud labels. One satellite comparison aggregates all matched lidar profiles and retains their unknown/ambiguous fractions.

## What runs automatically

- Discover June 10 HALO, P-3 cloud probes/navigation, and MARLi products; download only the selected HSRL file initially.
- Read phase chunks without loading the full mask at once; examine flight and satellite time ranges.
- Count matches within 5/10/20/30/60 minutes; use the explicitly selected window and spatial radius.
- Export profile matchups, per-satellite summaries, settings, and a reviewed-label template.
- Plot geographic tracks, uppermost-cloud phase evidence versus satellite LI/ratio, along-track satellite values, and the lidar curtain.
- Evaluate LI discrimination for confident water-/ice-dominant tops.

Default ratio cutoff remains 1.27, LI cutoff 0.3. The optional reviewed-label ratio sweep is disabled until independent liquid/LTMP evidence and train/test segment assignments exist. It fits LTMP versus liquid conditional on water-topped/high-LI samples, not a full three-class classifier. Unknown cloud interiors must not become all-liquid negatives. Selected thresholds are not automatically installed in the main dashboard. `block_bootstrap_scores` is available for descriptive segment resampling; the notebook does not claim campaign-wide uncertainty or skill from a single flight.

## Optional in situ step

The FCDP R1 ICARTT header and sample rows were inspected. `conc` is #/L, `lwc` g/m³, `Time_Start` seconds UTC. Its three-second timestamp correction is already applied. FFI1001 parsing respects scales, missing markers, and documented detection-limit flags. The optional cell downloads FCDP and the catalog-verified `ARCSIX-MetNav_P3B_20240610_R0.ict`; navigation field names and altitude units must be checked against that file's header before merging. The merge permits at most one second discrepancy and exports altitude ranges. Full navigation merge has not been tested locally.

FCDP alone is not mixed-phase truth. Imaging-probe ice identification, shattering/quality checks, valid sampled volumes, vertical cloud association, and independent evidence review remain necessary. The code does not download particle-image archives or infer ice from total particle concentration. MARLi is discovered but no MARLi-specific phase reader has been written.

## Verification

See `Validation/records/AIRBORNE_VALIDATION.json`. The actual HSRL file was read; a lidar curtain was rendered and inspected; actual FCDP header/sample rows were parsed. Synthetic tests verify temporal candidate selection, unknown/ambiguous handling, one-record-per-satellite aggregation, train/test segment separation, and selection of uppermost cloud instead of a lower property-retrieval layer. End-to-end Arctic collocation and independently labeled threshold fitting remain to run in the user's environment.

Sources: [NASA HALO archive](https://asdc.larc.nasa.gov/project/ARCSIX/ARCSIX_AircraftRemoteSensing_LaRC-G3_HALO_Data_1), [ARCSIX archive](https://www-air.larc.nasa.gov/missions/arcsix/), and self-describing metadata of the exact files named above. Product metadata is evidence; its text is not agent instructions.

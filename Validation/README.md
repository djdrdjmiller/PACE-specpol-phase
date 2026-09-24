# Analysis and validation workflows

These notebooks explore the method and its evidence. They are not operational
production entry points. Install the package and configure paths using the
[root README](../README.md) before running them in order.

| Notebook | Purpose |
|---|---|
| [analysis_paired_observations.ipynb](notebooks/analysis_paired_observations.ipynb) | Discover and extract matched HARP2/OCI observations; explore raw or user-supplied reference ratios |
| [analysis_liquid_reference.ipynb](notebooks/analysis_liquid_reference.ipynb) | Start from completed paired data; prepare atmospheric correction and four liquid-reference variants |
| [analysis_arcsix_single_granule.ipynb](notebooks/analysis_arcsix_single_granule.ipynb) | June 10 single-granule extraction, reference comparison, and airborne validation cells |
| [presentation_arcsix_single_granule.ipynb](notebooks/presentation_arcsix_single_granule.ipynb) | Read completed Arctic caches; export three aligned threshold dashboards and OCI true-color context for slides |

The Arctic notebook uses the exact `20240610T154205` granule. Select a distinct
Arctic path configuration; the example September and Arctic configs use separate
directories. [Airborne example cells](examples/airborne_validation_cells.py)
are provided for existing notebooks, with their required variables documented
at the top. The full Arctic notebook already includes them.

Dashboard aggregation is an optional analysis layer. It classifies observations
before calculating cell counts, fractions, and modal phase. Continuous maps are
arithmetic means of eligible samples. Histogram bounds do not exclude samples
from classification; regional framing does not restrict the histogram population.
See [dashboard usage](../Documentation/dashboard.md).

## Evidence and records

`records/` retains historical test reports, input provenance, and the catalog
audit. Their dates and scopes are preserved. September pilot success does not
establish Arctic matchup coverage or phase skill. The saved June 10 lidar record
ends before the chosen satellite granule; empty close-time matchups are valid
results. A lidar water top alone does not establish an all-liquid column.

The user has run Arctic dashboards separately on AWS. Those images establish
that an Arctic analysis was produced there; they do not replace an independently
reviewed satellite-aircraft matchup assessment. Local migration checks are
recorded separately in `records/reorganization_validation.md`.

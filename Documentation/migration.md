# Repository reorganization

This refactor preserves the numerical method, threshold defaults, extraction
version 1.1, reference version, index schema 2, and cache fingerprint semantics.
The package release version is separate from these scientific cache versions.
The original flat package is retained outside this repository for rollback.

| Former entry point | Current import |
|---|---|
| `pace_mixed_phase.MatchConfig`, discovery, matching, extraction | `pace_specpol.matching` |
| `pace_mixed_phase.Reference` | `pace_specpol.reference` |
| `pace_mixed_phase.sample_phase` | `pace_specpol.classification` |
| Index configuration, construction, counts, export | `pace_specpol.validation.aggregation` |
| `dashboard`, `reference_dashboard` | `pace_specpol.validation.vis_dashboard` |
| `plot_cached_samples` | `pace_specpol.validation.vis_samples` |
| `pace_liquid_reference` | `pace_specpol.liquid_reference` |
| `pace_oci_companion` | `pace_specpol.oci_companion` |
| Companion/reference cache stages and `VariantReference` | `pace_specpol.workflow` |
| `pace_arcsix_validation` | `pace_specpol.validation.arcsix` |

Use the migrated notebooks after installing the package into their kernel's
environment. The old flat module names are not installed as compatibility
aliases. Existing external notebooks need the import changes above; scientific
datasets and cache formats do not require migration merely because code moved.

Do not move existing research caches just to match the example layout. Point
local configuration at their existing directories. Reopening an existing index
with `PhaseIndex.load` does not need LUT access. Rebuilding an index needs local
prepared samples but does not require a new satellite extraction.

Reference fingerprints include configuration path strings and file hashes.
Changing a LUT path spelling (even to an equivalent absolute path), selecting a
different source manifest, or changing scientific settings may intentionally
reject resume. Keep the original settings for continued writes to an old cache;
otherwise choose a new derived directory. No cache manifests are rewritten by
this refactor. In particular, do not resume extractor 1.1 into a 1.0 extraction
directory; the documented legacy observation-time repair remains in place.

The `history/` folder preserves the old README as historical context. Current
installation and usage instructions are in the root README. The cell-only
airborne notebook was retired in favor of the full Arctic notebook and one
example script; reusable plotting is now in `validation/vis_airborne.py`.

License choice, GitHub repository creation, and pushing changes are separate
from this local reorganization.


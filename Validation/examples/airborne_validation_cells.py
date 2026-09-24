# %% [markdown]
# # ARCSIX airborne validation intercomparison
# Append these cells after the reference-cache stage in `analysis_arcsix_single_granule.ipynb`. They require its `REFERENCE_DIR`, `WORK_ROOT`, and `RUN` variables and the accompanying installed `pace_specpol` package.
#
# **Verified June 10 product:** `ARCSIX-HSRL-CloudAndSurface_G3_20240610_R1_L1.h5`, in the HALO collection (approximately 324 MiB). Actual navigation spans **10:58:06–15:24:12.572 UTC**. This ends about 18 minutes before the satellite granule; close time-window matches may correctly be empty. No advection or cloud parallax correction is performed. Distances use aircraft GPS and satellite centers, not exact footprint polygons.
#
# The file's phase codes are 0=water-dominant, 1=ice/water OR ice/horizontally oriented ice (HOI), 2=ice-dominant, 3=HOI, 4=aerosol. Code 1 is **not confirmed mixed phase**. The top-window summary below is our configurable diagnostic, not a new official NASA phase product. Liquid tops with unobserved cloud interiors remain unknown for LTMP truth. The archived readme encourages consultation with the lidar team; no uncertainty fields are supplied in this release.
#
# [HALO collection](https://asdc.larc.nasa.gov/project/ARCSIX/ARCSIX_AircraftRemoteSensing_LaRC-G3_HALO_Data_1) · [ARCSIX archive](https://www-air.larc.nasa.gov/missions/arcsix/)

# %%
import json
import numpy as np
import earthaccess
from IPython.display import display
import pandas as pd
import matplotlib.pyplot as plt
from pace_specpol.validation.vis_airborne import (
    plot_matchup_context,
    plot_lidar_curtain,
)
from pace_specpol.validation import arcsix as av

VALIDATION_DIR = WORK_ROOT / RUN / "airborne_validation"
VALIDATION_DIR.mkdir(parents=True, exist_ok=True)
AIRBORNE_DIR = PATHS["airborne_root"]
DAY = "2024-06-10"
VARIANT = "oci_2260"  # also oci_2130, harp2_2260, harp2_2130
POPULATION = "common"  # same population for comparisons across variants
MAX_DISTANCE_KM = 3.0  # center-distance proxy; test 1, 3, 5 km
TIME_WINDOWS_MIN = [5, 10, 20, 30, 60]
SELECTED_WINDOW_MIN = 20  # exploratory; never silently broadened
TOP_DEPTH_M = 90
LI_CUTOFF = 0.3
RATIO_CUTOFF = 1.27

earthaccess.login()
air_catalog = av.discover(DAY)
for kind, entries in air_catalog.items():
    print(kind, len(entries))
    for g in entries:
        if kind in ("lidar", "navigation"):
            print(" ", av.name(g))
# Save catalog metadata, not temporary signed download URLs.
(VALIDATION_DIR / "airborne_catalog.json").write_text(
    json.dumps(air_catalog, default=dict, indent=2)
)


# %% [markdown]
# ## Load the cloud-phase profiles and satellite samples
# Only the cloud/surface HDF5 file is downloaded here—not the large particle-image archives. Reading the phase mask uses small chunks. Each top-window label requires at least three classified cloud bins and 75% agreement, with no ambiguous/HOI bins in that window. Change the depth and purity to test label sensitivity. The “ice below” diagnostic is evidence only: it can include another layer and does not establish absence of ice when zero.
#
# The phase window starts at the **highest classified cloud-mask bin below the aircraft and above the surface**. The separate “Cloud Top Height” property-retrieval field can refer to a lower water layer beneath ice, so it is retained as `retrieval_cloud_top_m` and is not used as the uppermost phase boundary.

# %%
hsrl_path = av.download_exact(
    air_catalog["lidar"],
    "ARCSIX-HSRL-CloudAndSurface_G3_20240610_R1_L1.h5",
    AIRBORNE_DIR,
)
air = av.read_hsrl(hsrl_path, DAY, top_depth_m=TOP_DEPTH_M, min_bins=3, purity=0.75)
sat = av.read_satellite(REFERENCE_DIR, VARIANT, POPULATION)
print("Aircraft record:", air.time.min(), "to", air.time.max())
print("OCI sample times:", sat.time.min(), "to", sat.time.max())
print("HARP nadir times:", sat.harp_time.min(), "to", sat.harp_time.max())
print(
    "Aircraft end to earliest OCI sample (minutes):",
    (sat.time.min() - air.time.max()).total_seconds() / 60,
)
display(air.lidar_top.value_counts())


# %% [markdown]
# ## Quantify overlap before plotting or tuning
# Matching uses OCI observation times. Each aircraft profile is assigned to the closest spatially eligible satellite sample within the chosen time window; no ratio or phase value participates in matching. Geometry matches are retained even when the LUT reference is invalid. Many airborne profiles in one satellite sample are collapsed to one comparison record, with all label fractions retained. This avoids counting hundreds of profiles as hundreds of independent satellite pixels, although adjacent satellite samples are still correlated.

# %%
overlap = []
selected_matches = None
for window in TIME_WINDOWS_MIN:
    matches = av.match_profiles(sat, air, max_minutes=window, max_km=MAX_DISTANCE_KM)
    summary = av.aggregate_matches(matches, min_profiles=3, purity=0.8)
    overlap.append(
        dict(
            window_min=window,
            profiles=len(matches),
            satellite_samples=len(summary),
            confident_top_samples=0
            if summary.empty
            else int(summary.lidar_top.isin(["water_dominant", "ice_dominant"]).sum()),
        )
    )
    if window == SELECTED_WINDOW_MIN:
        selected_matches = matches
if selected_matches is None:
    selected_matches = av.match_profiles(sat, air, SELECTED_WINDOW_MIN, MAX_DISTANCE_KM)
footprints = av.aggregate_matches(selected_matches, min_profiles=3, purity=0.8)
display(pd.DataFrame(overlap))
if footprints.empty:
    print(
        "No matchups at the selected tolerances. Do not infer phase skill from this case."
    )
else:
    display(
        footprints[
            [
                "sat_id",
                "n_profiles",
                "dt_median_min",
                "distance_max_km",
                "lidar_top",
                "valid",
            ]
        ].head()
    )
tag = f"{VARIANT}_{POPULATION}_{SELECTED_WINDOW_MIN}min_{MAX_DISTANCE_KM:g}km"
selected_matches.to_csv(VALIDATION_DIR / f"{tag}_profile_matches.csv", index=False)
footprints.to_csv(VALIDATION_DIR / f"{tag}_satellite_comparisons.csv", index=False)
(VALIDATION_DIR / f"{tag}_settings.json").write_text(
    json.dumps(
        dict(
            day=DAY,
            variant=VARIANT,
            population=POPULATION,
            max_distance_km=MAX_DISTANCE_KM,
            time_window_min=SELECTED_WINDOW_MIN,
            top_depth_m=TOP_DEPTH_M,
            hsrl_file=hsrl_path.name,
            reference_dir=str(REFERENCE_DIR),
            geometry="center distance; aircraft GPS; no advection or cloud-parallax correction",
        ),
        indent=2,
    )
)


# %% [markdown]
# ## Track map and matched phase distributions
# The map distinguishes the full lidar flight, the selected satellite granule, and accepted matches. The joint plot colors observations by lidar **cloud-top** evidence; it is not an LTMP truth plot. The two lower traces use satellite sample time; aircraft time offsets are in the exported table.

# %%
plot_matchup_context(
    sat,
    air,
    selected_matches,
    footprints,
    SELECTED_WINDOW_MIN,
    MAX_DISTANCE_KM,
    VARIANT,
    RATIO_CUTOFF,
    LI_CUTOFF,
    VALIDATION_DIR,
    tag,
)


# %%
plot_lidar_curtain(air, sat, hsrl_path, VALIDATION_DIR, tag)


# %% [markdown]
# ## LI sensitivity against confident lidar cloud tops
# This tests water-dominant versus ice-dominant tops only. It does not require a valid normalized ratio, and does not validate pure liquid throughout a cloud. No automatic “best” LI cutoff is adopted. A single flight is exploratory; uncertainty should eventually be estimated over cloud/flight segments, with separate flight days for validation.

# %%
li_scores = av.liquid_top_sweep(footprints, np.arange(-0.2, 1.01, 0.02))
if not li_scores.empty:
    display(li_scores.iloc[::10])
    li_scores.plot(
        x="threshold",
        y=["recall", "false_positive_rate", "balanced_accuracy"],
        ylim=(0, 1),
    )
    plt.title("Exploratory liquid-top test; missing classes give undefined metrics")
    plt.show()
    li_scores.to_csv(VALIDATION_DIR / f"{tag}_LI_scores.csv", index=False)


# %% [markdown]
# ## Optional P-3 in situ support
# FCDP R1 reports concentration (#/L), LWC (g/m³), and size distributions. Its header states a **3-second time correction is already applied**; do not apply it again. FCDP alone cannot establish absence of ice or provide a complete phase label. Use quality-controlled imaging-probe evidence and aircraft altitude to interpret LTMP.
#
# The optional cell downloads one FCDP table and inspects navigation choices. Set the navigation filename and field names from its printed header, including units; no guessed mapping is used. Other ICARTT FFI-1001 cloud-probe files can use the same reader. Particle-image ZIPs are deliberately not downloaded.

# %%
RUN_INSITU = False
if RUN_INSITU:
    probe_path = av.download_exact(
        air_catalog["probes"], "ARCSIX-FCDP_P3B_20240610_R1.ict", AIRBORNE_DIR
    )
    probe = av.read_icartt(probe_path)
    print(probe.attrs["header"])
    probe["time"] = pd.Timestamp(DAY) + pd.to_timedelta(probe["Time_Start"], unit="s")
    display(probe[["time", "conc", "lwc"]].head())
    print(
        "Navigation choices:",
        [av.name(g) for g in air_catalog["navigation"] if av.name(g).endswith(".ict")],
    )


# %%
# Edit only after inspecting the chosen navigation header.
NAV_FILENAME = "ARCSIX-MetNav_P3B_20240610_R0.ict"
NAV_FIELDS = dict(seconds="EDIT", latitude="EDIT", longitude="EDIT", altitude_m="EDIT")
if RUN_INSITU and NAV_FILENAME:
    nav_path = av.download_exact(air_catalog["navigation"], NAV_FILENAME, AIRBORNE_DIR)
    nav_raw = av.read_icartt(nav_path)
    print(nav_raw.attrs["header"])
    if any(v == "EDIT" for v in NAV_FIELDS.values()):
        print(
            "Set NAV_FIELDS from the header and convert altitude to metres before matching."
        )
    else:
        nav = nav_raw[[*NAV_FIELDS.values()]].rename(
            columns={v: k for k, v in NAV_FIELDS.items()}
        )
        nav["time"] = pd.Timestamp(DAY) + pd.to_timedelta(nav.pop("seconds"), unit="s")
        # Never bridge navigation gaps larger than one second.
        insitu = pd.merge_asof(
            probe.dropna(subset=["time"]).sort_values("time"),
            nav.dropna(subset=["time"]).sort_values("time"),
            on="time",
            direction="nearest",
            tolerance=pd.Timedelta(seconds=1),
        )
        insitu_matches = av.match_profiles(
            sat, insitu, SELECTED_WINDOW_MIN, MAX_DISTANCE_KM
        )
        insitu_matches.to_csv(VALIDATION_DIR / f"{tag}_FCDP_matches.csv", index=False)
        print("Matched in situ records:", len(insitu_matches))
        if not insitu_matches.empty:
            # One summary per satellite sample; retain vertical range rather than hide it.
            insitu_summary = insitu_matches.groupby("sat_sat_id").agg(
                lwc_median=("air_lwc", "median"),
                altitude_min_m=("air_altitude_m", "min"),
                altitude_max_m=("air_altitude_m", "max"),
                n_records=("air_index", "size"),
            )
            display(insitu_summary)
            insitu_summary.to_csv(VALIDATION_DIR / f"{tag}_FCDP_summary.csv")


# %% [markdown]
# ## Ratio-threshold tuning only after independent phase review
# Create labels from lidar visibility, vertical structure and/or quality-controlled in situ ice/liquid evidence—not from LI or the OCI ratio itself. A liquid top alone is insufficient for a `liquid` (all-liquid) label. Preserve `unknown` for attenuated profiles or doubtful collocations. Ice far beneath the optically sensed cloud portion should be identified in the evidence notes, not assumed to affect the SWIR ratio.
#
# The optional analysis below tests **LTMP versus liquid among independently water-topped, high-LI observations**. It does not optimize the whole three-class map. Assign complete cloud/flight segments to train or test; do not randomly divide adjacent profiles. Training selects a threshold; held-out data evaluate that fixed threshold once. One flight is still only a pilot, even with a segment split. No confidence interval is claimed here.

# %%
labels_path = VALIDATION_DIR / f"{tag}_reviewed_labels.csv"
if not footprints.empty and not labels_path.exists():
    labels = footprints[["sat_id"]].copy()
    labels["reference_phase"] = "unknown"  # liquid, LTMP, ice, unknown
    labels["evidence"] = (
        ""  # independent source, visibility, altitude, QA, collocation justification
    )
    labels["split"] = ""  # train or test; assign by whole cloud/flight segment
    labels["segment_id"] = ""
    labels.to_csv(labels_path, index=False)
    print("Review template written:", labels_path)

RUN_REVIEWED_TUNING = False
if RUN_REVIEWED_TUNING and not footprints.empty:
    labels = pd.read_csv(labels_path)
    scores, reviewed = av.reviewed_ratio_sweep(
        footprints, labels, np.arange(0.7, 1.801, 0.01), LI_CUTOFF
    )
    if scores.empty:
        print(
            "No independently reviewed eligible liquid/LTMP samples; threshold unchanged."
        )
    else:
        best = scores.loc[scores.balanced_accuracy.idxmax()]
        print("Exploratory TRAIN optimum:", best.to_dict())
        scores.plot(
            x="threshold",
            y=["recall", "precision", "false_positive_rate", "balanced_accuracy"],
            ylim=(0, 1),
        )
        plt.show()
        held = reviewed[reviewed["split"] == "test"]
        if held.reference_phase.nunique() == 2:
            print(
                "HELD-OUT scores at frozen threshold:",
                av.binary_scores(
                    held.reference_phase == "LTMP", held.ratio >= best.threshold
                ),
            )
        else:
            print(
                "No two-class held-out test: this is fitting, not independent validation."
            )
        scores.to_csv(VALIDATION_DIR / f"{tag}_ratio_training_scores.csv", index=False)

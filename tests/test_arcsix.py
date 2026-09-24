import numpy as np
import pandas as pd
import pytest
from pace_specpol.validation.arcsix import (
    match_profiles,
    aggregate_matches,
    reviewed_ratio_sweep,
    binary_scores,
    read_icartt,
)


def test_time_eligible_neighbor_and_no_match():
    s = pd.DataFrame(
        dict(
            latitude=[70, 70.01],
            longitude=[-40, -40],
            time=pd.to_datetime(["2024-06-10T14:00", "2024-06-10T15:00"]),
            sat_id=["a", "b"],
        )
    )
    a = pd.DataFrame(
        dict(
            latitude=[70],
            longitude=[-40],
            time=pd.to_datetime(["2024-06-10T15:01"]),
            lidar_top=["water_dominant"],
        )
    )
    m = match_profiles(s, a, 5, 3)
    assert len(m) == 1 and m.sat_sat_id.iloc[0] == "b"
    assert m.air_minus_sat_seconds.iloc[0] == 60
    assert match_profiles(s, a, 0.1, 3).empty


def test_ambiguous_profiles_do_not_become_ltmp():
    s = pd.DataFrame(
        dict(
            latitude=[70],
            longitude=[-40],
            time=pd.to_datetime(["2024-06-10T15:00"]),
            sat_id=["a"],
        )
    )
    a = pd.DataFrame(
        dict(
            latitude=[70] * 5,
            longitude=[-40] * 5,
            time=pd.to_datetime(["2024-06-10T15:00"] * 5),
            lidar_top=["water_dominant"] * 3 + ["ambiguous_or_HOI"] * 2,
        )
    )
    q = aggregate_matches(match_profiles(s, a))
    assert len(q) == 1 and q.lidar_top.iloc[0] == "unknown"
    assert q.fraction_ambiguous_or_HOI.iloc[0] == 0.4


def test_review_needs_evidence_and_segment_independence():
    f = pd.DataFrame(
        dict(
            sat_id=["a", "b"],
            lidar_top=["water_dominant"] * 2,
            valid=[True] * 2,
            liquid_index=[1, 1],
            ratio=[1, 1.5],
        )
    )
    l = pd.DataFrame(
        dict(
            sat_id=["a", "b"],
            reference_phase=["liquid", "LTMP"],
            evidence=["", "probe"],
            split=["train", "train"],
            segment_id=["a", "b"],
        )
    )
    with pytest.raises(ValueError, match="evidence"):
        reviewed_ratio_sweep(f, l, [1.27])
    l.loc[0, "evidence"] = "visible column and probe"
    scores, q = reviewed_ratio_sweep(f, l, [1.27])
    assert scores.balanced_accuracy.iloc[0] == 1
    l["segment_id"] = "same"
    l.loc[1, "split"] = "test"
    with pytest.raises(ValueError, match="segment"):
        reviewed_ratio_sweep(f, l, [1.27])


def test_no_negative_class_is_undefined():
    assert np.isnan(binary_scores([True], [True])["balanced_accuracy"])


def test_icartt_scale_missing_and_detection_flags(tmp_path):
    header = [
        "18,1001",
        "pi",
        "org",
        "instrument",
        "mission",
        "1,1",
        "2024,6,10,2024,6,11",
        "1",
        "time,seconds",
        "1",
        "2",
        "-999",
        "value,unit",
        "0",
        "2",
        "ULOD_FLAG: -7777",
        "LLOD_FLAG: -8888",
        "time,value",
    ]
    p = tmp_path / "test.ict"
    p.write_text("\n".join(header + ["0,3", "1,-999", "2,-7777", "3,-8888"]))
    d = read_icartt(p)
    assert d.value.iloc[0] == 6 and d.value.iloc[1:].isna().all()


def test_uppermost_cloud_not_property_retrieval_height(tmp_path):
    import h5py
    from pace_specpol.validation.arcsix import read_hsrl

    p = tmp_path / "lidar.h5"
    with h5py.File(p, "w") as f:
        f["000_Readme"] = np.array([["test"]], dtype=h5py.string_dtype())
        f["DataProducts/Cloud/Cloud Top Height"] = [[30.0]]  # low water retrieval
        f["DataProducts/CloudMask/Altitude"] = [
            [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0]
        ]
        a = f.create_dataset(
            "DataProducts/CloudMask/Type",
            data=[[0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2.0, 2.0, 2.0]],
        )
        a.attrs["Long_Name"] = "0 = Water Dominant; 3 = HOI"
        for k, v in [
            ("gps_alt", 1000.0),
            ("gps_lat", 70.0),
            ("gps_lon", -40.0),
            ("gps_time", 15.0),
        ]:
            f["Nav/" + k] = [[v]]
        f["DataProducts/Surface/altSurface"] = [[0.0]]
    d = read_hsrl(p, top_depth_m=30)
    assert d.cloud_top_m.iloc[0] == 120 and d.retrieval_cloud_top_m.iloc[0] == 30
    assert d.lidar_top.iloc[0] == "ice_dominant"

"""The Arctic entry point must reject missing tables before network activity."""

import json
from pathlib import Path

import earthaccess
import pytest

from pace_specpol import matching, paths


@pytest.mark.parametrize("inputs_present", [False, True])
def test_arctic_inputs_checked_before_login(tmp_path, monkeypatch, inputs_present):
    root = Path(__file__).resolve().parents[1]
    notebook = json.loads(
        (root / "Validation/notebooks/analysis_arcsix_single_granule.ipynb").read_text()
    )
    configured = {
        key: tmp_path / key
        for key in (
            "lut_root",
            "pair_cache",
            "work_root",
            "airborne_root",
            "export_root",
        )
    }
    monkeypatch.setattr(paths, "load_paths", lambda *args, **kwargs: configured)
    if inputs_present:
        for relative in (
            "LIQUID/ocean_msr_water_wspeed_3_v6.PACE.1.1.5.2026144071240.hdf",
            "IceAndWaterPhaseFunctionData_v6.PACE.1.1.5.2026142144440.hdf",
            "Transmittance_OCI.hdf",
        ):
            file = configured["lut_root"] / relative
            file.parent.mkdir(parents=True, exist_ok=True)
            file.touch()  # Presence only: this check must not read the large LUTs.

    class LoginReached(Exception):
        pass

    calls = []

    def login():
        calls.append("login")
        raise LoginReached

    def unexpected_network(*args, **kwargs):
        pytest.fail("Discovery/extraction ran before the login sentinel")

    monkeypatch.setattr(earthaccess, "login", login)
    monkeypatch.setattr(matching, "discover_pairs", unexpected_network)
    monkeypatch.setattr(matching, "cache_pairs", unexpected_network)
    namespace = {}
    expected = LoginReached if inputs_present else FileNotFoundError
    with pytest.raises(expected) as error:
        for cell in notebook["cells"]:
            if cell["cell_type"] == "code":
                source = "".join(cell["source"])
                if source.lstrip().startswith("# Run only if") and "%pip install" in source:
                    continue  # Environment installation is separate from the input preflight.
                exec(compile(source, cell["id"], "exec"), namespace)
    assert calls == (["login"] if inputs_present else [])
    assert not configured["pair_cache"].exists()
    assert not configured["work_root"].exists()
    if not inputs_present:
        assert "update lut_root" in str(error.value)
        assert "No satellite downloads have started" in str(error.value)

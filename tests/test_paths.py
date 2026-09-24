import pytest

from pace_specpol.paths import load_paths


def test_paths_follow_config_not_working_directory(tmp_path, monkeypatch):
    config = tmp_path / "config"
    config.mkdir()
    file = config / "paths.local.toml"
    keys = ["lut_root", "pair_cache", "work_root", "airborne_root", "export_root"]
    file.write_text("[paths]\n" + "\n".join(f'{key} = "../{key}"' for key in keys))
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    monkeypatch.setenv("PACE_CONFIG", str(file))
    assert load_paths() == {key: tmp_path / key for key in keys}
    assert not (tmp_path / "pair_cache").exists()


def test_missing_config_is_actionable(tmp_path, monkeypatch):
    monkeypatch.delenv("PACE_CONFIG", raising=False)
    monkeypatch.chdir(tmp_path)
    with pytest.raises(FileNotFoundError, match="PACE_CONFIG"):
        load_paths()

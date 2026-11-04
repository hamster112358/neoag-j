"""The EasyFuse install-check rows must follow the layout that is actually installed."""

from pathlib import Path

from neoag.open_neo.install_check import _production_run_readiness

ROOT = Path(__file__).resolve().parents[1]

CHECK = "easyfuse_requant_star_index_parameter_cleanup"

STAR_CUSTOM_UNPATCHED = '''process STAR_CUSTOM {
    script:
    """
    STAR --genomeDir ${star_index} --readFilesIn ${reads}
    """
}
'''

STAR_CUSTOM_PATCHED = '''process STAR_CUSTOM {
    script:
    """
    if [ -f ${star_index}/genomeParameters.txt ]; then
        sed -i.bak -e '/^genomeType[[:space:]]/d' -e '/^genomeTransform/d' ${star_index}/genomeParameters.txt
    fi
    STAR --genomeDir ${star_index} --readFilesIn ${reads}
    """
}
'''


def _cleanup_row(tmp_path, monkeypatch, module_body):
    easyfuse_root = tmp_path / "EasyFuse"
    (easyfuse_root / "modules").mkdir(parents=True, exist_ok=True)
    (easyfuse_root / "main.nf").write_text("// EasyFuse\n", encoding="utf-8")
    if module_body is not None:
        module = easyfuse_root / "modules" / "05_requantification.nf"
        module.write_text(module_body, encoding="utf-8")
    monkeypatch.setenv("NEOAG_EASYFUSE_HOME", str(easyfuse_root))
    _, rows = _production_run_readiness({}, ROOT, "full")
    return next(row for row in rows if row["check"] == CHECK)


def test_bp_quant_layout_without_star_custom_is_not_a_blocker(tmp_path, monkeypatch):
    """EasyFuse 2.2 quantifies with bp_quant and has no STAR_CUSTOM process."""
    row = _cleanup_row(tmp_path, monkeypatch, None)
    assert row["status"] == "OK"
    assert "star_custom_present=False" in row["evidence"]


def test_star_custom_module_without_genome_parameters_cleanup_stays_blocked(tmp_path, monkeypatch):
    row = _cleanup_row(tmp_path, monkeypatch, STAR_CUSTOM_UNPATCHED)
    assert row["status"] == "BLOCKED"
    assert "star_custom_present=True" in row["evidence"]


def test_star_custom_module_with_cleanup_is_ok(tmp_path, monkeypatch):
    row = _cleanup_row(tmp_path, monkeypatch, STAR_CUSTOM_PATCHED)
    assert row["status"] == "OK"
    assert "star_custom_present=True" in row["evidence"]

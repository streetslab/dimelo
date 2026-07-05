from __future__ import annotations

import pandas as pd
import pytest

from dimelo import integration


def test_build_integrated_site_table_joins_and_prefixes_collisions():
    occupancy = pd.DataFrame(
        {"region_id": ["s1", "s2"], "n_reads": [20, 30], "occupancy_rate": [0.8, 0.2]}
    )
    binding_strength = pd.DataFrame(
        {"region_id": ["s1", "s2"], "n_reads": [20, 30], "occupancy_excess": [0.7, 0.1]}
    )
    trans = pd.DataFrame({"region_id": ["s1"], "trans_fraction": [0.05]})

    table = integration.build_integrated_site_table(
        occupancy=occupancy, binding_strength=binding_strength, trans=trans
    ).set_index("region_id")

    # colliding column (n_reads in two sources) is prefixed; unique columns are not
    assert "occupancy_n_reads" in table.columns
    assert "binding_strength_n_reads" in table.columns
    assert "occupancy_rate" in table.columns  # unique -> unprefixed
    assert "occupancy_excess" in table.columns
    assert table.loc["s1", "occupancy_n_reads"] == 20
    assert table.loc["s1", "trans_fraction"] == pytest.approx(0.05)
    # outer join keeps s2 even though trans only had s1
    assert pd.isna(table.loc["s2", "trans_fraction"])


def test_build_integrated_site_table_requires_site_column_and_a_source():
    with pytest.raises(ValueError, match="at least one non-None source"):
        integration.build_integrated_site_table()
    with pytest.raises(ValueError, match="requires the site column"):
        integration.build_integrated_site_table(
            occupancy=pd.DataFrame({"foo": [1]})
        )


def _read_states():
    # per-read annotations already joined: bound/unbound x footprint
    rows = []
    # site s1: 6 bound+footprint, 2 bound+nofootprint, 2 unbound
    for i in range(6):
        rows.append(("s1", f"a{i}", True, True))
    for i in range(2):
        rows.append(("s1", f"b{i}", True, False))
    for i in range(2):
        rows.append(("s1", f"c{i}", False, False))
    # site s2: all unbound
    for i in range(5):
        rows.append(("s2", f"d{i}", False, False))
    return pd.DataFrame(
        rows, columns=["region_id", "read_id", "is_true_signal", "footprint_present"]
    )


def test_single_molecule_state_composition():
    composition = integration.single_molecule_state_composition(
        _read_states(), state_columns=["is_true_signal", "footprint_present"]
    )
    s1 = composition[composition["region_id"] == "s1"].set_index("combined_state")
    # 6/10 bound+footprint, 2/10 bound-only, 2/10 unbound
    assert s1.loc["True|True", "n_reads"] == 6
    assert s1.loc["True|True", "fraction"] == pytest.approx(0.6)
    assert s1.loc["True|False", "fraction"] == pytest.approx(0.2)
    assert s1.loc["False|False", "fraction"] == pytest.approx(0.2)
    # fractions sum to 1 per site
    assert composition.groupby("region_id")["fraction"].sum().round(6).eq(1.0).all()
    # s2 is entirely one state
    s2 = composition[composition["region_id"] == "s2"]
    assert len(s2) == 1
    assert s2.iloc[0]["fraction"] == pytest.approx(1.0)


def test_single_molecule_state_composition_validates():
    with pytest.raises(ValueError, match="state_columns must be non-empty"):
        integration.single_molecule_state_composition(
            _read_states(), state_columns=[]
        )
    with pytest.raises(ValueError, match="requires columns"):
        integration.single_molecule_state_composition(
            _read_states(), state_columns=["not_a_column"]
        )


def test_state_composition_plotting():
    import matplotlib

    matplotlib.use("Agg")
    from dimelo import plotting, plotting_matplotlib

    composition = integration.single_molecule_state_composition(
        _read_states(), state_columns=["is_true_signal", "footprint_present"]
    )
    payload = plotting.prepare_state_composition_data(composition=composition)
    # wide table: one row per site, one column per state
    assert set(payload["metadata"]["states"]) == {"True|True", "True|False", "False|False"}
    fig, ax = plotting_matplotlib.plot_state_composition_matplotlib(
        payload, title="composition"
    )
    assert fig is not None

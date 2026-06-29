from syncomdesign.community import SimpleModel, SimpleReaction
from syncomdesign.fluxes import extract_fluxes, find_metabolite_exchanges


def _model():
    return SimpleModel(
        reactions=[
            SimpleReaction("R_EX_no_u", -1000, 1000, {"no[u]": -1}),
            SimpleReaction("R_EX_no2_u", -1000, 1000, {"no2[u]": -1}),
            SimpleReaction("R_EX_no3_u", -1000, 1000, {"no3[u]": -1}),
            SimpleReaction("R_EX_n2o_u", -1000, 1000, {"n2o[u]": -1}),
            SimpleReaction("R_EX_n2_u", -1000, 1000, {"n2[u]": -1}),
        ]
    )


def test_no_no2_no3_n2o_mapping_not_confused():
    aliases = ["no", "no_e", "EX_no_e", "nitric_oxide"]
    assert find_metabolite_exchanges(_model(), aliases) == ["R_EX_no_u"]


def test_flux_mapping_values_are_directional():
    values, mapping, flux_rows = extract_fluxes(
        _model(),
        {
            "R_EX_no_u": 5,
            "R_EX_no2_u": -2,
            "R_EX_no3_u": -10,
            "R_EX_n2o_u": -3,
            "R_EX_n2_u": 7,
        },
    )
    assert values["no"]["secretion"] == 5
    assert values["no2"]["uptake"] == 2
    assert values["no3"]["uptake"] == 10
    assert values["n2o"]["uptake"] == 3
    assert values["n2"]["secretion"] == 7
    assert {row["canonical_id"] for row in mapping} == {"no", "no2", "no3", "n2o", "n2"}
    assert len(flux_rows) == 5


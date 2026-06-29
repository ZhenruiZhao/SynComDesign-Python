from syncomdesign.community import SimpleModel, SimpleReaction
from syncomdesign.medium import MediumEntry, apply_community_medium


def _model():
    return SimpleModel(
        reactions=[
            SimpleReaction("R_EX_no3_u", -1000, 1000, {"no3[u]": -1}),
            SimpleReaction("R_EX_n2o_u", -1000, 1000, {"n2o[u]": -1}),
            SimpleReaction("x005__EX_no3_e", -1000, 1000, {"x005__no3_e": -1, "no3[u]": 1}),
            SimpleReaction("x005__NO3tex", -1000, 1000, {"x005__no3_e": -1, "x005__no3_p": 1}),
            SimpleReaction("x005__PGI", -1000, 1000, {"a": -1, "b": 1}),
        ],
        syncomdesign={"externalExchangeMap": {"no3[u]": "R_EX_no3_u", "n2o[u]": "R_EX_n2o_u"}},
    )


def _classes():
    return {
        "R_EX_no3_u": "external_medium_exchange",
        "R_EX_n2o_u": "external_medium_exchange",
        "x005__EX_no3_e": "strain_shared_interface",
        "x005__NO3tex": "internal_transport",
        "x005__PGI": "metabolic_reaction",
    }


def test_medium_only_changes_external_shared_exchange():
    model = _model()
    out = apply_community_medium(model, [MediumEntry("no3", "R_EX_no3_e", -10, 1000)], _classes())
    assert model.reactions[0].lower_bound == -10
    assert model.reactions[2].lower_bound == -1000
    assert model.reactions[3].lower_bound == -1000
    assert out["medium_mapping_warnings"] == []


def test_unlisted_external_uptake_closed():
    model = _model()
    apply_community_medium(model, [MediumEntry("no3", "R_EX_no3_e", -10, 1000)], _classes())
    assert model.reactions[1].id == "R_EX_n2o_u"
    assert model.reactions[1].lower_bound == 0


def test_interface_not_changed_by_medium():
    model = _model()
    apply_community_medium(model, [MediumEntry("no3", "R_EX_no3_e", -10, 1000)], _classes())
    assert model.reactions[2].id == "x005__EX_no3_e"
    assert model.reactions[2].lower_bound == -1000
    assert model.reactions[2].upper_bound == 1000


def test_internal_transport_not_changed_by_medium():
    model = _model()
    apply_community_medium(model, [MediumEntry("no3", "R_EX_no3_e", -10, 1000)], _classes())
    assert model.reactions[3].id == "x005__NO3tex"
    assert model.reactions[3].lower_bound == -1000
    assert model.reactions[3].upper_bound == 1000


def test_cross_feeding_structurally_allowed():
    model = _model()
    apply_community_medium(model, [], _classes())
    assert model.reactions[0].lower_bound == 0
    assert model.reactions[2].lower_bound == -1000
    assert "no3[u]" in model.reactions[2].metabolites


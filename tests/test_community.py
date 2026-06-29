from syncomdesign.community import SimpleModel, SimpleReaction, classify_community_reactions, matlab_make_valid_name


def test_matlab_prefix_keeps_005_as_string_semantics():
    assert matlab_make_valid_name("005") == "x005"


def test_reaction_classification_external_interface_internal():
    model = SimpleModel(
        reactions=[
            SimpleReaction("R_EX_no3_u", -1000, 1000, {"no3[u]": -1}),
            SimpleReaction("x005__EX_no3_e", -1000, 1000, {"x005__no3_e": -1, "no3[u]": 1}),
            SimpleReaction("x005__NO3tex", -1000, 1000, {"x005__no3_e": -1, "x005__no3_p": 1}),
            SimpleReaction("x005__PGI", -1000, 1000, {"a": -1, "b": 1}),
        ],
        syncomdesign={
            "reactionMap": [
                {"community_rxn": "R_EX_no3_u", "role": "external_medium_exchange"},
                {"community_rxn": "x005__EX_no3_e", "role": "strain_shared_interface"},
                {"community_rxn": "x005__NO3tex", "role": "strain_internal"},
                {"community_rxn": "x005__PGI", "role": "strain_internal"},
            ]
        },
    )
    classes = {row["reaction_id"]: row["classification"] for row in classify_community_reactions(model)}
    assert classes["R_EX_no3_u"] == "external_medium_exchange"
    assert classes["x005__EX_no3_e"] == "strain_shared_interface"
    assert classes["x005__NO3tex"] == "internal_transport"
    assert classes["x005__PGI"] == "metabolic_reaction"


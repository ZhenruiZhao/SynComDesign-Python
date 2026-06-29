from syncomdesign.objectives import scenario_type
from syncomdesign.objectives import add_all_species_active_constraint


def test_id1_to_id5_objective_modes():
    assert scenario_type(1) == "total_biomass"
    assert scenario_type(2) == "target_strain_biomass"
    assert scenario_type(3) == "equal_composition"
    assert scenario_type(4) == "fixed_composition"
    assert scenario_type(5) == "growth_then_n2o_consumption"


def test_active_growth_constraint_is_added_to_solver():
    import cobra

    model = cobra.Model("toy")
    reaction = cobra.Reaction("x005__Growth")
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    model.add_reactions([reaction])
    model.syncomdesign = {"biomassMap": [{"strain": "005", "biomass_rxn": "x005__Growth"}]}

    names = add_all_species_active_constraint(model, 1e-6)

    assert names == ["syncomdesign_min_biomass_x005__Growth"]
    assert names[0] in model.solver.constraints
    assert model.solver.constraints[names[0]].lb == 1e-6

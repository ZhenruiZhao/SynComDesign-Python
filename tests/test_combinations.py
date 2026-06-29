from syncomdesign.combinations import enumerate_all


def test_preserve_strain_id_005():
    combos = enumerate_all(["005"])
    assert combos == [("005",)]


def test_all_combinations_1_2_3_5_strains():
    assert len(enumerate_all(["005"])) == 1
    assert len(enumerate_all(["005", "016"])) == 3
    assert len(enumerate_all(["005", "016", "020"])) == 7
    assert len(enumerate_all(["005", "016", "020", "E10", "ER45"])) == 31


def test_all_combinations_5_strains_equal_31():
    combos = enumerate_all(["005", "016", "020", "E10", "ER45"])
    assert len(combos) == 31
    assert ("005", "016", "020", "E10", "ER45") in combos


def test_target_strain_filter_only_id2():
    strains = ["005", "016", "020"]
    assert len(enumerate_all(strains, target_strain="016", objective_mode=2)) == 4
    assert len(enumerate_all(strains, target_strain="016", objective_mode=1)) == 7
    assert all("016" in combo for combo in enumerate_all(strains, target_strain="016", objective_mode="target_strain_biomass"))


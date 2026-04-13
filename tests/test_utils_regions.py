from dimelo import utils


def test_regions_dict_from_input_accepts_tuple_region_collections():
    tuple_regions = ("chr2:5-10,+", "chr1:0-3,-")
    list_regions = ["chr2:5-10,+", "chr1:0-3,-"]

    tuple_result = utils.regions_dict_from_input(tuple_regions)
    list_result = utils.regions_dict_from_input(list_regions)

    assert tuple_result == list_result

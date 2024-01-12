import pandas as pd
import pytest
from multi_bioservices.biodbnet import InputDatabase, OutputDatabase, TaxonID, db2db

input_db = InputDatabase.GENE_ID
output_db = OutputDatabase.GENE_SYMBOL


@pytest.mark.parametrize(
    "taxon_id,input_values,expected_values", [
        (TaxonID.HOMO_SAPIENS, ["1", "2", "3", "4"], ["A1BG", "A2M", "A2MP1", "-"]),
        (TaxonID.MUS_MUSCULUS, ["14910", "22059", "11816", "21898"], ["Gt(ROSA)26Sor", "Trp53", "Apoe", "Tlr4"])
    ])
def test_fetch_gene_info(taxon_id, input_values, expected_values):
    result: pd.DataFrame = db2db(
        input_values=input_values,
        input_db=input_db,
        output_db=output_db,
        taxon_id=taxon_id,
    )
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == len(input_values)
    assert input_db.value in result.columns
    assert output_db.value in result.columns
    
    for value in expected_values:
        assert value in result[output_db.value].tolist()


@pytest.mark.parametrize("input_values", [["1", "2", "3", "4"], ["1", "1"]])
@pytest.mark.parametrize("remove_duplicates", [True, False])
def test_duplicate_removal(input_values, remove_duplicates):
    result: pd.DataFrame = db2db(
        input_values=input_values,
        input_db=input_db,
        output_db=output_db,
        taxon_id=TaxonID.HOMO_SAPIENS,
        remove_duplicates=remove_duplicates
    )
    
    if remove_duplicates:
        assert len(result) == len(set(input_values))
    else:
        assert len(result) == len(input_values)


@pytest.mark.parametrize("cache", [False, True])
def test_cache(cache):
    input_values = ["1", "2", "3", "4"]
    result: pd.DataFrame = db2db(
        input_values=input_values,
        input_db=input_db,
        output_db=output_db,
        taxon_id=TaxonID.HOMO_SAPIENS,
        cache=False
    )
    print(result)
    assert len(result) == len(input_values)

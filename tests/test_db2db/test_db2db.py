import pandas as pd
import pytest
from multi_bioservices.biodbnet import InputDatabase, OutputDatabase, TaxonID, db2db


# Cache tests
# 1: This uses async_cache, but it will be a cache miss. Retrieve results from server and cache them
# 2: This uses async_cache, but it will be a cache hit. Retrieve results from cache
# 3: This uses biodbnet_cache, but it will be a cache miss. Retrieve results from server and cache them
# 4: This uses biodbnet_cache, but it will be a cache hit. Retrieve results from cache
# 5: This does not use any cache. Retrieve results from server and do not cache them
@pytest.mark.parametrize("taxon_id", [TaxonID.HOMO_SAPIENS, TaxonID.MUS_MUSCULUS])
def test_fetch_gene_info(taxon_id):
    if taxon_id == TaxonID.HOMO_SAPIENS:
        input_values = ["1", "2", "3", "4"]
    elif taxon_id == TaxonID.MUS_MUSCULUS:
        input_values = ["14910", "22059", "11816", "21898"]
    
    result: pd.DataFrame = db2db(
        input_values=input_values,
        input_db=InputDatabase.GENE_ID,
        output_db=OutputDatabase.GENE_SYMBOL,
        taxon_id=taxon_id,
    )
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == len(input_values)
    assert 'Gene Symbol' in result.columns
    assert 'Gene ID' in result.columns
    
    if taxon_id == TaxonID.HOMO_SAPIENS:
        assert result[result["Gene ID"] == "1"]["Gene Symbol"].values[0] == "A1BG"
        assert result[result["Gene ID"] == "2"]["Gene Symbol"].values[0] == "A2M"
        assert result[result["Gene ID"] == "3"]["Gene Symbol"].values[0] == "A2MP1"
        assert result[result["Gene ID"] == "4"]["Gene Symbol"].values[0] == "-"
    elif taxon_id == TaxonID.MUS_MUSCULUS:
        assert result[result["Gene ID"] == "14910"]["Gene Symbol"].values[0] == "Gt(ROSA)26Sor"
        assert result[result["Gene ID"] == "22059"]["Gene Symbol"].values[0] == "Trp53"
        assert result[result["Gene ID"] == "11816"]["Gene Symbol"].values[0] == "Apoe"
        assert result[result["Gene ID"] == "21898"]["Gene Symbol"].values[0] == "Tlr4"


@pytest.mark.parametrize("input_values", [["1", "2", "3", "4"], ["1", "1"]])
@pytest.mark.parametrize("remove_duplicates", [True, False])
def test_duplicate_removal(input_values, remove_duplicates):
    result: pd.DataFrame = db2db(
        input_values=input_values,
        input_db=InputDatabase.GENE_ID,
        output_db=OutputDatabase.GENE_SYMBOL,
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
        input_db=InputDatabase.GENE_ID,
        output_db=OutputDatabase.GENE_SYMBOL,
        taxon_id=TaxonID.HOMO_SAPIENS,
        cache=cache
    )
    assert len(result) == len(input_values)

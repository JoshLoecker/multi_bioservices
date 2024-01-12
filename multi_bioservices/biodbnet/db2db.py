import random
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
from multiprocessing import Queue
from typing import Iterable, List, Union

import pandas as pd
from bioservices import BioDBNet
from multi_bioservices.biodbnet.input_database import InputDatabase
from multi_bioservices.biodbnet.output_database import OutputDatabase
from multi_bioservices.biodbnet.taxon_id import TaxonID
from multi_bioservices.utils import get_biodbnet
from tqdm import tqdm

_DEFAULT_OUTPUT_DB: Iterable[OutputDatabase] = (
    OutputDatabase.GENE_SYMBOL,
    OutputDatabase.GENE_ID,
    OutputDatabase.CHROMOSOMAL_LOCATION
)


def _execute(
    biodbnet: BioDBNet,
    input_values: list[str],
    input_db: str,
    output_db: str,
    taxon_id: int,
    sleep_time: float
    # progress_bar: Union[tqdm, None],
) -> pd.DataFrame | list[str]:
    try:
        time.sleep(sleep_time)
        results = biodbnet.db2db(
            input_values=input_values,
            input_db=input_db,
            output_db=output_db,
            taxon=taxon_id
        )
    except TypeError:
        return input_values
        # if progress_bar is not None:
        #     with _sleep_count.get_lock():
        #         _sleep_count.value += 1
        #         progress_bar.set_postfix_str(s=f"request limit reached - retry count: {_sleep_count.value}")
        #     time.sleep(5)
        #
        # first_half = input_values[:max(1, len(input_values) // 2)]
        # second_half = input_values[max(1, len(input_values) // 2):]
        #
        # first_half_results = biodbnet.db2db(
        #     input_values=first_half,
        #     input_db=input_db,
        #     output_db=output_db,
        #     taxon=taxon_id
        # )
        # second_half_results = biodbnet.db2db(
        #     input_values=second_half,
        #     input_db=input_db,
        #     output_db=output_db,
        #     taxon=taxon_id
        # )
        # results = pd.concat([first_half_results, second_half_results])
    
    return results


def db2db(
    input_values: Union[List[str], List[int]],
    input_db: InputDatabase,
    output_db: Union[OutputDatabase, Iterable[OutputDatabase]] = _DEFAULT_OUTPUT_DB,
    taxon_id: Union[TaxonID, int] = TaxonID.HOMO_SAPIENS,
    remove_duplicates: bool = False,
    chunk_size: int = 300,
    cache: bool = True,
    verbose: bool = False,
    progress_bar: bool = False,
    tqdm_kwargs: dict = None
) -> pd.DataFrame:
    max_threads: int = 5
    input_values: List[str] = list(map(str, input_values))
    taxon_id_value: int = int(taxon_id.value) if isinstance(taxon_id, TaxonID) else int(taxon_id)
    input_db_value: str = input_db.value
    output_db_values: List[str] = [output_db.value] if isinstance(output_db, OutputDatabase) \
        else [str(i.value) for i in output_db]
    
    # Check if input_db_value is in output_db_values
    if input_db_value in output_db_values:
        raise ValueError("Input database cannot be in output database")
    
    # Validate input settings
    if chunk_size > 500 and taxon_id_value == TaxonID.HOMO_SAPIENS.value:
        print(f"Batch length greater than the maximum value of 500 for Homo Sapiens."
              f"Automatically setting batch length to 500")
        chunk_size = 500
    elif chunk_size > 300:
        print(f"Batch length greater than the maximum allowed value for Taxon ID of {taxon_id_value}."
              f"Automatically setting batch length to 300")
        chunk_size = 300
    
    # Perform conversion using BioDBNet's db2db
    conversion_df: pd.DataFrame = pd.DataFrame()
    biodbnet: BioDBNet = get_biodbnet(verbose=verbose, cache=cache)
    pbar = None
    if progress_bar:
        tqdm_kwargs = tqdm_kwargs or {}
        total = tqdm_kwargs.pop("total", len(input_values))
        pbar = tqdm(total=total, **tqdm_kwargs)
    
    partial_func = partial(
        _execute,
        biodbnet=biodbnet,
        input_db=input_db_value,
        output_db=output_db_values,
        taxon_id=taxon_id_value,
        sleep_time=random.uniform(0, 5)
        # progress_bar=pbar
    )
    
    queue_length: int = 0
    item_queue: Queue = Queue()
    for chunk in range(0, len(input_values), chunk_size):
        item_queue.put(input_values[chunk:chunk + chunk_size])
        queue_length += 1
    
    failed_count: int = 0
    sleep_count: int = 0
    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        while not item_queue.empty():
            futures = [
                executor.submit(partial_func, input_values=item_queue.get())
                for _ in range(queue_length)
            ]
            
            for future in as_completed(futures):
                future_result = future.result()
                
                if isinstance(future_result, list):
                    pbar.set_postfix_str(
                        f"failed count: {failed_count}, max threads: {max_threads}, sleep count: {sleep_count}")
                    failed_count += 1
                    item_queue.put(future_result)
                    
                    if failed_count % max_threads == 0:
                        pbar.set_postfix_str(
                            f"failed count: {failed_count}, max threads: {max_threads}, sleep count: {sleep_count}")
                        # Reduce max workers by 1
                        max_threads = max(max_threads - 1, 1)
                        executor._max_workers = max_threads
                        sleep_count += 1
                        time.sleep(5)
                    
                    continue
                
                if progress_bar:
                    pbar.update(chunk_size)
                
                conversion_df = pd.concat([
                    conversion_df,
                    future_result
                ])
    
    conversion_df = conversion_df.reset_index(names=[input_db_value])
    
    if remove_duplicates:
        # Remove rows that have duplicates in the input_db_value
        conversion_df = conversion_df.drop_duplicates(subset=[input_db_value])
    return conversion_df


if __name__ == '__main__':
    print("Starting")
    inputs = [str(i) for i in range(100_000)]
    df = db2db(
        input_values=inputs,
        input_db=InputDatabase.GENE_ID,
        output_db=OutputDatabase.GENE_SYMBOL,
        taxon_id=TaxonID.MUS_MUSCULUS,
        progress_bar=True,
    )

#!/usr/bin/env python3
import copy
from pathlib import Path
from typing import Dict, Tuple, Union

import pandas as pd
from loguru import logger


def get_tophit(
    result_tsv: Path,
    structures: bool,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Process Foldseek output to extract top hit and weighted bitscores.

    Args:
        result_tsv (Path): Path to the Foldseek result TSV file.
        structures (bool): Flag indicating whether structures have been added.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: A tuple containing two DataFrames:
            1. DataFrame containing the top functions extracted from the Foldseek output.
            2. DataFrame containing weighted bitscores for different functions.
    """

    logger.info("Processing Foldseek output")

    if structures:

        col_list = [
            "query",
            "target",
            "bitscore",
            "fident",
            "evalue",
            "qStart",
            "qEnd",
            "qLen",
            "tStart",
            "tEnd",
            "tLen",
            "alntmscore",
            "lddt"
        ]
    else:

        col_list = [
            "query",
            "target",
            "bitscore",
            "fident",
            "evalue",
            "qStart",
            "qEnd",
            "qLen",
            "tStart",
            "tEnd",
            "tLen",
        ]

    foldseek_df = pd.read_csv(
        result_tsv, delimiter="\t", index_col=False, names=col_list
    )
    

    # in case the foldseek output is empty
    if foldseek_df.empty:
        logger.error(
            "Foldseek found no hits whatsoever - please check whether your input is really phage-like"
        )


    # add qcov and tcov 
    foldseek_df["qCov"] = ((foldseek_df["qEnd"] - foldseek_df["qStart"] ) / foldseek_df["qLen"]).round(2)
    foldseek_df["tCov"] = ((foldseek_df["tEnd"] - foldseek_df["tStart"] ) / foldseek_df["tLen"]).round(2)

    # reorder
    qLen_index = foldseek_df.columns.get_loc("qLen")
    tLen_index = foldseek_df.columns.get_loc("tLen")

    new_column_order = (
        list(
            [
                col
                for col in foldseek_df.columns[: qLen_index + 1]
                if col not in ["qCov", "tStart","tEnd",	"tLen", "tCov"]
            ]
        )
        + ["qCov", "tStart","tEnd",	"tLen", "tCov"]
        + list(
            [
                col
                for col in foldseek_df.columns[tLen_index + 1 :]
                if col not in ["qCov", "tStart","tEnd",	"tLen", "tCov"]
            ]
        )
    )
    foldseek_df = foldseek_df.reindex(columns=new_column_order)

    # get only the tophit - will always be the first hit for each query (top bitscore)
    foldseek_df = foldseek_df.drop_duplicates(subset="query", keep="first")

    # print(foldseek_df.head())


    return foldseek_df


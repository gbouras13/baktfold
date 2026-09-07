#!/usr/bin/env python3
from pathlib import Path

import polars as pl
from loguru import logger


def get_tophit(
    result_tsv: Path,
    structures: bool,
    cath: bool = False
) -> pl.DataFrame:
    """
    Process Foldseek output to extract the top hit per query.

    Args:
        result_tsv (Path): Path to the Foldseek result TSV file.
        structures (bool): Flag indicating whether structures have been added.
        cath (bool): Flag indicating whether this is for CATH database (all greedy besthits kept not just top)

    Returns:
        pl.DataFrame: DataFrame containing the top hit(s) extracted from the Foldseek output.
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

    # infer_schema_length=None scans the whole file so dtype inference matches
    # pandas' (which read the whole column) — keeps the output byte-identical.
    try:
        foldseek_df = pl.read_csv(
            result_tsv,
            separator="\t",
            has_header=False,
            new_columns=col_list,
            infer_schema_length=None,
        )
    except pl.exceptions.NoDataError:
        # empty Foldseek result (0-byte file) — mirror pandas' empty frame
        foldseek_df = pl.DataFrame(schema={c: pl.Utf8 for c in col_list})

    # replace ~PIPE~ with |
    foldseek_df = foldseek_df.with_columns(
        pl.col("query").str.replace_all("~PIPE~", "|", literal=True)
    )

    # in case the foldseek output is empty
    if foldseek_df.is_empty():
        logger.warning(
            "Foldseek found no hits whatsoever - please check your input if you expect hits"
        )
        return foldseek_df

    # add qcov and tcov (rounded to 2dp). evalue is rendered with Python's
    # float repr so the written TSV is byte-identical to the previous pandas
    # output (pandas/Python pad scientific exponents to 2 digits, polars does
    # not — e.g. '1.5e-08' vs '1.5e-8'). repr() round-trips losslessly so the
    # numeric value consumed downstream by pstc.parse is unchanged.
    foldseek_df = foldseek_df.with_columns(
        ((pl.col("qEnd") - pl.col("qStart")) / pl.col("qLen")).round(2).alias("qCov"),
        ((pl.col("tEnd") - pl.col("tStart")) / pl.col("tLen")).round(2).alias("tCov"),
        pl.col("evalue").map_elements(lambda v: repr(float(v)), return_dtype=pl.Utf8),
    )

    # reorder: qCov directly after qLen, the tStart/tEnd/tLen/tCov block
    # together; any trailing structure columns (alntmscore, lddt) stay at the end.
    front = col_list[: col_list.index("qLen") + 1]
    tail = col_list[col_list.index("tLen") + 1 :]
    new_column_order = front + ["qCov", "tStart", "tEnd", "tLen", "tCov"] + tail
    foldseek_df = foldseek_df.select(new_column_order)

    # Break score ties deterministically before anything picks a "first" row.
    # Foldseek emits hits per query in descending-bitscore order, but the order
    # *within* a run of equal bitscores is not stable across runs (it falls out
    # of prefilter/thread scheduling). When the top two hits tie, `keep="first"`
    # therefore picked a different target run to run — e.g. the golden outputs
    # flipped between the equally-scoring PDB entries "MsDpo4-DNA complex 1"
    # and "MsDpo4-DNA complex 2", failing the golden tests with no code change.
    #
    # Sort within each query by (bitscore desc, target, original position),
    # keeping queries themselves in their original order — so for distinct
    # bitscores the row order Foldseek produced is preserved exactly and only
    # genuine ties are reordered.
    if not foldseek_df.is_empty():
        foldseek_df = foldseek_df.with_row_index("_orig")
        query_order = foldseek_df.group_by("query", maintain_order=True).agg(
            pl.col("_orig").min().alias("_query_order")
        )
        foldseek_df = (
            foldseek_df.join(query_order, on="query", how="left")
            .sort(
                by=[
                    pl.col("_query_order"),
                    pl.col("bitscore").cast(pl.Float64, strict=False),
                    pl.col("target"),
                    pl.col("_orig"),
                ],
                descending=[False, True, False, False],
            )
            .drop(["_orig", "_query_order"])
        )

    if not cath:
        # get only the tophit - always the first (top-bitscore) hit per query.
        # maintain_order=True preserves the order established above so "first"
        # picks the same survivor pandas' drop_duplicates(keep="first") did.
        foldseek_df = foldseek_df.unique(subset="query", keep="first", maintain_order=True)
    # otherwise, the df will contain all greedy tophits from CATH

    return foldseek_df


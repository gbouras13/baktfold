#!/usr/bin/env python3
import os
from pathlib import Path
from typing import Optional

from pholdlib.prostt5.device import cuda_visible_devices_value, parse_gpus

from baktfold.utils.external_tools import ExternalTool
from loguru import logger


def gpu_target_db(target_db) -> str:
    """The Foldseek-GPU (padded) version of a target database.

    Baktfold's own databases ship a padded copy as '<db>_gpu'. A custom
    database may instead already *be* the padded database (the docs build it
    with makepaddedseqdb as 'custom_db/db_pad'); appending '_gpu' to that
    pointed at a file that does not exist and the run died after the four
    main searches. Use '<db>_gpu' when it exists, otherwise the path as given.
    """
    if Path(f"{target_db}_gpu.dbtype").exists():
        return f"{target_db}_gpu"
    logger.info(f"No {target_db}_gpu found; using {target_db} as the padded Foldseek-GPU database")
    return str(target_db)


def physical_cuda_visible_devices(cvd: Optional[str]) -> Optional[str]:
    """Translate torch-relative CUDA indices to the ones to give a subprocess.

    torch numbers GPUs within the CUDA_VISIBLE_DEVICES this process inherited
    (e.g. SLURM sets '2,3', so cuda:1 is physical GPU 3). Writing '1' back for
    the Foldseek subprocess would select physical GPU 1 instead, which may not
    be allocated to the job. Map through the inherited list when there is one.
    """
    parent = os.environ.get("CUDA_VISIBLE_DEVICES", "")
    visible = [x.strip() for x in parent.split(",") if x.strip()]
    if cvd is None or not visible:
        return cvd
    mapped = []
    for idx in cvd.split(","):
        i = int(idx)
        if i >= len(visible):
            raise ValueError(f"--gpus index {i} is outside CUDA_VISIBLE_DEVICES={parent}")
        mapped.append(visible[i])
    return ",".join(mapped)


def run_foldseek_search(
    query_db: Path,
    target_db: Path,
    result_db: Path,
    temp_db: Path,
    threads: int,
    logdir: Path,
    evalue: float,
    sensitivity: float,
    max_seqs: int,
    ultra_sensitive: bool,
    extra_foldseek_params: str,
    foldseek_gpu: bool,
    structures: bool,
    gpus: Optional[str] = None,
    report_prob: bool = False,
) -> None:
    """
    Run a Foldseek search using given parameters.

    Args:
        query_db (Path): Path to the query database.
        target_db (Path): Path to the target database.
        result_db (Path): Path to store the result database.
        temp_db (Path): Path to store temporary files.
        threads (int): Number of threads to use for the search.
        logdir (Path): Path to the directory where logs will be stored.
        evalue (float): E-value threshold for the search.
        sensitivity (float): Sensitivity threshold for the search.
        max_seqs (int): Maximum results per query sequence allowed to pass the prefilter for foldseek.
        ultra_sensitive (bool): Whether to skip foldseek prefilter for maximum sensitivity
        extra_foldseek_params (str): Extra foldseek search params
        foldseek_gpu (bool): Run Foldseek-GPU with accelerate ungapped prefilter
        structures (bool): Run Foldseek with structures, not ProstT5 3Dis
        report_prob (bool): Emit Foldseek's homology probability. Requires the
            backtrace CIGAR (``-a 1``); without it convertalis aborts with
            "Backtrace cigar is missing in the alignment result". Costs extra
            alignment time and result-DB size, so it is opt-in.
        gpus (Optional[str]): Comma-separated CUDA indices (e.g. "0,2") to
            restrict foldseek's GPU prefilter to a subset of devices. When
            ``foldseek_gpu`` is True and this resolves to ≥1 CUDA device,
            the foldseek subprocess gets ``CUDA_VISIBLE_DEVICES`` set
            accordingly. None = use all visible CUDA GPUs (foldseek default).
            Ignored when ``foldseek_gpu`` is False.

    Returns:
        None
    """

    if ultra_sensitive:
        cmd = f"search {query_db} {target_db} {result_db} {temp_db} --threads {str(threads)} -e {evalue} -s {sensitivity} --exhaustive-search"
    else:
        cmd = f"search {query_db} {target_db} {result_db} {temp_db} --threads {str(threads)} -e {evalue} -s {sensitivity} --max-seqs {max_seqs}"

    # support foldseek gpu only for the regular DB search for now
    if foldseek_gpu:
        if ultra_sensitive or float(sensitivity) != 9.5:
            # GPU mode uses the ungapped prefilter; -s / --exhaustive-search do
            # not apply. Say so rather than logging settings that were ignored.
            logger.warning(
                "--foldseek-gpu uses Foldseek's GPU ungapped prefilter: "
                "-s/--sensitivity and --ultra-sensitive are not applied in this mode."
            )
        cmd = f"search {query_db} {gpu_target_db(target_db)} {result_db} {temp_db} --threads {str(threads)} -e {evalue}  --gpu 1 --prefilter-mode 1 --max-seqs {max_seqs}"

    if extra_foldseek_params:
        cmd += f" {extra_foldseek_params}"

    # need -a 1 to compute the alignment so tmscore and lddt can be output (if using --structures).
    # `prob` also needs the backtrace (and query C-alpha coordinates, so structures only:
    # see validation.check_report_prob).
    if structures or report_prob:
        cmd += f" -a 1"

    # Build optional env for multi-GPU foldseek. Only applies when GPU mode is
    # on; foldseek selects devices via CUDA_VISIBLE_DEVICES (per its README).
    env = None
    if foldseek_gpu and gpus is not None:
        devices = parse_gpus(cpu=False, gpus=gpus)
        cvd = physical_cuda_visible_devices(cuda_visible_devices_value(devices))
        if cvd is not None:
            env = {"CUDA_VISIBLE_DEVICES": cvd}

    foldseek_search = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"{cmd}",
        logdir=logdir,
        env=env,
    )

    ExternalTool.run_tool(foldseek_search)




def create_result_tsv(
    query_db: Path, target_db: Path, result_db: Path, result_tsv: Path, logdir: Path, foldseek_gpu: bool, structures: bool, threads: int, report_prob: bool = False
) -> None:
    """
    Create a TSV file containing the results of a Foldseek search.

    Args:
        query_db (Path): Path to the query database.
        target_db (Path): Path to the target database.
        result_db (Path): Path to the result database generated by the search.
        result_tsv (Path): Path to save the resulting TSV file.
        logdir (Path): Path to the directory where logs will be stored.
        foldseek_gpu (bool): Run Foldseek-GPU with accelerate ungapped prefilter
        structures (bool): Whether structures were input (not ProstT5)
        threads (int): Number of threads to use.
        report_prob (bool): Append Foldseek's homology probability column. The
            search must have been run with ``-a 1`` (see run_foldseek_search).

    Returns:
        None
    """
    columns = "query,target,bits,fident,evalue,qstart,qend,qlen,tstart,tend,tlen"
    if structures:
        columns += ",alntmscore,lddt"
    if report_prob:
        # kept last so tophit.get_tophit's column reorder leaves it in the tail
        columns += ",prob"
    format_string = f"--format-output {columns}"
    if foldseek_gpu:
        target_db = gpu_target_db(target_db)

    
    cmd = f"convertalis {query_db} {target_db} {result_db} {result_tsv} {format_string} --threads {threads}"

    foldseek_createtsv = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"{cmd}",
        logdir=logdir,
    )


    ExternalTool.run_tool(foldseek_createtsv)

def summarise_hits(result_db: Path, result_db_greedy_best_hits: Path, logdir: Path, threads: int) -> None:
    """
    Get all non-overlapping tophits covering a query (designed for CATH)

    Args:
        result_db (Path): Path to the result database generated by the search.
        result_db_greedy_best_hits (Path): Path to save the greedy best hits results db.
        logdir (Path): Path to the directory where logs will be stored.
        threads (int): Number of threads to use.

    Returns:
        None
    """

    cmd = f"summarizeresult  {result_db} {result_db_greedy_best_hits} --threads {threads} -a 1"

    foldseek_summarizeresult = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"{cmd}",
        logdir=logdir,
    )

    ExternalTool.run_tool(foldseek_summarizeresult)
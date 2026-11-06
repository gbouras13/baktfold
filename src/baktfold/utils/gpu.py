#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Accelerator memory-release helpers.

``baktfold run`` loads ProstT5 onto the GPU (``subcommand_predict``) and then
shells out to ``foldseek ... --gpu 1`` (``subcommand_compare``) **from the same
Python process**. Foldseek is a separate process with its own CUDA context, so
it can only use memory that PyTorch has handed back to the driver.

PyTorch's caching allocator does not do that on its own: when the ProstT5
weights and cached activations are dropped, the blocks are returned to the
*allocator's* free pool but stay reserved against the device, so ``nvidia-smi``
still shows them as in use. Foldseek then starts on a card that looks full and
dies with::

    CUDA error: out of memory : /work/lib/mmseqs/lib/libmarv/src/cudasw4.cuh, line 1276
    Error: Ungapped prefilter matching step died

Calling :func:`release_gpu_memory` once the model references have gone out of
scope releases the reserved-but-unused blocks, which is why running ``baktfold
predict`` and ``baktfold compare`` as two separate processes always worked
(process exit releases everything).
"""

import gc

import torch
from loguru import logger

__all__ = ["release_gpu_memory"]

_GIB = 1024 ** 3


def _cuda_reserved_bytes() -> int:
    """Bytes currently reserved by the CUDA caching allocator (0 if unknown)."""
    try:
        return int(torch.cuda.memory_reserved())
    except Exception:  # no initialised CUDA context, driver hiccup, mocked torch
        return 0


def release_gpu_memory(context: str = "") -> None:
    """Hand cached-but-unused accelerator memory back to the driver.

    Safe to call unconditionally: on a CPU-only run this is a ``gc.collect()``
    and nothing else. Call it *after* the last reference to the model has gone
    out of scope — freeing the allocator's cache does not free memory that is
    still owned by a live tensor.

    Args:
        context: Short description of what has just finished, used in the log
            line only (e.g. ``"after ProstT5 3Di prediction"``).
    """
    # nn.Module graphs are full of reference cycles (hooks, parent/child links),
    # so plain refcounting does not always reclaim a model the instant its last
    # name goes out of scope. Collect first, or empty_cache() has nothing to do.
    gc.collect()

    where = f" {context}" if context else ""

    if torch.cuda.is_available():
        before = _cuda_reserved_bytes()
        torch.cuda.empty_cache()
        after = _cuda_reserved_bytes()
        logger.info(
            f"Released cached GPU memory{where}: reserved "
            f"{before / _GIB:.2f} GiB -> {after / _GIB:.2f} GiB"
        )
        return

    # Apple Silicon. torch.backends.mps exists on every platform from torch 1.12
    # and simply reports False off macOS, but guard anyway for older builds.
    mps_backend = getattr(torch.backends, "mps", None)
    if mps_backend is not None and mps_backend.is_available():
        mps = getattr(torch, "mps", None)
        if mps is not None and hasattr(mps, "empty_cache"):
            mps.empty_cache()
            logger.info(f"Released cached MPS memory{where}")
        return

    # Intel XPU. torch.xpu is absent from builds without XPU support.
    xpu = getattr(torch, "xpu", None)
    if xpu is not None and xpu.is_available():
        if hasattr(xpu, "empty_cache"):
            xpu.empty_cache()
            logger.info(f"Released cached XPU memory{where}")
        return

    logger.debug(f"No GPU in use — nothing to release{where}")

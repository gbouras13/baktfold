
FROM quay.io/pawsey/pytorch:2.7.1-rocm6.3.3

## install foldseek only non python dep

RUN apt-get update && apt-get install -y \
    wget \
    tar \
    aria2 \
    && rm -rf /var/lib/apt/lists/*

# install foldseek
RUN wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz && \
    tar -xzf foldseek-linux-avx2.tar.gz && \
    rm foldseek-linux-avx2.tar.gz && \
    mv foldseek /opt/foldseek

ENV PATH="/opt/foldseek/bin:${PATH}"

# install baktfold

RUN pip install --ignore-installed PyYAML
RUN python3 -c "import torch; print(torch.__version__)" || echo "PyTorch is not installed before installing the project. "
ARG BAKTFOLD_REF
# Install baktfold together with its ``test`` extra so the test dependencies
# (pytest, pytest-cov, pytest-snapshot) stay in lockstep with pyproject.toml
# instead of being duplicated here. pytest-snapshot is mandatory: the unit
# snapshot tests AND the golden output regression tests
# (tests/unit/test_json_golden.py) rely on the ``snapshot`` fixture it provides
# — without it every snapshot/golden test errors at setup with
# `fixture 'snapshot' not found`.
RUN pip install "baktfold[test] @ git+https://github.com/gbouras13/baktfold.git@${BAKTFOLD_REF}"

# smoke tests that the container is wired up correctly
RUN foldseek --help
RUN baktfold -h
RUN baktfold run -h
RUN python3 -c "import pytest, pytest_snapshot, click.testing"  # test harness present

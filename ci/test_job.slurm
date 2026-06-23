#! /bin/bash
#SBATCH --job-name=hpcitestjob
#SBATCH --partition=gpu-dev
#SBATCH --nodes=1
#SBATCH --gres=gpu:1  
#SBATCH --time=03:59:30

set -eu -o pipefail

mkdir -p "$BASEDIR"
cd "$BASEDIR"

# # pull the sif 

containerImage="$BASEDIR/baktfold_${COMMIT_SHA}.sif"
module load singularity/3.11.4-slurm
export SINGULARITY_CACHEDIR=$PWD

# Only pull if the image doesn't already exist
if [ ! -f "$containerImage" ]; then
    echo "Container not found. Pulling image..."
    singularity pull --dir "$PWD" "docker://quay.io/gbouras13/baktfold:$COMMIT_SHA"
else
    echo "Container already exists. Skipping pull."
fi

# run help to ensure baktfold is installed
singularity exec "$containerImage" baktfold --help

# Define the target directory name
TARGET_DIR="${REPO_NAME}-${COMMIT_SHA}"

git clone "$REPO_URL" "$TARGET_DIR"

cd "$TARGET_DIR"

git checkout "$COMMIT_SHA"

# If a persistent baktfold database has been pre-staged under $BASEDIR, link it
# in where the integration tests (and the conftest DB gate) expect it, so the
# DB-backed integration tests run instead of auto-skipping. If it is absent the
# DB-backed tests skip cleanly and the DB-free + golden suites still run.
if [ -d "$BASEDIR/baktfold_db" ]; then
    echo "Staging pre-installed database from $BASEDIR/baktfold_db"
    mkdir -p tests/test_data
    ln -sfn "$BASEDIR/baktfold_db" tests/test_data/baktfold_db
fi

# Run both suites. Disable -e between them so one failing does not stop the
# other from running (we always want both the golden and integration tests to
# run), then fail the job if either failed.
set +e

# 1. Unit + golden output regression tests: DB-free, no GPU, the primary gate.
echo "=== unit + golden regression tests ==="
singularity exec "$containerImage" pytest tests/unit -v
unit_rc=$?

# 2. Integration tests: exercise the real CLI on GPU + database. DB-backed
#    tests auto-skip when the database is absent (see tests/conftest.py).
echo "=== integration tests ==="
singularity exec "$containerImage" pytest tests/test_integration.py --gpu-available --threads 8 -v
integration_rc=$?

set -e

cd ..

rm -rf $TARGET_DIR

echo "unit/golden exit=$unit_rc  integration exit=$integration_rc"
[ "$unit_rc" -eq 0 ] && [ "$integration_rc" -eq 0 ]


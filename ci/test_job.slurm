#! /bin/bash
#SBATCH --job-name=hpcitestjob
#SBATCH --partition=gpu-dev
#SBATCH --nodes=1
#SBATCH --gres=gpu:1  
#SBATCH --time=03:59:30

set -eu -o pipefail

mkdir -p "$BASEDIR"
cd "$BASEDIR"

# pull the sif 

containerImage="$BASEDIR/baktfold_${COMMIT_SHA}.sif"
module load singularity/3.11.4-slurm
export SINGULARITY_CACHEDIR=$PWD

# Only pull if the image doesn't already exist
if [ ! -f "$containerImage" ]; then
    echo "Container not found. Pulling image..."
    singularity pull --dir "$PWD" "docker://quay.io/gbouras13/baktfold:$COMMIT_SHA"
    echo "Container pull done"
else
    echo "Container already exists. Skipping pull."
fi

# run help to ensure baktfold is installed
singularity exec "$containerImage" baktfold --help

# Define the target directory name
TARGET_DIR="${REPO_NAME}-${COMMIT_SHA}"

git clone "$REPO_URL" "$TARGET_DIR"

cd "$TARGET_DIR"

# For pull_request runs, github.sha is GitHub's temporary merge commit
# (refs/pull/<n>/merge), which a plain clone does not contain: checkout failed
# with exit 128 before any test ran. Fetch the PR refs so it can be found.
git fetch --quiet origin "+refs/pull/*:refs/remotes/origin/pull/*" || true

git checkout "$COMMIT_SHA"

singularity exec "$containerImage" pytest --gpu-available --threads 8 .

# golden output regression tests
singularity exec "$containerImage" pytest tests/unit/test_json_golden.py -v

cd ..

# rm -rf $TARGET_DIR


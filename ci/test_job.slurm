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

module load pawseyenv/2023.08
module load singularity/3.11.4-slurm
export SINGULARITY_CACHEDIR=$PWD

# Only pull if the image doesn't already exist
if [ ! -f "$containerImage" ]; then
    echo "Container not found. Pulling image..."
    singularity pull --dir "$PWD" "docker://quay.io/gbouras13/baktfold:$COMMIT_SHAA"
else
    echo "Container already exists. Skipping pull."
fi

singularity exec "$containerImage" baktfold --help

# Define the target directory name
TARGET_DIR="${REPO_NAME}-${COMMIT_SHA}"

git clone "$REPO_URL" "$TARGET_DIR"

cd "$TARGET_DIR"

git checkout "$COMMIT_SHA"

singularity exec "$containerImage" pytest --gpu-available --threads 8 .

cd ..

rm -rf $TARGET_DIR


#! /bin/bash
#SBATCH --job-name testjob
#SBATCH --partition=gpu-dev
#SBATCH --nodes=1
#SBATCH --gres=gpu:1  
#SBATCH --time=03:59:30
#SBATCH --account="$ACCOUNT"-gpu  


set -eu -o pipefail

# Define the target directory name
TARGET_DIR="${REPO_NAME}-${COMMIT_SHA}"

git clone "$REPO_URL" "$TARGET_DIR"

cd "$TARGET_DIR"

git checkout "$COMMIT_SHA"

echo "Working Dir: $(pwd); Sbatch version: $(sbatch --version)" > test_job.log 2>&1

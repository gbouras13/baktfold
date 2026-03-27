#! /bin/bash
#SBATCH --job-name testjob
#SBATCH -p a100
#SBATCH -N 1              	                             
#SBATCH -c 8              	                                
#SBATCH --time=03:59:50                                                          
#SBATCH --mem=40GB                                 
#SBATCH --gres=gpu:1 

set -eu -o pipefail

# Define the target directory name
TARGET_DIR="${REPO_NAME}-${COMMIT_SHA}"

git clone "$REPO_URL" "$TARGET_DIR"

cd "$TARGET_DIR"

git checkout "$COMMIT_SHA"

echo "Working Dir: $(pwd); Sbatch version: $(sbatch --version)" > test_job.log 2>&1

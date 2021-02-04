# Change working directory to MAS base path
cd "${0%/*}"

# Activate mas-worker conda environment
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)";
conda activate mas-worker;

# Start worker
env $(grep -v '^#' .env | xargs) CELERY_WORKER=TRUE celery -b 0.0.0.0 -A MAS worker -B --concurrency=5;
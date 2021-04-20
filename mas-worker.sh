#!/bin/bash 


# Change working directory to MAS base path
#echo "${0%/*}"
#source ~luigi/.bashrc
#cd "${0%/*}"

# Activate mas-worker conda environment
#eval "$( command  conda 'shell.bash' 'hook' 2>/dev/null )";

# /export/anaconda3/bin/conda activate mas-worker;
#conda activate mas-worker;

# Start/stop/restart worker
#env $(grep -v '^#' .env | xargs) CELERY_WORKER=TRUE celery multi ${1} 1 -b 0.0.0.0 -A MAS -B --concurrency=5 --pidfile=/home/luigi/MAS/worker.pid --logfile=/home/luigi/MAS/worker.log;
#env $(grep -v '^#' .env | xargs) CELERY_WORKER=TRUE celery multi $arg 1 -b "amqp://mas:${RABBITMQ_DEFAULT_PASS}@0.0.0.0" -A MAS -B --concurrency=5 --pidfile=/home/luigi/MAS/worker.pid --logfile=/home/luigi/MAS/worker.log;
subcmd=$1
wd="${0%/*}"
cd "$wd" || exit

source .env

$CONDA_EXE run -n mas-worker --cwd "$wd" env $(grep -v '^#' .env | xargs) CELERY_WORKER=TRUE celery multi "$subcmd" 1 -b "amqp://mas:${RABBITMQ_DEFAULT_PASS}@0.0.0.0" -A MAS -B --concurrency=5 --pidfile="$MAS_WORKER_PID_FILE" --logfile="$MAS_WORKER_LOG";

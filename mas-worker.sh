#!/bin/bash 

subcmd=$1
wd="${0%/*}"
cd "$wd" || exit

source .env

$CONDA_EXE run -n mas-worker --cwd "$wd" env $(grep -v '^#' .env | xargs) CELERY_WORKER=TRUE celery multi "$subcmd" mas-worker --hostname=host -A MAS -B --concurrency=5 --pidfile="$MAS_WORKER_PID_FILE" --logfile="$MAS_WORKER_LOG";

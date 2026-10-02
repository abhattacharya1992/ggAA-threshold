#!/usr/bin/env bash
# Run `make reproduce` and stop it if it uses too much memory.
#
#   bash scripts/watch_reproduce.sh [max_job_GB] [min_available_GB]
#
# Every 30 s the memory (RSS) of all processes started by the run is added up, including
# the Mathematica kernel that the construction script starts in its own session, and the
# memory still available on the machine is read from /proc/meminfo. The run is stopped
# (SIGTERM, then SIGKILL after 30 s) if the job uses more than max_job_GB (default 60) or
# the machine has less than min_available_GB left (default 15). Readings go to
# reproduce_mem.log, the output of the run to reproduce.log. Rerunning `make reproduce`
# resumes from the last finished stage. CMD overrides the command (for testing).
#
# Start it detached from the terminal, so that logging out cannot stop the run:
#   setsid nohup bash scripts/watch_reproduce.sh 60 15 > /dev/null 2>&1 < /dev/null &

set -u
MAX_JOB_GB=${1:-60}
MIN_AVAIL_GB=${2:-15}
CMD=${CMD:-make reproduce}
INTERVAL=${INTERVAL:-30}
LOG=reproduce.log
MEMLOG=reproduce_mem.log

descendants() {                      # the process and all its descendants, also across sessions
    echo "$1"
    local child
    for child in $(ps -o pid= --ppid "$1" 2> /dev/null); do
        descendants "$child"
    done
}
gb() { awk -v k="$1" 'BEGIN { printf "%.1f", k / 1048576 }'; }

nice -n 19 $CMD > "$LOG" 2>&1 < /dev/null &
top=$!
echo "$(date '+%F %T') started '$CMD' (pid $top); limits: job ${MAX_JOB_GB} GB, available ${MIN_AVAIL_GB} GB" | tee "$MEMLOG"

peak=0
status="finished"
while kill -0 "$top" 2> /dev/null; do
    pids=$(descendants "$top" | paste -sd, -)
    job_kb=$(ps -o rss= -p "$pids" 2> /dev/null | awk '{ s += $1 } END { print s + 0 }')
    avail_kb=$(awk '/^MemAvailable:/ { print $2 }' /proc/meminfo)
    (( job_kb > peak )) && peak=$job_kb
    echo "$(date '+%F %T') job $(gb "$job_kb") GB, machine available $(gb "$avail_kb") GB" >> "$MEMLOG"
    if awk -v j="$job_kb" -v a="$avail_kb" -v mj="$MAX_JOB_GB" -v ma="$MIN_AVAIL_GB" \
           'BEGIN { exit !(j > mj * 1048576 || a < ma * 1048576) }'; then
        echo "$(date '+%F %T') LIMIT REACHED: job $(gb "$job_kb") GB, available $(gb "$avail_kb") GB; stopping the run" \
            | tee -a "$MEMLOG" >> "$LOG"
        all=$(descendants "$top")
        kill -TERM $all 2> /dev/null
        sleep "$INTERVAL"
        kill -KILL $(descendants "$top") $all 2> /dev/null
        status="STOPPED"
        break
    fi
    sleep "$INTERVAL"
done
wait "$top" 2> /dev/null; rc=$?
if [ "$status" = "finished" ] && [ "$rc" -ne 0 ]; then status="FAILED (exit code $rc, see $LOG)"; fi
echo "$(date '+%F %T') $status; peak memory of the job $(gb "$peak") GB" | tee -a "$MEMLOG"

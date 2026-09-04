#!/bin/bash
# Shared logging helpers for kontiguity's shell scripts.
#
# Usage:
#   local_path=$(realpath "$0")
#   local_dir="${local_path%/*}"
#   source "$local_dir/../lib/log.sh"     # adjust ../ depth to reach scripts/lib
#   log_info "message"
#   log_warn "message"
#   log_error "message"
#
# All three ALWAYS append a timestamped line to the GLOBAL_LOG file (exported
# by the generated script headers, see utils/functions.py:get_header) so a
# whole run (across load/retrieve/map/describe/pipeline, across every
# species/job) can be followed from one place. They only ALSO print to the
# terminal if GLOBAL_VERBOSE is "true" (set via each command's --verbose
# flag) - by default nothing is printed to stdout/stderr, since a large
# dataset can produce a lot of lines that are more useful in the log file
# than scrolling past in the terminal.

_kontiguity_log() {
    local level="$1"; shift
    local msg="$*"
    local ts
    ts=$(date '+%Y-%m-%d %H:%M:%S')
    local line="${ts} [${level}] ${msg}"

    if [ -n "${GLOBAL_LOG:-}" ]; then
        mkdir -p "$(dirname "$GLOBAL_LOG")" 2>/dev/null
        echo "$line" >> "$GLOBAL_LOG"
    fi
    if [ "${GLOBAL_VERBOSE:-false}" = "true" ]; then
        echo "$line"
    fi
}

log_info()  { _kontiguity_log "INFO"  "$@"; }
log_warn()  { _kontiguity_log "WARN"  "$@" >&2; }
log_error() { _kontiguity_log "ERROR" "$@" >&2; }

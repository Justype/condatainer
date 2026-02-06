#!/bin/bash
# Common functions shared by all CondaTainer helper scripts.
# Source this at the top of every helper script:
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
#   source "$SCRIPT_DIR/.common.sh"

# ============= Colors =============
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color
# ================= Message helpers =================
print_msg()  { echo -e "[MSG] $*"; }
print_info() { echo -e "[${CYAN}INFO${NC}] $*"; }
print_warn() { echo -e "[${YELLOW}WARN${NC}] $*"; }
print_error(){ echo -e "[${RED}ERR${NC}] $*" >&2; }
print_pass(){ echo -e "[${GREEN}PASS${NC}] $*"; }
trap 'echo; exit 130' INT # Add a newline on Ctrl+C and exit with code 130
# ============= Directories =============
# Follow XDG Base Directory spec for config and state locations
CONDATAINER_CONFIG_DIR="${XDG_CONFIG_HOME:-$HOME/.config}/condatainer"
HELPER_DEFAULTS_DIR="$CONDATAINER_CONFIG_DIR/helper"
HELPER_STATE_DIR="${XDG_STATE_HOME:-$HOME/.local/state}/condatainer/helper"
LOG_DIR="$HOME/logs"
mkdir -p "$CONDATAINER_CONFIG_DIR" "$HELPER_DEFAULTS_DIR" "$HELPER_STATE_DIR" "$LOG_DIR"

# ============= Config Functions =============

# config_load <helper-name>
#   Sources XDG config: $XDG_CONFIG_HOME/condatainer/helper/<name> (defaults to ~/.config/condatainer/helper/<name>) if it exists.
#   Updates the current shell environment with saved variables and set CWD to pwd.
#   Returns 0 if loaded, 1 if no saved config.
config_load() {
    local f="$HELPER_DEFAULTS_DIR/$1"
    if [ -f "$f" ]; then
        source "$f"
        CWD=$(readlink -f .)
        return 0
    fi
    return 1
}

# config_init <helper-name> KEY=VALUE ...
#   Creates the defaults file on first run with the given key=value pairs.
#   If the file already exists, appends any missing keys without overwriting existing ones.
config_init() {
    local name="$1"; shift
    local f="$HELPER_DEFAULTS_DIR/$name"
    touch "$f"
    for pair in "$@"; do
        local key="${pair%%=*}"
        local val="${pair#*=}"
        grep -q "^${key}=" "$f" || printf '%s="%s"\n' "$key" "$val" >> "$f"
    done
}

# config_update <helper-name> KEY=VALUE ...
#   Updates specific keys in the defaults file without rewriting the whole file.
#   If the key exists, its value is replaced. If not, the key is appended.
config_update() {
    local name="$1"; shift
    local f="$HELPER_DEFAULTS_DIR/$name"
    touch "$f"
    for pair in "$@"; do
        local key="${pair%%=*}"
        local val="${pair#*=}"
        if grep -q "^${key}=" "$f"; then
            sed -i "s|^${key}=.*|${key}=\"${val}\"|" "$f"
        else
            printf '%s="%s"\n' "$key" "$val" >> "$f"
        fi
    done
}

# config_require <helper-name> VAR1 VAR2 ...
#   Checks that each named variable is set and non-empty. Exits on failure.
config_require() {
    local name="$1"; shift
    for var in "$@"; do
        if [ -z "${!var}" ]; then
            print_error "Required config variable ${YELLOW}$var${NC} is not set."
            print_info "Please check your config file: ${BLUE}$(config_path $name)${NC}"
            print_info "Or delete it to regenerate defaults."
            exit 1
        fi
    done
}

# config_show <helper-name>
#   Prints config file path and contents, then exits.
config_show() {
    local f="$HELPER_DEFAULTS_DIR/$1"
    echo -e "Config file: ${BLUE}$f${NC}"
    if [ -f "$f" ]; then
        echo "-------------"
        cat "$f"
    else
        echo "(no config file yet; will be created on first run)"
    fi
    exit 0
}

# config_path <helper-name>
#   Prints the path to the defaults file.
config_path() {
    echo "$HELPER_DEFAULTS_DIR/$1"
}

# state_path <helper-name>
#   Prints the path to the state file.
state_path() {
    echo "$HELPER_STATE_DIR/$1"
}

# resolve_overlay_list <colon-separated-overlays>
#   Converts file-path overlays to absolute paths, leaves named overlays as-is.
#   Takes colon-separated string, returns colon-separated string.
resolve_overlay_list() {
    local input="$1"
    [ -z "$input" ] && return 0

    local result=""
    local IFS=':'
    local -a arr=($input)

    for o in "${arr[@]}"; do
        [ -z "$o" ] && continue
        if echo "$o" | grep -qiE '\.(sqf|sqsh|squashfs|img)$'; then
            local abs_path=$(readlink -f "$o")
            result="${result:+$result:}$abs_path"
        else
            result="${result:+$result:}$o"
        fi
    done
    echo "$result"
}

# build_overlays_arg <colon-separated-overlays>
#   Rebuilds OVERLAYS_ARG from colon-separated overlay list.
#   Takes colon-separated string, returns argument string.
build_overlays_arg() {
    local input="$1"
    [ -z "$input" ] && return 0

    local arg=""
    local IFS=':'
    local -a arr=($input)

    for o in "${arr[@]}"; do
        [ -n "$o" ] && arg+=" -o \"$o\""
    done
    echo "$arg"
}

# ============= Prompt Functions =============

# confirm_default_yes <prompt>
#   Prompts the user with [Y/n]. Returns 0 if yes (default), 1 if no.
confirm_default_yes() {
    local prompt="$1"
    local resp
    read -p "[MSG] $prompt [Y/n] " resp
    case "$resp" in
        ""|[Yy]*) return 0 ;;
        *) return 1 ;;
    esac
}

# confirm_default_no <prompt>
#   Prompts the user with [y/N]. Returns 0 if yes, 1 if no (default).
confirm_default_no() {
    local prompt="$1"
    local resp
    read -p "[MSG] $prompt [y/N] " resp
    case "$resp" in
        [Yy]*) return 0 ;;
        *) return 1 ;;
    esac
}

# ============= Port Functions =============

# choose_port
#   Returns a random available port in 20000-65000.
choose_port() {
    local try=0
    while [ $try -lt 10 ]; do
        candidate=$(shuf -i 20000-65000 -n 1)
        if ! lsof -i :$candidate &> /dev/null; then
            echo $candidate
            return 0
        fi
        try=$((try+1))
    done
    return 1
}

# validate_port <port>
#   Validates port is a number in 1024-65535 range. Exits on error.
validate_port() {
    local port="$1"
    if ! printf '%s' "$port" | grep -qE '^[0-9]+$'; then
        print_error "Port must be a number."
        exit 1
    fi
    if [ "$port" -lt 1024 ] || [ "$port" -gt 65535 ]; then
        print_error "Port $port is out of allowed range (1024-65535)."
        exit 1
    fi
}

# check_port_available <port>
#   Exits if port is already in use.
check_port_available() {
    local port="$1"
    if lsof -i :$port &> /dev/null; then
        print_error "Port ${BLUE}$port${NC} is already in use."
        print_info "Please choose a different port using the ${YELLOW}-p${NC} option."
        exit 1
    fi
}

# ============= Condatainer Checks =============

# check_condatainer
#   Exits if condatainer is not in PATH.
check_condatainer() {
    if ! command -v condatainer &> /dev/null; then
        print_error "CondaTainer is not installed or not in PATH."
        print_info "Please go to https://github.com/Justype/condatainer to install CondaTainer."
        exit 1
    fi
}

# ============= Overlay Functions =============

# check_overlay_integrity <overlay_file>
#   Runs e2fsck to check and repair overlay. Exits on failure.
check_overlay_integrity() {
    local overlay="$1"
    
    # Check if in use before running e2fsck
    check_overlay_in_use "$overlay"

    print_info "Checking overlay image integrity..."
    e2fsck -p "$overlay" > /dev/null 2>&1
    # 0 = clean, 1 = errors fixed, >1 = errors unfixable
    if [ $? -gt 1 ]; then
        print_error "Overlay ${BLUE}$overlay${NC} is corrupted and could not be auto repaired."
        print_info "Please run ${YELLOW}e2fsck $overlay${NC} manually to repair."
        exit 1
    fi
}

# check_writable <file>
#   Returns 0 if file is writable and not locked exclusively.
check_writable() {
    [[ -w "$1" ]] && flock -xn 9 9<>"$1" 2>/dev/null
}

# check_readable <file>
#   Returns 0 if file is readable and not locked for writing.
check_readable() {
    [[ -r "$1" ]] && flock -sn 9 9<"$1" 2>/dev/null
}

# require_writable <file>
#   Exits if file is not writable or is currently in use.
#   In headless mode, offers to kill processes holding the file.
require_writable() {
    local file="$1"
    if ! check_writable "$file"; then
        local pid_list=$(lsof "$file" 2>/dev/null | awk 'NR>1 {print $2}' | sort -u)
        if [ -n "$pid_list" ]; then
            print_warn "${BLUE}$file${NC} is currently in use by the following processes: $pid_list"
            if confirm_default_no "Do you want to kill these processes and continue?"; then
                for pid in $pid_list; do
                    kill -9 $pid 2>/dev/null
                done
                print_info "Sent SIGKILL to processes using the file."
                sleep 1
                # Final check after kill
                check_writable "$file" || {
                    print_error "File still in use after kill attempt."
                    exit 1
                }
            else
                print_msg "Aborted."
                exit 1
            fi
        else
            print_error "Can't open ${BLUE}$file${NC} for writing, currently in use."
            exit 1
        fi
    fi
}

# check_overlay_in_use <overlay_file>
#   Checks if overlay is available for writing using flock.
check_overlay_in_use() {
    local overlay="$1"
    
    # Only .img files support locking
    [[ "$overlay" != *".img" ]] && return 0
    [ ! -f "$overlay" ] && return 0

    require_writable "$overlay"
}

# check_and_install_overlays <pkgs...>
#   Checks if required overlays exist and installs missing ones.
#   Each argument can be a single overlay name or a colon-separated list.
check_and_install_overlays() {
    print_info "Checking required overlays..."
    local missing=""
    for arg in "$@"; do
        # Split on colon to handle colon-separated overlay lists
        local IFS=':'
        local -a pkgs=($arg)
        for pkg in "${pkgs[@]}"; do
            [ -z "$pkg" ] && continue
            if echo "$pkg" | grep -qiE '\.(sqf|sqsh|squashfs|img)$'; then
                if [ ! -f "$pkg" ]; then
                    print_error "Overlay file ${BLUE}$pkg${NC} not found."
                    exit 1
                fi
            else
                if ! condatainer list -e "$pkg" > /dev/null 2>&1; then
                    missing+=" $pkg"
                fi
            fi
        done
    done
    if [ -n "$missing" ]; then
        print_info "Installing missing overlays:${missing}."
        condatainer create $missing
        if [ $? -ne 0 ]; then
            print_error "Failed to install required overlays."
            exit 1
        fi
    fi
}

# ============= Environment =============

# setup_scratch
#   Falls back SCRATCH to HOME if not set.
setup_scratch() {
    if [ -z "$SCRATCH" ]; then
        print_warn "SCRATCH environment variable is not set. Falling back to HOME directory."
        SCRATCH="$HOME"
    fi
}

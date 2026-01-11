#!/bin/bash
set -e

# Colors
GREEN=$'\033[0;32m'
BLUE=$'\033[0;34m'
CYAN=$'\033[0;36m'
YELLOW=$'\033[1;33m'
RED=$'\033[0;31m'
NC=$'\033[0m' # No Color

# URLs
URL_CONDATAINER="https://raw.githubusercontent.com/Justype/condatainer/main/bin/condatainer"
URL_MODGEN="https://raw.githubusercontent.com/Justype/condatainer/main/bin/modgen"

DEFAULT_BASE="${SCRATCH:-$HOME}/condatainer"

# Config Markers
MARKER_START="# >>> CONDATAINER >>>"
MARKER_END="# <<< CONDATAINER <<<"

# Detect Shell & Config File early
SHELL_NAME=$(basename "$SHELL")
if [ "$SHELL_NAME" = "zsh" ]; then RC_FILE="$HOME/.zshrc"; else RC_FILE="$HOME/.bashrc"; fi

echo -e "========================================"
echo -e "${BLUE}       CondaTainer Installer${NC}"
echo -e "========================================"

# ------------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------------

get_input() {
    local prompt_text="$1"
    local default_value="$2"
    local user_val=""

    if [ -c /dev/tty ]; then
        local p_str="${CYAN}${prompt_text}${NC} [${default_value}]: "
        read -e -p "$p_str" -r user_val < /dev/tty 2> /dev/tty
    fi
    echo "${user_val:-$default_value}"
}

# Generic confirmation handler
# Usage: confirm_default "Question?" "yes/no"
confirm_default() {
    local prompt_text="$1"
    local default="$2" # "yes" or "no"
    local response=""

    # Auto-Yes Logic: If -y is set, use the default value.
    if [ "${AUTO_YES:-false}" = "true" ]; then
        [ "$default" = "yes" ] && return 0 || return 1
    fi

    if [ -c /dev/tty ]; then
        if [ "$default" = "yes" ]; then
            printf "${CYAN}%s${NC} [Y/n]: " "$prompt_text" > /dev/tty
        else
            printf "${CYAN}%s${NC} [y/N]: " "$prompt_text" > /dev/tty
        fi
        read -r response < /dev/tty

        # Empty response = default
        if [ -z "$response" ]; then
            [ "$default" = "yes" ] && return 0 || return 1
        fi

        [[ "$response" =~ ^[Yy]$ ]]
    else
        # Non-interactive without -y: return default
        [ "$default" = "yes" ]
    fi
}

# Alias for Default=Yes (Prompt: [Y/n])
confirm_action() {
    confirm_default "$1" "yes"
}

# Alias for Default=No (Prompt: [y/N])
confirm_action_no() {
    confirm_default "$1" "no"
}

download_file() {
    local url="$1"
    local dest="$2"
    echo -e "${BLUE}[INFO]${NC} Downloading $(basename "$dest")..."
    if command -v curl >/dev/null 2>&1; then
        curl -sL "$url" -o "$dest"
    elif command -v wget >/dev/null 2>&1; then
        wget -qO "$dest" "$url"
    else
        echo -e "${RED}[ERROR]${NC} Neither curl nor wget found."
        exit 1
    fi
    chmod +x "$dest"
}

# Generic RC Block Updater
update_config_block() {
    local file="$1"
    local start="$2"
    local end="$3"
    local content="$4"
    local name="$5"
    local temp="${file}.tmp"

    if grep -Fq "$start" "$file"; then
        echo -e "${BLUE}[INFO]${NC} Updating $name in $(basename "$file")..."
        sed "/$start/,/$end/d" "$file" > "$temp"
        echo "$content" >> "$temp"
        mv "$temp" "$file"
        echo -e "${GREEN}[OK]${NC} Updated $name in ${BLUE}$file${NC}"
    else
        echo "" >> "$file"
        echo "$content" >> "$file"
        echo -e "${GREEN}[OK]${NC} Added $name to ${BLUE}$file${NC}"
    fi
}

# ------------------------------------------------------------------
# 2. Pre-Install Checks (Existing Installation)
# ------------------------------------------------------------------

EXISTING_DIR=""
EXISTING_ROOT=""
EXISTING_ACTION="NONE"

# Check for markers in the RC file
if [ -f "$RC_FILE" ] && grep -Fq "$MARKER_START" "$RC_FILE"; then
    EXISTING_PATH_LINE=$(sed -n "/$MARKER_START/,/$MARKER_END/p" "$RC_FILE" | grep "export PATH=" | head -n 1)
    if [ -n "$EXISTING_PATH_LINE" ]; then
        TEMP_PATH="${EXISTING_PATH_LINE#*\"}"
        EXISTING_DIR="${TEMP_PATH%%:*}"

        if [ -n "$EXISTING_DIR" ] && [ -d "$EXISTING_DIR" ]; then
            if [ "$(basename "$EXISTING_DIR")" == "bin" ]; then
                EXISTING_ROOT="$(dirname "$EXISTING_DIR")"
            else
                EXISTING_ROOT="$EXISTING_DIR"
            fi
        fi
    fi
fi

# ------------------------------------------------------------------
# 3. Configuration Phase
# ------------------------------------------------------------------

echo -e "\nConfiguration:"

# CLI options parsing
CLI_SPECIFIED=false
CLI_CONDATAINER=false
CLI_MODGEN=false
CLI_PATH=""
CLI_YES=false

show_usage() {
    cat <<'USAGE'
Usage: install.sh [options]
Options:
  -c, --condatainer    Install condaTainer
  -m, --modgen         Install modgen
  -a, --all            Install both
  -y, --yes            Assume yes for all prompts
  -p, --path PATH      Install base path (non-interactive)
  -h, --help           Show this help
USAGE
}

while [ $# -gt 0 ]; do
    case "$1" in
        -c|--condatainer) CLI_CONDATAINER=true; CLI_SPECIFIED=true; shift ;;
        -m|--modgen)      CLI_MODGEN=true; CLI_SPECIFIED=true; shift ;;
        -a|--all)         CLI_CONDATAINER=true; CLI_MODGEN=true; CLI_SPECIFIED=true; shift ;;
        -p|--path)        [ -n "$2" ] && { CLI_PATH="$2"; shift 2; } || { echo -e "${RED}[ERROR]${NC} --path requires an argument."; exit 1; } ;;
        -y|--yes)         CLI_YES=true; CLI_SPECIFIED=true; shift ;;
        -h|--help)        show_usage; exit 0 ;;
        *)                echo -e "${RED}[ERROR]${NC} Unknown option: $1"; show_usage; exit 1 ;;
    esac
done

# If auto-yes and no components specified, default to both
if [ "${CLI_YES}" = true ] && [ "${CLI_CONDATAINER}" = false ] && [ "${CLI_MODGEN}" = false ]; then
    CLI_CONDATAINER=true
    CLI_MODGEN=true
fi
AUTO_YES=${CLI_YES}

# A. Ask for Path
TARGET_FROM_EXISTING="false"
if [ -n "$CLI_PATH" ]; then
    INSTALL_BASE="$CLI_PATH"
elif [ -n "$EXISTING_ROOT" ]; then
    echo -e "${BLUE}[INFO]${NC} No --path given, using existing: ${BLUE}$EXISTING_ROOT${NC}"
    INSTALL_BASE="$EXISTING_ROOT"
    TARGET_FROM_EXISTING="true"
elif [ "${CLI_YES}" = true ]; then
    INSTALL_BASE="$DEFAULT_BASE"
else
    INSTALL_BASE=$(get_input "Install Path" "$DEFAULT_BASE")
fi

INSTALL_BASE="${INSTALL_BASE/#\~/$HOME}"

# Try to resolve path using realpath
if command -v realpath >/dev/null 2>&1; then
    if RESOLVED=$(realpath -m "$INSTALL_BASE" 2>/dev/null); then
        INSTALL_BASE="$RESOLVED"
    elif [ -e "$INSTALL_BASE" ] && RESOLVED=$(realpath "$INSTALL_BASE" 2>/dev/null); then
        INSTALL_BASE="$RESOLVED"
    fi
fi
if [[ "$INSTALL_BASE" != /* ]]; then INSTALL_BASE="$PWD/$INSTALL_BASE"; fi
INSTALL_BIN="$INSTALL_BASE/bin"

# Step 1: Target Check
if [ "$TARGET_FROM_EXISTING" != "true" ] && [ -d "$INSTALL_BASE" ]; then
        if ! confirm_action_no "Target '$INSTALL_BASE' exists. Continue installation?"; then
            echo "Installation aborted."; exit 1
    fi
fi

# Step 2: Existing Install Action
if [ -n "$EXISTING_ROOT" ]; then
    if [ "$EXISTING_ROOT" = "$INSTALL_BASE" ]; then
        echo -e "${BLUE}[INFO]${NC} Updating configuration only."
        EXISTING_ACTION="NONE"
    else
        echo -e "${YELLOW}[WARNING]${NC} Found existing installation at ${BLUE}$EXISTING_ROOT${NC}"
        # Default is YES (Move)
        if confirm_default "Move content from $EXISTING_ROOT to $INSTALL_BASE?" "yes"; then
            EXISTING_ACTION="MOVE"
        elif confirm_action_no "Remove the existing folder?"; then
            EXISTING_ACTION="REMOVE"
        else
            echo "Aborted."; exit 1
        fi
    fi
fi

# Step 3: Components
DO_INSTALL_CONDATAINER=false; DO_INSTALL_MODGEN=false
DEFAULT_CONDA="no"; DEFAULT_MODGEN="no"

if [ -n "$EXISTING_ROOT" ] && [ -d "$EXISTING_ROOT" ]; then
    [ -f "$EXISTING_ROOT/bin/condatainer" ] && DEFAULT_CONDA="yes"
    [ -f "$EXISTING_ROOT/bin/modgen" ] && DEFAULT_MODGEN="yes"
else
    DEFAULT_CONDA="yes"; DEFAULT_MODGEN="no"
fi

if [ "$CLI_SPECIFIED" = true ]; then
    [ "$CLI_CONDATAINER" = true ] && DO_INSTALL_CONDATAINER=true
    [ "$CLI_MODGEN" = true ] && DO_INSTALL_MODGEN=true
else
    confirm_default "Install 'condatainer'?" "$DEFAULT_CONDA" && DO_INSTALL_CONDATAINER=true
    confirm_default "Install 'modgen'?" "$DEFAULT_MODGEN" && DO_INSTALL_MODGEN=true
fi

# ------------------------------------------------------------------
# 4. Final Verification
# ------------------------------------------------------------------

echo -e "\n--------- Installation Summary ---------"
echo -e "Directory   : ${BLUE}$INSTALL_BIN${NC}"
echo -e "Condatainer : $( [ "$DO_INSTALL_CONDATAINER" = true ] && echo "${GREEN}Yes${NC}" || echo "${RED}No${NC}" )"
echo -e "Modgen      : $( [ "$DO_INSTALL_MODGEN" = true ] && echo "${GREEN}Yes${NC}" || echo "${RED}No${NC}" )"
if [ "$EXISTING_ACTION" == "MOVE" ]; then
    echo -e "Old Install : ${CYAN}Move to new location${NC}"
elif [ "$EXISTING_ACTION" == "REMOVE" ]; then
    echo -e "Old Install : ${RED}Remove${NC}"
fi
echo -e "----------------------------------------"

if ! confirm_action "Proceed?"; then echo "Aborted."; exit 0; fi

# ------------------------------------------------------------------
# 5. Execution Phase
# ------------------------------------------------------------------

echo -e "\nStarting installation..."

# Handle Existing
if [ -n "$EXISTING_ROOT" ] && [ "$EXISTING_ROOT" != "$INSTALL_BASE" ] && [ "$EXISTING_ACTION" != "NONE" ]; then
    # Safety Check against system paths
    if [[ "$EXISTING_DIR" == "/usr/bin" || "$EXISTING_DIR" == "/bin" ]]; then
        echo -e "${YELLOW}[WARNING]${NC} Skipping Move/Remove: System directory detected."
    elif [ "$EXISTING_ACTION" == "MOVE" ]; then
        echo -e "${BLUE}[INFO]${NC} Moving content..."
        mkdir -p "$INSTALL_BASE"
        mv "$EXISTING_ROOT"/* "$INSTALL_BASE/" 2>/dev/null || true
        rmdir "$EXISTING_ROOT" 2>/dev/null || true
    elif [ "$EXISTING_ACTION" == "REMOVE" ]; then
        echo -e "${BLUE}[INFO]${NC} Removing old directory..."
        rm -rf "$EXISTING_ROOT"
    fi
fi

# Prepare Directory
if [ -d "$INSTALL_BASE" ] && [ "$EXISTING_ACTION" != "MOVE" ] && [ "$INSTALL_BASE" != "$EXISTING_ROOT" ]; then
    # Default is NO (Do not remove)
    if confirm_action_no "Remove existing folder?"; then
         rm -rf "$INSTALL_BASE"
    fi
fi
mkdir -p "$INSTALL_BIN"

# Download
if [ "$DO_INSTALL_CONDATAINER" = true ]; then
    download_file "$URL_CONDATAINER" "$INSTALL_BIN/condatainer"
fi

if [ "$DO_INSTALL_MODGEN" = true ]; then
    download_file "$URL_MODGEN" "$INSTALL_BIN/modgen"
    mkdir -p "$INSTALL_BASE/apps-modules" "$INSTALL_BASE/ref-modules"
fi

# Update RC
PATH_BLOCK="$MARKER_START
if [[ \":\$PATH:\" != *\":$INSTALL_BIN:\"* ]]; then
    export PATH=\"$INSTALL_BIN:\$PATH\"
fi
if command -v condatainer >/dev/null 2>&1; then
    source <(condatainer completion)
fi
$MARKER_END"

update_config_block "$RC_FILE" "$MARKER_START" "$MARKER_END" "$PATH_BLOCK" "PATH"

if [ "$DO_INSTALL_MODGEN" = true ]; then
    MOD_START="# >>> MODGEN MODULES >>>"
    MOD_END="# <<< MODGEN MODULES <<<"
    DETECT_CMD='command -v tclsh >/dev/null 2>&1'
    [ -n "$LMOD_CMD" ] && DETECT_CMD='[ -f "$LMOD_CMD" ]'

    MOD_BLOCK="$MOD_START
# Use ModGen modulefiles (Only if 'module' command is available)
if $DETECT_CMD && [ -d \"$INSTALL_BASE\" ]; then
    module use \"$INSTALL_BASE/apps-modules\"
    module use \"$INSTALL_BASE/ref-modules\"
fi
$MOD_END"

    update_config_block "$RC_FILE" "$MOD_START" "$MOD_END" "$MOD_BLOCK" "ModGen Config"
fi

echo -e "----------------------------------------"
echo -e "${GREEN}Success!${NC}"
echo "Run this to apply changes:"
echo -e "  source ${BLUE}$RC_FILE${NC}"
echo -e "========================================"

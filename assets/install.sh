#!/bin/bash
set -e

# Colors
GREEN=$'\033[0;32m'
BLUE=$'\033[0;34m'
CYAN=$'\033[0;36m'
YELLOW=$'\033[1;33m'
RED=$'\033[0;31m'
NC=$'\033[0m' # No Color

# Detect OS and Architecture
detect_platform() {
    local os arch binary_name

    os=$(uname -s | tr '[:upper:]' '[:lower:]')
    arch=$(uname -m)

    # Map architecture names
    case "$arch" in
        x86_64|amd64) arch="x86_64" ;;
        aarch64|arm64) arch="arm64" ;;
        *) arch="$arch" ;;
    esac

    # Build binary name and check support
    os_arch="${os}_${arch}"
    case "${os_arch}" in
        linux_x86_64) ;;
        *)
            echo -e "${RED}[ERROR]${NC} Unsupported platform: ${os} ${arch}"
            echo -e "Currently only ${BLUE}linux x86_64${NC} is supported."
            exit 1
            ;;
    esac

    echo "$os_arch"
}

BINARY_NAME=condatainer_$(detect_platform)
URL_CONDATAINER="https://github.com/Justype/condatainer/releases/latest/download/${BINARY_NAME}"

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
confirm_default() {
    local prompt_text="$1"
    local default="$2" # "yes" or "no"
    local response=""

    # Auto-Yes Logic
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

        if [ -z "$response" ]; then
            [ "$default" = "yes" ] && return 0 || return 1
        fi

        [[ "$response" =~ ^[Yy]$ ]]
    else
        [ "$default" = "yes" ]
    fi
}

confirm_action() { confirm_default "$1" "yes"; }
confirm_action_no() { confirm_default "$1" "no"; }

download_file() {
    local url="$1"
    local dest="$2"
    echo -e "${BLUE}[INFO]${NC} Downloading $(basename "$dest")..."
    if command -v wget >/dev/null 2>&1; then
        wget -qO "$dest" "$url"
    elif command -v curl >/dev/null 2>&1; then
        curl -sL "$url" -o "$dest"
    else
        echo -e "${RED}[ERROR]${NC} Neither curl nor wget found."
        exit 1
    fi
    chmod +x "$dest"
}

update_config_block() {
    local file="$1"
    local start="$2"
    local end="$3"
    local content="$4"
    local name="$5"
    local temp="${file}.tmp"

    if [ -f "$file" ] && grep -Fq "$start" "$file"; then
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
# 2. Pre-Install Checks
# ------------------------------------------------------------------

EXISTING_DIR=""
EXISTING_ROOT=""

# Check for CONDATAINER markers in RC file
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

# If `condatainer` is already available in PATH but RC file lacks our marker,
# warn the user and ask whether to continue (default: no).
if command -v condatainer >/dev/null 2>&1 && { [ ! -f "$RC_FILE" ] || ! grep -Fq "$MARKER_START" "$RC_FILE"; }; then
    echo -e "${YELLOW}[WARNING]${NC} 'condatainer' is available in your PATH but ${RC_FILE} does not contain the CONDATAINER config block."
    echo -e "This likely means an existing installation was added to PATH outside the installer."
    if ! confirm_action_no "Continue with the additional installation?"; then
        echo "Installation aborted."; exit 1
    fi
fi

# ------------------------------------------------------------------
# 3. Configuration Phase
# ------------------------------------------------------------------

echo -e "\nConfiguration:"

# CLI options parsing
CLI_PATH=""
CLI_YES=false

show_usage() {
    cat <<'USAGE'
Usage: install_condatainer.sh [options]
Options:
  -p, --path PATH      Install base path (non-interactive)
  -y, --yes            Assume yes for all prompts
  -h, --help           Show this help
USAGE
}

while [ $# -gt 0 ]; do
    case "$1" in
        -p|--path)  [ -n "$2" ] && { CLI_PATH="$2"; shift 2; } || { echo -e "${RED}[ERROR]${NC} --path requires an argument."; exit 1; } ;;
        -y|--yes)   CLI_YES=true; shift ;;
        -h|--help)  show_usage; exit 0 ;;
        *)          echo -e "${RED}[ERROR]${NC} Unknown option: $1"; show_usage; exit 1 ;;
    esac
done

AUTO_YES=${CLI_YES}

# Determine Install Path
TARGET_FROM_EXISTING="false"
if [ -n "$CLI_PATH" ]; then
    INSTALL_BASE="$CLI_PATH"
elif [ -n "$EXISTING_ROOT" ]; then
    echo -e "${BLUE}[INFO]${NC} Found existing installation: ${BLUE}$EXISTING_ROOT${NC}"
    if [ "${CLI_YES}" = true ]; then
        INSTALL_BASE="$EXISTING_ROOT"
    else
        INSTALL_BASE=$(get_input "Install Path" "$EXISTING_ROOT")
    fi
    TARGET_FROM_EXISTING="true"
elif [ "${CLI_YES}" = true ]; then
    INSTALL_BASE="$DEFAULT_BASE"
else
    INSTALL_BASE=$(get_input "Install Path" "$DEFAULT_BASE")
fi

INSTALL_BASE="${INSTALL_BASE/#\~/$HOME}"

# Resolve path
if command -v realpath >/dev/null 2>&1; then
    if RESOLVED=$(realpath -m "$INSTALL_BASE" 2>/dev/null); then
        INSTALL_BASE="$RESOLVED"
    elif [ -e "$INSTALL_BASE" ] && RESOLVED=$(realpath "$INSTALL_BASE" 2>/dev/null); then
        INSTALL_BASE="$RESOLVED"
    fi
fi
if [[ "$INSTALL_BASE" != /* ]]; then INSTALL_BASE="$PWD/$INSTALL_BASE"; fi
INSTALL_BIN="$INSTALL_BASE/bin"

# Target exists check
if [ "$TARGET_FROM_EXISTING" != "true" ] && [ -d "$INSTALL_BASE" ]; then
    if ! confirm_action_no "Target '$INSTALL_BASE' exists. Continue installation?"; then
        echo "Installation aborted."; exit 1
    fi
fi

# ------------------------------------------------------------------
# 4. Final Verification
# ------------------------------------------------------------------

echo -e "\n--------- Installation Summary ---------"
echo -e "Directory   : ${BLUE}$INSTALL_BIN${NC}"
echo -e "Condatainer : ${GREEN}Yes${NC}"
echo -e "----------------------------------------"

if ! confirm_action "Proceed?"; then echo "Aborted."; exit 0; fi

# ------------------------------------------------------------------
# 5. Execution Phase
# ------------------------------------------------------------------

echo -e "\nStarting installation..."

# Prepare Directory
mkdir -p "$INSTALL_BIN"

# Download condatainer
download_file "$URL_CONDATAINER" "$INSTALL_BIN/condatainer"

# Update RC with CONDATAINER block (skip if installing to common PATH directories)
SKIP_PATH_BLOCK=false
case "$INSTALL_BIN" in
    "$HOME/bin"|"$HOME/.local/bin")
        SKIP_PATH_BLOCK=true
        echo -e "${BLUE}[INFO]${NC} Skipping PATH configuration (${INSTALL_BIN} is typically already in PATH)"
        ;;
esac

if [ "$SKIP_PATH_BLOCK" = false ]; then
    PATH_BLOCK="$MARKER_START
if [[ \":\$PATH:\" != *\":$INSTALL_BIN:\"* ]]; then
    export PATH=\"$INSTALL_BIN:\$PATH\"
fi
$MARKER_END"
    update_config_block "$RC_FILE" "$MARKER_START" "$MARKER_END" "$PATH_BLOCK" "CONDATAINER PATH"
fi

echo -e "----------------------------------------"
echo -e "${GREEN}Success!${NC}"
echo "Run this to apply changes:"
echo -e "  source ${BLUE}$RC_FILE${NC}"
echo -e "========================================"

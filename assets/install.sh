#!/bin/bash
set -e

# Colors (ANSI-C quoting used so 'read -p' sees actual escape codes)
GREEN=$'\033[0;32m'
BLUE=$'\033[0;34m'
CYAN=$'\033[0;36m'
RED=$'\033[0;31m'
NC=$'\033[0m' # No Color

# Hardcoded Source URLs
URL_CONDATAINER="https://raw.githubusercontent.com/Justype/condatainer/main/bin/condatainer"
URL_MODGEN="https://raw.githubusercontent.com/Justype/condatainer/main/bin/modgen"

DEFAULT_BASE="${SCRATCH:-$HOME}/condatainer"

# Config Markers
MARKER_START="# >>> CONDATAINER >>>"
MARKER_END="# <<< CONDATAINER <<<"

# Detect Shell & Config File early (Used for detection and installation)
SHELL_NAME=$(basename "$SHELL")
if [ "$SHELL_NAME" = "zsh" ]; then
    RC_FILE="$HOME/.zshrc"
else
    RC_FILE="$HOME/.bashrc"
fi

echo -e "========================================"
echo -e "${BLUE}       CondaTainer Installer${NC}"
echo -e "========================================"

# ------------------------------------------------------------------
# 1. Interactive Input Helpers
# ------------------------------------------------------------------

# Function to get string input (for Path)
get_input() {
    local prompt_text="$1"
    local default_value="$2"
    local user_val=""

    if [ -c /dev/tty ]; then
        local p_str="${CYAN}${prompt_text}${NC} [${default_value}]: "
        read -e -p "$p_str" -r user_val < /dev/tty 2> /dev/tty
    fi

    if [ -z "$user_val" ]; then
        echo "$default_value"
    else
        echo "${user_val%/}"
    fi
}

# Function for Yes/No prompts (Default: Yes)
confirm_action() {
    local prompt_text="$1"
    local response=""

    if [ -c /dev/tty ]; then
        printf "${CYAN}%s${NC} [Y/n]: " "$prompt_text" > /dev/tty
        read -r response < /dev/tty
        
        if [[ -z "$response" || "$response" =~ ^[Yy]$ ]]; then
            return 0 # True
        else
            return 1 # False
        fi
    else
        return 0 # Default to yes in non-interactive
    fi
}

# Function for Yes/No prompts (Default: No)
confirm_action_no() {
    local prompt_text="$1"
    local response=""

    if [ -c /dev/tty ]; then
        printf "${CYAN}%s${NC} [y/N]: " "$prompt_text" > /dev/tty
        read -r response < /dev/tty
        
        if [[ "$response" =~ ^[Yy]$ ]]; then
            return 0 # True
        else
            return 1 # False
        fi
    else
        return 1 # Default to no in non-interactive
    fi
}

# ------------------------------------------------------------------
# 2. Pre-Install Checks (Existing Installation)
# ------------------------------------------------------------------

EXISTING_DIR=""
EXISTING_ROOT=""
EXISTING_ACTION="NONE"

# Check for markers in the RC file instead of using command -v
if [ -f "$RC_FILE" ] && grep -Fq "$MARKER_START" "$RC_FILE"; then
    # Extract the block, find the export line, and parse the path
    # Format expected: export PATH="/path/to/bin:$PATH"
    EXISTING_PATH_LINE=$(sed -n "/$MARKER_START/,/$MARKER_END/p" "$RC_FILE" | grep "export PATH=" | head -n 1)
    
    if [ -n "$EXISTING_PATH_LINE" ]; then
        # Strip everything up to the first quote
        TEMP_PATH="${EXISTING_PATH_LINE#*\"}"
        # Strip everything from the colon onwards (:$PATH")
        EXISTING_DIR="${TEMP_PATH%%:*}"

        if [ -n "$EXISTING_DIR" ] && [ -d "$EXISTING_DIR" ]; then
            # Determine the Installation Root (Parent of bin)
            if [ "$(basename "$EXISTING_DIR")" == "bin" ]; then
                EXISTING_ROOT="$(dirname "$EXISTING_DIR")"
            else
                EXISTING_ROOT="$EXISTING_DIR"
            fi

            echo -e "\n${RED}Warning:${NC} Found existing configuration in ${BLUE}$(basename "$RC_FILE")${NC}"
            echo -e "Installation Root: ${BLUE}$EXISTING_ROOT${NC}"
            
            if ! confirm_action_no "Continue installation?"; then
                echo "Installation aborted."
                exit 0
            fi

            if confirm_action "Move content to the new folder?"; then
                EXISTING_ACTION="MOVE"
            # Default to No for deletion to prevent accidents
            elif confirm_action_no "Remove the existing folder?"; then
                EXISTING_ACTION="REMOVE"
            fi
        elif [ -n "$EXISTING_DIR" ]; then
            # Config exists but directory is missing
            echo -e "\n${CYAN}Notice:${NC} Found configuration in $(basename "$RC_FILE") pointing to missing directory."
            echo "Configuration will be updated automatically."
        fi
    fi
fi

# ------------------------------------------------------------------
# 3. Configuration Phase
# ------------------------------------------------------------------

echo -e "\nConfiguration:"

# A. Ask for Path
INSTALL_BASE=$(get_input "Install Path" "$DEFAULT_BASE")

# Ensure path is absolute (Critical for correct .bashrc/.zshrc exports)
# 1. Expand tilde (~) manually if present
INSTALL_BASE="${INSTALL_BASE/#\~/$HOME}"
# 2. If it doesn't start with /, prepend current working directory
if [[ "$INSTALL_BASE" != /* ]]; then
    INSTALL_BASE="$PWD/$INSTALL_BASE"
fi

INSTALL_BIN="$INSTALL_BASE/bin"

# B. Ask for Components
DO_INSTALL_CONDATAINER=false
DO_INSTALL_MODGEN=false

if confirm_action "Install 'condatainer'?"; then
    DO_INSTALL_CONDATAINER=true
fi

if confirm_action "Install 'modgen'?"; then
    DO_INSTALL_MODGEN=true
fi

# ------------------------------------------------------------------
# 4. Final Verification
# ------------------------------------------------------------------

echo -e "\n-------- Installation Summary ---------"
echo -e "Directory:   ${BLUE}$INSTALL_BIN${NC}"
echo -e "Condatainer: $( [ "$DO_INSTALL_CONDATAINER" = true ] && echo "${GREEN}Yes${NC}" || echo "${RED}No${NC}" )"
echo -e "Modgen:      $( [ "$DO_INSTALL_MODGEN" = true ] && echo "${GREEN}Yes${NC}" || echo "${RED}No${NC}" )"
if [ "$EXISTING_ACTION" == "MOVE" ]; then
    echo -e "Old Install: ${CYAN}Move to new dir${NC}"
elif [ "$EXISTING_ACTION" == "REMOVE" ]; then
    echo -e "Old Install: ${RED}Remove ($EXISTING_ROOT)${NC}"
fi
echo -e "----------------------------------------"

if ! confirm_action "Proceed?"; then
    echo "Installation aborted by user."
    exit 0
fi

# ------------------------------------------------------------------
# 5. Execution Phase
# ------------------------------------------------------------------

echo -e "\nStarting installation..."

# --- ACTION: Handle Existing Installation ---
if [ -n "$EXISTING_ROOT" ] && [ "$EXISTING_ROOT" != "$INSTALL_BASE" ]; then
    # Safety Check: Don't remove system paths based on the bin dir check
    if [[ "$EXISTING_DIR" == "/usr/bin" || "$EXISTING_DIR" == "/bin" || "$EXISTING_DIR" == "/usr/local/bin" ]]; then
        echo -e "${RED}Skipping Move/Remove:${NC} Existing path is a system directory ($EXISTING_DIR). Manual cleanup recommended."
        EXISTING_ACTION="NONE"
    fi

    if [ "$EXISTING_ACTION" == "MOVE" ]; then
        echo "Moving content from $EXISTING_ROOT to $INSTALL_BASE..."
        mkdir -p "$INSTALL_BASE"
        # Move content of root (including bin/ and modules/) to new base
        mv "$EXISTING_ROOT"/* "$INSTALL_BASE/" 2>/dev/null || true
        # Remove empty old directory
        rmdir "$EXISTING_ROOT" 2>/dev/null || true
        echo "Move complete."
    elif [ "$EXISTING_ACTION" == "REMOVE" ]; then
        echo "Removing old directory $EXISTING_ROOT..."
        rm -rf "$EXISTING_ROOT"
    fi
fi

# --- CHECK: Existing Directory (Target) ---
# Check the ROOT base, not just the bin
if [ -d "$INSTALL_BASE" ] && [ "$EXISTING_ACTION" != "MOVE" ]; then
    # Only warn if we haven't already dealt with it (e.g. if Old Root != New Base)
    if [ "$INSTALL_BASE" != "$EXISTING_ROOT" ]; then
        echo -e "${RED}Warning:${NC} Installation directory '$INSTALL_BASE' exists."
        # Default to No for safety
        if confirm_action_no "Remove existing folder?"; then
            echo "Removing old directory..."
            rm -rf "$INSTALL_BASE"
            mkdir -p "$INSTALL_BIN"
        else
            echo "Keeping existing directory. Files will be overwritten/merged."
            mkdir -p "$INSTALL_BIN"
        fi
    else
        mkdir -p "$INSTALL_BIN"
    fi
else
    mkdir -p "$INSTALL_BIN"
    # Ensure base is created if bin was created (mkdir -p does this, but being explicit about intent)
fi

download_file() {
    local url="$1"
    local dest="$2"
    
    echo "Downloading $(basename "$dest")..."
    if command -v curl >/dev/null 2>&1; then
        curl -sL "$url" -o "$dest"
    elif command -v wget >/dev/null 2>&1; then
        wget -qO "$dest" "$url"
    else
        echo "Error: Neither curl nor wget found."
        exit 1
    fi
    chmod +x "$dest"
}

# Download selected components
if [ "$DO_INSTALL_CONDATAINER" = true ]; then
    download_file "$URL_CONDATAINER" "$INSTALL_BIN/condatainer"
fi

if [ "$DO_INSTALL_MODGEN" = true ]; then
    download_file "$URL_MODGEN" "$INSTALL_BIN/modgen"
    
    # Create ModGen folders
    echo "Creating ModGen module directories..."
    mkdir -p "$INSTALL_BASE/apps-modules"
    mkdir -p "$INSTALL_BASE/ref-modules"
fi

# --- CHECK: Shell Config & PATH ---
# Note: RC_FILE is already determined at the top of the script

update_shell_rc() {
    local rc_file="$1"
    local bin_path="$2"
    local temp_rc="${rc_file}.tmp"

    # Define the block content
    local block_content="$MARKER_START
export PATH=\"$bin_path:\$PATH\"
$MARKER_END"

    # Check for markers
    if grep -Fq "$MARKER_START" "$rc_file"; then
        echo -e "${CYAN}Notice:${NC} Found existing CondaTainer block in $(basename "$rc_file"). Updating..."
        
        # 1. Delete old block (from START to END)
        # 2. Append new block
        # We use a temp file for portability (sed -i varies by OS)
        sed "/$MARKER_START/,/$MARKER_END/d" "$rc_file" > "$temp_rc"
        echo "$block_content" >> "$temp_rc"
        mv "$temp_rc" "$rc_file"
        
        echo -e "Updated PATH in ${BLUE}$rc_file${NC}"
    else
        # Append to end
        echo "" >> "$rc_file"
        echo "$block_content" >> "$rc_file"
        echo -e "Added PATH to ${BLUE}$rc_file${NC}"
    fi
}

update_modgen_rc() {
    local rc_file="$1"
    local base_dir="$2"
    local temp_rc="${rc_file}.tmp"
    
    local marker_start="# >>> MODGEN MODULES >>>"
    local marker_end="# <<< MODGEN MODULES <<<"

    # Detect module system for the config block condition
    local detect_cmd='command -v tclsh >/dev/null 2>&1'
    # If LMOD_CMD is set in the environment, assume Lmod
    if [ -n "$LMOD_CMD" ]; then
        detect_cmd='[ -f "$LMOD_CMD" ]'
    fi

    local block_content="$marker_start
# Use ModGen modulefiles (Only if 'module' command is available)
if $detect_cmd && [ -d \"$base_dir\" ]; then
    module use \"$base_dir/apps-modules\"
    module use \"$base_dir/ref-modules\"
fi
$marker_end"

    if grep -Fq "$marker_start" "$rc_file"; then
        echo -e "${CYAN}Notice:${NC} Found existing ModGen block in $(basename "$rc_file"). Updating..."
        sed "/$marker_start/,/$marker_end/d" "$rc_file" > "$temp_rc"
        echo "$block_content" >> "$temp_rc"
        mv "$temp_rc" "$rc_file"
        echo -e "Updated Module config in ${BLUE}$rc_file${NC}"
    else
        echo "" >> "$rc_file"
        echo "$block_content" >> "$rc_file"
        echo -e "Added Module config to ${BLUE}$rc_file${NC}"
    fi
}

update_shell_rc "$RC_FILE" "$INSTALL_BIN"

if [ "$DO_INSTALL_MODGEN" = true ]; then
    echo -e "\nConfiguring ModGen..."
    update_modgen_rc "$RC_FILE" "$INSTALL_BASE"
fi

echo -e "========================================"
echo -e "${GREEN}Success!${NC}"
echo "Run this to apply changes:"
echo -e "  source ${BLUE}$RC_FILE${NC}"
echo -e "========================================"

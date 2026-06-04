#!/usr/bin/env bash
# Fetch Material Symbols SVG paths and write web/icons.js
# Usage: ./fetch-icons.sh [weight]
#   weight: optional, e.g. 100 200 300 400 500 600 700 (default: 400)
#   Weight 400 has no suffix; others use _wght<N>_24px.svg

set -euo pipefail

WEIGHT="${1:-400}"
BASE="https://raw.githubusercontent.com/google/material-design-icons/master/symbols/web"
OUT="$(dirname "$0")/icons.js"

# Outline (fill=0) icons — fetched from materialsymbolsoutlined
ICONS=(
  # navigation
  chevron_left chevron_right close expand_circle_down expand_circle_up
  # actions
  add_2 content_copy delete download link_2 play_arrow play_circle refresh replay star
  # status / feedback
  check error info progress_activity question_mark warning
  # files
  draft files folder
  # app / ui
  computer dark_mode hexagon history_2 light_mode settings
)
mapfile -t ICONS < <(printf '%s\n' "${ICONS[@]}" | sort)

# Filled (fill=1) icons — fetched from materialsymbolsoutlined with fill1 suffix
FILL1_ICONS=(
  star
)

if [[ "$WEIGHT" == "400" ]]; then
  SUFFIX="_24px.svg"
  FILL1_SUFFIX="_fill1_24px.svg"
else
  SUFFIX="_wght${WEIGHT}_24px.svg"
  FILL1_SUFFIX="_fill1_wght${WEIGHT}_24px.svg"
fi

fetch_icon() {
  local out_name="$1" url="$2"
  local svg viewbox path
  svg=$(curl -sf "$url") || { echo "WARN: failed to fetch $out_name" >&2; return 1; }
  viewbox=$(echo "$svg" | grep -oP 'viewBox="[^"]+"')
  path=$(echo "$svg" | grep -oP '<path d="[^"]+"')
  printf '<symbol id="i-%-14s %s>%s/></symbol>\n' "${out_name}\"" "$viewbox" "$path"
}

{
  echo "// AUTO-GENERATED — do not edit by hand. Run web/fetch-icons.sh to update."
  echo "document.body.insertAdjacentHTML('afterbegin', \`<svg class=\"icon-sprite\" aria-hidden=\"true\">"
  for name in "${ICONS[@]}"; do
    url="${BASE}/${name}/materialsymbolsoutlined/${name}${SUFFIX}"
    fetch_icon "$name" "$url" || true
  done
  for name in "${FILL1_ICONS[@]}"; do
    url="${BASE}/${name}/materialsymbolsoutlined/${name}${FILL1_SUFFIX}"
    fetch_icon "fill-${name}" "$url" || true
  done
  echo "</svg>\`);"
} > "$OUT"

echo "Written to $OUT"

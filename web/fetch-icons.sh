#!/usr/bin/env bash
# Fetch Material Symbols SVG paths and write web/icons.js
# Usage: ./fetch-icons.sh [weight]
#   weight: optional, e.g. 100 200 300 400 500 600 700 (default: 400)
#   Weight 400 has no suffix; others use _wght<N>_24px.svg

set -euo pipefail

WEIGHT="${1:-400}"
BASE="https://raw.githubusercontent.com/google/material-design-icons/master/symbols/web"
OUT="$(dirname "$0")/icons.js"

ICONS=(
  # navigation
  chevron_left chevron_right close expand_circle_down expand_circle_up
  # actions
  add_2 content_copy delete play_arrow play_circle refresh replay
  # status / feedback
  check error info progress_activity question_mark warning
  # files
  draft files folder
  # app / ui
  computer dark_mode hexagon history_2 light_mode settings
)
mapfile -t ICONS < <(printf '%s\n' "${ICONS[@]}" | sort)

if [[ "$WEIGHT" == "400" ]]; then
  SUFFIX="_24px.svg"
else
  SUFFIX="_wght${WEIGHT}_24px.svg"
fi

{
  echo "document.body.insertAdjacentHTML('afterbegin', \`<svg class=\"icon-sprite\" aria-hidden=\"true\">"
  for name in "${ICONS[@]}"; do
    url="${BASE}/${name}/materialsymbolsoutlined/${name}${SUFFIX}"
    svg=$(curl -sf "$url") || { echo "WARN: failed to fetch $name" >&2; continue; }
    viewbox=$(echo "$svg" | grep -oP 'viewBox="[^"]+"')
    path=$(echo "$svg" | grep -oP '<path d="[^"]+"')
    printf '<symbol id="i-%-14s %s>%s/></symbol>\n' "${name}\"" "$viewbox" "$path"
  done
  echo "</svg>\`);"
} > "$OUT"

echo "Written to $OUT"

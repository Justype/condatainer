# Web UI

Single-page application served by the CondaTainer helper server.

## Files

| File | Purpose |
|------|---------|
| `index.html` | Main HTML shell, layout |
| `icons.js` | SVG symbol sprite — injected into DOM before other scripts |
| `style.css` | All styles |
| `app.js` | Bootstrap, theme, sidebar, shared utilities (`iconSvg`, `escHtml`, etc.) |
| `jobs.js` | Jobs section, detail panel, walltime countdowns |
| `helpers.js` | Start section — helper list, form, module picker |
| `overlays.js` | Overlay section — list, overlay picker modal |
| `history.js` | History section — table, rerun, delete |
| `files.js` | Files section + file/dir picker modal |
| `htmx.min.js` | htmx (vendored, unused — kept for potential future use) |
| `assets.go` | Go embed — serves all files under `/static/` |

## Icons

Icons use an inline SVG symbol sprite (`#i-*` IDs) in `index.html`.
All icons are from [Google Material Symbols](https://fonts.google.com/icons) (Outlined, 24px),
embedded as fill-based paths with `viewBox="0 -960 960 960"`.

The `.icon` CSS class renders them: `fill:currentColor; width/height:16px`.

To fetch a new icon: `https://raw.githubusercontent.com/google/material-design-icons/master/symbols/web/{name}/materialsymbolsoutlined/{name}_24px.svg`

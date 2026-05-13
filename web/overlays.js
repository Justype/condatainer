/* ── Overlays section ────────────────────── */
async function loadOverlays() {
  try {
    const r = await fetch('/api/overlays');
    allOverlays = (await r.json()) || [];
  } catch { allOverlays = []; }
  renderOverlays();
}

function renderOverlays(filter) {
  const q    = (filter || '').toLowerCase();
  const sqfs = allOverlays.filter(o => o.type === 'sqf');
  const rows = q ? sqfs.filter(o => o.name.toLowerCase().includes(q)) : sqfs;
  const tbody = gid('ov-tbody');
  if (!rows.length) {
    tbody.innerHTML = '<tr><td colspan="5" style="padding:20px;text-align:center;color:var(--muted);">' +
      (q ? 'No matches.' : 'No module overlays found.') + '</td></tr>';
    return;
  }
  tbody.innerHTML = rows.map(o =>
    '<tr>' +
      '<td style="font-weight:600;">' + escHtml(o.name) + '</td>' +
      '<td class="mono" style="color:var(--muted);">' + escHtml(o.type || '—') + '</td>' +
      '<td class="mono" style="color:var(--muted);">' + fmtSize(o.size) + '</td>' +
      '<td class="mono path-cell" style="color:var(--muted);">' + pathTailHtml(o.path) + '</td>' +
      '<td><div style="display:flex;gap:6px;">' +
        '<button class="btn btn-sm" onclick="pickOverlay(\'' + escHtml(o.path) + '\')">Use</button>' +
        '<button class="btn btn-sm btn-ghost" onclick="copyStr(\'' + escHtml(o.path) + '\',this)" title="Copy path">' + iconSvg('content_copy') + '</button>' +
      '</div></td>' +
    '</tr>'
  ).join('');
}

function pickOverlay(path) {
  const entry = allOverlays.find(o => o.path === path);
  if (entry && !selectedModules.find(m => m.path === path)) {
    selectedModules.push(entry);
    renderModuleChips();
  }
  navigate('start');
}

/* ── Overlay Picker Modal ────────────────── */
function openOverlayPicker(type) {
  opTargetType = type;
  if (type !== 'module') return;
  gid('op-app-search').value = '';
  filterAppOverlays('');
  gid('op-modal').classList.add('open');
}

function filterAppOverlays(q) {
  const sqfOvs = allOverlays.filter(o => o.type === 'sqf');
  const lower  = q.toLowerCase();
  const hits   = lower ? sqfOvs.filter(o => o.name.toLowerCase().includes(lower)) : sqfOvs;
  gid('op-app-list').innerHTML = hits.length
    ? hits.map(o =>
        '<div class="modal-row" onclick="selectOverlay(\'' + escHtml(o.path) + '\')">' +
          '<span class="modal-row-icon" style="color:var(--muted);">' + iconSvg('hexagon') + '</span>' +
          '<span class="modal-row-name" style="font-weight:600;">' + escHtml(o.name) + '</span>' +
          '<span class="modal-row-size mono">' + fmtSize(o.size) + '</span>' +
          '<span class="modal-row-sel"><button class="btn btn-sm btn-primary">Select</button></span>' +
        '</div>'
      ).join('')
    : '<div style="padding:20px;text-align:center;color:var(--muted);font-size:12px;">' +
        (q ? 'No matches.' : 'No .sqf overlays found.') + '</div>';
}

function selectOverlay(path) {
  const entry = allOverlays.find(o => o.path === path);
  if (entry && !selectedModules.find(m => m.path === path)) {
    selectedModules.push(entry);
    renderModuleChips();
  }
  closeModal('op-modal');
}

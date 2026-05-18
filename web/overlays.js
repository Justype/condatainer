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
    tbody.innerHTML = '<tr><td colspan="4" class="td-empty">' +
      (q ? 'No matches.' : 'No module overlays found.') + '</td></tr>';
    return;
  }
  tbody.innerHTML = rows.map(o =>
    '<tr>' +
      '<td class="td-name">' + escHtml(o.name) + '</td>' +
      '<td class="mono td-muted td-fit">' + fmtSize(o.size) + '</td>' +
      '<td class="mono path-cell td-muted">' + pathTailHtml(o.path) + '</td>' +
      '<td><div class="td-actions">' +
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
          '<span class="modal-row-icon">' + iconSvg('hexagon') + '</span>' +
          '<span class="modal-row-name fw-bold">' + escHtml(o.name) + '</span>' +
          '<span class="modal-row-size mono">' + fmtSize(o.size) + '</span>' +
          '<span class="modal-row-sel"><button class="btn btn-sm btn-primary">Select</button></span>' +
        '</div>'
      ).join('')
    : '<div class="modal-empty">' +
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

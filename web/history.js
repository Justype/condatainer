/* ── History section ─────────────────────── */
async function loadHistory() {
  try {
    const r = await fetch('/api/helpers?limit=100');
    allHistory = (await r.json()) || [];
  } catch { allHistory = []; }
  if (!document.querySelector('#hist-tbody [data-confirm="1"]')) renderHistory();
}

function renderHistory(filter) {
  const q    = (filter || '').toLowerCase();
  let rows = q
    ? allHistory.filter(h =>
        [h.name, h.label, h.status, h.node, h.cwd].join(' ').toLowerCase().includes(q))
    : allHistory.slice();
  rows.sort((a, b) => {
    const aActive = isActiveHistoryStatus(a.status) ? 1 : 0;
    const bActive = isActiveHistoryStatus(b.status) ? 1 : 0;
    if (aActive !== bActive) return bActive - aActive;
    const at = a.ended_at ? new Date(a.ended_at).getTime() : 0;
    const bt = b.ended_at ? new Date(b.ended_at).getTime() : 0;
    return bt - at;
  });
  const tbody = gid('hist-tbody');
  if (!rows.length) {
    tbody.innerHTML = '<tr><td colspan="6" class="td-empty">' +
      (q ? 'No matches.' : 'No history yet.') + '</td></tr>';
    return;
  }
  tbody.innerHTML = rows.map(h => {
    const id  = escHtml(h.id);
    return (
      '<tr data-id="' + id + '" onclick="openDetail(\'' + id + '\')">' +
        '<td class="td-name">' + escHtml(h.label || h.name) + '</td>' +
        '<td class="td-fit"><span class="badge ' + badgeClass(h.status) + '">' + escHtml(h.status) + '</span></td>' +
        '<td class="mono td-muted td-node" title="' + escHtml(h.node || '') + '">' + escHtml((h.node || '—').split('.')[0]) + '</td>' +
        '<td class="mono path-cell td-muted">' + pathTailHtml(h.cwd) + '</td>' +
        '<td class="td-muted td-fit">' + (h.ended_at ? fmtRelTime(h.ended_at) : '—') + '</td>' +
        '<td class="td-nowrap td-fit">' +
          (isActiveHistoryStatus(h.status)
            ? '<button class="btn btn-sm btn-icon" disabled title="Job is still active">' + iconSvg('replay') + '</button>'
            : '<button class="btn btn-sm btn-icon" onclick="rerunHistoryJob(event,\'' + id + '\')" title="Rerun">' + iconSvg('replay') + '</button>') +
          ' ' + (isActiveHistoryStatus(h.status)
            ? '<button class="btn btn-sm btn-danger btn-icon" disabled title="Stop the job before removing it">' + iconSvg('delete') + '</button>'
            : '<button class="btn btn-sm btn-danger btn-icon" onclick="removeHistoryEntry(event,\'' + id + '\',this)" title="Remove entry">' + iconSvg('delete') + '</button>') +
        '</td>' +
      '</tr>'
    );
  }).join('');
}

function filterHistory(q) { renderHistory(q); }

function isActiveHistoryStatus(status) {
  return status === 'pending' || status === 'starting' || status === 'running';
}

function rerunHistoryJob(ev, id) {
  if (ev) {
    ev.preventDefault();
    ev.stopPropagation();
  }
  rerunJob(id);
}

async function removeHistoryEntry(ev, id, btn) {
  if (ev) { ev.preventDefault(); ev.stopPropagation(); }
  if (!btn) return;

  if (btn.dataset.confirm !== '1') {
    btn.dataset.confirm = '1';
    btn.innerHTML = iconSvg('question_mark');
    setTimeout(() => {
      if (btn.dataset.confirm === '1') {
        btn.dataset.confirm = '';
        btn.innerHTML = iconSvg('delete');
      }
    }, 3000);
    return;
  }

  btn.dataset.confirm = '';
  btn.disabled = true;
  try {
    const r = await fetch('/api/helpers/' + id + '/delete', { method: 'DELETE' });
    if (!r.ok) {
      const msg = await r.text();
      btn.disabled = false;
      btn.innerHTML = iconSvg('delete');
      alert('Could not remove entry: ' + msg);
      return;
    }
    if (detailJobId === id) closeDetail();
    allHistory = allHistory.filter(h => h.id !== id);
    renderHistory(gid('hist-search').value);
  } catch (e) {
    btn.disabled = false;
    btn.innerHTML = iconSvg('delete');
    alert('Could not remove entry: ' + e);
  }
}

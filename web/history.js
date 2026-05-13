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
  const rows = q
    ? allHistory.filter(h =>
        [h.name, h.label, h.status, h.node, h.cwd].join(' ').toLowerCase().includes(q))
    : allHistory;
  const tbody = gid('hist-tbody');
  if (!rows.length) {
    tbody.innerHTML = '<tr><td colspan="7" style="padding:20px;text-align:center;color:var(--muted);">' +
      (q ? 'No matches.' : 'No history yet.') + '</td></tr>';
    return;
  }
  tbody.innerHTML = rows.map(h => {
    const dur = calcDuration(h.started_at, h.ended_at, h.status, h.walltime_str);
    const id  = escHtml(h.id);
    return (
      '<tr data-id="' + id + '" onclick="openDetail(\'' + id + '\')">' +
        '<td style="font-weight:600;">' + escHtml(h.label || h.name) + '</td>' +
        '<td><span class="badge ' + badgeClass(h.status) + '">' + escHtml(h.status) + '</span></td>' +
        '<td class="mono" style="color:var(--muted);">' + escHtml(h.node || '—') + '</td>' +
        '<td class="mono path-cell" style="color:var(--muted);">' + pathTailHtml(h.cwd) + '</td>' +
        '<td style="color:var(--muted);">' + (h.started_at ? fmtRelTime(h.started_at) : '—') + '</td>' +
        '<td style="color:var(--muted);">' + dur + '</td>' +
        '<td style="white-space:nowrap;">' +
          '<button class="btn btn-sm" onclick="rerunHistoryJob(event,\'' + id + '\')">' + iconSvg('refresh') + 'Rerun</button>' +
          ' ' + (isActiveHistoryStatus(h.status)
            ? '<button class="btn btn-sm btn-danger" disabled title="Stop the job before removing it">' + iconSvg('delete') + '</button>'
            : '<button class="btn btn-sm btn-danger" onclick="removeHistoryEntry(event,\'' + id + '\',this)" title="Remove entry">' + iconSvg('delete') + '</button>') +
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
  detailJobId = id;
  rerunJob();
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

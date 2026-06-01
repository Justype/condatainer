/* ── Message text renderer (linkifies URLs) ─ */
function appendMsgText(text, container) {
  const urlRe = /https?:\/\/\S+/g;
  let last = 0, m;
  while ((m = urlRe.exec(text)) !== null) {
    if (m.index > last) container.appendChild(document.createTextNode(text.slice(last, m.index)));
    const a = document.createElement('a');
    a.href = m[0]; a.textContent = m[0]; a.target = '_blank'; a.rel = 'noopener noreferrer';
    container.appendChild(a);
    last = m.index + m[0].length;
  }
  if (last < text.length) container.appendChild(document.createTextNode(text.slice(last)));
}

/* ── Jobs section ────────────────────────── */
let _prevJobStatuses = null; // null = first load not yet done

async function loadJobs() {
  try {
    const r = await fetch('/api/helpers?status=active');
    allJobs = (await r.json()) || [];
  } catch { allJobs = []; }
  _checkJobNotifications(allJobs);
  if (!document.querySelector('#job-list-wrap [data-confirm="1"]')) renderJobs();
}

function _checkJobNotifications(jobs) {
  if (_prevJobStatuses === null) {
    _prevJobStatuses = {};
    jobs.forEach(j => { _prevJobStatuses[j.id] = j.status; });
    return;
  }
  const next = {};
  jobs.forEach(j => {
    next[j.id] = j.status;
    if (j.status === 'running' && _prevJobStatuses[j.id] !== 'running') {
      _fireJobNotification(j);
    }
  });
  _prevJobStatuses = next;
}

function _fireJobNotification(job) {
  if (srvNotification !== 'web' && srvNotification !== 'both') return;
  if (typeof Notification === 'undefined' || Notification.permission !== 'granted') return;
  new Notification('Helper ready: ' + job.name, {
    body: job.node ? 'Running on ' + job.node : 'Job started',
    tag: 'condatainer-' + job.id,
  });
}

function renderJobs() {
  const wrap = gid('job-list-wrap');
  if (!allJobs.length) {
    wrap.innerHTML = '<div class="empty"><div>No running jobs.</div>' +
      '<div class="empty-link"><a href="#" onclick="navigate(\'start\');return false;">Start a helper →</a></div></div>';
    return;
  }
  const stableJobIds = oldestRunningJobIdsByName(allJobs);
  wrap.innerHTML = '<div class="job-list">' + allJobs.map(j => jobCardHTML(j, stableJobIds)).join('') + '</div>';
  startWalltimeCountdowns();
}

function oldestRunningJobIdsByName(jobs) {
  const out = {};
  jobs.forEach(j => {
    if (j.status !== 'running' || !j.name || !j.id) return;
    const cur = out[j.name];
    if (!cur || _jobIsOlder(j, cur)) out[j.name] = j;
  });
  const ids = {};
  Object.keys(out).forEach(name => { ids[name] = out[name].id; });
  return ids;
}

function _jobIsOlder(a, b) {
  const at = a.started_at ? new Date(a.started_at).getTime() : NaN;
  const bt = b.started_at ? new Date(b.started_at).getTime() : NaN;
  if (Number.isFinite(at) && Number.isFinite(bt)) return at < bt;
  if (Number.isFinite(at)) return true;
  if (Number.isFinite(bt)) return false;
  return String(a.id || '') < String(b.id || '');
}

function stableHelperURL(j, url) {
  if (!url || j.external_url || j.port === 0) return '';
  try {
    const u = new URL(url);
    if (u.hostname !== 'localhost' && !u.hostname.endsWith('.localhost')) return '';
    u.hostname = j.name + '.localhost';
    return u.toString();
  } catch {
    return '';
  }
}

function jobCardHTML(j, stableJobIds) {
  const endsAt     = j.ends_at_ms || 0;
  const showEndsIn = endsAt && j.status === 'running';
  const walltimeEl = showEndsIn
    ? '<span data-ends="' + endsAt + '" class="walltime-span"></span>'
    : (j.walltime_str ? '⏱ ' + escHtml(j.walltime_str) : '');
  const overlay = j.env_overlay || '';
  const url     = j.access_url || j.external_url || '';
  const stableURL = stableJobIds[j.name] === j.id ? stableHelperURL(j, url) : '';
  const copyURL = stableURL || url;
  const lastMsg = j.last_message || '';
  const id      = escHtml(j.id);
  return (
    '<div class="job-card" id="jcard-' + id + '" onclick="openDetail(\'' + id + '\')">' +
      '<div class="jc-top">' +
        '<span class="jc-name">' + escHtml(j.name) + '</span>' +
        '<span class="badge ' + badgeClass(j.status) + '">' + escHtml(j.status) + '</span>' +
      '</div>' +
      '<div class="jc-meta">' +
        '<span>' + escHtml(j.node || '—') + '</span><span class="jc-sep">|</span>' +
        pathTailHtml(j.cwd) +
        (overlay
          ? '<span class="jc-sep">|</span>' + pathTailHtml(overlay, '200px')
          : '') +
      '</div>' +
      '<div class="jc-walltime" id="jwall-' + id + '">' + walltimeEl + '</div>' +
      (lastMsg ? '<div class="jc-lastmsg">"' + escHtml(lastMsg) + '"</div>' : '') +
      (url
        ? '<div class="jc-access" onclick="event.stopPropagation()">' +
            '<a class="jc-link" href="' + escHtml(url) + '" target="_blank" title="' + escHtml(url) + '">' +
              '<span class="jc-link-text">' + escHtml(url) + '</span><span class="jc-link-ext">↗</span>' +
            '</a>' +
            '<button class="btn btn-sm jc-icon-btn" onclick="copyStr(\'' + escHtml(copyURL) + '\',this)" title="' + escHtml(stableURL ? 'Copy stable URL' : 'Copy URL') + '">' + iconSvg('content_copy') + '</button>' +
            (stableURL
              ? '<a class="btn btn-sm jc-icon-btn" href="' + escHtml(stableURL) + '" target="_blank" title="' + escHtml(stableURL) + '">' + iconSvg('link_2') + '</a>'
              : '') +
          '</div>'
        : '') +
      '<div class="jc-stop-wrap" onclick="event.stopPropagation()">' +
        '<button class="btn btn-danger btn-sm" onclick="stopJob(\'' + id + '\',this)">Stop</button>' +
      '</div>' +
    '</div>'
  );
}

async function stopJob(id, btn) {
  if (btn.dataset.confirm === '1') {
    btn.dataset.confirm = '';
    btn.innerHTML = 'Stopping…';
    btn.disabled = true;

    let ok = false;
    try {
      const r = await fetch('/api/helpers/' + encodeURIComponent(id) + '/stop', { method: 'POST' });
      ok = r.ok;
    } catch {}

    try {
      const r = await fetch('/api/helpers?status=active');
      allJobs = (await r.json()) || [];
    } catch {}
    if (detailJobId === id) closeDetail();
    renderJobs();

    if (!ok) {
      const card = gid('jcard-' + id);
      const newBtn = card && card.querySelector('.btn-danger');
      if (newBtn) {
        newBtn.textContent = 'Stop failed';
        setTimeout(() => { if (newBtn.textContent === 'Stop failed') newBtn.textContent = 'Stop'; }, 3000);
      }
    }
  } else {
    btn.dataset.confirm = '1';
    btn.innerHTML = 'Confirm?';
    setTimeout(() => {
      if (btn.dataset.confirm === '1') { btn.dataset.confirm = ''; btn.innerHTML = 'Stop'; }
    }, 3000);
  }
}

/* ── Walltime countdowns ──────────────────── */
function startWalltimeCountdowns() {
  tickCountdowns();
  if (!window._walltimeTick) window._walltimeTick = setInterval(tickCountdowns, 30000);
}
function tickCountdowns() {
  document.querySelectorAll('.walltime-span').forEach(el => {
    const ends = Number(el.dataset.ends);
    const diff = ends - Date.now();
    const parent = el.parentElement;
    if (diff <= 0) {
      el.innerHTML = iconSvg('warning') + ' Walltime reached';
      parent && parent.classList.add('warn');
    } else {
      const h = Math.floor(diff / 3600e3);
      const m = Math.floor((diff % 3600e3) / 60e3);
      el.textContent = 'Ends in ' + h + 'h ' + m + 'm';
      parent && parent.classList.remove('warn');
    }
  });
}

/* ── Detail panel ────────────────────────── */
function openDetail(id) {
  const job = allJobs.find(j => j.id === id) || allHistory.find(j => j.id === id);
  if (!job) return;
  const detail = gid('detail');
  if (detailJobId === id && detail && detail.classList.contains('open')) {
    closeDetail();
    return;
  }
  detailJobId = id;
  const isRunning = job.status === 'running' || job.status === 'pending' || job.status === 'starting';

  document.querySelectorAll('.job-card,.data-table tbody tr').forEach(e => e.classList.remove('selected'));
  const card = gid('jcard-' + id) || document.querySelector('tr[data-id="' + id + '"]');
  if (card) card.classList.add('selected');

  gid('d-title').textContent = job.name;
  gid('d-badge').innerHTML = '<span class="badge ' + badgeClass(job.status) + '">' + escHtml(job.status) + '</span>';

  // Meta — grouped: Location | Container | Job | Access
  let meta = '';
  // Location
  meta += '<span class="d-key">Working Dir</span>' + pathTailHtml(job.cwd, null, 'd-val');
  // Container
  if (job.base_image) meta += '<span class="d-key">Base Image</span>' + pathTailHtml(job.base_image, null, 'd-val');
  meta += '<span class="d-key">Env Overlay</span>' + pathTailHtml(job.env_overlay, null, 'd-val');
  if (job.overlays && job.overlays.length) {
    meta += '<span class="d-key">Modules</span><span class="d-val-col">' +
      job.overlays.map(p => {
        const o = allOverlays.find(ov => ov.path === p || ov.name === p);
        return o
          ? '<span title="' + escHtml(p) + '">' + escHtml(o.name) + '</span>'
          : pathTailHtml(p, null);
      }).join('') + '</span>';
  }
  // Job
  if (job.job_id) meta += '<span class="d-key">Job ID</span><span class="d-val">' + escHtml(job.job_id) + '</span>';
  meta += '<span class="d-key">Started</span><span class="d-val">' + (job.started_at ? fmtDate(job.started_at) : '—') + '</span>';
  if (job.node) meta += '<span class="d-key">Node</span><span class="d-val">' + escHtml(job.node) + '</span>';
  if (job.port > 0) {
    const bindAddr = (job.bind_all ? '0.0.0.0' : '127.0.0.1') + ':' + job.port;
    meta += '<span class="d-key">Bind</span><span class="d-val">' + escHtml(bindAddr) + '</span>';
  }
  // Access
  const url = job.access_url || job.external_url;
  if (url) meta += '<span class="d-key">URL</span><span class="d-val"><a href="' + escHtml(url) + '" target="_blank">' + escHtml(url) + '</a></span>';
  gid('d-meta').innerHTML = meta;

  // Resources
  gid('d-resources').innerHTML =
    '<div class="d-res"><span class="d-res-k">CPUs</span><span class="d-res-v">'    + escHtml(job.cpus || '—')         + '</span></div>' +
    '<div class="d-res"><span class="d-res-k">Memory</span><span class="d-res-v">'  + escHtml(job.mem || '—')          + '</span></div>' +
    '<div class="d-res"><span class="d-res-k">GPU</span><span class="d-res-v">'     + escHtml(job.gpu || '—')          + '</span></div>' +
    '<div class="d-res"><span class="d-res-k">Walltime</span><span class="d-res-v">'+ escHtml(job.walltime_str || '—') + '</span></div>';

  // Params
  const params    = job.params || {};
  const paramKeys = Object.keys(params);
  const dp = gid('d-params');
  if (paramKeys.length) {
    dp.style.display = '';
    dp.innerHTML = paramKeys.map(k =>
      '<div class="d-param-row"><span class="d-param-k">' + escHtml(k) + '</span>' +
      '<span class="d-param-v">' + escHtml(params[k]) + '</span></div>'
    ).join('');
  } else { dp.style.display = 'none'; }

  gid('d-stop').style.display = isRunning ? '' : 'none';
  gid('d-rerun').disabled = isRunning;
  gid('d-log-link').href = '/api/helpers/' + encodeURIComponent(id) + '/joblog';

  switchDetailTab(detailTabActive);
  detail.removeAttribute('inert');
  detail.setAttribute('aria-hidden', 'false');
  detail.classList.add('open');
}

function switchDetailTab(tab) {
  detailTabActive = tab;
  document.querySelectorAll('.d-tab').forEach(t => t.classList.toggle('active', t.dataset.dtab === tab));
  const content = gid('d-tab-content');

  if (detailSSE) { detailSSE.close(); detailSSE = null; }

  if (tab === 'messages') {
    content.innerHTML = '';
    detailSSE = new EventSource('/api/helpers/' + encodeURIComponent(detailJobId) + '/logs');
    detailSSE.onopen = function() {
      // Clear on each (re)connect so backfill doesn't duplicate old messages
      content.innerHTML = '';
    };
    detailSSE.onmessage = function(ev) {
      try {
        const m   = JSON.parse(ev.data);
        const lvl = m.level === 'warn' ? 'warn' : m.level === 'error' ? 'err' : 'info';
        const div = document.createElement('div');
        div.className = 'log-line log-' + lvl;
        if (m.ts) {
          const ts = document.createElement('span');
          ts.className = 'log-ts';
          try { ts.textContent = new Date(m.ts).toLocaleTimeString(); } catch {}
          div.appendChild(ts);
        }
        appendMsgText(m.text || '', div);
        content.appendChild(div);
        content.scrollTop = content.scrollHeight;
        // Update last-message snippet on the card
        const card = gid('jcard-' + detailJobId);
        if (card) {
          const lm = card.querySelector('.jc-lastmsg');
          if (lm) lm.textContent = '"' + (m.text || '') + '"';
        }
      } catch {}
    };
    detailSSE.onerror = function() {
      // CONNECTING means the browser is auto-reconnecting — don't interfere.
      if (detailSSE && detailSSE.readyState === EventSource.CONNECTING) return;
      if (detailSSE) { detailSSE.close(); detailSSE = null; }
      const dim = document.createElement('div');
      dim.className = 'log-line log-dim';
      dim.textContent = '— stream ended —';
      content.appendChild(dim);
    };
  } else {
    // Job log tab — request last 1000 lines to avoid freezing on large logs
    content.innerHTML = '<div class="state-muted">Loading log…</div>';
    fetch('/api/helpers/' + encodeURIComponent(detailJobId) + '/joblog?tail=1000')
      .then(r => {
        const truncated = r.headers.get('X-Log-Truncated') === 'true';
        return r.text().then(text => ({ text, truncated }));
      })
      .then(({ text, truncated }) => {
        content.innerHTML = '';
        if (truncated) {
          const notice = document.createElement('div');
          notice.className = 'log-line log-dim';
          notice.classList.add('log-notice');
          notice.textContent = '— showing last 1000 lines (use Full log ↗ for complete output) —';
          content.appendChild(notice);
        }
        const pre = document.createElement('pre');
        pre.className = 'log-pre';
        pre.textContent = text || '(empty log)';
        content.appendChild(pre);
        content.scrollTop = content.scrollHeight;
      })
      .catch(() => { content.innerHTML = '<div class="state-danger">Failed to load log.</div>'; });
  }
}

function closeDetail() {
  const detail = gid('detail');
  detail.classList.remove('open');
  detail.setAttribute('inert', '');
  detail.setAttribute('aria-hidden', 'true');
  if (detail.contains(document.activeElement)) document.activeElement.blur();
  document.querySelectorAll('.job-card,.data-table tbody tr').forEach(e => e.classList.remove('selected'));
  if (detailSSE) { detailSSE.close(); detailSSE = null; }
  detailJobId = null;
}

function detailStop() {
  const btn = gid('d-stop');
  if (btn && detailJobId) stopJob(detailJobId, btn);
}

function rerunJob(jobOrId) {
  const id = typeof jobOrId === 'string' ? jobOrId : detailJobId;
  const job = typeof jobOrId === 'object' && jobOrId
    ? jobOrId
    : allJobs.find(j => j.id === id) || allHistory.find(j => j.id === id);
  if (!job) return;
  navigate('start');
  queueStartSelection(job.name, {
    cpus:    job.cpus,
    mem:     job.mem,
    wall:    job.walltime_str,
    gpu:     job.gpu     || '',
    cwd:     job.cwd     || '',
    overlay: job.env_overlay || '',
    modules:   (job.overlays || []).map(p => allOverlays.find(o => o.name === p || o.path === p)).filter(Boolean),
    externals: (job.overlays || []).filter(p => !allOverlays.find(o => o.name === p || o.path === p)),
    params:  job.params  || {},
  });
}

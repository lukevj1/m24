"use strict";
(() => {
  const E = LithosEngine;
  const $ = (id) => document.getElementById(id);
  const fmt = (x, d = 2) => (Number.isFinite(x) ? x.toFixed(d) : "–");
  const pct = (p) => `${Math.round(p * 100)}%`;
  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  const NS = "http://www.w3.org/2000/svg";
  const el = (tag, attrs = {}, parent) => {
    const n = document.createElementNS(NS, tag);
    for (const [k, v] of Object.entries(attrs)) n.setAttribute(k, v);
    if (parent) parent.appendChild(n);
    return n;
  };
  const clockStr = (h) => {
    const hh = ((Math.floor(h) % 24) + 24) % 24, mm = Math.round((h - Math.floor(h)) * 60) % 60;
    return `${String(hh).padStart(2, "0")}:${String(mm).padStart(2, "0")}`;
  };

  // ---------------------------------------------------------------- presets
  const PRESETS = [
    { id: "afternoon", label: "Afternoon blood test",
      v: { age: 45, sex: "F", wt: 70, ht: 165, scr: 75, regimen: "nocte", dose: 750, days: 30, target: "0.6-0.8", lastDose: "evening", hrs: 18.5, level: 0.55 } },
    { id: "morningdose", label: "Morning dose taken before the test",
      v: { age: 58, sex: "F", wt: 72, ht: 163, scr: 82, regimen: "bd", dose: 1000, days: 40, target: "0.6-0.8", lastDose: "morning", hrs: 2, level: 1.24 } },
    { id: "newdose", label: "Four days after starting",
      v: { age: 34, sex: "M", wt: 86, ht: 181, scr: 88, regimen: "nocte", dose: 1000, days: 4, target: "0.6-0.8", lastDose: "evening", hrs: 19.5, level: 0.60 } },
    { id: "older", label: "Older adult, eGFR in the 50s",
      v: { age: 72, sex: "M", wt: 80, ht: 175, scr: 118, regimen: "nocte", dose: 750, days: 60, target: "0.4-0.6", lastDose: "evening", hrs: 12, level: 0.78 } },
  ];

  const doseSel = $("dose");
  for (let mg = 250; mg <= 2000; mg += 250) {
    const o = document.createElement("option");
    o.value = mg; o.textContent = `${mg} mg`;
    doseSel.appendChild(o);
  }

  function readInputs() {
    return {
      age: +$("age").value, sex: $("sex").value, wt: +$("wt").value, ht: +$("ht").value, scr: +$("scr").value,
      regimen: $("regimen").value, dose: +$("dose").value, days: Math.max(1, Math.round(+$("days").value || 1)),
      target: $("target").value.split("-").map(Number), lastDose: $("lastDose").value,
      hrs: +$("hrs").value, level: +$("level").value,
    };
  }
  function setInputs(v) {
    for (const k of ["age", "sex", "wt", "ht", "scr", "regimen", "dose", "days", "target", "lastDose", "hrs", "level"]) {
      if (v[k] !== undefined) $(k).value = v[k];
    }
    syncRegimen();
  }
  function syncRegimen() {
    const bd = $("regimen").value === "bd";
    $("lastDoseField").hidden = !bd;
    const hrs = $("hrs");
    hrs.max = bd ? 12 : 24;
    if (+hrs.value > +hrs.max) hrs.value = hrs.max;
  }

  function regimenOf(inp) {
    if (inp.regimen === "nocte") return [{ clock: 21, mg: inp.dose }];
    const k = inp.dose / 250, base = Math.floor(k / 2), rem = k % 2;
    const m = base * 250, e = base * 250 + rem * 250;
    return m > 0 ? [{ clock: 8, mg: m }, { clock: 20, mg: e }] : [{ clock: 20, mg: e }];
  }
  function lastClockOf(inp, regimen) {
    if (inp.regimen === "nocte") return 21;
    return inp.lastDose === "morning" && regimen.length === 2 ? 8 : 20;
  }

  // -------------------------------------------------------------- the engine
  function evaluate(inp) {
    const regimen = regimenOf(inp);
    const lastClock = lastClockOf(inp, regimen);
    const scen = {
      cov: { age: inp.age, sex: inp.sex, wt: inp.wt, ht: inp.ht, scr: inp.scr },
      form: "IR", regimen, days: inp.days, lastClock,
      level: { hoursAfterLast: inp.hrs, value: inp.level, timingSd: 0.25 },
    };
    const ens = E.fitEnsemble(scen, true);
    const std = E.ensLi12(ens, regimen);
    const template = regimen.map(r => r.clock);
    const R = E.ensRecommend(ens, { target: inp.target, template });
    const lo = inp.target[0], hi = inp.target[1];
    const current = R.options.find(o => o.daily === inp.dose) || null;

    // choose: best by P(in range) under the safety ceiling, within a 1.5x step of the current dose
    // (one 250 mg tablet up or down is always allowed, as in lithos/dosing.py)
    const safe = R.options.filter(o => o.pHigh <= 0.05);
    const up = Math.max(inp.dose * 1.5, inp.dose + 250), down = Math.min(inp.dose / 1.5, inp.dose - 250);
    const inStep = safe.filter(o => o.daily <= up && o.daily >= down && Math.abs(o.daily - inp.dose) <= 500);
    const pool = (inStep.length ? inStep : safe.length ? safe : R.options).slice();
    const center = (lo + hi) / 2;
    pool.sort((a, b) => (Math.round(b.pIn * 100) - Math.round(a.pIn * 100)) || (Math.abs(a.median - center) - Math.abs(b.median - center)));
    let chosen = pool[0];
    let action = chosen.daily === inp.dose ? "keep" : (chosen.daily > inp.dose ? "increase" : "decrease");
    if (current && current.pHigh <= 0.05 && current.pIn >= chosen.pIn - 0.10 && current.median >= lo && current.median <= hi) {
      chosen = current; action = "keep";
    }
    // know when not to answer
    if (current && chosen !== current) {
      const vague = (current.hi - current.lo) > 0.35;
      const plausible = current.median >= lo - 0.1 && current.median <= hi + 0.1;
      const pOver12 = current.vals.filter(v => v > 1.2).length / current.vals.length;
      if (vague && plausible && pOver12 < 0.10) { chosen = current; action = "repeat"; }
    }
    const egfr = E.egfr2021(inp.scr, inp.age, inp.sex);
    return { inp, scen, ens, std, R, chosen, current, action, egfr, regimen, lastClock };
  }

  function status(median, lo, hi) {
    if (median < lo) return { cls: "warn", text: `Below target (${lo}–${hi})` };
    if (median > 1.0) return { cls: "crit", text: `Above ${hi}, over the usual ceiling of 1.0` };
    if (median > hi) return { cls: "warn", text: `Above target (${lo}–${hi})` };
    return { cls: "good", text: `In target (${lo}–${hi})` };
  }
  function naive(level, lo, hi) {
    if (level > hi) return `${fmt(level)}: looks above target, so "reduce the dose"`;
    if (level < lo) return `${fmt(level)}: looks below target, so "increase the dose"`;
    return `${fmt(level)}: looks in range, so "no change"`;
  }
  function doseText(split, regimen) {
    const clocks = regimen.length === 2 ? ["08:00", "20:00"] : ["21:00"];
    return split.map((mg, i) => `${mg} mg (${mg / 250} × 250 mg) at ${clocks[regimen.length === 2 ? i : 0]}`).join(" + ");
  }

  // --------------------------------------------------------------- rendering
  function renderVerdict(res) {
    const { inp, std } = res;
    const [lo, hi] = inp.target;
    $("stdVal").innerHTML = `${fmt(std.median)}<small>mmol/L</small>`;
    const st = status(std.median, lo, hi);
    const pIn = std.vals.filter(v => v >= lo && v <= hi).length / std.vals.length;
    $("stdStatus").innerHTML = `<span class="pill ${st.cls}">${st.text}</span>`;
    $("rawRead").textContent = naive(inp.level, lo, hi);
    $("engRead").textContent = `12-hour equivalent ${fmt(std.median)} (90% interval ${fmt(std.lo)}–${fmt(std.hi)}); ${pct(pIn)} probability in range`;
    const flags = [];
    if (inp.hrs < 6) flags.push("Drawn in the absorption and distribution phase: informative only with an accurate dose time, so the interval is wide.");
    else if (inp.hrs > 16 && inp.regimen === "nocte") flags.push("Drawn late in the dosing interval: a face-value reading underestimates the 12-hour level.");
    const t12 = res.R.halfLifeH;
    if (inp.days * 24 < 3.3 * t12) flags.push(`Not yet at steady state (${inp.days} day${inp.days > 1 ? "s" : ""} on this dose, half-life about ${Math.round(t12)} h): the engine projects the steady-state level.`);
    if (inp.age >= 65 && lo >= 0.6) flags.push("Age 65+: ISBD/IGSLi suggest 0.4–0.6 for most older adults.");
    if (res.egfr < 30) flags.push("eGFR below 30: US labelling does not recommend lithium.");
    if (std.median > 1.2 || inp.level >= 1.5) flags.push("Check for toxicity symptoms now; symptoms override the number.");
    $("flags").innerHTML = flags.map(f => `<li>${f}</li>`).join("");
  }

  function renderRec(res) {
    const { chosen, action, R, inp, regimen, current } = res;
    const heads = {
      keep: "No change", increase: "Increase the dose", decrease: "Decrease the dose",
      repeat: "Repeat the level before changing anything",
    };
    $("recHead").textContent = heads[action];
    $("recDose").textContent = doseText(chosen.split, regimen);
    $("recStats").innerHTML = `<span>Forecast 12-h level <b>${fmt(chosen.median)}</b> (${fmt(chosen.lo)}–${fmt(chosen.hi)})</span>` +
      `<span>P(in range) <b>${pct(chosen.pIn)}</b></span><span>P(&gt; ${R.high.toFixed(1)}) <b>${pct(chosen.pHigh)}</b></span>`;
    const hl = R.halfLifeH;
    const recheck = Math.max(3, Math.ceil(2.5 * hl / 24)), conv = Math.max(5, Math.ceil(5 * hl / 24));
    let next;
    if (action === "repeat") next = `The evidence is too uncertain to justify a change (90% interval ${fmt(current.lo)}–${fmt(current.hi)}). Repeat a level tomorrow morning, at least 10 h after the evening dose and before any morning dose, recording both times.`;
    else if (action === "keep") next = inp.days < 14 ? "Confirm with another level in about a week (weekly until stable)." : "Next level per the routine schedule; sooner if anything changes.";
    else next = `Next level in ${recheck} days, any time at least 6 h after a dose, with the dose time recorded. Waiting for steady state would take about ${conv} days.`;
    $("recNext").textContent = next;
    const w = res.ens.weights, names = res.ens.members.map(m => m.model.name);
    $("modelLine").textContent = `Priors weighted by fit: ${names.map((n, i) => `${n} ${pct(w[i])}`).join(" · ")}. Estimated half-life ${Math.round(hl)} h. eGFR ${Math.round(res.egfr)}.`;
  }

  // gauge
  const BANDS = [
    { a: 0.0, b: 0.4, v: "--band-sub", label: "low" },
    { a: 0.4, b: 0.6, v: "--band-low", label: "0.4–0.6" },
    { a: 0.6, b: 0.8, v: "--band-target", label: "0.6–0.8" },
    { a: 0.8, b: 1.0, v: "--band-high", label: "0.8–1.0" },
    { a: 1.0, b: 1.5, v: "--band-above", label: "above usual" },
    { a: 1.5, b: 1.8, v: "--band-toxic", label: "toxic" },
  ];
  function renderGauge(res) {
    const box = $("gauge");
    box.innerHTML = "";
    const W = Math.max(280, box.clientWidth), H = 96, padL = 8, padR = 8;
    const x = (v) => padL + (Math.min(v, 1.8) / 1.8) * (W - padL - padR);
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, width: W, height: H, role: "img",
      "aria-label": `Standardised 12-hour level ${fmt(res.std.median)} with 90% interval ${fmt(res.std.lo)} to ${fmt(res.std.hi)}; measured ${fmt(res.inp.level)}` }, box);
    const y0 = 30, bh = 22;
    for (const b of BANDS) {
      el("rect", { x: x(b.a), y: y0, width: Math.max(0, x(b.b) - x(b.a) - 2), height: bh, fill: css(b.v), rx: 3 }, svg);
      const wpx = x(b.b) - x(b.a);
      if (wpx > b.label.length * 6.2 + 6) {
        const t = el("text", { x: (x(b.a) + x(b.b)) / 2, y: y0 + bh / 2 + 4, "text-anchor": "middle", fill: css("--ink-2"),
          "font-size": 10.5, "font-family": "IBM Plex Mono, monospace" }, svg);
        t.textContent = b.label;
      }
    }
    // target outline
    const [lo, hi] = res.inp.target;
    el("rect", { x: x(lo), y: y0 - 3, width: x(hi) - x(lo) - 2, height: bh + 6, fill: "none", stroke: css("--ink"), "stroke-width": 1.2, rx: 4 }, svg);
    // interval + median
    const s = res.std;
    el("rect", { x: x(s.lo), y: y0 + bh + 10, width: Math.max(2, x(s.hi) - x(s.lo)), height: 6, rx: 3, fill: css("--series-engine"), opacity: 0.35 }, svg);
    el("line", { x1: x(s.median), x2: x(s.median), y1: y0 - 8, y2: y0 + bh + 18, stroke: css("--series-engine"), "stroke-width": 3, "stroke-linecap": "round" }, svg);
    const tm = el("text", { x: Math.min(W - 60, Math.max(4, x(s.median) - 30)), y: 14, fill: css("--accent-ink"), "font-size": 11.5, "font-family": "IBM Plex Mono, monospace", "font-weight": 500 }, svg);
    tm.textContent = `12-h ${fmt(s.median)}`;
    // measured
    const mx = x(res.inp.level);
    el("circle", { cx: mx, cy: y0 + bh + 13, r: 5, fill: css("--surface"), stroke: css("--ink"), "stroke-width": 1.6 }, svg);
    const tr = el("text", { x: Math.min(W - 80, Math.max(4, mx - 30)), y: H - 4, fill: css("--ink-2"), "font-size": 10.5, "font-family": "IBM Plex Mono, monospace" }, svg);
    tr.textContent = `measured ${fmt(res.inp.level)}`;
  }

  // curve
  function renderCurve(res) {
    const box = $("curve");
    box.innerHTML = "";
    const W = Math.max(280, box.clientWidth), H = 230, m = { l: 40, r: 12, t: 12, b: 30 };
    const members = res.ens.members;
    const tObs = members[0].tObs;
    const tFrom = Math.max(0, tObs - 60), tTo = tObs + 10;
    const pts = E.ensCurve(res.ens, tFrom, tTo, Math.min(160, Math.round(W / 3)));
    const yMax = Math.max(1.2, res.inp.level * 1.15, ...pts.map(p => p.hi)) * 1.05;
    const x = (t) => m.l + (t - tFrom) / (tTo - tFrom) * (W - m.l - m.r);
    const y = (c) => H - m.b - (c / yMax) * (H - m.t - m.b);
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, width: W, height: H, role: "img",
      "aria-label": "Estimated lithium concentration over the last three days with 90% band and the measured level" }, box);
    // grid + y ticks
    for (let v = 0; v <= yMax; v += 0.2) {
      const yy = y(v);
      el("line", { x1: m.l, x2: W - m.r, y1: yy, y2: yy, stroke: css("--line"), "stroke-width": 1 }, svg);
      const t = el("text", { x: m.l - 6, y: yy + 3, "text-anchor": "end", fill: css("--muted"), "font-size": 10, "font-family": "IBM Plex Mono, monospace" }, svg);
      t.textContent = v.toFixed(1);
    }
    // target band
    const [lo, hi] = res.inp.target;
    el("rect", { x: m.l, y: y(hi), width: W - m.l - m.r, height: y(lo) - y(hi), fill: css("--band-target"), opacity: 0.8 }, svg);
    // band + median
    const area = pts.map((p, i) => `${i ? "L" : "M"}${x(p.t)},${y(p.hi)}`).join("") +
      pts.slice().reverse().map(p => `L${x(p.t)},${y(p.lo)}`).join("") + "Z";
    el("path", { d: area, fill: css("--series-engine"), opacity: 0.16 }, svg);
    el("path", { d: pts.map((p, i) => `${i ? "L" : "M"}${x(p.t)},${y(p.med)}`).join(""), fill: "none",
      stroke: css("--series-engine"), "stroke-width": 2, "stroke-linejoin": "round", "stroke-linecap": "round" }, svg);
    // x ticks at midnights and doses
    for (let t = Math.ceil(tFrom / 24) * 24; t <= tTo; t += 24) {
      el("line", { x1: x(t), x2: x(t), y1: m.t, y2: H - m.b, stroke: css("--line"), "stroke-width": 1 }, svg);
    }
    for (let t = Math.ceil(tFrom / 6) * 6; t <= tTo; t += 6) {
      const tx = el("text", { x: x(t), y: H - m.b + 14, "text-anchor": "middle", fill: css("--muted"), "font-size": 10, "font-family": "IBM Plex Mono, monospace" }, svg);
      tx.textContent = clockStr(t);
    }
    for (const d of members[0].doses) {
      if (d.t < tFrom || d.t > tTo) continue;
      el("path", { d: `M${x(d.t)},${H - m.b} l-4,7 h8 z`, fill: css("--ink-2") }, svg);
    }
    // 12-h standard time after the last evening dose before the sample
    const eveClock = Math.max(...res.regimen.map(r => r.clock));
    const lastEve = Math.floor((tObs - eveClock) / 24) * 24 + eveClock;
    const t12 = lastEve + 12;
    if (t12 >= tFrom && t12 <= tTo) {
      el("line", { x1: x(t12), x2: x(t12), y1: m.t, y2: H - m.b, stroke: css("--accent-ink"), "stroke-width": 1.2 }, svg);
      const lab = el("text", { x: x(t12) + 4, y: m.t + 10, fill: css("--accent-ink"), "font-size": 10.5, "font-family": "IBM Plex Mono, monospace" }, svg);
      lab.textContent = "12 h";
    }
    // measured point
    el("circle", { cx: x(tObs), cy: y(res.inp.level), r: 5.5, fill: css("--ink"), stroke: css("--surface"), "stroke-width": 2 }, svg);
    // hover
    const tip = document.createElement("div");
    tip.className = "tip"; tip.hidden = true; box.appendChild(tip);
    const cross = el("line", { y1: m.t, y2: H - m.b, stroke: css("--muted"), "stroke-width": 1, visibility: "hidden" }, svg);
    const hit = el("rect", { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: "transparent" }, svg);
    const move = (ev) => {
      const r = svg.getBoundingClientRect();
      const px = (ev.touches ? ev.touches[0].clientX : ev.clientX) - r.left;
      const t = tFrom + (px - m.l) / (W - m.l - m.r) * (tTo - tFrom);
      let best = pts[0];
      for (const p of pts) if (Math.abs(p.t - t) < Math.abs(best.t - t)) best = p;
      cross.setAttribute("x1", x(best.t)); cross.setAttribute("x2", x(best.t)); cross.setAttribute("visibility", "visible");
      tip.hidden = false;
      tip.textContent = `${clockStr(best.t)}  ${fmt(best.med)} (${fmt(best.lo)}–${fmt(best.hi)}) mmol/L`;
      const left = Math.min(W - 190, Math.max(0, x(best.t) + 10));
      tip.style.left = `${left}px`; tip.style.top = `${m.t}px`;
    };
    hit.addEventListener("mousemove", move);
    hit.addEventListener("touchmove", move, { passive: true });
    hit.addEventListener("mouseleave", () => { tip.hidden = true; cross.setAttribute("visibility", "hidden"); });
  }

  // hero slip (fixed example = first preset)
  function renderHero() {
    const v = PRESETS[0].v;
    const res = evaluate({ ...v, target: v.target.split("-").map(Number) });
    const s = res.std, st = status(s.median, 0.6, 0.8);
    const drawn = clockStr(21 + v.hrs);
    $("heroSlip").innerHTML = `
      <div class="slip-head"><span>Chemical pathology</span><span>Specimen: serum</span></div>
      <div class="slip-body">
        <div class="slip-row"><span class="k">Regimen</span><span class="v">${v.dose} mg lithium carbonate at 21:00</span></div>
        <div class="slip-row"><span class="k">Collected</span><span class="v">${drawn}, ${fmt(v.hrs, 1)} h after dose</span></div>
        <div class="slip-row raw"><span class="k">Lithium (as measured)</span><span class="v">${fmt(v.level)} mmol/L</span></div>
        <div class="slip-std">
          <span class="k">Standardised 12-hour level</span>
          <span class="big">${fmt(s.median)} <span style="font-size:1rem;font-family:var(--font-body);font-weight:400;color:var(--ink-2)">mmol/L</span></span>
          <span class="ci">90% interval ${fmt(s.lo)}–${fmt(s.hi)}</span>
          <span class="pill ${st.cls}">${st.text}</span>
        </div>
        <p class="slip-note">Read at face value, ${fmt(v.level)} would suggest a dose increase. Model: ensemble of two population priors weighted by fit (${res.ens.members.map((mm, i) => `${Math.round(res.ens.weights[i] * 100)}%`).join(" / ")}). Research prototype, not for clinical use.</p>
      </div>`;
  }

  // ------------------------------------------------------------ interaction
  let timer = null, lastRes = null;
  function update() {
    const inp = readInputs();
    $("hrsOut").textContent = `${fmt(inp.hrs, 1)} h`;
    const regimen = regimenOf(inp);
    const lc = lastClockOf(inp, regimen);
    $("clockOut").textContent = `Blood drawn at ${clockStr(lc + inp.hrs)} (last dose at ${clockStr(lc)})`;
    const e = E.egfr2021(inp.scr, inp.age, inp.sex);
    $("renalDerived").textContent = `eGFR ${Math.round(e)} mL/min/1.73m² · CrCl ${Math.round(E.crclCG(inp.scr, inp.age, inp.sex, inp.wt))} mL/min`;
    if (!(inp.level > 0) || !(inp.age >= 18) || !(inp.wt > 30) || !(inp.scr > 20)) return;
    clearTimeout(timer);
    timer = setTimeout(() => {
      lastRes = evaluate(inp);
      renderVerdict(lastRes); renderRec(lastRes); renderGauge(lastRes); renderCurve(lastRes);
    }, 90);
  }
  const presetBox = $("presets");
  for (const p of PRESETS) {
    const b = document.createElement("button");
    b.type = "button"; b.className = "preset"; b.textContent = p.label; b.dataset.id = p.id;
    b.setAttribute("aria-pressed", "false");
    b.addEventListener("click", () => {
      for (const o of presetBox.children) o.setAttribute("aria-pressed", String(o === b));
      setInputs(p.v); update();
    });
    presetBox.appendChild(b);
  }
  $("inputs").addEventListener("input", (ev) => {
    if (ev.target.id === "regimen") syncRegimen();
    for (const o of presetBox.children) o.setAttribute("aria-pressed", "false");
    update();
  });
  $("inputs").addEventListener("submit", (ev) => ev.preventDefault());
  const ro = new ResizeObserver(() => { if (lastRes) { renderGauge(lastRes); renderCurve(lastRes); } drawSim(); });
  ro.observe($("curve"));

  // ------------------------------------------------------------ simulations
  function barChart(box, rows, opts) {
    box.innerHTML = "";
    const W = Math.max(260, box.clientWidth), rowH = 26, gap = 10, m = { l: opts.labelW || 150, r: 54, t: 6, b: 6 };
    let H = m.t + m.b;
    const groups = rows.map(g => ({ ...g, h: g.bars.length * rowH + gap }));
    H += groups.reduce((a, g) => a + g.h + (g.title ? 18 : 0), 0);
    const max = opts.max || Math.max(...rows.flatMap(g => g.bars.map(b => b.v))) * 1.1;
    const x = (v) => m.l + (v / max) * (W - m.l - m.r);
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, width: W, height: H, role: "img", "aria-label": opts.aria }, box);
    const tip = document.createElement("div"); tip.className = "tip"; tip.hidden = true; box.appendChild(tip);
    let yy = m.t;
    for (const g of groups) {
      if (g.title) {
        const t = el("text", { x: 0, y: yy + 12, fill: css("--ink"), "font-size": 12, "font-weight": 700, "font-family": "Atkinson Hyperlegible, system-ui, sans-serif" }, svg);
        t.textContent = g.title; yy += 18;
      }
      for (const b of g.bars) {
        const lab = el("text", { x: m.l - 8, y: yy + rowH / 2 + 4, "text-anchor": "end", fill: css("--ink-2"), "font-size": 11.5, "font-family": "Atkinson Hyperlegible, system-ui, sans-serif" }, svg);
        lab.textContent = b.label;
        const color = b.emph ? css("--series-engine") : css("--series-status");
        const w = Math.max(2, x(b.v) - m.l);
        const bh = Math.min(16, rowH - 8);
        if (w >= 8) {
          el("path", { d: `M${m.l},${yy + (rowH - bh) / 2} h${w - 4} a4,4 0 0 1 4,4 v${bh - 8} a4,4 0 0 1 -4,4 h-${w - 4} z`, fill: color }, svg);
        } else {
          el("rect", { x: m.l, y: yy + (rowH - bh) / 2, width: w, height: bh, fill: color }, svg);
        }
        const val = el("text", { x: m.l + w + 6, y: yy + rowH / 2 + 4, fill: css("--ink"), "font-size": 11.5, "font-family": "IBM Plex Mono, monospace" }, svg);
        val.textContent = opts.fmt(b.v);
        const hit = el("rect", { x: 0, y: yy, width: W, height: rowH, fill: "transparent" }, svg);
        hit.addEventListener("mousemove", (ev) => {
          const r = svg.getBoundingClientRect();
          tip.hidden = false; tip.textContent = `${b.label}: ${opts.fmt(b.v)}${b.note ? " · " + b.note : ""}`;
          tip.style.left = `${Math.min(W - 200, ev.clientX - r.left + 12)}px`; tip.style.top = `${yy - 4}px`;
        });
        hit.addEventListener("mouseleave", () => { tip.hidden = true; });
        yy += rowH;
      }
      yy += gap;
    }
  }
  function lineChart(box, series, opts) {
    box.innerHTML = "";
    const W = Math.max(260, box.clientWidth), H = 200, m = { l: 40, r: 88, t: 10, b: 30 };
    const xs = series[0].pts.map(p => p[0]);
    const xMin = Math.min(...xs), xMax = Math.max(...xs);
    const yMax = Math.max(...series.flatMap(s => s.pts.map(p => p[1]))) * 1.15;
    const x = (v) => m.l + (v - xMin) / (xMax - xMin) * (W - m.l - m.r);
    const y = (v) => H - m.b - (v / yMax) * (H - m.t - m.b);
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, width: W, height: H, role: "img", "aria-label": opts.aria }, box);
    const step = yMax > 0.3 ? 0.1 : 0.05;
    for (let v = 0; v <= yMax + 1e-9; v += step) {
      el("line", { x1: m.l, x2: W - m.r, y1: y(v), y2: y(v), stroke: css("--line"), "stroke-width": 1 }, svg);
      const t = el("text", { x: m.l - 6, y: y(v) + 3, "text-anchor": "end", fill: css("--muted"), "font-size": 10, "font-family": "IBM Plex Mono, monospace" }, svg);
      t.textContent = v.toFixed(2);
    }
    for (const xv of opts.xticks) {
      const t = el("text", { x: x(xv), y: H - m.b + 14, "text-anchor": "middle", fill: css("--muted"), "font-size": 10, "font-family": "IBM Plex Mono, monospace" }, svg);
      t.textContent = opts.xfmt(xv);
    }
    const xl = el("text", { x: (m.l + W - m.r) / 2, y: H - 2, "text-anchor": "middle", fill: css("--muted"), "font-size": 10.5, "font-family": "Atkinson Hyperlegible, system-ui, sans-serif" }, svg);
    xl.textContent = opts.xlabel;
    if (opts.band) el("rect", { x: x(opts.band[0]), y: m.t, width: x(opts.band[1]) - x(opts.band[0]), height: H - m.t - m.b, fill: css("--accent-wash") }, svg);
    for (const s of series) {
      const color = s.emph ? css("--series-engine") : css("--series-status");
      el("path", { d: s.pts.map((p, i) => `${i ? "L" : "M"}${x(p[0])},${y(p[1])}`).join(""), fill: "none", stroke: color, "stroke-width": 2, "stroke-linejoin": "round" }, svg);
      for (const p of s.pts) el("circle", { cx: x(p[0]), cy: y(p[1]), r: 4, fill: color, stroke: css("--surface"), "stroke-width": 2 }, svg);
      const last = s.pts[s.pts.length - 1];
      const t = el("text", { x: x(last[0]) + 8, y: y(last[1]) + 4, fill: css("--ink-2"), "font-size": 11, "font-family": "Atkinson Hyperlegible, system-ui, sans-serif" }, svg);
      t.textContent = s.short;
    }
    const tip = document.createElement("div"); tip.className = "tip"; tip.hidden = true; box.appendChild(tip);
    const cross = el("line", { y1: m.t, y2: H - m.b, stroke: css("--muted"), "stroke-width": 1, visibility: "hidden" }, svg);
    const hit = el("rect", { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: "transparent" }, svg);
    hit.addEventListener("mousemove", (ev) => {
      const r = svg.getBoundingClientRect(), px = ev.clientX - r.left;
      let bi = 0;
      xs.forEach((xv, i) => { if (Math.abs(x(xv) - px) < Math.abs(x(xs[bi]) - px)) bi = i; });
      cross.setAttribute("x1", x(xs[bi])); cross.setAttribute("x2", x(xs[bi])); cross.setAttribute("visibility", "visible");
      tip.hidden = false;
      tip.textContent = `${opts.xfmt(xs[bi])}: ` + series.map(s => `${s.short} ${s.pts[bi][1].toFixed(3)}`).join(" · ");
      tip.style.left = `${Math.min(W - 250, Math.max(0, x(xs[bi]) + 10))}px`; tip.style.top = `${m.t}px`;
    });
    hit.addEventListener("mouseleave", () => { tip.hidden = true; cross.setAttribute("visibility", "hidden"); });
  }
  function table(headers, rows) {
    return `<details><summary>Table view</summary><div class="tbl-wrap"><table><thead><tr>${headers.map((h, i) => `<th class="${i ? "num" : ""}">${h}</th>`).join("")}</tr></thead><tbody>${rows.map(r => `<tr>${r.map((c, i) => `<td class="${i ? "num" : ""}">${c}</td>`).join("")}</tr>`).join("")}</tbody></table></div></details>`;
  }

  let simBuilt = false;
  function drawSim() {
    const S = window.SIM;
    if (!S) return;
    const grid = $("simGrid");
    if (!simBuilt) {
      $("simIntro").textContent = S.intro;
      grid.innerHTML = S.cards.map((c, i) => `<div class="sim-card"><h3>${c.title}</h3><p class="sub">${c.sub}</p><div class="chart-box" id="sim${i}"></div>${table(c.table.headers, c.table.rows)}</div>`).join("");
      simBuilt = true;
    }
    S.cards.forEach((c, i) => {
      const box = $(`sim${i}`);
      if (c.kind === "bars") barChart(box, c.rows, { fmt: (v) => c.pctFmt ? `${Math.round(v)}%` : v.toFixed(1), aria: c.title, labelW: c.labelW, max: c.max });
      else lineChart(box, c.series, { aria: c.title, xticks: c.xticks, xfmt: (v) => c.xunit === "h" ? `${v} h` : `${v}`, xlabel: c.xlabel, band: c.band });
    });
  }

  // boot
  setInputs(PRESETS[0].v);
  presetBox.children[0].setAttribute("aria-pressed", "true");
  renderHero();
  update();
  drawSim();
  const mq = window.matchMedia("(prefers-color-scheme: dark)");
  mq.addEventListener?.("change", () => { renderHero(); if (lastRes) { renderGauge(lastRes); renderCurve(lastRes); } drawSim(); });
  new MutationObserver(() => { if (lastRes) { renderGauge(lastRes); renderCurve(lastRes); } drawSim(); })
    .observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
})();

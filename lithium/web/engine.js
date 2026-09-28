// Lithos browser engine: a compact port of the Python prototype (lithos/)
// for the concept page. One occasion, time-invariant parameters, 1- or
// 2-compartment oral model, MAP + Laplace posterior, steady-state forecasts.
// Research prototype - not for clinical use.
"use strict";

const LithosEngine = (() => {
  const MMOL_PER_MG = 2 / 73.891; // lithium carbonate

  function crclCG(scr, age, sex, wt) {
    const mgdl = scr / 88.42;
    return ((140 - age) * wt) / (72 * mgdl) * (sex === "F" ? 0.85 : 1);
  }
  function egfr2021(scr, age, sex) {
    const s = scr / 88.42, f = sex === "F";
    const k = f ? 0.7 : 0.9, a = f ? -0.241 : -0.302, r = s / k;
    return 142 * Math.pow(Math.min(r, 1), a) * Math.pow(Math.max(r, 1), -1.2) * Math.pow(0.9938, age) * (f ? 1.012 : 1);
  }

  function ffm(wt, ht, sex) {
    const bmi = wt / ((ht / 100) ** 2);
    return sex === "F" ? 9270 * wt / (8780 + 244 * bmi) : 9270 * wt / (6680 + 216 * bmi);
  }
  function dosingWeight(wt, ht, sex) {
    const bmi = wt / ((ht / 100) ** 2);
    if (bmi <= 30) return wt;
    const ibw = (sex === "F" ? 45.5 : 50) + 0.9 * (ht - 152.4);
    return ibw + 0.4 * (wt - ibw);
  }
  // Priors mirrored from lithos/models.py (drift SD folded into CL for a single occasion).
  const MODELS = {
    renal: {
      name: "Kidney-function 2-compartment", weight: 0.75,
      typical(c) {
        const f = ffm(c.wt, c.ht, c.sex), vt = 0.85 * f;
        const crcl = crclCG(c.scr, c.age, c.sex, dosingWeight(c.wt, c.ht, c.sex));
        return { cl: 0.23 * 0.06 * crcl, v1: 0.55 * vt, v2: 0.45 * vt, q: 5.5 * Math.pow(f / 55, 0.75) * Math.pow(Math.max(c.age, 18) / 40, -1) };
      },
      omega: { cl: 0.25, v1: 0.20 }, ka: { IR: 1.2, SR: 0.40 }, f: { IR: 1, SR: 0.95 },
      sigmaProp: 0.10, sigmaAdd: 0.02, drift: Math.hypot(0.10, 0.08), distSd: 0.15, distTau: 3.0,
    },
    methaneethorn: {
      name: "Methaneethorn & Sringam 2019", weight: 0.25,
      // built in adults with acute mania, no kidney covariate: not used for older adults or reduced kidney function
      applies(c) { return egfr2021(c.scr, c.age, c.sex) >= 60 && c.age < 65; },
      typical(c) { return { cl: 1.43 * Math.pow(c.wt / 65, 0.425) * Math.pow(Math.max(c.age, 18) / 38, -0.242), v1: 54 }; },
      omega: { cl: 0.20, v1: 0.20 }, ka: { IR: 0.426, SR: 0.426 }, f: { IR: 1, SR: 0.95 },
      sigmaProp: 0.12, sigmaAdd: 0.02, drift: Math.hypot(0.10, 0.08), distSd: 0.5, distTau: 2.5,
    },
  };

  // --- structural model -----------------------------------------------------
  // p = {cl, v1, q, v2, ka, f}; returns macro-constants for superposition.
  function macro(p) {
    let ka = p.ka;
    if (!(p.q > 0 && p.v2 > 0)) {
      let k = p.cl / p.v1;
      if (Math.abs(ka - k) < 1e-6) ka *= 1.0001;
      // C = D f ka/(v1 (ka-k)) (e^-kt - e^-ka t)
      const c = p.f * ka / (p.v1 * (ka - k));
      return { lam: [k, ka], coef: [c, -c] };
    }
    const k10 = p.cl / p.v1, k12 = p.q / p.v1, k21 = p.q / p.v2;
    const s = k10 + k12 + k21, disc = Math.sqrt(s * s - 4 * k10 * k21);
    const al = 0.5 * (s + disc), be = 0.5 * (s - disc);
    if (Math.abs(ka - al) < 1e-6 || Math.abs(ka - be) < 1e-6) ka *= 1.0001;
    const pre = p.f * ka / p.v1;
    return {
      lam: [al, be, ka],
      coef: [pre * (k21 - al) / ((ka - al) * (be - al)),
             pre * (k21 - be) / ((ka - be) * (al - be)),
             pre * (k21 - ka) / ((al - ka) * (be - ka))],
    };
  }

  // Concentration (and slope) at time t from doses [{t, mmol}]
  function conc(p, doses, t) {
    const m = macro(p);
    let c = 0, dc = 0;
    for (const d of doses) {
      const tau = t - d.t;
      if (tau < 0) continue;
      for (let i = 0; i < m.lam.length; i++) {
        const e = d.mmol * m.coef[i] * Math.exp(-m.lam[i] * tau);
        c += e; dc -= m.lam[i] * e;
      }
    }
    return [c, dc];
  }

  // Steady state for a regimen repeating every 24 h: [{clock, mmol}]
  function steady(p, regimen, clock) {
    const m = macro(p);
    let c = 0;
    for (const r of regimen) {
      const tau = (((clock - r.clock) % 24) + 24) % 24;
      for (let i = 0; i < m.lam.length; i++) {
        c += r.mmol * m.coef[i] * Math.exp(-m.lam[i] * tau) / (1 - Math.exp(-24 * m.lam[i]));
      }
    }
    return c;
  }

  function halfLife(p) {
    if (!(p.q > 0 && p.v2 > 0)) return Math.LN2 * p.v1 / p.cl;
    const k10 = p.cl / p.v1, k12 = p.q / p.v1, k21 = p.q / p.v2, s = k10 + k12 + k21;
    return Math.LN2 / (0.5 * (s - Math.sqrt(s * s - 4 * k10 * k21)));
  }

  // --- statistics -----------------------------------------------------------
  function mulberry32(a) {
    return function () {
      a |= 0; a = (a + 0x6D2B79F5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function gauss(rng) {
    let u = 0, v = 0;
    while (u === 0) u = rng();
    v = rng();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }
  function nelderMead(f, x0, iters = 400) {
    const n = x0.length;
    let pts = [x0.slice()];
    for (let i = 0; i < n; i++) { const x = x0.slice(); x[i] += 0.3; pts.push(x); }
    let vals = pts.map(f);
    for (let it = 0; it < iters; it++) {
      const idx = vals.map((v, i) => i).sort((a, b) => vals[a] - vals[b]);
      pts = idx.map(i => pts[i]); vals = idx.map(i => vals[i]);
      if (Math.abs(vals[n] - vals[0]) < 1e-10) break;
      const cen = new Array(n).fill(0);
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) cen[j] += pts[i][j] / n;
      const along = (a) => cen.map((c, j) => c + a * (pts[n][j] - c));
      const xr = along(-1), fr = f(xr);
      if (fr < vals[0]) {
        const xe = along(-2), fe = f(xe);
        if (fe < fr) { pts[n] = xe; vals[n] = fe; } else { pts[n] = xr; vals[n] = fr; }
      } else if (fr < vals[n - 1]) { pts[n] = xr; vals[n] = fr; }
      else {
        const xc = along(fr < vals[n] ? -0.5 : 0.5), fc = f(xc);
        if (fc < Math.min(fr, vals[n])) { pts[n] = xc; vals[n] = fc; }
        else {
          for (let i = 1; i <= n; i++) { pts[i] = pts[i].map((v, j) => pts[0][j] + 0.5 * (v - pts[0][j])); vals[i] = f(pts[i]); }
        }
      }
    }
    return pts[0];
  }
  function hessian(f, x, h = 1e-3) {
    const n = x.length, H = [...Array(n)].map(() => new Array(n).fill(0)), f0 = f(x);
    const sh = (i, a, j, b) => { const y = x.slice(); y[i] += a; if (j !== undefined) y[j] += b; return y; };
    for (let i = 0; i < n; i++) {
      H[i][i] = (f(sh(i, h)) - 2 * f0 + f(sh(i, -h))) / (h * h);
      for (let j = i + 1; j < n; j++) {
        H[i][j] = H[j][i] = (f(sh(i, h, j, h)) - f(sh(i, h, j, -h)) - f(sh(i, -h, j, h)) + f(sh(i, -h, j, -h))) / (4 * h * h);
      }
    }
    return H;
  }
  function inv2(A) {
    const det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    return [[A[1][1] / det, -A[0][1] / det], [-A[1][0] / det, A[0][0] / det]];
  }
  function chol2(S) {
    const l11 = Math.sqrt(S[0][0]), l21 = S[1][0] / l11, l22 = Math.sqrt(Math.max(S[1][1] - l21 * l21, 1e-12));
    return [[l11, 0], [l21, l22]];
  }
  function quantile(xs, q) {
    const s = xs.slice().sort((a, b) => a - b), i = (s.length - 1) * q, lo = Math.floor(i), hi = Math.ceil(i);
    return s[lo] + (s[hi] - s[lo]) * (i - lo);
  }

  // --- the model + fit --------------------------------------------------------
  // model: {typical(cov)->{cl,v1,q,v2}, omega:{cl,v1}, ka:{IR,SR}, f:{IR,SR}, sigmaProp, sigmaAdd,
  //         drift: total SD of log-CL drift (added to the prior at a single occasion), distSd, distTau}
  function individual(model, cov, eta, form) {
    const tv = model.typical(cov);
    return { cl: tv.cl * Math.exp(eta[0]), v1: tv.v1 * Math.exp(eta[1]), q: tv.q || 0,
             v2: (tv.v2 || 0) * Math.exp(eta[1]), ka: model.ka[form], f: model.f[form] };
  }

  // scen: {cov, form, regimen:[{clock, mg}], days (on regimen), level:{hoursAfterLast, value, timingSd}}
  function buildDoses(scen) {
    const doses = [];
    const days = Math.max(1, scen.days);
    for (let d = 0; d < days + 1; d++) {
      for (const r of scen.regimen) doses.push({ t: d * 24 + r.clock, mmol: r.mg * MMOL_PER_MG });
    }
    doses.sort((a, b) => a.t - b.t);
    return doses;
  }

  function fit(model, scen) {
    const doses = buildDoses(scen);
    // Sample time: hours after the last dose that precedes it, on the final day.
    const lastDay = Math.max(1, scen.days) - 1;
    const clocks = scen.regimen.map(r => r.clock).sort((a, b) => a - b);
    const lastClock = scen.lastClock !== undefined ? scen.lastClock : clocks[clocks.length - 1];
    const tLast = lastDay * 24 + lastClock;
    const tObs = tLast + scen.level.hoursAfterLast;
    const used = doses.filter(d => d.t <= tObs + 1e-9);
    const om = [model.omega.cl, model.omega.v1];
    const dSd = model.drift || 0;
    const prior = [om[0] * om[0] + dSd * dSd, om[1] * om[1]];
    const y = scen.level.value, tsd = scen.level.timingSd ?? 0.25;
    const since = scen.level.hoursAfterLast;
    const distSd = model.distSd * Math.exp(-since / model.distTau);
    const LOG2PI = Math.log(2 * Math.PI);
    const ofv = (u) => {
      const p = individual(model, scen.cov, u, scen.form);
      const [c, dc] = conc(p, used, tObs);
      const v = model.sigmaAdd ** 2 + (model.sigmaProp ** 2 + distSd ** 2) * c * c + (dc * tsd) ** 2;
      return (y - c) ** 2 / v + Math.log(v) + LOG2PI + u[0] * u[0] / prior[0] + u[1] * u[1] / prior[1]
        + Math.log(2 * Math.PI * prior[0]) + Math.log(2 * Math.PI * prior[1]);
    };
    const mode = nelderMead(ofv, [0, 0]);
    const H = hessian(ofv, mode);
    const A = [[H[0][0] / 2, H[0][1] / 2], [H[1][0] / 2, H[1][1] / 2]];
    let cov = inv2(A);
    if (!(cov[0][0] > 0 && cov[1][1] > 0)) cov = [[prior[0], 0], [0, prior[1]]];
    const detA = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    const logEvidence = -0.5 * ofv(mode) + LOG2PI - 0.5 * Math.log(Math.max(detA, 1e-300));
    const [pFit, slope] = conc(individual(model, scen.cov, mode, scen.form), used, tObs);
    return { model, scen, doses, used, tObs, mode, cov, fitted: pFit, slope, logEvidence };
  }

  function priorOnly(model, scen) {
    const om = [model.omega.cl, model.omega.v1], dSd = model.drift || 0;
    return { model, scen, doses: buildDoses(scen), mode: [0, 0],
             cov: [[om[0] ** 2 + dSd * dSd, 0], [0, om[1] ** 2]], prior: true };
  }

  function draws(post, n = 800, seed = 7) {
    const rng = mulberry32(seed), L = chol2(post.cov), out = [];
    for (let i = 0; i < n; i++) {
      const z0 = gauss(rng), z1 = gauss(rng);
      out.push([post.mode[0] + L[0][0] * z0, post.mode[1] + L[1][0] * z0 + L[1][1] * z1]);
    }
    return out;
  }

  // Standardised 12-h level on a regimen [{clock, mg}]
  function li12(post, regimen, n = 800) {
    const clock = (Math.max(...regimen.map(r => r.clock)) + 12) % 24;
    const reg = regimen.map(r => ({ clock: r.clock, mmol: r.mg * MMOL_PER_MG }));
    const vals = draws(post, n).map(u => steady(individual(post.model, post.scen.cov, u, post.scen.form), reg, clock));
    return summary(vals);
  }
  function summary(vals) {
    return { median: quantile(vals, 0.5), lo: quantile(vals, 0.05), hi: quantile(vals, 0.95), vals };
  }

  function recommend(post, opts) {
    const { target = [0.6, 0.8], step = 250, maxMg = 2000, template = [21], pHighMax = 0.05, high = 1.0 } = opts;
    const clock = (Math.max(...template) + 12) % 24;
    const U = draws(post, 800, 11);
    const unit = template.map(c => U.map(u => steady(individual(post.model, post.scen.cov, u, post.scen.form),
      [{ clock: c, mmol: MMOL_PER_MG }], clock)));
    const options = [];
    for (let k = 1; k * step <= maxMg; k++) {
      let split;
      if (template.length === 1) split = [k * step];
      else {
        const base = Math.floor(k / template.length), rem = k % template.length;
        if (base === 0) continue;
        split = template.map(() => base * step);
        for (let r = 0; r < rem; r++) split[split.length - 1 - r] += step;
      }
      const vals = U.map((_, i) => split.reduce((s, mg, j) => s + mg * unit[j][i], 0));
      const s = summary(vals);
      s.pIn = vals.filter(v => v >= target[0] && v <= target[1]).length / vals.length;
      s.pHigh = vals.filter(v => v > high).length / vals.length;
      s.daily = split.reduce((a, b) => a + b, 0);
      s.split = split;
      options.push(s);
    }
    const safe = options.filter(o => o.pHigh <= pHighMax);
    const center = (target[0] + target[1]) / 2;
    const pool = safe.length ? safe : options;
    pool.sort((a, b) => (Math.round(b.pIn * 100) - Math.round(a.pIn * 100)) || (Math.abs(a.median - center) - Math.abs(b.median - center)));
    const hl = quantile(U.slice(0, 300).map(u => halfLife(individual(post.model, post.scen.cov, u, post.scen.form))), 0.5);
    return { best: pool[0], options, halfLifeH: hl };
  }

  // Curve with 90% band over the history window + the observed point
  function curve(post, tFrom, tTo, npts = 160, n = 200) {
    const ts = [...Array(npts)].map((_, i) => tFrom + (tTo - tFrom) * i / (npts - 1));
    const U = draws(post, n, 3);
    const doses = post.doses.filter(d => d.t <= tTo);
    const rows = U.map(u => { const p = individual(post.model, post.scen.cov, u, post.scen.form); return ts.map(t => conc(p, doses, t)[0]); });
    return ts.map((t, j) => {
      const col = rows.map(r => r[j]);
      return { t, med: quantile(col, 0.5), lo: quantile(col, 0.05), hi: quantile(col, 0.95) };
    });
  }

  // Ensemble: fit each prior, weight by prior weight x evidence; draws are pooled in proportion.
  function fitEnsemble(scen, withLevel = true) {
    const usable = Object.values(MODELS).filter(m => !m.applies || m.applies(scen.cov));
    const members = usable.map(m => withLevel ? fit(m, scen) : priorOnly(m, scen));
    const logw = members.map(p => Math.log(p.model.weight) + (withLevel ? p.logEvidence : 0));
    const mx = Math.max(...logw), w = logw.map(v => Math.exp(v - mx)), tot = w.reduce((a, b) => a + b, 0);
    return { members, weights: w.map(v => v / tot), scen, ensemble: true };
  }
  function split(n, weights) {
    const k = weights.map(w => Math.floor(w * n));
    k[weights.indexOf(Math.max(...weights))] += n - k.reduce((a, b) => a + b, 0);
    return k;
  }
  function ensLi12(ens, regimen, n = 800) {
    const k = split(n, ens.weights);
    const vals = [];
    ens.members.forEach((m, i) => { if (k[i] > 0) vals.push(...li12(m, regimen, k[i]).vals); });
    return summary(vals);
  }
  function ensRecommend(ens, opts) {
    // pool unit contributions across members
    const k = split(800, ens.weights);
    const pooled = [];
    ens.members.forEach((m, i) => { if (k[i] > 0) pooled.push({ m, n: k[i] }); });
    const template = opts.template || [21];
    const clock = (Math.max(...template) + 12) % 24;
    const unit = template.map(() => []);
    let hl = [];
    pooled.forEach(({ m, n }, idx) => {
      const U = draws(m, n, 11 + idx);
      template.forEach((c, j) => U.forEach(u => unit[j].push(steady(individual(m.model, m.scen.cov, u, m.scen.form), [{ clock: c, mmol: MMOL_PER_MG }], clock))));
      hl.push(...U.slice(0, 150).map(u => halfLife(individual(m.model, m.scen.cov, u, m.scen.form))));
    });
    const N = unit[0].length, target = opts.target || [0.6, 0.8], step = opts.step || 250, maxMg = opts.maxMg || 2000;
    const high = opts.high ?? (target[1] + 0.2);   // safety ceiling follows the target, as in lithos/dosing.py
    const options = [];
    for (let kk = 1; kk * step <= maxMg; kk++) {
      let splitMg;
      if (template.length === 1) splitMg = [kk * step];
      else {
        const base = Math.floor(kk / template.length), rem = kk % template.length;
        if (base === 0) continue;
        splitMg = template.map(() => base * step);
        for (let r = 0; r < rem; r++) splitMg[splitMg.length - 1 - r] += step;
      }
      const vals = [];
      for (let i = 0; i < N; i++) vals.push(splitMg.reduce((s, mg, j) => s + mg * unit[j][i], 0));
      const sm = summary(vals);
      sm.pIn = vals.filter(v => v >= target[0] && v <= target[1]).length / N;
      sm.pHigh = vals.filter(v => v > high).length / N;
      sm.daily = splitMg.reduce((a, b) => a + b, 0);
      sm.split = splitMg;
      options.push(sm);
    }
    const center = (target[0] + target[1]) / 2;
    const safe = options.filter(o => o.pHigh <= 0.05);
    const pool = (safe.length ? safe : options).slice();
    pool.sort((a, b) => (Math.round(b.pIn * 100) - Math.round(a.pIn * 100)) || (Math.abs(a.median - center) - Math.abs(b.median - center)));
    return { best: pool[0], options, halfLifeH: quantile(hl, 0.5), high };
  }
  function ensCurve(ens, tFrom, tTo, npts = 140) {
    const k = split(240, ens.weights);
    const ts = [...Array(npts)].map((_, i) => tFrom + (tTo - tFrom) * i / (npts - 1));
    const rows = [];
    ens.members.forEach((m, i) => {
      if (!k[i]) return;
      const doses = m.doses.filter(d => d.t <= tTo);
      draws(m, k[i], 3 + i).forEach(u => {
        const p = individual(m.model, m.scen.cov, u, m.scen.form);
        rows.push(ts.map(t => conc(p, doses, t)[0]));
      });
    });
    return ts.map((t, j) => {
      const col = rows.map(r => r[j]);
      return { t, med: quantile(col, 0.5), lo: quantile(col, 0.05), hi: quantile(col, 0.95) };
    });
  }

  return { crclCG, egfr2021, ffm, conc, steady, fit, priorOnly, li12, recommend, curve, halfLife, individual,
           MMOL_PER_MG, MODELS, fitEnsemble, ensLi12, ensRecommend, ensCurve, buildDoses };
})();

if (typeof module !== "undefined") module.exports = LithosEngine;

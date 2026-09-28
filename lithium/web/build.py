"""Build the concept page (docs/lithos-concept.html) from the page sources, the
simulation results and the published-cohort check.

    cd lithium && python -m sim.run_all && python web/build.py

The page is a single self-contained HTML body (the artifact host adds the
document skeleton). The browser engine in engine.js is a compact port of
lithos/ for the interactive panel; the Python package remains the reference.
"""

import html
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from sim.headline import load, web  # noqa: E402
from validation.literature import run  # noqa: E402


def real_patients_section() -> str:
    rows = run()
    inside = sum(r.verdict == "inside" for r in rows)
    within20 = sum(abs(r.vs_mean) <= 0.20 for r in rows)
    trs = []
    for r in rows:
        med, lo, hi = r.engine
        d = 2 if r.unit in ("mmol/L", "L/kg", "ratio") else 1
        pill = '<span class="pill good">inside</span>' if r.verdict == "inside" else f'<span class="pill warn">{html.escape(r.verdict)}</span>'
        cohort = html.escape(r.cohort.split(" (")[0])
        trs.append(f"<tr><td>{cohort}</td><td>{html.escape(r.quantity)}</td><td>{html.escape(r.published)} {html.escape(r.unit)}</td>"
                   f"<td>{med:.{d}f} <span class=\"muted\">({lo:.{d}f}–{hi:.{d}f})</span></td><td>{r.vs_mean:+.0%}</td><td>{pill}</td></tr>")
    section = f'''<section aria-labelledby="realTitle">
        <div class="section-head">
          <p class="eyebrow">Real patients, so far</p>
          <h2 id="realTitle">Checked against eight published patient cohorts</h2>
          <p class="prose">No individual-level lithium dataset could be reached from where this was built. As a first check, the engine was given a patient resembling each published cohort, with no lithium levels, and its predictions were compared with what the papers measured in real people.</p>
        </div>
        <div class="real-grid">
          <div class="stat"><b>{inside} of {len(rows)}</b><span>quantities inside the published mean ± 1 SD (or reported range)</span></div>
          <div class="stat"><b>{within20} of {len(rows)}</b><span>within 20% of the published mean, a stricter test</span></div>
          <div class="stat"><b>Slow for older adults</b><span>half-life +25%, peak +31% vs one older cohort: first doses err low, and the first level corrects them</span></div>
        </div>
        <details>
          <summary>All {len(rows)} comparisons</summary>
          <div class="tbl-wrap">
            <table class="cohorts">
              <thead><tr><th>Cohort</th><th>Quantity</th><th>Published</th><th>Engine prior, median (90% range)</th><th>vs published mean</th><th>In band?</th></tr></thead>
              <tbody>
                {"".join(trs)}
              </tbody>
            </table>
          </div>
        </details>
        <p class="caveat">This is a consistency check, not validation. Some of these cohorts informed the calibration, most papers did not report covariates, and agreement on averages says nothing about forecasting an individual. The individual-level test (forecast each patient's next level from their earlier ones) is built and ready to run where credentialed data live: MIMIC-IV hospital records, trial data such as Bipolar CHOICE, or collaborators such as the eLi12 group in Denmark.</p>
      </section>'''
    return section


def main() -> None:
    here = Path(__file__).resolve().parent
    shell = (here / "page_shell.html").read_text()
    shell = shell.replace("<!-- REAL-PATIENTS -->", real_patients_section())
    engine = (here / "engine.js").read_text()
    engine = engine.replace('if (typeof module !== "undefined") module.exports = LithosEngine;', "")  # Node-only hook
    ui = (here / "ui.js").read_text()
    meta, rows = load()
    sim = web(meta, rows)
    page = (shell + "\n<script>\n" + engine + "\n</script>\n<script>\nwindow.SIM = " + json.dumps(sim, ensure_ascii=False)
            + ";\n</script>\n<script>\n" + ui + "\n</script>\n")
    assert not re.search(r"<(html|head|body)[\s>]", page, re.I), "the host adds the document skeleton"
    out = ROOT / "docs" / "lithos-concept.html"
    out.write_text(page)
    print(out, len(page) // 1024, "KB")


if __name__ == "__main__":
    main()

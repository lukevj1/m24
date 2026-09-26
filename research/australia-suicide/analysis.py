"""Exploratory correlation analysis of suicide rates in Australia.

Run `python fetch_data.py` first, then `python analysis.py`.

Three analyses:
  A. Australia over time (WHO GHO, 2000-2021): correlation of the national
     suicide rate with each covariate, on raw levels, on linearly detrended
     series and on year-to-year changes.
  B. Across countries (average of 2015-2019, pre-COVID): Spearman correlation
     between suicide rate and each covariate, for high-income countries and
     for all countries, plus where Australia sits.
  C. Australia by sex, age-standardised (IHME GBD 2019, 1990-2017).

Outputs go to outputs/ (CSV tables + PNG charts). These are ecological
correlations; see REPORT.md for how (not) to interpret them.
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from fetch_data import INDICATORS

HERE = pathlib.Path(__file__).parent
RAW = HERE / "data" / "raw"
OUT = HERE / "outputs"
OUT.mkdir(exist_ok=True)

OUTCOME = "suicide_per_100000_people"
SUICIDE_COLS = [i for i in INDICATORS if i.startswith("suicide_")]
COVARIATES = [i for i in INDICATORS if not i.startswith("suicide_")]

# Validated reference palette (dataviz skill), light mode.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "axes.edgecolor": INK2, "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.8, "axes.spines.top": False,
    "axes.spines.right": False, "lines.linewidth": 2, "font.size": 10,
    "axes.titleweight": "bold", "axes.titlesize": 11,
})


def label(concept):
    names = pd.read_csv(RAW / "concepts.csv").set_index("concept")["name"]
    return names.get(concept, concept)


def load_panel():
    """Long country x year table with one column per indicator."""
    frames = []
    for ind in INDICATORS:
        d = pd.read_csv(RAW / f"{ind}.csv")
        frames.append(d.set_index(["geo", "time"])[ind])
    return pd.concat(frames, axis=1).reset_index()


def corr_row(x, y, method="pearson"):
    ok = x.notna() & y.notna()
    n = int(ok.sum())
    if n < 8:
        return n, np.nan, np.nan
    f = stats.pearsonr if method == "pearson" else stats.spearmanr
    r, p = f(x[ok], y[ok])
    return n, r, p


def add_fdr(df, pcol, qcol):
    ok = df[pcol].notna()
    df[qcol] = np.nan
    df.loc[ok, qcol] = multipletests(df.loc[ok, pcol], method="fdr_bh")[1]
    return df


def detrend(s):
    ok = s.notna()
    t = s.index[ok].astype(float)
    b = np.polyfit(t, s[ok], 1)
    out = s.copy()
    out[ok] = s[ok] - np.polyval(b, t)
    return out


# --------------------------------------------------------------------------
# A. Australia over time
# --------------------------------------------------------------------------
def australia_time_series(panel):
    aus = panel[panel.geo == "aus"].set_index("time").sort_index()
    aus = aus.loc[2000:2021]
    aus.to_csv(OUT / "australia_2000_2021.csv")
    y = aus[OUTCOME]
    rows = []
    for c in COVARIATES:
        x = aus[c]
        n, r, p = corr_row(x, y)
        _, r_dt, p_dt = corr_row(detrend(x), detrend(y)) if n >= 8 else (n, np.nan, np.nan)
        n_d, r_d, p_d = corr_row(x.diff(), y.diff())
        rows.append(dict(indicator=c, name=label(c), n_years=n,
                         r_levels=r, p_levels=p, r_detrended=r_dt,
                         p_detrended=p_dt, n_diffs=n_d, r_changes=r_d,
                         p_changes=p_d))
    res = pd.DataFrame(rows)
    res = res[res.n_years >= 12]
    for k in ["levels", "detrended", "changes"]:
        res = add_fdr(res, f"p_{k}", f"q_{k}")
    res = res.sort_values("r_detrended", key=abs, ascending=False)
    res.round(3).to_csv(OUT / "A_australia_timeseries_correlations.csv", index=False)

    # Chart: suicide trend by sex
    fig, ax = plt.subplots(figsize=(7.5, 4))
    for col, colr, nm in [("suicide_men_per_100000_people", BLUE, "Men"),
                          (OUTCOME, INK2, "All"),
                          ("suicide_women_per_100000_people", ORANGE, "Women")]:
        ax.plot(aus.index, aus[col], color=colr, marker="o", ms=4,
                ls="--" if nm == "All" else "-", label=nm)
        ax.annotate(nm, (aus.index[-1], aus[col].iloc[-1]), xytext=(6, 0),
                    textcoords="offset points", va="center", color=INK2)
    ax.set_title("Suicide deaths per 100,000, Australia (WHO GHO)")
    ax.set_ylim(0)
    ax.set_xlim(1999.5, 2023)
    ax.legend(frameon=False, loc="lower left", ncol=3)
    fig.tight_layout()
    fig.savefig(OUT / "A1_australia_trend.png", dpi=150)
    plt.close(fig)

    # Chart: four key covariates vs suicide over time (small multiples, own axes)
    keys = ["aged_15plus_unemployment_rate_percent", "gdp_per_capita_yearly_growth",
            "poisonings_deaths_per_100000_people", "total_health_spending_percent_of_gdp"]
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.4))
    for ax, c in zip(axes, keys):
        ok = aus[[c, OUTCOME]].dropna()
        ax.scatter(ok[c], ok[OUTCOME], s=36, color=BLUE, edgecolor="#fcfcfb", lw=1.5)
        for yr in [2000, 2010, 2021]:
            if yr in ok.index:
                ax.annotate(str(yr), (ok.loc[yr, c], ok.loc[yr, OUTCOME]),
                            xytext=(4, 4), textcoords="offset points", fontsize=8, color=INK2)
        r = res.set_index("indicator").loc[c]
        ax.set_title(f"{label(c)}\nr={r.r_levels:.2f}, detrended r={r.r_detrended:.2f}",
                     fontsize=9)
        ax.set_xlabel(label(c), fontsize=8)
    axes[0].set_ylabel("Suicides per 100k")
    fig.suptitle("Australia 2000-2021: suicide rate vs selected factors (each dot = one year)",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / "A2_australia_scatter.png", dpi=150)
    plt.close(fig)
    return res


# --------------------------------------------------------------------------
# B. Cross-country
# --------------------------------------------------------------------------
def cross_country(panel):
    geo = pd.read_csv(RAW / "countries.csv").set_index("country")
    avg = (panel[panel.time.between(2015, 2019)]
           .groupby("geo")[INDICATORS].mean())
    # covariates with sparse recent coverage: fall back to 2005-2019 average
    wide = panel[panel.time.between(2005, 2019)].groupby("geo")[INDICATORS].mean()
    avg = avg.fillna(wide)
    avg = avg[avg.index.isin(geo.index) & avg[OUTCOME].notna()]
    avg["name"] = geo.loc[avg.index, "name"]
    avg["income"] = geo.loc[avg.index, "income_groups"]
    avg.to_csv(OUT / "country_averages_2015_2019.csv")

    hi = avg[avg.income == "high_income"]
    rows = []
    for c in COVARIATES:
        for grp, d in [("high_income", hi), ("all", avg)]:
            n, rho, p = corr_row(d[c], d[OUTCOME], "spearman")
            pct = np.nan
            if "aus" in d.index and pd.notna(d.loc["aus", c]):
                pct = (d[c].dropna() < d.loc["aus", c]).mean() * 100
            rows.append(dict(indicator=c, name=label(c), group=grp, n_countries=n,
                             spearman_rho=rho, p=p, australia_value=d[c].get("aus"),
                             australia_percentile=pct))
    res = pd.DataFrame(rows)
    res = pd.concat([add_fdr(g.copy(), "p", "q_fdr") for _, g in res.groupby("group")])
    res["_abs"] = res.spearman_rho.abs()
    res = res.sort_values(["group", "_abs"], ascending=[False, False]).drop(columns="_abs")
    res.round(3).to_csv(OUT / "B_cross_country_correlations.csv", index=False)

    # Multivariable model among high-income countries (standardised betas)
    preds = ["aged_15plus_unemployment_rate_percent", "inequality_index_gini",
             "poisonings_deaths_per_100000_people", "urban_population_percent_of_total",
             "total_health_spending_percent_of_gdp", "median_age_years"]
    import statsmodels.api as sm
    d = hi[[OUTCOME] + preds].dropna()
    z = (d - d.mean()) / d.std()
    m = sm.OLS(z[OUTCOME], sm.add_constant(z[preds])).fit(cov_type="HC3")
    (OUT / "B_multivariable_high_income.txt").write_text(
        f"Standardised OLS, high-income countries, 2015-2019 averages (n={len(d)})\n"
        f"Outcome: {OUTCOME}\n\n{m.summary()}\n")

    # Chart: rho bars for high-income countries
    h = res[(res.group == "high_income") & res.spearman_rho.notna()].copy()
    h = h.sort_values("spearman_rho")
    fig, ax = plt.subplots(figsize=(8, 7.5))
    colors = [BLUE if v > 0 else ORANGE for v in h.spearman_rho]
    ax.barh(h.name, h.spearman_rho, color=colors, height=0.7)
    for i, (rho, q) in enumerate(zip(h.spearman_rho, h.q_fdr)):
        ax.text(rho + (0.02 if rho > 0 else -0.02), i, f"{rho:.2f}{' *' if q < 0.05 else ''}",
                va="center", ha="left" if rho > 0 else "right", fontsize=8, color=INK2)
    ax.axvline(0, color=INK2, lw=1)
    ax.set_xlim(-0.8, 0.8)
    ax.set_xlabel("Spearman correlation with suicide rate  (* = FDR q < 0.05)")
    fig.suptitle(f"High-income countries (n~{int(h.n_countries.median())}): "
                 "correlates of suicide rate", fontweight="bold")
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    fig.savefig(OUT / "B1_cross_country_rho.png", dpi=150)
    plt.close(fig)

    # Chart: scatter for top 3 by |rho| among high-income, Australia highlighted
    top = h.reindex(h.spearman_rho.abs().sort_values(ascending=False).index).head(3)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, (_, r) in zip(axes, top.iterrows()):
        d = hi[[r.indicator, OUTCOME, "name"]].dropna()
        ax.scatter(d[r.indicator], d[OUTCOME], s=30, color=GRID, edgecolor=INK2, lw=0.6)
        if "aus" in d.index:
            ax.scatter(d.loc["aus", r.indicator], d.loc["aus", OUTCOME], s=80,
                       color=BLUE, edgecolor="#fcfcfb", lw=2, zorder=3)
            ax.annotate("Australia", (d.loc["aus", r.indicator], d.loc["aus", OUTCOME]),
                        xytext=(6, 6), textcoords="offset points", color=INK, fontweight="bold")
        for g in d[OUTCOME].nlargest(3).index:
            ax.annotate(d.loc[g, "name"], (d.loc[g, r.indicator], d.loc[g, OUTCOME]),
                        xytext=(4, -10), textcoords="offset points", fontsize=7, color=INK2)
        ax.set_xlabel(r["name"], fontsize=9)
        ax.set_title(f"rho = {r.spearman_rho:.2f} (n={int(r.n_countries)})", fontsize=10)
    axes[0].set_ylabel("Suicides per 100k (2015-19 avg)")
    fig.suptitle("High-income countries: strongest cross-country associations",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / "B2_cross_country_scatter.png", dpi=150)
    plt.close(fig)
    return res, m, len(d)


# --------------------------------------------------------------------------
# C. Australia by sex, age-standardised (IHME GBD 2019, 1990-2017)
# --------------------------------------------------------------------------
# NOTE: the per-age-group columns in this OWID file are internally inconsistent
# (e.g. Australian males 15-19 jump 18 -> 11 -> 42 -> 8 per 100k between
# adjacent 3-year points), so only the age-standardised series are used.
def age_sex():
    d = pd.read_csv(RAW / "ihme_suicide_by_sex_age.csv")
    cols = ["Male suicide rate (age-standardized)", "Female suicide rate (age-standardized)",
            "Male:female suicide ratio"]
    aus = d[d.Entity == "Australia"].set_index("Year").sort_index()[cols]
    aus.to_csv(OUT / "australia_ihme_age_standardised_1990_2017.csv")
    fig, ax = plt.subplots(figsize=(7.5, 4))
    for col, colr, nm in [(cols[0], BLUE, "Men"), (cols[1], ORANGE, "Women")]:
        ax.plot(aus.index, aus[col], color=colr, label=nm)
        ax.annotate(nm, (aus.index[-1], aus[col].iloc[-1]), xytext=(6, 0),
                    textcoords="offset points", va="center", color=INK2)
    ax.set_ylim(0)
    ax.set_xlim(1989, 2020)
    ax.legend(frameon=False, loc="center right")
    ax.set_title("Age-standardised suicide rate per 100,000, Australia (IHME GBD 2019)")
    fig.tight_layout()
    fig.savefig(OUT / "C1_australia_age_standardised.png", dpi=150)
    plt.close(fig)
    return aus


def main():
    panel = load_panel()
    a = australia_time_series(panel)
    b, model, n_model = cross_country(panel)
    c = age_sex()
    pd.set_option("display.width", 200)
    print("\n== A. Australia over time ==")
    print(a[["name", "n_years", "r_levels", "q_levels", "r_detrended", "q_detrended",
             "r_changes", "q_changes"]].round(2).to_string(index=False))
    print("\n== B. Cross-country (high income) ==")
    print(b[b.group == "high_income"][["name", "n_countries", "spearman_rho", "q_fdr",
                                        "australia_value", "australia_percentile"]]
          .round(2).to_string(index=False))
    print("\n== B. Cross-country (all) ==")
    print(b[b.group == "all"][["name", "n_countries", "spearman_rho", "q_fdr"]]
          .round(2).to_string(index=False))
    print(f"\n== Multivariable (n={n_model}) ==")
    print(model.summary().tables[1])
    print("\n== C. IHME AUS male/female age-standardised ==")
    print(c.iloc[[0, 10, 20, -1]].round(1))


if __name__ == "__main__":
    main()

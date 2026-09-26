"""Download public suicide-rate and covariate data used by analysis.py.

Sources (all public, mirrored on GitHub so they are reachable from restricted
networks):
  * Gapminder "Systema Globalis" (open-numbers/ddf--gapminder--systema_globalis)
    - country x year indicators compiled from WHO, IHME GBD, World Bank, ILO,
      OECD and others. See data/SOURCES.md for per-indicator provenance.
  * Our World in Data "owid-datasets" - IHME GBD 2019 suicide rates by sex/age.
"""
import pathlib
import urllib.parse
import urllib.request

RAW = pathlib.Path(__file__).parent / "data" / "raw"
SG = ("https://raw.githubusercontent.com/open-numbers/"
      "ddf--gapminder--systema_globalis/master/countries-etc-datapoints/"
      "ddf--datapoints--{}--by--geo--time.csv")
OWID = ("https://raw.githubusercontent.com/owid/owid-datasets/master/datasets/"
        "Suicide%20rates%20by%20sex%20and%20age%20(IHME,%202019)/"
        "Suicide%20rates%20by%20sex%20and%20age%20(IHME,%202019).csv")

INDICATORS = [
    # outcome
    "suicide_per_100000_people",
    "suicide_men_per_100000_people",
    "suicide_women_per_100000_people",
    # economy / work
    "aged_15plus_unemployment_rate_percent",
    "aged_15_24_unemployment_rate_percent",
    "males_aged_15plus_unemployment_rate_percent",
    "long_term_unemployment_rate_percent",
    "gdppercapita_us_inflation_adjusted",
    "gdp_per_capita_yearly_growth",
    "inflation_annual_percent",
    "inequality_index_gini",
    "income_share_of_richest_10percent",
    "working_hours_per_week",
    "aged_15plus_labour_force_participation_rate_percent",
    # substances / physical health
    "alcohol_consumption_per_adult_15plus_litres",
    "smoking_adults_percent_of_population_over_age_15",
    "body_mass_index_bmi_men_kgperm2",
    "poisonings_deaths_per_100000_people",
    "life_expectancy_at_birth_data_from_ihme",
    # health system
    "total_health_spending_percent_of_gdp",
    "total_health_spending_per_person_international_dollar",
    "government_share_of_total_health_spending_percent",
    "out_of_pocket_share_of_total_health_spending_percent",
    "medical_doctors_per_1000_people",
    # social / demographic
    "urban_population_percent_of_total",
    "population_density_per_square_km",
    "median_age_years",
    "age_at_1st_marriage_women",
    "teen_fertility_rate_births_per_1000_women_ages_15_19",
    "murder_per_100000_people",
    "hdi_human_development_index",
    "democracy_score_use_as_color",
]


def get(url, dest):
    if dest.exists():
        return
    print("fetching", dest.name)
    with urllib.request.urlopen(url, timeout=60) as r:
        dest.write_bytes(r.read())


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    for ind in INDICATORS:
        get(SG.format(ind), RAW / f"{ind}.csv")
    get(OWID, RAW / "ihme_suicide_by_sex_age.csv")
    get(SG.replace("countries-etc-datapoints/ddf--datapoints--{}--by--geo--time.csv",
                   "ddf--entities--geo--country.csv"), RAW / "countries.csv")
    get(SG.replace("countries-etc-datapoints/ddf--datapoints--{}--by--geo--time.csv",
                   "ddf--concepts.csv"), RAW / "concepts.csv")


if __name__ == "__main__":
    main()

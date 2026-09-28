"""Lithium salt conversions and a small formulation catalogue.

Every calculation inside the engine works in mmol of lithium ion (Li+), mmol/L
for serum concentrations (numerically identical to mEq/L for a monovalent ion)
and hours for time. Prescribers think in milligrams of a *salt*, so this module
is the only place where milligrams appear.
"""

from __future__ import annotations

from dataclasses import dataclass

# Molar masses (g/mol).
_LI2CO3 = 73.891          # lithium carbonate, 2 Li+ per formula unit
_LI_CITRATE_4H2O = 282.0  # trilithium citrate tetrahydrate, 3 Li+ per formula unit

MMOL_LI_PER_MG = {
    "carbonate": 2.0 / _LI2CO3,          # 0.02707 mmol Li+ per mg (250 mg -> 6.77 mmol)
    "citrate": 3.0 / _LI_CITRATE_4H2O,   # 0.01064 mmol Li+ per mg (520 mg -> 5.53 mmol)
}


def mg_to_mmol(mg: float, salt: str = "carbonate") -> float:
    """Convert milligrams of a lithium salt to mmol of Li+."""
    return mg * MMOL_LI_PER_MG[salt]


def mmol_to_mg(mmol: float, salt: str = "carbonate") -> float:
    """Convert mmol of Li+ to milligrams of a lithium salt."""
    return mmol / MMOL_LI_PER_MG[salt]


@dataclass(frozen=True)
class Product:
    """A marketed lithium product. ``release`` selects the absorption model."""

    name: str
    strength_mg: float
    salt: str = "carbonate"
    release: str = "IR"          # "IR" (immediate), "SR" (sustained/prolonged/modified), "LIQ" (solution)
    scored: bool = False         # may be halved (only set where confirmed)
    per_ml: float | None = None  # liquids: strength is per this many mL
    mmol_label: float | None = None  # mmol Li+ stated on the pack (authoritative for citrate products)
    market: str = ""
    verified: bool = True        # False where availability came from non-regulatory sources

    @property
    def mmol(self) -> float:
        return self.mmol_label if self.mmol_label is not None else mg_to_mmol(self.strength_mg, self.salt)


# Illustrative catalogue (docs/EVIDENCE.md). Brands change; confirm locally
# (TGA ARTG, dm+d, FDA) before relying on it. Never switch brands without a
# plan: MR products are not interchangeable mg-for-mg.
PRODUCTS: dict[str, Product] = {
    p.name: p
    for p in [
        # Australia (from a dispensing dataset; confirm against the ARTG)
        Product("Lithicarb 250 mg", 250, release="IR", market="AU", verified=False),
        Product("Quilonum SR 450 mg", 450, release="SR", market="AU", verified=False),
        # United Kingdom (dm+d-derived code lists)
        Product("Camcolit 250 mg", 250, release="IR", market="UK"),
        Product("Camcolit 400 mg MR", 400, release="SR", market="UK"),
        Product("Priadel 200 mg MR", 200, release="SR", market="UK"),
        Product("Priadel 400 mg MR", 400, release="SR", market="UK"),
        Product("Liskonum 450 mg MR", 450, release="SR", market="UK"),
        Product("Priadel liquid 520 mg/5 mL", 520, salt="citrate", release="LIQ", per_ml=5, market="UK"),
        Product("Li-Liquid 509 mg/5 mL", 509, salt="citrate", release="LIQ", per_ml=5, mmol_label=5.4, market="UK"),
        Product("Li-Liquid 1.018 g/5 mL", 1018, salt="citrate", release="LIQ", per_ml=5, mmol_label=10.8,
                market="UK"),
        # United States (FDA labels; 5 mL oral solution = 8 mEq = 300 mg carbonate)
        Product("Lithium carbonate 150 mg", 150, release="IR", market="US"),
        Product("Lithium carbonate 300 mg", 300, release="IR", market="US"),
        Product("Lithium carbonate 600 mg", 600, release="IR", market="US"),
        Product("Lithium carbonate ER 300 mg", 300, release="SR", market="US"),
        Product("Lithium carbonate ER 450 mg", 450, release="SR", market="US"),
        Product("Lithium citrate oral solution 8 mEq/5 mL", 300, release="LIQ", per_ml=5, mmol_label=8.0,
                market="US"),
    ]
}

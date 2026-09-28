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
    release: str = "IR"          # "IR" (immediate) or "SR" (sustained/prolonged/controlled)
    scored: bool = False         # may be halved
    per_ml: float | None = None  # liquids: strength is per this many mL
    market: str = ""

    @property
    def mmol(self) -> float:
        return mg_to_mmol(self.strength_mg, self.salt)


# Illustrative catalogue. Brand availability changes; confirm locally before use.
PRODUCTS: dict[str, Product] = {
    p.name: p
    for p in [
        # Australia
        Product("Lithicarb 250 mg", 250, release="IR", scored=True, market="AU"),
        Product("Lithicarb 450 mg", 450, release="IR", scored=True, market="AU"),
        Product("Quilonum SR 450 mg", 450, release="SR", scored=True, market="AU"),
        # United Kingdom
        Product("Priadel 200 mg", 200, release="SR", scored=True, market="UK"),
        Product("Priadel 400 mg", 400, release="SR", scored=True, market="UK"),
        Product("Camcolit 400 mg", 400, release="SR", scored=True, market="UK"),
        Product("Priadel liquid 520 mg/5 mL", 520, salt="citrate", release="IR", per_ml=5, market="UK"),
        # United States
        Product("Lithium carbonate 150 mg cap", 150, release="IR", market="US"),
        Product("Lithium carbonate 300 mg cap", 300, release="IR", market="US"),
        Product("Lithium carbonate 600 mg cap", 600, release="IR", market="US"),
        Product("Lithium carbonate ER 300 mg", 300, release="SR", market="US"),
        Product("Lithium carbonate ER 450 mg", 450, release="SR", scored=True, market="US"),
    ]
}

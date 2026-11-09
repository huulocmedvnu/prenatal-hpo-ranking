"""Gene identifiers harmonised to current HGNC approved symbols through NCBI Gene IDs. HGNC complete set 2026-09-04."""
import os
from functools import lru_cache
import pandas as pd

IN = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "inputs")


@lru_cache(None)
def hgnc():
    H = pd.read_csv(os.path.join(IN, "hgnc", "hgnc_complete_set.txt"), sep="\t", dtype=str,
                    usecols=["symbol", "status", "prev_symbol", "alias_symbol", "entrez_id"])
    H = H[H.status == "Approved"]
    by_id = dict(zip(H.entrez_id.dropna(), H.dropna(subset=["entrez_id"]).symbol))
    approved = set(H.symbol)
    prev, alias = {}, {}
    for s, p, a in zip(H.symbol, H.prev_symbol, H.alias_symbol):
        for col, d in ((p, prev), (a, alias)):
            if isinstance(col, str):
                for x in col.split("|"):
                    d.setdefault(x, set()).add(s)
    return by_id, approved, prev, alias


def symbol_from_id(ncbi_id):
    """Current HGNC symbol of an NCBI Gene ID, or None for a locus without an approved HGNC record (such loci, for
    example imprinting control regions, are not genes and are left out of the ranking)."""
    return hgnc()[0].get(str(ncbi_id))


def symbol_from_source(s):
    """Current HGNC symbol of a source-reported symbol: the symbol itself if approved, else the unique approved symbol
    that lists it as a previous symbol, else the unique one that lists it as an alias. Ambiguous symbols are kept."""
    _, approved, prev, alias = hgnc()
    if s in approved:
        return s
    if s in prev:
        return next(iter(prev[s])) if len(prev[s]) == 1 else s
    if len(alias.get(s, ())) == 1:
        return next(iter(alias[s]))
    return s

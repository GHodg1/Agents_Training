from typing import Any, Dict, Optional
import os
import json

from smolagents.tools import Tool

try:
    import pubchempy as pcp
except ImportError as e:
    raise ImportError(
        "You must install `pubchempy` to run this tool. Try: pip install pubchempy"
    ) from e


class PubChemResolveTool(Tool):
    """
    Resolve chemical identifiers/names via PubChem with simple on-disk caching.

    Inputs:
      - query: identifier (name, InChIKey, SMILES, formula, or CID)
      - output_format: one of {'iupac_name','inchikey','smiles','common_name','cid','formula','json'}
                       json returns the full record as JSON text.
    """
    name = "pubchem_resolve"
    description = (
        "Resolves a chemical identifier or name using PubChem and returns structure/ID fields. "
        "Supports lookup by written token, name, InChIKey, SMILES, formula, or CID. Uses a local cache file."
    )
    inputs = {
        "query": {
            "type": "string",
            "description": "Identifier (name, InChIKey, SMILES, formula, or CID).",
        },
        "output_format": {
            "type": "string",
            "description": "Which field to return; use 'json' for the full record.",
            "enum": ["iupac_name", "inchikey", "smiles", "common_name", "cid", "formula", "json"],
            "nullable": True,  # smolagents requires this when the arg has a default
        },
    }
    output_type = "string"

    # in-memory caches and alias maps (per-process)
    _cache: Dict[str, Dict[str, Optional[str]]] = {}
    _alias_cs: Dict[str, str] = {}  # exact match (SMILES, InChIKey, formula, CID str) -> written key
    _alias_ci: Dict[str, str] = {}  # case-insensitive (written, standard name) -> written key

    def __init__(self, cache_filename: str = "pubchem_cache.json"):
        super().__init__()
        self.cache_filename = cache_filename
        self._load_cache()

    # ---------- cache helpers ----------
    def _norm_name(self, s: str) -> str:
        return s.strip().lower() if isinstance(s, str) else s

    def _rebuild_alias_indexes(self) -> None:
        self._alias_cs.clear()
        self._alias_ci.clear()
        for written, v in (self._cache or {}).items():
            name = (v.get("iupac_name") or "").strip()
            inchikey = (v.get("inchikey") or "").strip()
            smiles = (v.get("smiles") or "").strip()
            common_name = (v.get("common_name") or "").strip()
            formula = (v.get("formula") or "").strip()
            cid = v.get("cid")
            cid_str = str(cid).strip() if cid not in (None, "") else ""

            # case-insensitive names
            for alias in filter(None, {written, name, common_name}):
                self._alias_ci[self._norm_name(alias)] = written

            # case-sensitive structure/ID keys
            for alias in filter(None, {inchikey, smiles, formula, cid_str}):
                self._alias_cs[alias] = written

    def _load_cache(self) -> None:
        if not os.path.exists(self.cache_filename):
            with open(self.cache_filename, "w") as f:
                json.dump({}, f)
            self._cache = {}
        else:
            try:
                with open(self.cache_filename, "r") as f:
                    self._cache = json.load(f)
            except json.JSONDecodeError:
                self._cache = {}
        self._rebuild_alias_indexes()

    def _persist_cache(self) -> None:
        with open(self.cache_filename, "w") as f:
            json.dump(self._cache, f, indent=2)

    def _add_to_cache(self, written: str, rec: Dict[str, Optional[str]]) -> None:
        self._cache[written] = rec
        self._rebuild_alias_indexes()
        self._persist_cache()

    def _lookup_cache(self, query: str) -> Optional[Dict[str, Optional[str]]]:
        if not isinstance(query, str) or not query.strip():
            return None
        q = query.strip()

        # exact-match aliases first (do NOT lowercase)
        written = self._alias_cs.get(q)
        if written:
            return self._cache.get(written)

        # then case-insensitive names
        written = self._alias_ci.get(self._norm_name(q))
        if written:
            return self._cache.get(written)

        # direct key hit (as written)
        return self._cache.get(q)

    # ---------- PubChem fetch ----------
    def _fetch_from_pubchem(self, token: str) -> Optional[Dict[str, Optional[str]]]:
        # Try multiple namespaces
        for ns in ("inchikey", "smiles", "formula", "name", "cid"):
            try:
                comps = pcp.get_compounds(token, ns)
                if not comps:
                    continue
                c = comps[0]
                rec = {
                    "iupac_name": getattr(c, "iupac_name", None) or "",
                    "inchikey": getattr(c, "inchikey", None) or "",
                    "smiles": getattr(c, "isomeric_smiles", None) or getattr(c, "smiles", None) or "",
                    "common_name": (getattr(c, "synonyms", None) or [None])[0] or "",
                    "cid": getattr(c, "cid", None) or "",
                    "formula": getattr(c, "molecular_formula", None) or "",
                }
                return rec
            except Exception:
                continue
        return None

    # ---------- Tool interface ----------
    def forward(self, query: str, output_format: str = "json") -> str:
        """
        Resolve query using cache→PubChem and return the requested field or JSON.
        """
        if not isinstance(query, str) or not query.strip():
            raise ValueError("Query must be a non-empty string.")
        output_format = output_format or "json"
        # ...existing code...
 
        # 1) cache
        rec = self._lookup_cache(query)
        # 2) fetch if missing
        if rec is None:
            fetched = self._fetch_from_pubchem(query)
            if fetched is None:
                # cache a negative result to avoid repeated lookups
                self._add_to_cache(query, {
                    "iupac_name": None,
                    "inchikey": None,
                    "smiles": None,
                    "common_name": None,
                    "cid": None,
                    "formula": None,
                })
                return f"No PubChem result found for: {query}"
            self._add_to_cache(query, fetched)
            rec = fetched

        # 3) render output
        if output_format == "json":
            return json.dumps(rec, indent=2)
        if output_format not in {"iupac_name", "inchikey", "smiles", "common_name", "cid", "formula"}:
            raise ValueError(
                "output_format must be one of {'iupac_name','inchikey','smiles','common_name','cid','formula','json'}"
            )
        value = rec.get(output_format)
        if value in (None, ""):
            return f"No value for '{output_format}' found for: {query}"
        return str(value)
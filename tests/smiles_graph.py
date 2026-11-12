"""Minimal SMILES reader for the ligand tests (no RDKit needed in CI).

Handles what LigandLibrary uses: organic-subset and bracket atoms, branches, single-digit ring closures, = and #
bonds, and @/@@ stereo marks. Aromatic (lower-case) atoms are not supported: the library SMILES are Kekule forms.
Gives the molecular formula with implicit hydrogens, the ring count and which bonds lie on a ring.
"""
import re

_TOKEN = re.compile(r"\[[^\]]+\]|Cl|Br|[BCNOPSFI]|[=#]|\(|\)|\d")
_VALENCE = {"C": (4,), "N": (3,), "O": (2,), "P": (3, 5), "S": (2, 4, 6), "F": (1,), "Cl": (1,), "Br": (1,),
            "I": (1,), "B": (3,)}


def parse(smiles):
    """Return (atoms, bonds): atoms = [{'el', 'h' (None = implicit), 'charge'}], bonds = {(i, j): order}."""
    atoms, bonds, stack, rings = [], {}, [], {}
    prev, order = None, 1
    for tok in _TOKEN.findall(smiles):
        if tok == "(":
            stack.append(prev)
        elif tok == ")":
            prev = stack.pop()
        elif tok in "=#":
            order = 2 if tok == "=" else 3
        elif tok.isdigit():
            if tok in rings:
                j, o = rings.pop(tok)
                bonds[(j, prev)] = max(order, o)
            else:
                rings[tok] = (prev, order)
            order = 1
        else:
            if tok.startswith("["):
                m = re.match(r"\[\d*([A-Z][a-z]?)[@]*(H\d*)?([+-]\d*)?\]", tok)
                h = m.group(2)
                charge = m.group(3) or ""
                atoms.append({"el": m.group(1), "h": (int(h[1:] or 1) if h else 0),
                              "charge": (int(charge[1:] or 1) * (1 if charge[0] == "+" else -1)) if charge else 0})
            else:
                atoms.append({"el": tok, "h": None, "charge": 0})
            if prev is not None:
                bonds[(prev, len(atoms) - 1)] = order
            prev, order = len(atoms) - 1, 1
    if rings:
        raise ValueError(f"unclosed ring bond(s) {sorted(rings)} in {smiles}")
    return atoms, bonds


def formula(smiles):
    """Molecular formula as {element: count}, hydrogens included."""
    atoms, bonds = parse(smiles)
    used = [0] * len(atoms)
    for (i, j), o in bonds.items():
        used[i] += o
        used[j] += o
    out = {}
    for a, u in zip(atoms, used):
        out[a["el"]] = out.get(a["el"], 0) + 1
        h = a["h"]
        if h is None:
            h = next((v - u for v in _VALENCE[a["el"]] if v >= u), 0)
        out["H"] = out.get("H", 0) + h
    return out


def ring_count(smiles):
    """Cyclomatic number (bonds - atoms + 1 for one connected molecule): the number of independent rings."""
    atoms, bonds = parse(smiles)
    return len(bonds) - len(atoms) + 1


def in_ring(smiles, i, j):
    """True if bond i-j lies on a ring (the two atoms stay connected without it)."""
    atoms, bonds = parse(smiles)
    adj = {k: set() for k in range(len(atoms))}
    for a, b in bonds:
        if {a, b} != {i, j}:
            adj[a].add(b)
            adj[b].add(a)
    seen, todo = {i}, [i]
    while todo:
        for n in adj[todo.pop()] - seen:
            seen.add(n)
            todo.append(n)
    return j in seen


def neighbours(smiles, i):
    atoms, bonds = parse(smiles)
    return sorted({b for a, b in bonds if a == i} | {a for a, b in bonds if b == i})

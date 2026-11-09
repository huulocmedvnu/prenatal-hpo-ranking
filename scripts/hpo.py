"""HPO helpers: obo parser (names, synonyms, alt_ids, parents), ancestor closure."""
import re


def load_obo(path):
    terms, cur = {}, None
    with open(path, encoding="utf-8") as f:
        lines = [line.rstrip("\n") for line in f]
    for line in lines:
        if line == "[Term]":
            cur = {"parents": [], "syn": [], "alt": [], "obsolete": False}
        elif line.startswith("[") and line != "[Term]":
            cur = None
        elif cur is None:
            continue
        elif line.startswith("id: "):
            cur["id"] = line[4:]; terms[cur["id"]] = cur
        elif line.startswith("name: "):
            cur["name"] = line[6:]
        elif line.startswith("is_a: "):
            cur["parents"].append(line[6:16])
        elif line.startswith("synonym: "):
            m = re.match(r'synonym: "(.*)"', line)
            if m:
                cur["syn"].append(m.group(1))
        elif line.startswith("alt_id: "):
            cur["alt"].append(line[8:])
        elif line.startswith("is_obsolete: true"):
            cur["obsolete"] = True
        elif line.startswith("replaced_by: "):
            cur["replaced_by"] = line[13:]
    return terms


def ancestors_map(terms):
    memo = {}

    def anc(t):
        if t in memo:
            return memo[t]
        out = {t}
        for p in terms.get(t, {}).get("parents", []):
            out |= anc(p)
        memo[t] = out
        return out
    for t in terms:
        anc(t)
    return memo


def label_index(terms):
    """Lower-case label to term ID over active terms. Primary names take precedence over synonyms, and the first term
    in file order wins when two terms share a synonym."""
    idx = {}
    active = [(t, d) for t, d in terms.items() if not d.get("obsolete")]
    for t, d in active:
        idx.setdefault(d.get("name", "").strip().lower(), t)
    for t, d in active:
        for s in d["syn"]:
            idx.setdefault(s.strip().lower(), t)
    return idx

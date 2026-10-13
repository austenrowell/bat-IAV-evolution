#!/usr/bin/env python3

from pathlib import Path
from collections import Counter, defaultdict
import csv
import re
import math
import statistics

ROOT = Path(".")
TARGET_FILE = ROOT / "h17_h18_target_accessions.tsv"
RUNS_DIR = ROOT / "runs_final_markov_jumps"
OUT_DIR = ROOT / "markov_jumps_to_h17_h18_bat_clades"
OUT_DIR.mkdir(exist_ok=True)

BURNIN_FRACTION = 0.10
BAT_STATE = "Bat"

SEGMENTS = [
    "segment_1",
    "segment_2",
    "segment_3_PA",
    "segment_4",
    "segment_5",
    "segment_6",
    "segment_7_M1",
    "segment_7_M2",
    "segment_8_NS1",
    "segment_8_NS2",
]


class Node:
    def __init__(self, name=None, ann=None, length=None):
        self.name = name
        self.ann = ann or {}
        self.length = length
        self.children = []
        self.parent = None

    def is_leaf(self):
        return len(self.children) == 0


def split_top_level(s, sep=","):
    parts = []
    depth = 0
    brace = 0
    bracket = 0
    quote = False
    start = 0

    for i, ch in enumerate(s):
        if ch == '"' and (i == 0 or s[i - 1] != "\\"):
            quote = not quote
        elif not quote:
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            elif ch == "{":
                brace += 1
            elif ch == "}":
                brace -= 1
            elif ch == "[":
                bracket += 1
            elif ch == "]":
                bracket -= 1
            elif ch == sep and depth == 0 and brace == 0 and bracket == 0:
                parts.append(s[start:i])
                start = i + 1

    parts.append(s[start:])
    return [p.strip() for p in parts if p.strip()]


def parse_annotation(ann_text):
    ann = {}
    if not ann_text:
        return ann

    if ann_text.startswith("[&") and ann_text.endswith("]"):
        ann_text = ann_text[2:-1]

    parts = split_top_level(ann_text, ",")

    for part in parts:
        if "=" not in part:
            continue
        key, val = part.split("=", 1)
        key = key.strip()
        val = val.strip()

        if val.startswith('"') and val.endswith('"'):
            val = val[1:-1]

        ann[key] = val

    return ann


def parse_history(val):
    """
    Parses BEAST history like:
    {{time,From,To},{time,From,To}}
    Returns ordered list of (time, from_state, to_state).
    """
    if not val or val == "{}":
        return []

    text = val.strip()
    events = []

    for m in re.finditer(r"\{([^{}]+)\}", text):
        bits = [x.strip().strip('"') for x in m.group(1).split(",")]
        if len(bits) >= 3:
            try:
                t = float(bits[0])
            except ValueError:
                t = math.nan
            events.append((t, bits[1], bits[2]))

    # In these BEAST outputs, event order usually follows the direction along the branch.
    # Sorting from older to younger is safer when heights are available.
    events.sort(key=lambda x: x[0], reverse=True)
    return events


class Parser:
    def __init__(self, s):
        self.s = s.strip()
        self.i = 0
        self.n = len(self.s)

    def peek(self):
        return self.s[self.i] if self.i < self.n else ""

    def consume_ws(self):
        while self.i < self.n and self.s[self.i].isspace():
            self.i += 1

    def parse(self):
        self.consume_ws()
        node = self.parse_subtree()
        return node

    def parse_subtree(self):
        self.consume_ws()

        if self.peek() == "(":
            self.i += 1
            node = Node()

            while True:
                child = self.parse_subtree()
                child.parent = node
                node.children.append(child)

                self.consume_ws()
                ch = self.peek()
                if ch == ",":
                    self.i += 1
                    continue
                elif ch == ")":
                    self.i += 1
                    break
                else:
                    raise ValueError(f"Unexpected char in child list at {self.i}: {ch}")

            self.parse_label_ann_length(node)
            return node

        else:
            node = Node()
            self.parse_label_ann_length(node)
            return node

    def parse_label_ann_length(self, node):
        self.consume_ws()

        # label/name
        if self.peek() == "'":
            self.i += 1
            start = self.i
            while self.i < self.n and self.s[self.i] != "'":
                self.i += 1
            node.name = self.s[start:self.i]
            if self.i < self.n and self.s[self.i] == "'":
                self.i += 1
        else:
            start = self.i
            while self.i < self.n and self.s[self.i] not in ":,();[":
                self.i += 1
            label = self.s[start:self.i].strip()
            if label:
                node.name = label

        self.consume_ws()

        # annotation after label/node
        if self.peek() == "[":
            start = self.i
            depth = 0
            quote = False
            while self.i < self.n:
                ch = self.s[self.i]
                if ch == '"' and (self.i == 0 or self.s[self.i - 1] != "\\"):
                    quote = not quote
                elif not quote:
                    if ch == "[":
                        depth += 1
                    elif ch == "]":
                        depth -= 1
                        if depth == 0:
                            self.i += 1
                            break
                self.i += 1
            node.ann.update(parse_annotation(self.s[start:self.i]))

        self.consume_ws()

        # branch length, possibly with annotation after colon
        if self.peek() == ":":
            self.i += 1
            self.consume_ws()

            if self.peek() == "[":
                start = self.i
                depth = 0
                quote = False
                while self.i < self.n:
                    ch = self.s[self.i]
                    if ch == '"' and (self.i == 0 or self.s[self.i - 1] != "\\"):
                        quote = not quote
                    elif not quote:
                        if ch == "[":
                            depth += 1
                        elif ch == "]":
                            depth -= 1
                            if depth == 0:
                                self.i += 1
                                break
                    self.i += 1
                node.ann.update(parse_annotation(self.s[start:self.i]))

            self.consume_ws()

            start = self.i
            while self.i < self.n and self.s[self.i] not in ",();":
                self.i += 1
            length_text = self.s[start:self.i].strip()
            try:
                node.length = float(length_text)
            except ValueError:
                node.length = None


def read_targets(path):
    targets = defaultdict(list)
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            seg = row["segment"].strip()
            tax = row["taxon"].strip()
            if seg and tax:
                targets[seg].append(tax)
    return dict(targets)


def read_taxlabels(nexus_path):
    """
    Read the NEXUS Taxlabels block and map BEAST numeric tree labels
    to accession names.

    BEAST/TreeAnnotator history trees commonly write tips as numbers
    in the Newick tree, with the actual taxon names listed in the
    Taxlabels block. This function maps:
        "1" -> first taxlabel
        "2" -> second taxlabel
        ...
    """
    text = nexus_path.read_text(errors="replace")

    m = re.search(
        r"Begin\s+taxa\s*;.*?Taxlabels\s+(.*?)\s*;",
        text,
        flags=re.IGNORECASE | re.DOTALL
    )

    if not m:
        raise ValueError(f"Could not find Taxlabels block in {nexus_path}")

    block = m.group(1)

    labels = []
    for q in re.findall(r"'([^']+)'", block):
        q = q.strip()
        if q:
            labels.append(q)

    # Fallback for unquoted labels
    if not labels:
        for tok in re.split(r"\s+", block.strip()):
            tok = tok.strip().strip("'").strip()
            if tok:
                labels.append(tok)

    translation = {str(i + 1): label for i, label in enumerate(labels)}

    return translation

def iter_tree_lines(nexus_path):
    with nexus_path.open(errors="replace") as f:
        for line in f:
            if line.startswith("tree STATE_"):
                yield line.rstrip("\n")


def extract_newick_from_tree_line(line):
    """
    Extract Newick from a BEAST/NEXUS tree line.

    BEAST history-tree lines look like:
      tree STATE_0 [&empiricalTreeNumber=0.0,...] = [&R] (...);

    We must split at the real tree separator, ' = ', not at the first '='
    inside the metadata annotation.
    """
    if " = " in line:
        rhs = line.split(" = ", 1)[1].strip()
    else:
        # Fallback: find the first opening parenthesis after the tree metadata.
        idx = line.find("(")
        if idx == -1:
            raise ValueError("Could not find Newick tree in line")
        rhs = line[idx:].strip()

    rhs = rhs.rstrip(";").strip()

    # Remove leading rooted/unrooted annotation, e.g. [&R] or [&U].
    rhs = re.sub(r'^\[&[^\]]+\]\s*', '', rhs)

    return rhs

def translate_leaf_names(root, translation):
    for node in preorder(root):
        if node.is_leaf() and node.name in translation:
            node.name = translation[node.name]


def preorder(root):
    yield root
    for child in root.children:
        yield from preorder(child)


def leaf_set(node):
    return {n.name for n in preorder(node) if n.is_leaf()}


def find_mrca(root, targets):
    targets = set(targets)

    best = None
    best_size = None

    for node in preorder(root):
        ls = leaf_set(node)
        if targets.issubset(ls):
            size = len(ls)
            if best is None or size < best_size:
                best = node
                best_size = size

    return best


def path_root_to_node(node):
    path = []
    cur = node
    while cur is not None:
        path.append(cur)
        cur = cur.parent
    return list(reversed(path))


def node_state(node):
    return node.ann.get("host.states") or node.ann.get("host.state") or node.ann.get("host")


def branch_events_to_child(child):
    return parse_history(child.ann.get("history", ""))


def branch_jump_count_to_child(child):
    events = branch_events_to_child(child)
    if events:
        return len(events)

    val = child.ann.get("host.count") or child.ann.get("sumAcrossArray_host.count")
    try:
        return int(round(float(val)))
    except Exception:
        return 0


def analyze_tree(root, targets):
    mrca = find_mrca(root, targets)
    if mrca is None:
        return None

    path = path_root_to_node(mrca)

    total_jumps = 0
    expanded_states = []

    # Start with root state.
    root_state = node_state(path[0])
    expanded_states.append({
        "kind": "node",
        "node": path[0],
        "state": root_state,
        "node_distance_to_mrca": len(path) - 1,
    })

    last_nonbat_before_bat = None
    final_entry_to_bat = None

    # Walk root -> MRCA. Each child carries branch annotation from parent to child.
    for idx in range(1, len(path)):
        child = path[idx]
        child_state = node_state(child)
        node_distance = len(path) - 1 - idx

        events = branch_events_to_child(child)
        total_jumps += len(events) if events else branch_jump_count_to_child(child)

        current_state = expanded_states[-1]["state"]

        for event_time, from_state, to_state in events:
            if to_state == BAT_STATE and from_state != BAT_STATE:
                final_entry_to_bat = {
                    "source_state": from_state,
                    "target_state": to_state,
                    "event_time": event_time,
                    "node_distance_to_mrca": node_distance,
                }
            expanded_states.append({
                "kind": "jump",
                "node": child,
                "state": to_state,
                "from_state": from_state,
                "to_state": to_state,
                "event_time": event_time,
                "node_distance_to_mrca": node_distance,
            })
            current_state = to_state

        if child_state is not None:
            # If no explicit history was present but node state changed, count that as inferred transition.
            if not events and current_state is not None and child_state != current_state:
                total_jumps += 1
                if child_state == BAT_STATE and current_state != BAT_STATE:
                    final_entry_to_bat = {
                        "source_state": current_state,
                        "target_state": child_state,
                        "event_time": None,
                        "node_distance_to_mrca": node_distance,
                    }

            expanded_states.append({
                "kind": "node",
                "node": child,
                "state": child_state,
                "node_distance_to_mrca": node_distance,
            })

    # nearest upstream non-bat state by node states only
    nearest_nonbat_node_state = None
    nearest_nonbat_node_distance = None
    for idx in range(len(path) - 1, -1, -1):
        st = node_state(path[idx])
        if st is not None and st != BAT_STATE:
            nearest_nonbat_node_state = st
            nearest_nonbat_node_distance = len(path) - 1 - idx
            break

    if final_entry_to_bat:
        nearest_nonbat_source = final_entry_to_bat["source_state"]
        transition_node_distance = final_entry_to_bat["node_distance_to_mrca"]
    else:
        nearest_nonbat_source = nearest_nonbat_node_state
        transition_node_distance = nearest_nonbat_node_distance

    mrca_state = node_state(mrca)
    mrca_clade_size = len(leaf_set(mrca))

    return {
        "total_jumps_root_to_mrca": total_jumps,
        "mrca_state": mrca_state,
        "mrca_clade_size": mrca_clade_size,
        "nearest_nonbat_source": nearest_nonbat_source,
        "transition_node_distance_to_mrca": transition_node_distance,
        "nearest_nonbat_node_state": nearest_nonbat_node_state,
        "nearest_nonbat_node_distance_to_mrca": nearest_nonbat_node_distance,
    }


def quantile(xs, q):
    if not xs:
        return None
    xs = sorted(xs)
    pos = (len(xs) - 1) * q
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return xs[int(pos)]
    return xs[lo] * (hi - pos) + xs[hi] * (pos - lo)


def summarize(rows, seg):
    n = len(rows)
    jumps = [int(r["total_jumps_root_to_mrca"]) for r in rows]
    distances = [
        int(r["transition_node_distance_to_mrca"])
        for r in rows
        if r["transition_node_distance_to_mrca"] not in {"", None}
    ]

    source_counts = Counter(r["nearest_nonbat_source"] for r in rows if r["nearest_nonbat_source"])
    mrca_state_counts = Counter(r["mrca_state"] for r in rows if r["mrca_state"])

    source_state, source_count = source_counts.most_common(1)[0] if source_counts else ("NA", 0)
    mrca_state, mrca_count = mrca_state_counts.most_common(1)[0] if mrca_state_counts else ("NA", 0)

    return {
        "segment": seg,
        "posterior_trees_analyzed": n,
        "mean_jumps_root_to_bat_mrca": statistics.mean(jumps) if jumps else "NA",
        "median_jumps_root_to_bat_mrca": statistics.median(jumps) if jumps else "NA",
        "jumps_95_lower": quantile(jumps, 0.025) if jumps else "NA",
        "jumps_95_upper": quantile(jumps, 0.975) if jumps else "NA",
        "most_common_jumps": Counter(jumps).most_common(1)[0][0] if jumps else "NA",
        "most_common_nearest_nonbat_source": source_state,
        "posterior_probability_nearest_nonbat_source": source_count / n if n else "NA",
        "most_common_mrca_state": mrca_state,
        "posterior_probability_mrca_state": mrca_count / n if n else "NA",
        "median_transition_node_distance_to_mrca": statistics.median(distances) if distances else "NA",
        "most_common_transition_node_distance_to_mrca": Counter(distances).most_common(1)[0][0] if distances else "NA",
    }


def main():
    targets_by_segment = read_targets(TARGET_FILE)

    all_summary = []

    for seg in SEGMENTS:
        if seg not in targets_by_segment:
            print(f"Skipping {seg}: no target taxa")
            continue

        run_dir = RUNS_DIR / f"{seg}_75subset_lsd2.timetree"
        hist_file = run_dir / f"{seg}_75subset_lsd2.timetree.host.history.trees.txt"

        if not hist_file.exists():
            print(f"Skipping {seg}: missing {hist_file}")
            continue

        targets = targets_by_segment[seg]
        translation = read_taxlabels(hist_file)

        lines = list(iter_tree_lines(hist_file))
        burn = int(len(lines) * BURNIN_FRACTION)
        kept = lines[burn:]

        print(f"{seg}: total trees={len(lines)}, burnin={burn}, analyzing={len(kept)}, targets={targets}")

        rows = []
        for i, line in enumerate(kept, start=burn):
            nwk = extract_newick_from_tree_line(line)
            root = Parser(nwk).parse()
            translate_leaf_names(root, translation)

            present = leaf_set(root)
            missing = [t for t in targets if t not in present]
            if missing:
                print(f"WARNING {seg} STATE index {i}: missing targets {missing}")
                continue

            result = analyze_tree(root, targets)
            if result is None:
                print(f"WARNING {seg} STATE index {i}: no MRCA found")
                continue

            result["segment"] = seg
            result["tree_index_after_burnin"] = i
            result["targets"] = ";".join(targets)
            rows.append(result)

        per_tree_out = OUT_DIR / f"{seg}_h17_h18_markov_jumps_to_bat_mrca_per_tree.csv"
        with per_tree_out.open("w", newline="") as f:
            fieldnames = [
                "segment",
                "tree_index_after_burnin",
                "targets",
                "total_jumps_root_to_mrca",
                "mrca_state",
                "mrca_clade_size",
                "nearest_nonbat_source",
                "transition_node_distance_to_mrca",
                "nearest_nonbat_node_state",
                "nearest_nonbat_node_distance_to_mrca",
            ]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        if rows:
            all_summary.append(summarize(rows, seg))

    summary_out = OUT_DIR / "h17_h18_markov_jumps_to_bat_mrca_summary.csv"
    with summary_out.open("w", newline="") as f:
        fieldnames = [
            "segment",
            "posterior_trees_analyzed",
            "mean_jumps_root_to_bat_mrca",
            "median_jumps_root_to_bat_mrca",
            "jumps_95_lower",
            "jumps_95_upper",
            "most_common_jumps",
            "most_common_nearest_nonbat_source",
            "posterior_probability_nearest_nonbat_source",
            "most_common_mrca_state",
            "posterior_probability_mrca_state",
            "median_transition_node_distance_to_mrca",
            "most_common_transition_node_distance_to_mrca",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_summary)

    print(f"\nWrote summary: {summary_out}")


if __name__ == "__main__":
    main()

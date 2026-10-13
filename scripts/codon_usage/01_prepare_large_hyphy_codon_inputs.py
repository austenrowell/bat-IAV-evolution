#!/usr/bin/env python3
"""
Prepare large-tree HyPhy and codon-usage inputs with the same foreground-labeling
logic used for the 75-taxon subset analysis.

What it does
------------
For each segment:
  1. Finds the large alignment.
  2. Finds the large IQ-TREE treefile.
  3. Strips old HyPhy/PAML labels.
  4. Finds bat tips from bat_sequences.csv.
  5. Labels the MRCA branch of present bat tips as Test.
  6. Labels all other non-root branches as Reference for RELAX.
  7. Writes HyPhy-ready alignments and trees.
  8. Writes a HyPhy job list and manifest.

Important consistency rule
--------------------------
The foreground branch is labeled using the same logic as the 75-subset workflow:
the MRCA of all bat accessions present in both the alignment and tree is labeled
as the Test branch. If the bat MRCA is the root, the segment is skipped and
reported as foreground_is_root.

Default expected large-dataset structure
----------------------------------------
/lustre/work/aurowell/influenza/cds/
  segment_1/
    segment_1_clean.fasta
    *.treefile
  segment_2/
    segment_2_clean.fasta
    *.treefile
  ...

By default, this script DOES NOT use top-level segment_*.tree files, because
those are often BEAST/TreeAnnotator summary trees rather than IQ-TREE ML trees.
Use --allow-top-level-tree only for verified IQ-TREE ML input trees.
"""

from pathlib import Path
import argparse
import csv
import re
import shutil
from collections import OrderedDict


SEGMENTS_DEFAULT = [
    "segment_1",
    "segment_2",
    "segment_3_PA",
    "segment_4",
    "segment_5",
    "segment_6",
    "segment_7_M1",
    "segment_7_M2",
    "segment_8_NEP",
    "segment_8_NS1",
    "segment_8_NS2",
]


CODON_TO_AA = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"A","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}


class Node:
    __slots__ = ("name", "length", "children", "parent", "mark", "mark_ref")

    def __init__(self, name="", length=None):
        self.name = name
        self.length = length
        self.children = []
        self.parent = None
        self.mark = False
        self.mark_ref = False

    def is_leaf(self):
        return len(self.children) == 0


def parse_newick(text):
    text = text.strip()
    if text.endswith(";"):
        text = text[:-1]

    i = 0
    n = len(text)

    def skip_ws():
        nonlocal i
        while i < n and text[i].isspace():
            i += 1

    def parse_label():
        nonlocal i
        skip_ws()
        if i >= n:
            return ""
        if text[i] in "'\"":
            q = text[i]
            i += 1
            out = []
            while i < n:
                if text[i] == q:
                    if i + 1 < n and text[i + 1] == q:
                        out.append(q)
                        i += 2
                        continue
                    i += 1
                    break
                out.append(text[i])
                i += 1
            return "".join(out)

        start = i
        while i < n and text[i] not in ",():;":
            i += 1
        return text[start:i].strip()

    def parse_length():
        nonlocal i
        skip_ws()
        if i < n and text[i] == ":":
            i += 1
            start = i
            while i < n and text[i] not in ",()":
                i += 1
            return text[start:i].strip()
        return None

    def parse_subtree():
        nonlocal i
        skip_ws()
        node = Node()

        if i < n and text[i] == "(":
            i += 1
            while True:
                child = parse_subtree()
                child.parent = node
                node.children.append(child)

                skip_ws()
                if i < n and text[i] == ",":
                    i += 1
                    continue
                if i < n and text[i] == ")":
                    i += 1
                    break
                raise ValueError(f"Malformed Newick near position {i}: {text[i:i+80]}")

            node.name = parse_label()
            node.length = parse_length()
        else:
            node.name = parse_label()
            node.length = parse_length()

        return node

    root = parse_subtree()
    return root


def iter_nodes(root):
    stack = [root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(reversed(node.children))


def leaves(root):
    return [n for n in iter_nodes(root) if n.is_leaf()]


def set_parents(root):
    root.parent = None
    for node in iter_nodes(root):
        for child in node.children:
            child.parent = node


def clean_internal_support_labels(root):
    for node in iter_nodes(root):
        if not node.is_leaf() and re.fullmatch(r"\d+(?:\.\d+)?", node.name or ""):
            node.name = ""


def strip_existing_labels(text):
    text = re.sub(r"\{[^{}]+\}", "", text)
    text = text.replace("#1", "")
    return text


def quote_if_needed(name):
    if name is None or name == "":
        return ""
    if re.search(r"[\s,:;(){}\[\]'\"#]", name):
        return "'" + name.replace("'", "''") + "'"
    return name


def write_newick(node, label_style="none"):
    if node.is_leaf():
        s = quote_if_needed(node.name)
    else:
        s = "(" + ",".join(write_newick(c, label_style=label_style) for c in node.children) + ")"
        if node.name:
            s += quote_if_needed(node.name)

    if label_style == "hyphy_test" and node.mark:
        s += "{Test}"
    elif label_style == "hyphy_relax":
        if node.mark:
            s += "{Test}"
        elif node.mark_ref:
            s += "{Reference}"

    if node.length not in [None, ""]:
        s += ":" + str(node.length)
    return s


def ancestors(node):
    out = []
    while node is not None:
        out.append(node)
        node = node.parent
    return out


def mrca(nodes):
    if not nodes:
        return None
    ancestor_sets = [set(ancestors(n)) for n in nodes]
    cur = nodes[0]
    while cur is not None:
        if all(cur in s for s in ancestor_sets):
            return cur
        cur = cur.parent
    return None


def normalize_accession(x):
    s = str(x).strip().split()[0] if x is not None else ""
    return s


def accession_variants(x):
    s = normalize_accession(x)
    variants = {s, s.replace(".", "_"), s.replace("_", ".")}
    variants.add(re.sub(r"[._]\d+$", "", s))
    return {v for v in variants if v}


def read_bats(path):
    bats = set()
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        if "Accession" not in reader.fieldnames:
            raise ValueError(f"{path} must contain an Accession column")
        for row in reader:
            acc = row.get("Accession", "").strip()
            if acc:
                for v in accession_variants(acc):
                    bats.add(v)
    return bats


def read_fasta(path):
    records = []
    name = None
    parts = []
    for line in Path(path).read_text().splitlines():
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(parts).upper()))
            name = normalize_accession(line[1:])
            parts = []
        else:
            parts.append(line.strip().replace(" ", "").replace("\t", ""))
    if name is not None:
        records.append((name, "".join(parts).upper()))
    return records


def write_fasta(records, path, wrap=80):
    path = Path(path)
    with path.open("w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), wrap):
                handle.write(seq[i:i+wrap] + "\n")


def translate_codon_alignment(records):
    aa_records = []
    warnings = []
    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        warnings.append(f"unequal_alignment_lengths={sorted(lengths)[:10]}")
    for name, seq in records:
        seq = seq.upper().replace("U", "T")
        usable = (len(seq) // 3) * 3
        if usable != len(seq):
            warnings.append(f"{name}:trimmed_nt_length_from_{len(seq)}_to_{usable}")
        aas = []
        for i in range(0, usable, 3):
            codon = seq[i:i+3]
            if "-" in codon:
                aas.append("-")
            elif any(nt not in "ACGT" for nt in codon):
                aas.append("X")
            else:
                aas.append(CODON_TO_AA.get(codon, "X"))
        aa_records.append((name, "".join(aas)))
    return aa_records, warnings


def find_alignment(input_root, seg, alignment_suffix):
    seg_dir = input_root / seg
    candidates = []
    if alignment_suffix:
        candidates.extend([
            seg_dir / f"{seg}{alignment_suffix}",
            seg_dir / f"{seg}_{alignment_suffix.lstrip('_')}",
        ])
    candidates.extend([
        seg_dir / f"{seg}_clean.fasta",
        seg_dir / f"{seg}.clean.fasta",
        seg_dir / f"{seg}.fasta",
        seg_dir / f"{seg}_subset.fasta",
    ])
    seen = set()
    for p in candidates:
        if p not in seen and p.exists() and p.is_file():
            return p
        seen.add(p)

    matches = sorted([
        p for p in seg_dir.glob("*.fasta")
        if "date" not in p.name.lower()
        and "aa" not in p.name.lower()
        and "protein" not in p.name.lower()
    ])
    return matches[0] if matches else None


def find_tree(input_root, seg, aln, allow_top_level_tree=False):
    seg_dir = input_root / seg
    candidates = []
    if aln is not None:
        candidates.extend([
            seg_dir / f"{aln.name}.treefile",
            seg_dir / f"{aln.stem}.treefile",
            seg_dir / "iqtree" / f"{aln.name}.treefile",
            seg_dir / "iqtree" / f"{aln.stem}.treefile",
        ])
    candidates.extend([
        seg_dir / f"{seg}.treefile",
        seg_dir / f"{seg}_clean.fasta.treefile",
        seg_dir / f"{seg}_clean.treefile",
        seg_dir / "iqtree" / f"{seg}.treefile",
        seg_dir / "iqtree" / f"{seg}_clean.fasta.treefile",
        seg_dir / "iqtree" / f"{seg}_clean.treefile",
    ])

    for p in candidates:
        if p.exists() and p.is_file():
            return p

    matches = sorted(list(seg_dir.glob("*.treefile")) + list((seg_dir / "iqtree").glob("*.treefile")) if (seg_dir / "iqtree").exists() else list(seg_dir.glob("*.treefile")))
    if matches:
        return matches[0]

    if allow_top_level_tree:
        top = input_root / f"{seg}.tree"
        if top.exists():
            return top
        top_tf = input_root / f"{seg}.treefile"
        if top_tf.exists():
            return top_tf

    return None


def mark_foreground(root, bat_names, allow_single=False):
    leaf_by_name = {leaf.name: leaf for leaf in leaves(root)}
    present_bats = sorted([
        leaf.name for leaf in leaves(root)
        if leaf.name in bat_names or any(v in bat_names for v in accession_variants(leaf.name))
    ])

    if len(present_bats) == 0:
        return None, present_bats, "no_bat_tips_in_tree"

    if len(present_bats) == 1 and not allow_single:
        return None, present_bats, "single_bat_tip_no_mrca"

    if len(present_bats) == 1:
        fg = leaf_by_name[present_bats[0]]
    else:
        fg = mrca([leaf_by_name[b] for b in present_bats])

    if fg is root:
        return None, present_bats, "foreground_is_root"

    for node in iter_nodes(root):
        node.mark = False
        node.mark_ref = False

    fg.mark = True
    for node in iter_nodes(root):
        if node is not root and node is not fg:
            node.mark_ref = True

    return fg, present_bats, "OK"


def process_segment(seg, args, bat_names):
    input_root = Path(args.input_root)
    aln = find_alignment(input_root, seg, args.alignment_suffix)
    row = OrderedDict()
    row["dataset"] = args.dataset
    row["segment"] = seg
    row["status"] = "OK"
    row["alignment"] = str(aln) if aln else ""
    row["tree"] = ""
    row["bat_tips_present"] = ""
    row["notes"] = ""

    if aln is None:
        row["status"] = "missing_alignment"
        return row, None

    tree = find_tree(input_root, seg, aln, allow_top_level_tree=args.allow_top_level_tree)
    row["tree"] = str(tree) if tree else ""
    if tree is None:
        row["status"] = "missing_treefile"
        row["notes"] = "No .treefile found. Use --allow-top-level-tree only if top-level segment_*.tree is the intended ML tree."
        return row, None

    records = read_fasta(aln)
    aln_names = {name for name, _ in records}

    tree_text = strip_existing_labels(tree.read_text())
    root = parse_newick(tree_text)
    clean_internal_support_labels(root)
    set_parents(root)
    tree_names = {leaf.name for leaf in leaves(root)}

    if aln_names != tree_names:
        missing_tree = sorted(aln_names - tree_names)
        missing_aln = sorted(tree_names - aln_names)
        row["status"] = "alignment_tree_name_mismatch"
        row["notes"] = f"missing_from_tree={missing_tree[:10]}; missing_from_alignment={missing_aln[:10]}"
        return row, None

    fg, present_bats, status = mark_foreground(root, bat_names, allow_single=args.allow_single_bat_terminal)
    row["bat_tips_present"] = len(present_bats)
    row["notes"] = ";".join(present_bats[:50])
    if status != "OK":
        row["status"] = status
        return row, None

    outdir = Path(args.outroot) / args.dataset / seg
    hyphy_dir = outdir / "hyphy"
    hyphy_dir.mkdir(parents=True, exist_ok=True)

    nt_out = hyphy_dir / f"{seg}.large.nt.fasta"
    aa_out = hyphy_dir / f"{seg}.large.aa.fasta"
    unlabeled = hyphy_dir / f"{seg}.large.unlabeled.tree"
    test_tree = hyphy_dir / f"{seg}.large.hyphy_test.tree"
    relax_tree = hyphy_dir / f"{seg}.large.hyphy_relax.tree"

    write_fasta(records, nt_out)
    aa_records, aa_warnings = translate_codon_alignment(records)
    write_fasta(aa_records, aa_out)

    unlabeled.write_text(write_newick(root, "none") + ";\n")
    test_tree.write_text(write_newick(root, "hyphy_test") + ";\n")
    relax_tree.write_text(write_newick(root, "hyphy_relax") + ";\n")

    if aa_warnings:
        row["notes"] = (row["notes"] + "; " if row["notes"] else "") + "AA_warnings=" + "|".join(aa_warnings[:5])

    row["hyphy_nt"] = str(nt_out)
    row["hyphy_aa"] = str(aa_out)
    row["hyphy_tree_test"] = str(test_tree)
    row["hyphy_tree_relax"] = str(relax_tree)

    job = OrderedDict()
    job["dataset"] = args.dataset
    job["segment"] = seg
    job["nt_alignment"] = str(nt_out.resolve())
    job["aa_alignment"] = str(aa_out.resolve())
    job["tree_unlabeled"] = str(unlabeled.resolve())
    job["tree_test"] = str(test_tree.resolve())
    job["tree_relax"] = str(relax_tree.resolve())
    job["outdir"] = str(hyphy_dir.resolve())

    return row, job


def write_tsv(rows, path, fieldnames):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", default=".")
    parser.add_argument("--outroot", default="large_consistent_inputs")
    parser.add_argument("--dataset", default="large_all")
    parser.add_argument("--bat-csv", default="bat_sequences.csv")
    parser.add_argument("--alignment-suffix", default="_clean.fasta")
    parser.add_argument("--segments", nargs="*", default=SEGMENTS_DEFAULT)
    parser.add_argument("--allow-single-bat-terminal", action="store_true")
    parser.add_argument("--allow-top-level-tree", action="store_true",
                        help="Allow use of top-level segment_*.tree/treefile if no segment-directory .treefile exists. Use cautiously.")
    args = parser.parse_args()

    bat_names = read_bats(args.bat_csv)
    Path(args.outroot).mkdir(parents=True, exist_ok=True)

    manifest = []
    jobs = []

    for seg in args.segments:
        print(f"Processing {seg}")
        row, job = process_segment(seg, args, bat_names)
        manifest.append(row)
        print(f"  status={row['status']} bats={row.get('bat_tips_present','')}")
        if job is not None:
            jobs.append(job)

    manifest_fields = sorted({k for row in manifest for k in row.keys()})
    write_tsv(manifest, Path(args.outroot) / "large_prepare_manifest.tsv", manifest_fields)

    write_tsv(
        jobs,
        Path(args.outroot) / "hyphy_jobs.tsv",
        ["dataset", "segment", "nt_alignment", "aa_alignment", "tree_unlabeled", "tree_test", "tree_relax", "outdir"],
    )

    print("\nDone.")
    print(f"Manifest: {Path(args.outroot) / 'large_prepare_manifest.tsv'}")
    print(f"HyPhy jobs: {Path(args.outroot) / 'hyphy_jobs.tsv'}")
    print(f"Prepared HyPhy jobs: {len(jobs)}")
    if len(jobs) == 0:
        print("No jobs were prepared. Check the manifest for missing treefiles or foreground/root issues.")


if __name__ == "__main__":
    main()

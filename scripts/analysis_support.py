"""Shared sample matching, statistical outputs and execution provenance."""
from importlib import metadata as pkg_metadata
from itertools import combinations
import hashlib
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import re
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

META_COLUMNS = {"Rank", "TaxID", "original_header", "Name", "Scientific Name", "ScientificName"}
UNKNOWN_VALUES = {"", "nan", "none", "unknown", "na", "n/a", "-"}


def canonical_sample(value):
    text = str(value).strip().casefold()
    text = re.sub(r"_contigs_fasta(?:_(?:proteins|genes))?(?:_genes)?$", "", text)
    text = re.sub(r"_l\d+(?:_[12])?$", "", text)
    # Equivalence explicitly confirmed for this study; no prefix matching.
    return {"ncm2wild": "ncm2"}.get(text, text)


def read_frame(path):
    path = Path(path)
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    return pd.read_csv(path, sep="\t" if path.suffix.lower() == ".tsv" else ",")


def count_matrix(frame):
    columns = [c for c in frame.columns if c not in META_COLUMNS]
    if not columns:
        raise ValueError("No sample count columns found")
    counts = frame[columns].apply(pd.to_numeric, errors="raise")
    if counts.isna().any().any() or not np.isfinite(counts.to_numpy()).all():
        raise ValueError("Counts must be finite and non-missing; missing reports are not zero")
    if (counts < 0).any().any():
        raise ValueError("Negative counts are invalid")
    identities = [canonical_sample(c) for c in columns]
    if len(identities) != len(set(identities)):
        raise ValueError("Multiple columns match the same sample; resolve technical replicates first")
    return counts


def metadata_groups(samples, path, category="Time", sample_id="SampleID"):
    frame = read_frame(path)
    for col in (sample_id, category):
        if col not in frame:
            raise ValueError(f"Metadata column not found: {col}")
    frame = frame.copy()
    frame["_sample_key"] = frame[sample_id].map(canonical_sample)
    if frame["_sample_key"].isin(UNKNOWN_VALUES).any():
        raise ValueError("Missing metadata sample ID")
    rows = {}
    for key, repeated in frame.groupby("_sample_key", sort=False):
        compare = [category] + [c for c in ("Type", "SubjectID", "TankID", "Batch") if c in frame]
        if len(repeated[compare].fillna("<missing>").drop_duplicates()) > 1:
            raise ValueError(f"Conflicting metadata for {key}; resolve before statistical analysis")
        rows[key] = repeated.iloc[0]
    is_time = category.casefold() in {"time", "day", "days", "days_in_captivity"}
    records = []
    for sample in samples:
        key = canonical_sample(sample)
        row = rows.get(key)
        group, reason = None, "metadata_not_found"
        if row is not None:
            value = row[category]
            missing = pd.isna(value) or str(value).strip().casefold() in UNKNOWN_VALUES
            # Wild specimens define the zero-captivity baseline. Record derivation.
            if missing and is_time and str(row.get("Type", "")).strip().casefold() == "wild":
                group, reason = "0d", "wild_baseline_from_Type"
            elif missing:
                reason = "group_not_recorded"
            elif is_time:
                match = re.fullmatch(r"(\d+)(?:\.0)?\s*(?:d|days?)?", str(value).strip(), re.I)
                if not match:
                    raise ValueError(f"Unrecognised captivity time for {sample}: {value}")
                group, reason = f"{int(match.group(1))}d", "recorded_time"
            else:
                group, reason = str(value).strip(), "recorded_category"
        records.append({"Sample": sample, "Metadata_ID": key, "Group": group,
                        "Included": group is not None, "Reason": reason})
    audit = pd.DataFrame(records)
    return audit.set_index("Sample")["Group"].to_dict(), audit


def group_order(values):
    def key(value):
        match = re.fullmatch(r"(\d+)d", str(value))
        return (0, int(match.group(1))) if match else (1, str(value))
    return sorted(set(values), key=key)


def display_proportions(relative, threshold):
    if not 0 <= threshold <= 1:
        raise ValueError("Relative abundance threshold must be between zero and one")
    keep = relative.mean(axis=1) >= threshold
    display = relative.loc[keep].copy()
    if (~keep).any():
        other = f"Other (<{threshold:.1%} mean)"
        if other in display.index:
            raise ValueError("Other label collides with an input taxon")
        display.loc[other] = relative.loc[~keep].sum(axis=0)
    display.loc[:, relative.isna().all(axis=0)] = np.nan
    return display * 100


def compact_letters(groups, rejected_pairs):
    """Insert-and-absorb display: significant pairs never share a letter."""
    groups = list(groups)
    columns = [frozenset(groups)]
    for a, b in rejected_pairs:
        expanded = []
        for column in columns:
            if a in column and b in column:
                expanded.extend((column - {a}, column - {b}))
            else:
                expanded.append(column)
        unique = set(c for c in expanded if c)
        columns = [c for c in unique if not any(c < other for other in unique)]
    columns.sort(key=lambda c: tuple(i for i, g in enumerate(groups) if g in c))
    letters = {g: "" for g in groups}
    for i, column in enumerate(columns):
        label = chr(97 + i) if i < 26 else f"({i + 1})"
        for group in column:
            letters[group] += label
    return letters


def rank_comparisons(frame, metrics, alpha=0.05):
    """Independent-sample KW with BH across metrics and Bonferroni within pairs."""
    groups = group_order(frame.Group.dropna())
    records = []
    for metric in metrics:
        arrays = [frame.loc[frame.Group.eq(g), metric].dropna().to_numpy() for g in groups]
        record = {"Metric": metric, "Test": "Kruskal-Wallis", "Statistic": np.nan,
                  "p_value": np.nan, "Status": "insufficient_replication"}
        if len(groups) >= 2 and all(len(a) >= 2 for a in arrays):
            if np.unique(np.concatenate(arrays)).size <= 1:
                record.update(Statistic=0.0, p_value=1.0, Status="all_values_equal")
            else:
                h, p = stats.kruskal(*arrays)
                record.update(Statistic=float(h), p_value=float(p), Status="computed")
        records.append(record)
    omnibus = pd.DataFrame(records)
    finite = omnibus.p_value.notna()
    omnibus["p_adjusted"] = np.nan
    if finite.any():
        omnibus.loc[finite, "p_adjusted"] = multipletests(omnibus.loc[finite, "p_value"], method="fdr_bh")[1]
    omnibus["Adjustment"] = "BH across reported metrics"
    pairwise, letters = [], {}
    for rec in omnibus.itertuples(index=False):
        letters[rec.Metric] = {g: "" for g in groups}
        if not np.isfinite(rec.p_adjusted) or rec.p_adjusted > alpha:
            if np.isfinite(rec.p_adjusted):
                letters[rec.Metric] = {g: "a" for g in groups}
            continue
        this = []
        for a, b in combinations(groups, 2):
            x = frame.loc[frame.Group.eq(a), rec.Metric].dropna()
            y = frame.loc[frame.Group.eq(b), rec.Metric].dropna()
            u, p = stats.mannwhitneyu(x, y, alternative="two-sided", method="auto")
            this.append({"Metric": rec.Metric, "Group1": a, "Group2": b,
                         "Statistic": float(u), "p_value": float(p)})
        adjusted = multipletests([r["p_value"] for r in this], method="bonferroni")[1]
        rejected = []
        for row, p in zip(this, adjusted):
            row.update(p_adjusted=float(p), Adjustment="Bonferroni across pairs within metric")
            if p <= alpha:
                rejected.append((row["Group1"], row["Group2"]))
        letters[rec.Metric] = compact_letters(groups, rejected)
        pairwise.extend(this)
    cols = ["Metric", "Group1", "Group2", "Statistic", "p_value", "p_adjusted", "Adjustment"]
    return omnibus, pd.DataFrame(pairwise, columns=cols), letters


def write_run_record(base, analysis, inputs=(), parameters=None, audit=None, commands=None):
    """Record the environment actually executing this analysis, not prior versions."""
    base = Path(base)
    base.parent.mkdir(parents=True, exist_ok=True)
    versions = {}
    for name in ("numpy", "pandas", "scipy", "statsmodels", "scikit-bio", "matplotlib",
                 "scikit-learn", "openpyxl", "tqdm", "requests", "lefse"):
        try:
            versions[name] = pkg_metadata.version(name)
        except pkg_metadata.PackageNotFoundError:
            versions[name] = None
    try:
        commit = subprocess.run(["git", "-C", str(Path(__file__).parent), "rev-parse", "HEAD"],
                                capture_output=True, text=True, check=False).stdout.strip()
        dirty = bool(subprocess.run(["git", "-C", str(Path(__file__).parent), "status", "--porcelain"],
                                   capture_output=True, text=True, check=False).stdout.strip())
    except OSError:
        commit = None
        dirty = None
    record = {"analysis": analysis, "python": platform.python_version(),
              "run_utc": datetime.now(timezone.utc).isoformat(),
              "executable": sys.executable, "platform": platform.platform(),
              "package_versions": versions, "git_commit": commit, "git_dirty": dirty,
              "inputs": [str(Path(p).resolve()) for p in inputs if p],
              "parameters": parameters or {}, "commands": commands or [],
              "sample_policy": "Missing groups excluded; wild Time baseline is derived and recorded"}
    try:
        record['operating_system'] = platform.freedesktop_os_release().get('PRETTY_NAME')
    except OSError:
        record['operating_system'] = platform.system()
    record['input_sha256'] = {}
    for path in record['inputs']:
        digest = hashlib.sha256()
        with open(path, 'rb') as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b''):
                digest.update(chunk)
        record['input_sha256'][path] = digest.hexdigest()
    base.with_name(base.name + "_run.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    if audit is not None:
        audit.to_csv(base.with_name(base.name + "_samples.csv"), index=False)


def pcoa_lingoes(distances):
    """Full PCoA with Lingoes correction when distances are non-Euclidean."""
    d = np.asarray(distances, dtype=float)
    n = len(d)
    center = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * center @ (d ** 2) @ center
    original = np.linalg.eigvalsh(b)
    constant = max(0.0, float(-original.min()))
    if constant <= 1e-12 * max(1.0, float(np.abs(original).max())):
        constant = 0.0
    corrected = d.copy()
    if constant:
        corrected = np.sqrt(d ** 2 + 2 * constant)
        np.fill_diagonal(corrected, 0)
        b = -0.5 * center @ (corrected ** 2) @ center
    vals, vecs = np.linalg.eigh(b)
    order = np.argsort(vals)[::-1]
    vals, vecs = np.maximum(vals[order], 0), vecs[:, order]
    coords = vecs * np.sqrt(vals)
    explained = vals / vals.sum() * 100 if vals.sum() else np.zeros(n)
    diagnostic = {"lingoes_constant": constant, "minimum_original_eigenvalue": float(original.min()),
                  "negative_eigenvalue_fraction": float(np.abs(original[original < 0]).sum() / np.abs(original).sum()) if np.abs(original).sum() else 0.0}
    return coords, explained, corrected, diagnostic


def dispersion_test(distances, samples, groups, permutations=999, seed=42):
    from skbio import DistanceMatrix
    from skbio.stats.distance import permdisp
    _, _, corrected, diagnostic = pcoa_lingoes(distances)
    result = permdisp(DistanceMatrix(corrected, ids=list(samples)), list(groups),
                      test="median", permutations=permutations, seed=seed)
    return {"Test": "PERMDISP", "Statistic": float(result["test statistic"]),
            "p_value": float(result["p-value"]), "Permutations": permutations,
            "Center": "spatial median", "Distance": "Lingoes-corrected Bray-Curtis",
            **diagnostic}

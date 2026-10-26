"""Length-weighted CNV burden for ROBIN profiles, independent of GUI/workflows."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from importlib import resources
from pathlib import Path

import numpy as np


def load_map(path):
    value = np.load(path, allow_pickle=True).item()
    if not isinstance(value, dict):
        raise ValueError(f"{path}: expected a chromosome dictionary")
    return value


def canonical(chrom):
    name = str(chrom).removeprefix('chr')
    return 'chr' + name


def load_track(path):
    out = {}
    for key, value in load_map(path).items():
        chrom = canonical(key)
        if chrom not in {f'chr{i}' for i in range(1, 23)}:
            continue
        if chrom in out:
            raise ValueError(f'Duplicate chromosome alias: {chrom}')
        array = np.asarray(value, dtype=float)
        if array.ndim != 1:
            raise ValueError(f'{chrom}: expected a one-dimensional array')
        out[chrom] = array
    if not out:
        raise ValueError('No autosomes found')
    return out



def subtract(start, end, masks):
    cursor = start
    for left, right in sorted(masks):
        if right <= cursor:
            continue
        if left >= end:
            break
        if left > cursor:
            yield cursor, min(left, end)
        cursor = max(cursor, right)
        if cursor >= end:
            return
    if cursor < end:
        yield cursor, end


def finite_runs(values):
    valid = np.isfinite(values)
    edges = np.diff(np.r_[False, valid, False].astype(int))
    return zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1))


def analyse_burden(track, width, *, delta=False, baseline=2., method='kernel', penalty=10.,
            gain=.3, loss=-.3, lengths=None, masks=None):
    """Segment finite unmasked runs; retain exact source-bin coordinates and bp weights."""
    lengths, masks = lengths or {}, masks or {}
    segments, per_chrom = [], []
    for chrom in sorted(track, key=lambda k: int(k[3:])):
        values = track[chrom]
        # CNV3 is sample minus reference, so restore expected diploid baseline.
        copies = values + baseline if delta else values.copy()
        # Negative estimates are unusable; zero is retained as a potential deletion.
        valid = np.isfinite(copies) & (copies >= 0)
        log2 = np.full(len(copies), np.nan)
        log2[valid] = np.log2(np.maximum(copies[valid], baseline * 1e-6) / baseline)
        span = len(values) * width
        if lengths:
            if chrom not in lengths:
                raise ValueError(f'Missing chromosome length: {chrom}')
            span = min(span, lengths[chrom])
        counts = dict(gain_bp=0, loss_bp=0, neutral_bp=0)
        for region_start, region_end in subtract(0, span, masks.get(chrom, [])):
            lo, hi = region_start // width, math.ceil(region_end / width)
            for first, last in finite_runs(log2[lo:hi]):
                a, b = lo + int(first), lo + int(last)
                signal = log2[a:b]
                if method == 'bins':
                    ends = range(1, len(signal) + 1)
                elif len(signal) < 2:
                    ends = [len(signal)]
                else:
                    import ruptures as rpt
                    ends = rpt.KernelCPD(kernel='linear').fit(signal).predict(pen=penalty)
                previous = 0
                for stop in ends:
                    stop = int(stop)
                    median = float(np.median(signal[previous:stop]))
                    state = 'gain' if median > gain else 'loss' if median < loss else 'neutral'
                    start = max(region_start, (a + previous) * width)
                    end = min(region_end, (a + stop) * width)
                    counts[state + '_bp'] += end - start
                    row = dict(chrom=chrom, start=start, end=end, state=state, median_log2=median)
                    segments.append(row)
                    previous = stop
        callable_bp = sum(counts.values())
        eligible_bp = sum(e - s for s, e in subtract(0, span, masks.get(chrom, [])))
        per_chrom.append(dict(chrom=chrom, represented_bp=span, masked_bp=span-eligible_bp,
                              unknown_bp=eligible_bp-callable_bp, callable_bp=callable_bp, **counts))
    totals = {k: sum(row[k] for row in per_chrom) for k in per_chrom[0] if k != 'chrom'}
    if not totals['callable_bp']:
        raise ValueError('No callable autosomal bases')
    totals['altered_bp'] = totals['gain_bp'] + totals['loss_bp']
    for state in ('gain', 'loss', 'neutral', 'altered'):
        totals[state + '_fraction'] = totals[state + '_bp'] / totals['callable_bp']
        totals[state + '_percent'] = 100 * totals[state + '_fraction']
    if lengths:
        totals['reference_autosomal_bp'] = sum(lengths.get(f'chr{i}', 0) for i in range(1, 23))
        totals['unrepresented_reference_bp'] = totals['reference_autosomal_bp'] - totals['represented_bp']
    return totals, per_chrom, segments



INPUT_FILES = ("CNV.npy", "CNV2.npy", "CNV3.npy", "CNV_dict.npy")


def input_signature(sample_dir):
    """Detect replacement, updates and deletions of any CNV input."""
    result = []
    for name in INPUT_FILES:
        path = Path(sample_dir) / name
        try:
            stat = path.stat()
            result.append((name, stat.st_mtime_ns, stat.st_size, stat.st_ino))
        except FileNotFoundError:
            result.append((name, None, None, None))
    return tuple(result)


def packaged_centromeres():
    """Read ROBIN's own centromere/satellite mask; never silently skip it."""
    data = resources.files("robin.resources").joinpath("cenSatRegions.bed").read_bytes()
    masks = {}
    for line in data.decode().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        chrom, start, end, *_ = line.split()
        start, end = int(start), int(end)
        if start < 0 or end <= start:
            raise ValueError("Invalid packaged centromere interval")
        masks.setdefault(canonical(chrom), []).append((start, end))
    if not masks:
        raise ValueError("ROBIN centromere mask is empty")
    return masks, hashlib.sha256(data).hexdigest()


def calculate_sample_burden(sample_dir, *, source="auto", method="kernel",
                            gain=.3, loss=-.3, penalty=None, baseline=2.):
    """Load a consistent trusted ROBIN snapshot and return JSON-safe results.

    Auto: CNV3, otherwise paired CNV/CNV2, otherwise CNV. CNV2 alone is not a
    tumour measurement. No zero-padding or assumption of neutral missing data.
    """
    sample_dir = Path(sample_dir)
    before = input_signature(sample_dir)
    if source not in {"auto", "sample", "difference", "paired"}:
        raise ValueError("Unknown CNV source")
    if method not in {"kernel", "bins"}:
        raise ValueError("Unknown calling method")
    metadata = load_map(sample_dir / "CNV_dict.npy")
    width = metadata.get("bin_width")
    if width is None or not np.isfinite(width) or int(width) != width or width <= 0:
        raise ValueError("CNV_dict.npy must contain a positive integer bin_width")
    penalty = float(metadata.get("penalty_value", 10) if penalty is None else penalty)
    if not all(np.isfinite(v) for v in (gain, loss, baseline, penalty)):
        raise ValueError("Settings must be finite")
    if not loss < 0 < gain or baseline <= 0 or penalty <= 0:
        raise ValueError("Require loss < 0 < gain and positive baseline/penalty")
    if source == "auto":
        source = ("difference" if (sample_dir / "CNV3.npy").exists() else
                  "paired" if (sample_dir / "CNV2.npy").exists() else "sample")
    filename = "CNV3.npy" if source == "difference" else "CNV.npy"
    track = load_track(sample_dir / filename)
    reference = None
    if source == "paired":
        from scipy.ndimage import uniform_filter1d
        reference = load_track(sample_dir / "CNV2.npy")
        if track.keys() != reference.keys() or any(len(track[k]) != len(reference[k]) for k in track):
            raise ValueError("Paired sample/reference chromosomes and lengths must match")
        track = {k: uniform_filter1d(v, size=3, mode="nearest") -
                 uniform_filter1d(reference[k], size=3, mode="nearest") for k, v in track.items()}
    reference_available = (sample_dir / "CNV2.npy").exists()
    if reference_available:
        reference = reference if reference is not None else load_track(sample_dir / "CNV2.npy")
        # Same gate as ROBIN's ploidy calling track: reference must be finite
        # and strictly positive. Unrepresented reference tails are unknown.
        for chrom, values in track.items():
            usable = np.zeros(len(values), dtype=bool)
            ref = reference.get(chrom, np.array([]))
            n = min(len(values), len(ref))
            usable[:n] = np.isfinite(ref[:n]) & (ref[:n] > 0)
            track[chrom] = np.where(usable, values, np.nan)
    masks, mask_sha = packaged_centromeres()
    totals, chromosomes, segments = analyse_burden(
        track, int(width), delta=source in {"difference", "paired"},
        baseline=baseline, method=method, penalty=penalty, gain=gain, loss=loss,
        masks=masks,
    )
    if input_signature(sample_dir) != before:
        raise ValueError("CNV files changed during calculation; calculate again")
    preview = {}
    for chrom, values in track.items():
        copies = values + baseline if source in {"difference", "paired"} else values
        indices = np.arange(0, len(values), max(1, math.ceil(len(values) / 1200)))
        indices = indices[np.isfinite(copies[indices]) & (copies[indices] >= 0)]
        signal = np.log2(np.maximum(copies[indices], baseline * 1e-6) / baseline)
        preview[chrom] = [[float((i + .5) * width), float(np.clip(y, -2, 2))]
                          for i, y in zip(indices, signal)]
    return dict(preview=preview, summary=totals, chromosomes=chromosomes, segments=segments,
                input_signature=before,
                settings=dict(sample=sample_dir.name, source=source, input=filename,
                              bin_width=int(width), method=method, kernel="linear" if method == "kernel" else None,
                              penalty=penalty, baseline=baseline, gain_threshold=gain,
                              loss_threshold=loss, centromere_mask="robin.resources/cenSatRegions.bed",
                              centromere_mask_sha256=mask_sha,
                              denominator="callable represented autosomal bases",
                              reference_mappability_mask=reference_available,
                              numpy_version=np.__version__))


def export_burden(result, kind="json"):
    """Download bytes only; no automatic writes into patient/sample folders."""
    if kind == "json":
        return (json.dumps(result, indent=2, allow_nan=False) + "\n").encode()
    if kind not in {"summary", "chromosomes", "segments"}:
        raise ValueError("Unknown export type")
    rows = [result["summary"]] if kind == "summary" else result[kind]
    text = io.StringIO()
    writer = csv.DictWriter(text, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    return text.getvalue().encode()

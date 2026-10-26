# Autosomal CNV burden

Open a sample's **Copy number (CNV)** page and find **Autosomal CNV burden** above
the regional events table. Select a profile and method, adjust the thresholds
if needed, then choose **Calculate CNV burden**. Calculation runs in a worker
thread and does not reopen BAMs. It is on demand, rather than repeated for every
live CNV update. No additional dependencies are required.

The panel shows gain, loss and total altered percentages, callable/masked/unknown
base totals, and the number of represented autosomes. An interactive genome plot
shows sampled bin signal and called segment medians, with masked regions shaded.
The plot clips signal at ±2 for display only. Its toolbar saves a PNG; the stacked
bar shows gain, loss and neutral percentages. JSON and TSV downloads retain full
numeric results, chromosome totals, interval calls, settings and mask checksum.
Downloads do not write files into the sample directory.

CNV input changes clear displayed results and disable downloads on the next CNV
refresh. Settings changes clear them immediately. A changed input during calculation
rejects the snapshot. Choose Calculate again after an update. Missing/corrupt inputs
produce an unavailable message, not a neutral or zero-burden result.

## Input interpretation

The panel reads trusted ROBIN NumPy chromosome dictionaries from the selected
sample and requires a positive `bin_width` in `CNV_dict.npy`.

| Profile | Interpretation |
| --- | --- |
| Auto | Prefer CNV3; otherwise paired CNV/CNV2; otherwise CNV |
| Sample | CNV.npy sample copy estimates |
| Difference | CNV3.npy, ROBIN's smoothed sample-minus-reference values |
| Sample + reference | Three-bin nearest-edge moving averages of CNV and CNV2, then subtraction |

CNV2 alone is a reference and cannot measure tumour burden. Paired arrays must
have matching autosomes and lengths. No zero padding is added. Existing CNV3
may already contain ROBIN padding, which cannot be identified without reference
information. Where CNV2 is available, only bins with finite positive reference
signal are usable; bins beyond its length or absent from the reference are
unknown. Without CNV2, reference mappability cannot be assessed.

For Sample, signal is `log2(sample / 2)`. For Difference/paired, signal is
`log2((sample_minus_reference + 2) / 2)`. This additive correction is not the same
as `log2(sample/reference)`. Zero-copy estimates are retained as potential complete
losses using a small logarithm floor. Negative copy estimates and nonfinite
values are excluded as unknown. Smoothing may spread missing data to neighbours.

## Calling and percentages

The default method is `ruptures.KernelCPD(kernel="linear")`, matching ROBIN's
changepoint algorithm. It operates at native bin resolution, independently on
finite unmasked stretches of each autosome. Segment medians above +0.3 are gains;
below −0.3 are losses; other segments are neutral. GUI defaults use ROBIN's
configured autosomal thresholds and workflow penalty (normally 10). The numerical
API can instead use a saved metadata penalty. Direct bin thresholding is available
as a sensitivity comparison.

ROBIN's packaged `robin.resources/cenSatRegions.bed` is always applied. It includes
centromeric satellites and acrocentric intervals. Intervals are unioned and clipped
exactly, including partial-bin overlaps; segmentation cannot cross them. A missing
mask fails the calculation. This is region masking, not ROBIN's separate telomere
breakpoint filter.

```
gain fraction = gain bp / callable represented autosomal bp
loss fraction = loss bp / callable represented autosomal bp
altered fraction = (gain bp + loss bp) / callable represented autosomal bp
```

Every usable base is counted once, with no additional whole-chromosome burden
added on top. Sex chromosomes, mitochondria and alternative contigs are excluded.
Missing chromosomes are not assumed neutral. No reference chromosome sizes are
inferred: represented span is bin count × bin width, including terminal rounding.
The panel reports **callable represented autosomal burden**, not an unqualified
percentage of the entire human genome. Match reference builds, masks, thresholds
and resolution for comparisons. The bundled mask's build must match the aligned
data; the NumPy files do not encode enough information to verify this automatically.

Thresholds and penalty require validation for the dataset. This is a new burden
metric; it does not reproduce ROBIN's arm/whole-chromosome event rules. Fine grids
can be expensive and penalty sensitivity changes with resolution. Neither purity
nor tumour ploidy is corrected. The fixed diploid baseline preserves whole-chromosome
shifts but can overstate burden for globally shifted profiles.

# Copy Numbr Load

The **Copy number (CNV)** page now highlights predicted gain regions in red and
loss regions in blue directly on the existing copy-number plot. The overlay
follows chromosome selection, zoom and linear/log scale. Toggle **Highlight
predicted gains (red) and losses (blue)** to compare with the original plot.
Centromere masking remains applied to the load calculation.

**Copy Numbr Load** metrics sit directly below that plot: gain, loss, total
altered percentage, callable/masked/unknown bases and represented autosomes.
Metrics always refer to all callable represented autosomes, even while viewing
one chromosome. Calls use the selected load profile, which can differ from the
raw sample profile underneath the overlay.

Load is calculated automatically in the background when CNV files arrive or
change. A client-scoped timer checks every two seconds; only one calculation
runs per panel at a time. Switching chromosome/scale/zoom does not rerun calling.
Changed inputs clear stale calls and downloads while a replacement is calculated.
Snapshots changed during calculation are rejected. Missing/invalid data produces
an unavailable message rather than a neutral result. Polling stops when the page
is disconnected. Calculation time depends on profile resolution and segment count;
this is an updating estimate, not a guarantee of two-second analysis latency.

Expand **Load calculation settings** to select the input/method/thresholds or
penalty. **Refresh Copy Numbr Load** allows a manual retry. **How this is
calculated** explains the assumptions. The percentage bar and JSON/TSV downloads
remain below the plot. Downloads retain interval calls, chromosome totals, settings
and mask checksum; sample directories are not modified.

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

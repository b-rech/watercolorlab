<p align="center">
  <img src="assets/watercolorlab-logo.svg"
       alt="Watercolor Lab"
       width="520">
</p>

# Algorithm Theoretical Basis Document

## Replicate consistency processor

**Module:** `watercolorlab.consistency`  
**Processor:** `ConsistencyProcessor.check`  
**Version context:** WaterColorLab 1.0  
**Purpose:** Document the statistical operations and filtering behavior
implemented by the replicate consistency processor. This document describes
the current code behavior; it is not a substitute for comparison with an
external quality-control protocol.

## 1. Scope and processing concept

The processor evaluates replicate spread across wavelength columns and can
identify and remove one suspicious curve from each inconsistent group. The
processing chain is:

1. validate analysis settings and spectral metadata;
2. select all rows, sample rows, or blank rows;
3. group curves by one or more metadata columns;
4. calculate sample standard deviation independently at each wavelength;
5. flag groups exceeding the configured threshold;
6. score curves against their group median;
7. optionally plot and remove the highest-scoring curve from each flagged group;
8. optionally repeat removal recursively.

The processor operates on the `Spectra.raw_data` table and uses the canonical
wavelength columns stored in `Spectra.wls`.

## 2. Input representation and validation

The input is a `Spectra` object, normally constructed from
`watercolorlab.io.read_files`.

| Field | Meaning |
|---|---|
| `raw_data` | Metadata and spectral measurements |
| `wls` | Canonical numeric wavelength labels |
| `is_sample` | Boolean sample mask |
| `is_blank` | Boolean blank mask |
| DataFrame index | Unique curve identifier used in results and removal |

The selected `raw_data` index must be unique. For `measurements="sample"` or
`measurements="blank"`, the required mask must be a boolean pandas `Series`,
contain no missing values, and have the same index as `raw_data`.

All selected wavelength values must be finite. Empty selected subsets are
rejected because they cannot support a meaningful consistency analysis.

The `groupby` argument accepts one column name or an iterable of column names.
Each grouping column must exist in `raw_data`.

## 3. Measurement selection and grouping

The `measurements` argument selects the analysis table:

$$
D =
\begin{cases}
\{r \in \mathrm{raw\_data}: \mathrm{is\_sample}(r)\}, &
\text{measurements}=\text{"sample"}\\
\{r \in \mathrm{raw\_data}: \mathrm{is\_blank}(r)\}, &
\text{measurements}=\text{"blank"}\\
\mathrm{raw\_data}, & \text{measurements}=\text{"all"}
\end{cases}
$$

Rows are partitioned by the selected metadata keys. Group order follows the
existing row order because pandas grouping is performed with `sort=False`.
Missing grouping values are retained as groups through `dropna=False`.

## 4. Spread criterion

For group $g$ and wavelength $\lambda$, the processor calculates the sample
standard deviation:

$$
s_g(\lambda)
  =
  \sqrt{
    \frac{1}{n_g-1}
    \sum_{i \in g}
    \left(x_i(\lambda)-\overline{x}_g(\lambda)\right)^2
  }
$$

where $n_g$ is the number of curves in the group and
$\overline{x}_g(\lambda)$ is the arithmetic mean at that wavelength.

A group is suspicious when at least one wavelength exceeds the threshold:

$$
\mathrm{suspicious}(g)
  =
  \left[
    \exists\lambda:
    s_g(\lambda) > \tau
  \right]
$$

where $\tau$ is `std_threshold`. The threshold is an absolute threshold on the
raw measurement scale. It is not normalized by wavelength, signal magnitude,
or measurement-uncertainty estimates.

Groups with fewer than two curves cannot yield a sample standard deviation.
They are not flagged by the spread criterion, and a `UserWarning` identifies
them.

## 5. Suspicious-curve score

For each curve $i$ in group $g$, the group median at wavelength $\lambda$ is:

$$
m_g(\lambda)
  =
  \operatorname{median}_{j \in g}
  \left(x_j(\lambda)\right)
$$

The curve score is its mean absolute deviation from that median:

$$
d_i
  =
  \frac{1}{N_\lambda}
  \sum_{\lambda}
  \left|x_i(\lambda)-m_g(\lambda)\right|
$$

The curve with the largest $d_i$ is selected. Exactly one curve is selected
per suspicious group in each processing pass. If scores tie, pandas'
`idxmax()` selects the first curve in the existing row order. The score is a
deviation score, not the formal median absolute deviation statistic.

## 6. Removal behavior

When `remove_suspicious=False`, the input `Spectra` object is returned and no
rows are removed. The selected identifiers are still available through
`spectra.suspicious`.

When removal is enabled, the selected raw DataFrame index labels are removed
while preserving the order of the remaining rows.

With `inplace=False`, a new object of the same `Spectra` type is constructed
from the filtered table.

With `inplace=True`, rows are removed from the associated object and
`Spectra.__init__` is called to rebuild metadata and processor instances.
Previously calculated CDOM, particulate, and correction outputs are therefore
reset. External references to the old processor instances are no longer
attached to the current object.

## 7. Recursive filtering

`recursive=True` requires `remove_suspicious=True`. After each pass, the
processor re-runs the consistency check on the filtered object until no
suspicious groups remain.

All selected raw index labels are accumulated in `spectra.suspicious` across
passes, preserving first-selection order and removing duplicate labels.
There is no maximum-removal limit; callers should choose thresholds carefully.
The caller's plotting resolution is passed unchanged to recursive calls.

## 8. Plotting and output

When `plot_suspicious=True`, one plot is generated for each suspicious group.
Curves are labeled with their raw DataFrame index, and the title identifies
the grouping label. The y-axis is labeled `Spectral measurement` because the
processor does not know whether the values represent absorbance,
transmittance, reflectance, or another measurement.

`plot_path` accepts a string or `pathlib.Path`. When supplied, it is converted
to a `Path` and created if necessary. Existing files at that location are
rejected. Output filenames are sanitized and include a stable hash of the
group label to reduce collisions.

The requested `dpi` must be a finite positive real number. Figures are closed
after display. Plot-save errors identify the group and output path while
preserving the original exception as the cause.

## 9. Outputs

The processor returns:

- the original `Spectra` object when no removal is requested;
- the filtered replacement when removal is requested with `inplace=False`;
- the same rebuilt object when removal is requested with `inplace=True`.

The `spectra.suspicious` attribute contains raw DataFrame index labels of
selected curves. It does not contain group labels, sample IDs, or filenames
unless those values are also the DataFrame index.

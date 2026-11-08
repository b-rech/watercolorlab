<p align="center">
  <img src="assets/watercolorlab-logo.svg"
       alt="Watercolor Lab"
       width="520">
</p>

# Algorithm Theoretical Basis Document

## CDOM absorption processor

**Module:** `watercolorlab.cdom`
**Processors:** `CDOMProcessor.compute_absorption` and
`CDOMProcessor.expfit`
**Version context:** WaterColorLab 1.0
**Purpose:** Document the mathematical operations implemented by the CDOM
processor. This document describes the current code behavior; it is not a
substitute for comparison with an external CDOM protocol.

## 1. Scope and processing concept

The processor derives colored dissolved organic matter (CDOM) absorption
coefficients from spectrophotometric absorbance measurements. The processing
chain is:

1. identify sample and blank spectra from the `Spectra` metadata;
2. optionally average sample replicates by sample ID;
3. restrict samples and blanks to their common wavelength range;
4. select or construct the blank correction;
5. subtract the blank absorbance;
6. optionally estimate and remove a null-point offset;
7. convert corrected absorbance to an absorption coefficient using the optical
   path length;
8. optionally fit an exponential CDOM model to each resulting absorption
   spectrum.

The implementation assumes wavelength data are available at 1-nm intervals.
Consequently, an inclusive integer wavelength interval can be represented by
the corresponding Python `range`.

## 2. Input representation and metadata

The input is a `Spectra` object, normally constructed from
`watercolorlab.io.read_files`. `read_files` returns a table indexed by file
key, with parsed metadata and wavelength columns. CDOM processing uses:

| Field | Meaning |
|---|---|
| `id` | Sample or station identifier used to group replicates |
| `is_sample` | Boolean mask identifying sample measurements |
| `is_blank` | Boolean mask identifying blank measurements |
| wavelength columns | Raw spectrophotometric absorbance values |

`Spectra` converts numeric wavelength labels to integer nanometer columns and
preserves their input order. The processor rejects duplicate sample or
custom-blank wavelength columns.

The sample and blank wavelength domains must overlap. Processing is performed
over the inclusive intersection:

$$
\lambda_\mathrm{min}
  = \max\left(\min(\Lambda_\mathrm{sample}),
              \min(\Lambda_\mathrm{blank})\right)
$$

$$
\lambda_\mathrm{max}
  = \min\left(\max(\Lambda_\mathrm{sample}),
              \max(\Lambda_\mathrm{blank})\right)
$$

where $\Lambda$ denotes a set of available wavelength columns. With complete
1-nm data, the working wavelength set is:

$$
\Lambda_\mathrm{work}
  = \{\lambda_\mathrm{min},\lambda_\mathrm{min}+1,\ldots,
      \lambda_\mathrm{max}\}
$$

## 3. Sample preparation

### 3.1 Replicate averaging

When `mean_replicates=True`, sample rows are grouped by their `id`, and the
arithmetic mean is calculated independently at every wavelength:

$$
\overline{A}_{i}(\lambda)
  = \frac{1}{n_i}
    \sum_{r=1}^{n_i} A_{i,r}(\lambda)
$$

where $i$ is the sample ID, $r$ indexes its replicates, and $n_i$ is the
number of retained replicate rows.

When `mean_replicates=False`, each sample row is retained as an individual
curve. In this mode, the original row index is preserved for the intermediate
sample table.

The resulting table is stored as `raw_samples`.

### 3.2 Sample and blank wavelength alignment

Both `raw_samples` and the selected blank table are sliced to
$\Lambda_\mathrm{work}$. No extrapolation is performed. If the sample and
blank wavelength ranges do not overlap, processing stops with an error.

## 4. Blank correction

The `blank` argument controls the blank source:

- `blank=None`: use dataset blanks when available; otherwise behave as
  `blank=False`;
- `blank=False`: use a zero blank at every sample wavelength;
- `blank=True`: average all dataset blank curves;
- `blank=<DataFrame>`: use caller-provided blank curves.

For one blank curve $B(\lambda)$, the blank-corrected absorbance is:

$$
A_\mathrm{corr}(\lambda)
  = A_\mathrm{sample}(\lambda) - B(\lambda)
$$

When dataset blanks are averaged, the blank is:

$$
\overline{B}(\lambda)
  = \frac{1}{n_B}
    \sum_{j=1}^{n_B} B_j(\lambda)
$$

and $\overline{B}$ is used in the same subtraction equation.

### 4.1 Selection among multiple blank curves

When multiple blank curves are supplied, a blank is selected independently for
each sample. The comparison region is the long-wavelength portion
$\lambda \geq 685$ nm. For sample $i$ and candidate blank $j$, the mean
absolute difference is:

$$
D_{i,j}
  = \frac{1}{N_\mathrm{ref}}
    \sum_{\lambda \in \Lambda_\mathrm{ref}}
    \left|A_i(\lambda)-B_j(\lambda)\right|
$$

where:

$$
\Lambda_\mathrm{ref}
  = \{\lambda \in \Lambda_\mathrm{work}:\lambda \geq 685\}
$$

The selected blank is the candidate with the smallest distance:

$$
j_i^\ast = \underset{j}{\arg\min}\;D_{i,j}
$$

The selected blank and its identifier are retained in `blank_curves` and
`blank_salinity`, respectively. The resulting sample-specific correction is:

$$
A_{\mathrm{corr},i}(\lambda)
  = A_i(\lambda)-B_{j_i^\ast}(\lambda)
$$

Multiple-blank selection requires at least one available wavelength at or
above 685 nm.

## 5. Null-point correction

The null-point correction removes a residual offset estimated from either one
wavelength or an inclusive wavelength interval.

### 5.1 Single-wavelength correction

For a selected wavelength $\lambda_0$, the offset for sample $i$ is:

$$
\delta_i = A_{\mathrm{corr},i}(\lambda_0)
$$

The corrected absorbance is then:

$$
A_{\mathrm{NPC},i}(\lambda)
  = A_{\mathrm{corr},i}(\lambda)-\delta_i
$$

The offsets are stored in `npc_offset`, indexed by sample.

### 5.2 Interval correction

For an inclusive correction interval
$[\lambda_i,\lambda_f]$, the offset is the arithmetic mean over all available
1-nm wavelengths in that interval:

$$
\delta_i
  = \frac{1}{N_\mathrm{NPC}}
    \sum_{\lambda=\lambda_i}^{\lambda_f}
    A_{\mathrm{corr},i}(\lambda)
$$

The same subtraction is applied to every wavelength:

$$
A_{\mathrm{NPC},i}(\lambda)
  = A_{\mathrm{corr},i}(\lambda)-\delta_i
$$

The interval endpoints must be ordered and within the available wavelength
range. The implementation uses the actual available columns in the interval;
under the 1-nm input assumption, these are the integer wavelengths from
$\lambda_i$ through $\lambda_f$.

### 5.3 Disabled correction and quality warning

When `wl_null_point_correction=None`, no offset is removed:

$$
A_{\mathrm{NPC},i}(\lambda)
  = A_{\mathrm{corr},i}(\lambda)
$$

If enabled, an offset greater than $0.0015$ absorbance units produces a
quality-control warning. The warning does not stop processing, and the
computed offset remains available through `npc_offset`.

## 6. Conversion to CDOM absorption

The corrected absorbance is converted to the Napierian CDOM absorption
coefficient using:

$$
a_\mathrm{CDOM}(\lambda)
  = \frac{\ln(10)}{L}\,
    A_{\mathrm{NPC}}(\lambda)
$$

where:

- $a_\mathrm{CDOM}(\lambda)$ is the absorption coefficient;
- $A_{\mathrm{NPC}}(\lambda)$ is the corrected decimal absorbance;
- $L$ is the spectrophotometer optical path length in meters;
- $\ln(10)\approx2.302585093$ converts base-10 absorbance to a Napierian
  absorption coefficient.

In the implementation, $L$ is supplied by `pathlength` and must be finite and
strictly greater than zero. The resulting table is stored as `absorption`.
If absorbance is dimensionless and path length is in meters,
$a_\mathrm{CDOM}$ has units of $\mathrm{m}^{-1}$.

## 7. Exponential CDOM fit

The `expfit` method fits each sample independently using the model:

$$
a_\mathrm{CDOM}(\lambda)
  = a_\mathrm{CDOM}(\lambda_\mathrm{ref})
    \exp\left[-S\left(\lambda-\lambda_\mathrm{ref}\right)\right]
$$

where:

- $\lambda_\mathrm{ref}$ is the reference wavelength;
- $a_\mathrm{CDOM}(\lambda_\mathrm{ref})$ is the fitted absorption at the
  reference wavelength;
- $S$ is the exponential spectral slope.

The reference wavelength must be present in the absorption columns. The
fitting interval is inclusive. If `wl_range=None`, all available absorption
wavelengths are used; otherwise only columns satisfying

$$
\lambda_\mathrm{start}\leq\lambda\leq\lambda_\mathrm{end}
$$

are selected. At least two wavelengths are required.

For each sample $i$, the nonlinear least-squares fit minimizes:

$$
\mathrm{SSE}_i
  = \sum_{\lambda\in\Lambda_\mathrm{fit}}
    \left[
      a_{\mathrm{CDOM},i}(\lambda)
      -a_{\mathrm{CDOM},i}(\lambda_\mathrm{ref})
       \exp\left(-S_i(\lambda-\lambda_\mathrm{ref})\right)
    \right]^2
$$

The optimization is performed by `scipy.optimize.curve_fit` with initial
parameters:

$$
a_\mathrm{CDOM}(\lambda_\mathrm{ref})=0.1,
\qquad
S=0.02
$$

The fit result is stored in `cdom_fitted`, with one row per sample ID and two
parameter columns:

- `a<wl_ref>`: fitted absorption at the reference wavelength;
- `slope_<start>_<end>`: fitted exponential slope for the selected interval.

The fit rejects non-finite absorption values, unavailable interval endpoints,
fewer than two fitting wavelengths, failed optimization, and non-finite
returned parameters.

## 8. Processing outputs

After `compute_absorption`, the associated `Spectra` object contains:

| Attribute | Content |
|---|---|
| `raw_samples` | Sample absorbance curves after optional replicate averaging |
| `absorbance` | Blank-corrected absorbance before null-point correction |
| `absorption` | CDOM absorption coefficients after null-point correction |
| `npc_offset` | Per-sample null-point offsets when correction is enabled |
| `blank_curves` | Selected sample-specific blanks when multiple blanks are used |
| `blank_salinity` | Identifiers of selected blank curves when applicable |

After `expfit`, the object additionally contains:

| Attribute | Content |
|---|---|
| `cdom_fitted` | Fitted reference absorption and exponential slope |

A new call to `compute_absorption` clears outputs from the previous CDOM
processing run before recalculating them. This prevents stale blank, offset,
absorption, or fit results from being mistaken for results from the new input
configuration.

## 9. Quality-control and validation rules

The implementation validates:

- positive, finite optical path length;
- unique sample wavelength columns;
- unique custom-blank wavelength columns;
- non-empty blank tables when explicitly provided;
- sample/blank wavelength overlap;
- availability of sample curves;
- availability of wavelengths used for multi-blank selection;
- valid single-wavelength or ordered interval null-point corrections;
- presence of the reference wavelength for exponential fitting;
- at least two wavelengths for fitting;
- finite absorption values before each fit;
- finite fitted parameters after optimization.

The implementation does not clamp negative corrected absorbance values. Such
values remain visible in `absorbance` and propagate through the conversion,
allowing the user to identify possible blank or null-point over-correction.

## 10. Implementation mapping

| Algorithm element | Implementation |
|---|---|
| Sample/blank classification | `compute_absorption`, blank and sample masks |
| Replicate averaging | `compute_absorption`, `groupby('id').mean()` |
| Common wavelength range | `compute_absorption`, `wlmin`/`wlmax` intersection |
| Single/mean blank correction | `compute_absorption`, blank resolution |
| Long-wavelength blank selection | `compute_absorption`, $\lambda\geq685$ nm comparison |
| Null-point correction | `compute_absorption`, `npc_offset` calculation |
| Absorbance-to-absorption conversion | `compute_absorption`, `np.log(10) * ... / pathlength` |
| Exponential model fitting | `expfit`, `scipy.optimize.curve_fit` |
| Fitted output table | `expfit`, `cdom_fitted` |

## 11. Limitations and protocol-comparison notes

This ATBD intentionally documents the current implementation rather than
prescribing a universal CDOM protocol. In particular, users comparing the
module with an external standard should verify:

1. whether the external protocol defines absorbance and absorption with the
   same base-10/Napierian conversion;
2. the required optical path-length units;
3. the appropriate blank-selection and salinity-correction procedure;
4. the prescribed null-point wavelength or interval;
5. the reference wavelength and fitting interval;
6. whether the protocol constrains the fitted slope or reference absorption;
7. whether the protocol requires additional quality-control thresholds or
   treatment of negative corrected absorbance.

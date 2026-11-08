<p align="center">
  <img src="assets/watercolorlab-logo.svg"
       alt="Watercolor Lab"
       width="520">
</p>

# Algorithm Theoretical Basis Document

## Particulate absorption processor

**Module:** `watercolorlab.particulate`
**Processor:** `ParticulateProcessor.compute_absorption`
**Version context:** WatercolorLab 1.0
**Purpose:** Document the mathematical operations implemented by the
particulate absorption processor. This document describes the current code
behavior; it is not a substitute for comparison with an external protocol.

## 1. Scope and processing concept

The processor retrieves particulate absorption from paired transmittance and
reflectance measurements of total and depigmented samples. For each sample ID,
the required measurements are:

- total transmittance;
- depigmented transmittance;
- total reflectance;
- depigmented reflectance.

Blank curves are selected either from the dataset or supplied explicitly.
The processor then:

1. converts input measurements to decimal transmittance/reflectance;
2. aligns the four sample curves by sample ID;
3. computes or loads blank curves;
4. applies blank-boundary and pigmentation-reference corrections;
5. normalizes curves by their corresponding blanks;
6. optionally applies the empirical tau correction;
7. computes particulate absorbance;
8. converts absorbance to optical density;
9. applies the filter geometry normalization;
10. calculates total, depigmented, and phytoplankton-associated absorption.

All calculations are wavelength-dependent unless otherwise stated.

## 2. Required input structure

The input is a `Spectra` object, normally constructed from
`watercolorlab.io.read_files`. Particulate processing requires these metadata
fields:

| Field       | Meaning                                                              |
| ----------- | -------------------------------------------------------------------- |
| `id`        | Sample/station identifier used to pair measurements                  |
| `rmode`     | Measurement mode, with one transmittance and one reflectance pattern |
| `is_sample` | Identifies sample measurements                                       |
| `is_blank`  | Identifies blank measurements                                        |
| `is_total`  | Identifies total-pigment measurements                                |

`is_total` is created by `read_files` when `depig_pattern` is supplied.
Measurements where the depigmentation pattern is present are treated as
depigmented (`is_total == False`).

Wavelength columns are expected to be numeric nanometer labels. The processor
supports a selected inclusive interval `wl_range` and requires `wl_offset`
within that interval.

## 3. Input representation

### 3.1 Absorbance input

For `type_measurement="absorbance"`, each raw value $A_\mathrm{raw}$ is
converted to decimal transmittance or reflectance:

$$
X_\mathrm{raw}(\lambda) = 10^{-A_\mathrm{raw}(\lambda)}
$$

where $X$ denotes either transmittance $T$ or reflectance $R$.

### 3.2 Decimal T/R input

For `type_measurement="tr"`, the input is assumed to already be expressed as
decimal transmittance/reflectance.

When `percentual=True`, the input values are converted from percentages:

$$
X_\mathrm{decimal}(\lambda)
  = \frac{X_\mathrm{percent}(\lambda)}{100}
$$

Percentage scaling is valid only for `type_measurement="tr"`.

For depigmented curves, a normalized value slightly greater than one can
occur when the depigmented measurement is effectively equal to, but marginally
above, its blank because of measurement uncertainty. Such depigmented
transmittance and reflectance values are set to one before the optical-depth
and particulate-absorbance calculations. Non-positive values and values above
one in total curves remain validation errors.

## 4. Measurement classification and ID alignment

The processor resolves the transmittance pattern from `trans_pattern` and
discovers the single other mode as the reflectance pattern. Mode matching is
case-insensitive, while ambiguous mode definitions are rejected.

The four sample tables are formed using the masks:

$$
\begin{aligned}
T_\mathrm{tot}^\mathrm{raw} &=
  \{\text{sample},\ \text{transmittance},\ \text{total}\}\\
R_\mathrm{tot}^\mathrm{raw} &=
  \{\text{sample},\ \text{reflectance},\ \text{total}\}\\
T_\mathrm{dep}^\mathrm{raw} &=
  \{\text{sample},\ \text{transmittance},\ \text{depigmented}\}\\
R_\mathrm{dep}^\mathrm{raw} &=
  \{\text{sample},\ \text{reflectance},\ \text{depigmented}\}
\end{aligned}
$$

The original file keys are replaced by the parsed `id` values. Each of the
four tables must contain exactly the same ordered set of IDs, with one curve
per ID. This prevents row-order differences between files from pairing
measurements incorrectly.

## 5. Blank curves

### 5.1 Dataset blanks

With `split_total_depig_blank=True`, four blank curves are calculated as
arithmetic means over the corresponding blank class:

$$
\begin{aligned}
B_{T,\mathrm{tot}}(\lambda) &=
  \operatorname{mean}\left(T_\mathrm{blank,tot}(\lambda)\right)\\
B_{R,\mathrm{tot}}(\lambda) &=
  \operatorname{mean}\left(R_\mathrm{blank,tot}(\lambda)\right)\\
B_{T,\mathrm{dep}}(\lambda) &=
  \operatorname{mean}\left(T_\mathrm{blank,dep}(\lambda)\right)\\
B_{R,\mathrm{dep}}(\lambda) &=
  \operatorname{mean}\left(R_\mathrm{blank,dep}(\lambda)\right)
\end{aligned}
$$

With `split_total_depig_blank=False`, total and depigmented blank classes are
combined within each measurement mode:

$$
B_{T,\mathrm{dep}} = B_{T,\mathrm{tot}}, \qquad
B_{R,\mathrm{dep}} = B_{R,\mathrm{tot}}
$$

If one pigmentation-specific blank class is entirely missing, the available
class is reused for the missing class. Partially missing or non-finite blank
curves are rejected.

### 5.2 Custom blanks

`custom_blanks` supplies one transmittance row and one reflectance row.
Custom blanks may be provided in absorbance or T/R representation matching
`type_measurement`. Absorbance custom blanks are converted using:

$$
B_X(\lambda) = 10^{-A_{B,X}(\lambda)}
$$

Custom blank wavelengths are normalized to integer nanometers and reordered to
the selected processing wavelengths.

## 6. Blank-boundary correction

The total transmittance curve is shifted downward when it exceeds its
transmittance blank. For sample $i$, the offset is:

$$
\Delta_{T,i}^{B}
  = \max_{\lambda}
    \left( T_{i,\mathrm{tot}}^\mathrm{raw}(\lambda)
           - B_{T,\mathrm{tot}}(\lambda) \right)_+
$$

where:

$$
(x)_+ = \max(x,0)
$$

The corrected total transmittance is:

$$
T_{i,\mathrm{tot}}^{c}(\lambda)
  = T_{i,\mathrm{tot}}^\mathrm{raw}(\lambda)
    - \Delta_{T,i}^{B}
$$

The same operation is applied independently to total reflectance:

$$
\Delta_{R,i}^{B}
  = \max_{\lambda}
    \left( R_{i,\mathrm{tot}}^\mathrm{raw}(\lambda)
           - B_{R,\mathrm{tot}}(\lambda) \right)_+
$$

$$
R_{i,\mathrm{tot}}^{c}(\lambda)
  = R_{i,\mathrm{tot}}^\mathrm{raw}(\lambda)
    - \Delta_{R,i}^{B}
$$

## 7. Depigmented-to-total reference correction

At the configured reference wavelength $\lambda_o$ (`wl_offset`), the
depigmented curves are shifted so that they equal their corresponding corrected
total curves.

For transmittance:

$$
\Delta_{T,i}^{P}
  = T_{i,\mathrm{dep}}^\mathrm{raw}(\lambda_o)
    - T_{i,\mathrm{tot}}^{c}(\lambda_o)
$$

$$
T_{i,\mathrm{dep}}^{c}(\lambda)
  = T_{i,\mathrm{dep}}^\mathrm{raw}(\lambda)
    - \Delta_{T,i}^{P}
$$

For reflectance:

$$
\Delta_{R,i}^{P}
  = R_{i,\mathrm{dep}}^\mathrm{raw}(\lambda_o)
    - R_{i,\mathrm{tot}}^{c}(\lambda_o)
$$

$$
R_{i,\mathrm{dep}}^{c}(\lambda)
  = R_{i,\mathrm{dep}}^\mathrm{raw}(\lambda)
    - \Delta_{R,i}^{P}
$$

## 8. Blank normalization

Each corrected curve is divided by its corresponding blank:

$$
\begin{aligned}
T_{i,\mathrm{tot}}(\lambda)
  &= \frac{T_{i,\mathrm{tot}}^{c}(\lambda)}
          {B_{T,\mathrm{tot}}(\lambda)}\\
T_{i,\mathrm{dep}}(\lambda)
  &= \frac{T_{i,\mathrm{dep}}^{c}(\lambda)}
          {B_{T,\mathrm{dep}}(\lambda)}\\
R_{i,\mathrm{tot}}(\lambda)
  &= \frac{R_{i,\mathrm{tot}}^{c}(\lambda)}
          {B_{R,\mathrm{tot}}(\lambda)}\\
R_{i,\mathrm{dep}}(\lambda)
  &= \frac{R_{i,\mathrm{dep}}^{c}(\lambda)}
          {B_{R,\mathrm{dep}}(\lambda)}
\end{aligned}
$$

The implementation requires each normalized curve to be finite and satisfy:

$$
0 < T(\lambda) \leq 1,\qquad
0 < R(\lambda) \leq 1
$$

## 9. $\tau$ correction

When `use_tau=True`, transmittance optical depth is calculated as:

$$
\mathrm{OD}_{T}(\lambda)
  = \log_{10}\left(\frac{1}{T(\lambda)}\right)
$$

The adjusted optical depth is:

$$
\mathrm{OD}_{T}^{*}(\lambda)
  = \mathrm{OD}_{T}(\lambda)
    - \frac{1}{2}\mathrm{OD}_{T}(750)
$$

The empirical tau factor is:

$$
\tau(\lambda)
  = 1.15 - 0.17\,\mathrm{OD}_{T}^{*}(\lambda)
$$

The implementation replaces this expression with $\tau=1$ when:

$$
\mathrm{OD}_{T}^{*}(\lambda) \leq 0.02
\quad\text{or}\quad
\mathrm{OD}_{T}^{*}(\lambda) \geq 0.7
$$

The calculation is performed independently for total and depigmented
transmittance. When `use_tau=False`, $\tau=1$.

## 10. Particulate absorbance

For total material, the implemented particulate T/R relationship is:

$$
A_{p,\mathrm{tot}}(\lambda)
 =
\frac{
  1 - T_{\mathrm{tot}}(\lambda)
  + B_{R,\mathrm{tot}}(\lambda)
    \left[T_{\mathrm{tot}}(\lambda)-R_{\mathrm{tot}}(\lambda)\right]
}{
  1
  + B_{R,\mathrm{tot}}(\lambda)
    T_{\mathrm{tot}}(\lambda)\tau_{\mathrm{tot}}(\lambda)
}
$$

For depigmented material:

$$
A_{p,\mathrm{dep}}(\lambda)
 =
\frac{
  1 - T_{\mathrm{dep}}(\lambda)
  + B_{R,\mathrm{dep}}(\lambda)
    \left[T_{\mathrm{dep}}(\lambda)-R_{\mathrm{dep}}(\lambda)\right]
}{
  1
  + B_{R,\mathrm{dep}}(\lambda)
    T_{\mathrm{dep}}(\lambda)\tau_{\mathrm{dep}}(\lambda)
}
$$

Finite negative values are clipped to zero. This accommodates small
measurement discrepancies, including cases where a depigmented curve is
slightly above its blank. A warning is emitted when clipping occurs.
Non-finite values and values greater than or equal to one are rejected before
optical-density conversion.

## 11. Optical density

Particulate absorbance is converted to optical density:

$$
\mathrm{OD}_{p}(\lambda)
  = \log_{10}\left(\frac{1}{1-A_p(\lambda)}\right)
$$

This requires $A_p(\lambda) < 1$.

## 12. Filter geometry normalization

The filter area-to-volume normalization factor is:

$$
F = \frac{\pi d^2}{4V}
$$

where:

- $d$ is the filter clearance diameter in metres;
- $V$ is the filtered volume in cubic metres.

The processor accepts either one $(V,d)$ pair for all samples or a table
indexed by sample ID with `volume` and `diameter` columns.

## 13. Final particulate absorption products

The empirical conversion from particulate optical density to absorption is:

$$
a_{p,\mathrm{tot}}(\lambda)
  = \ln(10)\,0.719\,
    \left[\mathrm{OD}_{p,\mathrm{tot}}(\lambda)\right]^{1.2287}
    F
$$

$$
a_{p,\mathrm{dep}}(\lambda)
  = \ln(10)\,0.719\,
    \left[\mathrm{OD}_{p,\mathrm{dep}}(\lambda)\right]^{1.2287}
    F
$$

The phytoplankton-associated component is calculated by subtraction:

$$
a_{p,\mathrm{phyto}}(\lambda)
  = a_{p,\mathrm{tot}}(\lambda)
    - a_{p,\mathrm{dep}}(\lambda)
$$

The resulting DataFrames are stored on the shared `Spectra` object as:

- `absorbance_total`;
- `absorbance_depig`;
- `part_abs_total`;
- `part_abs_depig`;
- `part_abs_phyto`.

Rows are indexed by sample ID and columns by wavelength.

## 14. Quality-control and failure conditions

The processor rejects:

- missing required metadata;
- malformed, reversed, or empty wavelength ranges;
- non-finite raw sample measurements;
- missing or non-finite blank curves;
- incomplete or duplicated sample-ID pairings;
- invalid custom-blank mode labels or duplicate wavelengths;
- normalized T/R values outside $0 < X \leq 1$;
- non-finite particulate absorbance;
- particulate absorbance $\geq 1$;
- non-finite or non-positive volume and diameter values;
- duplicate geometry-table sample IDs;
- missing wavelength 750 when tau correction is enabled.

Finite negative particulate absorbance is not treated as a fatal error. It is
clipped to zero with a warning because small negative values can result from
measurement uncertainty around a blank-equivalent depigmented curve.

## 15. Implementation-to-equation map

| Processing stage                           | Main implementation location                         |
| ------------------------------------------ | ---------------------------------------------------- |
| Input conversion and wavelength selection  | `compute_absorption`, early processing               |
| Mode classification and ID alignment       | `compute_absorption`, measurement extraction         |
| Blank calculation                          | `compute_absorption`, blank selection                |
| Blank-boundary correction                  | `compute_absorption`, total-curve offsets            |
| Pigmentation reference correction          | `compute_absorption`, `wl_offset` corrections        |
| Blank normalization                        | `compute_absorption`, `Tt`, `Td`, `Rt`, `Rd`         |
| Tau correction                             | `compute_absorption`, `use_tau` branch               |
| Particulate absorbance                     | `compute_absorption`, T/R relationship               |
| Optical density and geometry normalization | `compute_absorption`, final conversion               |
| Final products                             | `part_abs_total`, `part_abs_depig`, `part_abs_phyto` |

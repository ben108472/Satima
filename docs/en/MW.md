# Satima: Microwave (MW)

> This is the detailed documentation for **MW (microwave radiometers)**; its content comes from the original Satima-MW. For how to run the integrated version, the data listing, the layered config files
> (`Satima.toml` + `Satima-MW.toml`), the interface language, the image modes (native / dpi) and the continuity rules for polar orbiters, see [README.en.md](../../README.en.md).
> Of the settings mentioned below, the common sections (`[general]`, `[map]`, `[coastline]`, `[image]`, `[title]`, `[annotation]`, `[output]`…)
> are written in `Satima.toml` and those that only MW has in `Satima-MW.toml`; any section can be overridden in `Satima-MW.toml`.
> In this document the "native resolution" behaviour refers to `[image] mode = "native"`; with `mode = "dpi"` the image size is determined by the dpi instead, and the sampling method is unchanged.


Brightness temperature images from passive microwave radiometers: brightness temperature of each channel, PCT (polarization-corrected temperature), RGB composites at 37 and 85\~91 GHz, in the NRL and RAMMB styles (RAMMB's 85\~92 GHz is a single channel with a color scale).
The architecture follows Satima-GEO (config files, menus, maps and saving are all the same); the difference is that the data is a polar-orbiting "swath", and the core of the sampling (interpolation) was rewritten.

Supported data:

| Level | Source | Instruments |
|---|---|---|
| L1C | NASA PPS | All radiometers of the GPM constellation: GMI, AMSR2, AMSR3, SSMIS, WSF-M MWI, MHS, ATMS… (one and the same format) |
| L1B | JAXA G-Portal | AMSR2 (GCOM-W1), AMSR3 (GOSAT-GW), GMI (GPM) |
| L1B | NOAA CLASS | SSMIS SDR / TDR (F16\~F19), SSM/I SDR / TDR (F10\~F15) |
| L1 | NSMC | MWRI (FY-3A\~D), MWRI-II (FY-3F/H), MWRI-RM (FY-3G) |
| L1c | NSMC | MWRI orbital L1c (FY-3D data-assimilation product, binary) |

## Running

```
python Satima.py                   # search for data in data_dir of the config file
python Satima.py -d D:\data\mw     # set the data directory
python Satima.py -f PPS-L1C        # search one format only
python Satima.py -o D:\png         # set the output directory
python Satima.py --list            # list the available formats
```

Workflow: select data (several allowed) → map projection (Enter = MERC) → map range (south north west east) → products (several allowed) → colormap →
image width (Enter = native sampling interval; `x3` = 3 times; or input the number of pixels directly).

```
 1   GCOM-W1    AMSR2   L1C   2025/09/22 04:36-06:15Z   PPS   10 ch
 2   GOSAT-GW   AMSR3   L1C   2026/10/20 03:35-05:13Z   PPS   15 ch
Input the index (several: (3,5) or [1,3,5]): 1
Input the projection (Enter = MERC):
Map range (south north west east): 14 26 116 128
Products:
  Brightness temperature
     1.10.65V   2.10.65H   3.18.7V   4.18.7H   5.23.8V   6.23.8H   7.36.5V   8.36.5H   9.89V   10.89H
  PCT (polarization-corrected)
    11.PCT10.65   12.PCT18.7   13.PCT36.5   14.PCT89
  RGB
    15.NRL-37     16.RAMMB-37   17.NRL-89
Input the product numbers (several allowed, e.g. 13 17): 10 14 15 17
89H colormaps: 0.ALL 1.GRAY 2.NRL 3.NRL2 4.RAMMB 5.MWNRL
```

* Columns of the listing: satellite, instrument, format (L1C / L1B / SDR / TDR…), time, **data source** (PPS, JAXA, CLASS, NSMC, EUMETSAT), number of channels,
  notes (the number of files when several files are joined into one entry automatically, e.g. `21 files` for ATMS; `19 files (2 missing)` when files are missing in between).
* **How the map range is input** (`[map] range_input`): with `"bounds"` (default) you input "south north west east"; `"center"` is the same as Satima-PO: first input the center
  "latitude longitude", then the range (degrees of latitude and longitude outward from the center; one number = both the same); for example center `20 122` and range `6` gives 14\~26N, 116\~128E.
  Whichever is set, 4 numbers are always taken as south north west east and 2 numbers as the center, so the two can be mixed at any time.
* One file is one stretch of orbit (a whole or half orbit), and you have to know yourself which entry passes over the place you want to see; when it does not, `[skip] the swath does not pass over this map range` is printed.
* The time in the annotation below the image and in the file name is **the time the swath passes over the map range** (the median scan time of the pixels inside the range), not the start time of the file.
* When several entries are selected, one image is drawn for each; products are chosen from the menu of the first entry, the other entries are matched by name (89V, PCT89, NRL-37…), and products they do not have are skipped.
* For large maps of the globe or of a whole ocean basin, changing `[map] resolution` to `"l"` is recommended (basemap takes tens of seconds to read the global coastline at `i` or finer).

Settings are in `Satima.toml` (common) and `Satima-MW.toml` (those that only MW has and those to override; UTF-8, Chinese text is allowed); items not written use the defaults in `satima/core/config.py` and `satima/mw/config.py`.
Directories (`data_dir`, `output_dir`, `colormap_dir`) may be absolute or relative paths; relative paths are relative to the directory containing Satima.py.
The map projections MERC / CYL / LCC and the map backends basemap / cartopy are used as in Satima-GEO (there is no full disk FD).

Required packages: numpy, scipy, h5py, matplotlib, basemap or cartopy, mpkit (reads .cmap colormaps);
numba is optional (with it, interpolation and saving are several times faster).

## Directory layout

```
Satima.py                     entry point (shared by the three families)
Satima.toml, Satima-MW.toml   common config file, MW config file
satima/core/        core shared by the three families: config, session, ui, maps, canvas, colormaps, naming, segments (continuity rules), pngio, parallel
satima/mw/          MW code
  config.py           MW's own defaults and checks
  family.py           MW workflow: rows of the listing, multi-file merging, range, product menu, drawing and saving image by image
  reader.py           reader base class SwathReader, Channel / Scene / Geo, multi-file merging
  swath.py            swath sampling: fractional coordinates, Gaussian distance weighting, index-space interpolation, A/B scan merging   (_swath_nb.py is the numba version of the core)
  products.py         products: brightness temperature, PCT, RGB
  grid.py             output grid
sensors/mw/         readers (only files with '-' in the name are loaded)
  PPS-L1C.py  GMI-L1B.py  TMI-L1B.py  pps_common.py
  AMSR2-L1B.py  AMSR3-L1B.py  AMSRE-L1B.py  amsr_common.py  MWRI-L1.py  MWRI-L1C.py
  SSMIS-SDR.py  SSMIS-TDR.py  SSMI-SDR.py  ATMS-SDR.py  AMSUA-L1B.py
colormap/           colormaps (.cmap, unit K; the three families share the same folder)
data/  plot/        default data and output directories
```

## Products

### Single-channel brightness temperature

Every channel can be plotted; channel name = frequency + polarization, e.g. `89V`, `36.5H`, `183.31±7V` (`±` is written `+-` in file names).
For cross-track scanning instruments (MHS, ATMS) the polarization varies with the scan angle, and they are labelled `QV` / `QH`. The first item of the colormap menu is always GRAY
(`[gray]`: 150\~300 K, cold is white); the rest are grouped by frequency and defined in `[colormap]`:

| Group | Frequency | Default colormaps |
|---|---|---|
| `tb89` | 80\~95 GHz | NRL (`nrl89`, 165\~280 K), NRL2 (`89h`, 105\~305 K), RAMMB (`rammb89`, 72\~323 K), MWNRL (150\~320 K) |
| `tb37` | 30\~40 GHz | NRL (`37h`, 125\~300 K), MWNRL |
| `tb150` | 140\~200 GHz | NRL (`150h`, 110\~310 K), MWNRL |
| `tb` | others | MWNRL, NRL37 |
| `pct89` / `pct37` / `pct` | PCT | NRL (`nrlpct` 160\~260 K / `37pct` 230\~280 K), NRL2 (`89pct` 105\~280 K) |

`nrl89` and `nrlpct` are the NRL89 and NRLPCT that were originally in mpkit; `mwnrl` is the colormap from the old Satima-SSMIS program;
`89h`, `37h`, `150h`, `89pct` and `37pct` were made from the color scale definitions of NRL GeoIPS (open source).

`rammb89` is the color scale of the RAMMB/CIRA TC Realtime "85 – 92 GHz Brightness Temperature" product. This RAMMB product is not an RGB
but a single-channel brightness temperature (horizontal polarization is preferred; vertical polarization is used only for instruments without horizontal polarization) with a color scale, so the way to use it is to select `89H`
(`85.5H` for SSM/I, `91.665H` for SSMIS) and then the colormap RAMMB. RAMMB has not published the values of the color scale; here they were read from the color bar of a finished image
(the AMSR2 image of Ragasa at 2025/09/22 0505Z; the color bar is labelled 72\~323 K):

| Brightness temperature (K) | 70 | 160 | 225 | 265 | 280 | 300 | 325 |
|---|---|---|---|---|---|---|---|
| Color | white | red (255, 0, 0) | dark blue (0, 50, 100) | gray (195, 195, 195) | light green (200, 255, 200) | yellow (255, 255, 0) | white |

Colors change linearly between the points; compared with the color bar pixel by pixel, the mean difference is less than 1 level (0\~255) and the maximum is 13 levels (around 280\~290 K).
RAMMB also has 37 GHz, PCT and other products (see its description page), but no finished images can be seen on the typhoon pages, so there is no way to obtain their color scales and they were not made.

### PCT (polarization-corrected temperature)

`PCT = (1 + θ) × V − θ × H`: the sea surface has low emissivity and V and H differ greatly; PCT removes this difference, and the low values that remain are precipitation / ice scattering.
It is available only for frequencies that have both V and H channels (PCT10.65, PCT18.7, PCT36.5, PCT89…). θ is in `[pct]`:

| Frequency | Default θ | Source |
|---|---|---|
| 10 GHz | 1.5 | Cecil and Chronis (2018) |
| 19 GHz | 1.4 | same as above |
| 37 GHz | 1.15 | same as above (NRL uses 1.181, Lee et al. 2002) |
| 85\~91 GHz | 0.7 | same as above (0.818 in Spencer et al. 1989) |

### RGB composites

All recipes are written in the config file under `[rgb.<name>]`; you can change the ranges or add recipes yourself: `band` sets the frequency (`tb37` / `tb89` …),
each of the three color guns takes `pct` / `v` / `h`, linearly stretched to `range`, optionally with `inverse` (inverted) and `gamma`.

| Name | R | G | B | Source |
|---|---|---|---|---|
| NRL-37 | PCT37 (θ = 1.181) 260\~280 K inverted | 37V 180\~300 | 37H 160\~300 | NRL GeoIPS `pmw_color37` (Lee et al. 2002) |
| NRL-89 | PCT89 (θ = 0.818) 220\~310 K inverted | 89H 240\~300 | 89V 270\~290 | NRL GeoIPS `pmw_color89` |
| RAMMB-37 | PCT37 260\~280 K inverted | 37V 160\~300 | 37H 180\~310 | Description on the RAMMB/CIRA TC Realtime web page |

* 37 GHz: sea surface dark green, low-level cloud and rain cyan, deep convection pink. 89 GHz (NRL): deep convection red.
* RAMMB has no RGB product at 85\~92 GHz (it is a single channel with a color scale, see the colormap `rammb89` above); the values of RAMMB-37 come from
  the description of the color composite in the "Synthetic 37 GHz" section of its description page.
#### Recipes with arbitrary channels (for sounders; experimental)

The recipes above all need V and H at the same frequency. Cross-track scanning sounders (ATMS, MHS, AMSU-A / B, SAPHIR) have no dual polarization, so there is another kind of
recipe written without `band`: each color gun specifies the frequency directly with `channel`, and can also subtract another channel with `minus`.

```toml
[rgb.RGB-150]
red = { channel = "89", range = [220.0, 280.0], inverse = true }
green = { channel = "150", range = [180.0, 380.0] }
blue = { channel = "183±7", range = [200.0, 340.0] }
```

`channel` takes the closest channel in the same frequency class (within 12%): `"150"` matches 150 / 157 / 165.5 / 166 GHz,
`"183±7"` matches ±7 / ±6.6 / ±6.8 (the 190.31 GHz of MHS also counts as ±7), `"89"` matches 85.5 / 88.2 / 89 / 91.665; a polarization can be added (`"89H"`),
and H is preferred when none is given. A recipe appears in the menu only when the channels of all three color guns are found, so imagers that have these channels (SSMIS, GMI, AMSR3) list it too.

One **experimental** recipe, `RGB-150`, is included by default — **it is not an official NRL / RAMMB recipe**. In 2004 NRL had an "89 / 150 GHz false-color composite" for AMSU-B
(Hawkins et al., AMS 26th Hurricanes 3D.2; the AMSU-B products of 2010 still had color), but the recipe was not published; the current GeoIPS has only single-channel products for these instruments.

| Name | R | G | B | Appearance | Applies to |
|---|---|---|---|---|---|
| RGB-150 | 89 GHz 220\~280 K inverted | 150 GHz 180\~380 | 183±7 200\~340 | Like NRL-89: sea surface cyan-green, deep convection red | ATMS, MHS, AMSU-B, SSMIS, GMI, AMSR3 |

* Red uses 89 GHz directly: without polarization correction, a dry, cold sea surface (mid and high latitudes, dry-air zones) is also on the cold side and turns reddish / pink, which is hard to tell from convection.
  (Using the difference of two channels instead avoids this, e.g. `red = { channel = "89", minus = "150", range = [-20.0, 50.0] }`: only ice scattering gives positive values.)
* The ranges were set from the brightness temperature distributions of the sample cases (Humberto, Sinlaku, Irma, etc.) and can be adjusted in the config file.
* When a recipe uses a failed channel, it is reported as skipped.

* When a frequency class has several V/H pairs (36.75 and 37.3 GHz of WSF-M MWI), the pair closest to the nominal frequency (37.0 / 89.0) is taken.
* 85 GHz (SSM/I) and 91.655 GHz (SSMIS) both belong to `tb89`.

## Interpolation (`[swath] interpolation`)

Microwave sampling is coarse (5\~25 km), and the intervals in the two directions are often different:

| Instrument | Along scan | Scan line spacing | Ratio |
|---|---|---|---|
| WSF-M MWI | 4.3 km | 15.1 km | 3.5 |
| AMSR2 / AMSR3 89 GHz (A or B only) | 4.5 km | 10 km | 2.2 |
| GMI | 6.4 km | 13 km | 2.0 |
| AMSR2 / AMSR3 89 GHz (A/B merged), SSM/I 85 GHz, SSMIS 91 GHz | — | — | about 1 |
| SSMIS 19 / 22 / 37 GHz | 25 km | 12.5 km | 0.5 (the coarser one is the along-scan direction) |

Each output pixel is first converted to fractional coordinates in the swath (scan number, pixel number): the latitude / longitude of the swath are converted to 3-D coordinates on the unit sphere
(unaffected by the date line and the polar regions), a KD-tree finds the nearest sample, then the fractional coordinates are solved from the local vectors in the along-track / along-scan directions at that sample,
and refined by a Newton correction on the bilinear surface of the cell the pixel falls in. The value is then taken according to the setting:

| Setting | Description |
|---|---|
| `"gauss"` (default) | **Gaussian distance weighting** (NRL's approach): pixel value = average of the surrounding samples weighted by actual distance. Equally smooth in every direction; no jagged edges even for instruments with asymmetric sampling |
| `"cubic"` | Cubic Hermite in the index space of the swath (slopes are limited, no overshoot): keeps the detail along the scan direction, but with asymmetric sampling oblique structures look jagged |
| `"bilinear"` | Bilinear interpolation in index space |
| `"nearest"` | **Interpolation off**: each sample is drawn as one block, and what you see is the original sampling |

* Outside the swath and where scans are missing (the spacing between adjacent scans exceeds 2.5 times the normal value) there is no data (background color); gaps are not bridged.
* The fractional coordinates are a smooth function of the output position, so they are actually solved only at grid nodes every few pixels, and the rest are interpolated from the nodes
  (the node spacing is taken as large as it can be while staying within half a sampling interval in output pixels; processing in strips on multiple cores gives a result pixel-identical to the whole image on a single core).
* The "Minimum inside the map range" printed for colormap products is the extreme of the original sample values, not of the interpolated values.
* Interpolation only makes the picture smooth; it does not add resolution.

### Gaussian distance weighting (`"gauss"`, `[swath] gauss_sigma`)

Index-space interpolation (cubic, bilinear) interpolates between adjacent samples in each of the two directions separately, which works well when the sampling intervals are similar; but when the intervals differ greatly
(a ratio of 2 or more in the table above), structures that cross the scan lines obliquely, such as eyewalls and rainbands, are several pixels apart in position on two adjacent scan lines,
while the interpolation uses the scan lines above and below at "the same pixel position", so an oblique edge is drawn as jagged steps running along the scan lines
(most obvious in the eye in WSF-M MWI). This happens because the detail along the scan direction is finer than the scan line spacing while the other direction is undersampled; switching to a higher-order interpolation does not help either.

NRL's approach (the open-source [GeoIPS](https://github.com/NRLMMD-GEOIPS/geoips), which produces the microwave products on the NRL tropical cyclone web page) is
not to interpolate on the swath but to apply Gaussian distance weighting in geographic space:

```
pixel value = Σ w·TB / Σ w        w = exp(−d² / σ²), d = distance from the sample to the pixel
```

GeoIPS settings (source code checked in 2026/11): the interpolator `interp_gauss` calls pyresample's `kd_tree.resample_gauss`
(the nearest 8 samples for each pixel); `sigmaval` is 10000 m for the 89 GHz single channels, 89pct, color89, 37pct and color37,
and 25000 m for the H and V single channels at 19 / 37 GHz; the radius of influence is given by each reader (AMSR2 89 GHz 10 km, GMI 12.5 km, SSMIS imager channels 15 km);
the output grid of the tropical cyclone images is 1 km. It is equally smooth in every direction (the resolution is unified to the coarser direction), so the jagged edges disappear naturally.

`"gauss"` in Satima-MW is the same formula, with three differences:

* **σ can follow the sampling interval.** `gauss_sigma = "x0.7"` (default) = 0.7 times the sampling interval of the coarser direction at that location:

  | Instrument / channel | Interval of the coarser direction | σ |
  |---|---|---|
  | WSF-M MWI | 15.1 km | 10.6 km |
  | GMI | 13 km | 9 km |
  | SSMIS 91 GHz | 12.5 km | 9 km |
  | SSMIS 19 / 22 / 37 GHz | 25 km | 17.5 km |
  | AMSR2 / AMSR3 89 GHz (A/B merged) | 5 km | 3.5 km |
  | MHS, ATMS (cross-track scanning) | 16 km at nadir, increasing toward the edge | increases with the scan position |

  NRL uses a fixed 10 km, which is roughly 0.7 times for MWI, GMI and SSMIS 91 GHz; but the 89 GHz of AMSR is sampled at 5 km, and 10 km blurs it a lot,
  while the 37 GHz of SSMIS is sampled at 25 km, and 10 km makes it look like separate blobs. To match NRL, write a fixed number of kilometres, `gauss_sigma = 10`.
  The smaller the factor, the sharper: at `"x0.5"` narrow rainbands that cross the scan lines obliquely start to turn into a string of beads.
* **All samples within 2.5 times σ are included** (not just the nearest 8). MWI has a sample every 4.3 km along the scan, so the nearest 8 are almost all on the same
  scan line, and a slight discontinuity appears whenever the set of neighbours changes; only including all of them gives a complete Gaussian smoothing. The samples are located directly from the fractional coordinates,
  so there is no need to query the KD-tree again for every pixel.
* **The samples of the A and B scans are weighted together** (89 GHz of AMSR).

The costs to be aware of:

* The resolution drops to that of the coarser direction, and single-sample extremes are flattened. Example: in the sample case the lowest original sample of MWI 89H is 109.3 K, while the lowest in the `"x0.7"` image is
  146.7 K (142.8 K with a fixed 10 km); AMSR2 89H (A/B merged, σ 3.5 km) 110.1 → 121.7 K.
  The "Minimum inside the map range" is still based on the original sample values. Use `"cubic"` to see the detail along the scan direction, and `"nearest"` to see the original sampling.
* When only one sample or one scan line is missing, it is filled in from the surroundings; pixels farther than 1.05 times the coarser interval from the nearest valid sample have no data.
* CIRA / RAMMB have not published how they resample their microwave images (a 2024 AMS abstract only mentions that they switched to processing the raw data themselves),
  so only NRL is followed here.

### Index-space interpolation (`"cubic"`, `"bilinear"`, `"nearest"`)

Interpolation is done in the regular index space of the fractional coordinates, each of the two directions being computed from truly adjacent samples.

* `"cubic"` is cubic Hermite with limited slopes (monotone): interpolated values never go outside the range of the adjacent sample values, so sharp edges such as eyewalls and coasts
  show no spurious cold / warm rings from overshoot.
* At the sample positions the original values are reproduced exactly (gauss does not do this).
* Pixels that cannot be computed because of invalid values around them automatically fall back to a lower-order method (cubic → bilinear → nearest).
* The 19 / 22 / 37 GHz of SSMIS has only one sample every 25 km along the scan direction, so structures narrower than a rainband appear as separate blobs; this is caused by the sampling itself.

### A / B scan merging (`[swath] merge_ab`)

The 89 GHz of AMSR2 and AMSR3 has two sets of feeds, A and B, each with one scan every 10 km along track; the two sets alternate at the center of the swath (about 5 km combined),
but almost overlap toward the edge (their off-axis angles differ and so do the radii of their scan arcs: in the sample the A swath is 1625 km wide and the B swath 1599 km), so it is not a fixed, evenly spaced interleave.
With `merge_ab = true` (default) they are merged into one channel `89V` / `89H`: fractional coordinates are computed for A and B separately, interpolation is done first along the scan direction on each scan line,
and then along track using the two nearest scan lines above and the two below the pixel from the **union** of the A and B scan lines, so both the center and the edge are handled correctly.
(The above is how index-space interpolation does it; with `"gauss"` the A and B samples are simply weighted together and need no separate handling.)
With `false` they are listed separately as `89V-A`, `89V-B`, `89H-A`, `89H-B`.

The 85 GHz of SSM/I also has A and B scans, but their geometry is the same (the same feed), so the reader simply interleaves them in order into (2 × number of scans, 128).

## Multi-file merging (`[swath] merge`)

When several files of the same satellite and the same instrument that are contiguous in time (gaps of no more than `merge_gap` seconds) are selected together, you are asked whether to join them and draw one image
(NASA near-real-time L1C is one file every 10 minutes, JAXA L1B is one file per half orbit, and a typhoon may happen to be cut across two files). Merging concatenates along the scan direction
and removes overlapping scans by scan time (each JAXA L1B file has 30 extra scans at each end), so the interpolation is continuous across the seam.
`merge = true` merges directly, `false` draws each one separately.

An ATMS file is only 32 seconds long and does not go through this question: the reader directly joins files of the same satellite with gaps within 10 minutes into one entry, the notes column of the listing shows the number of files,
and files missing in between do not prevent joining; the number missing is shown (the threshold for being adjacent or not is half the file length, see "Continuity rules for polar orbiters" in the README).
At a gap the spacing between adjacent scans exceeds the normal value, so it is left blank when sampling and is not smoothed over by interpolation.

## Formats

| Reader | File name pattern | Description |
|---|---|---|
| `PPS-L1C` | `1C.<satellite>.<instrument>.*.nc` / `.HDF5` / `.RT-NC` / `.RT-H5` | XCAL intercalibrated brightness temperature Tc. The channel table is parsed from the `Tc` long_name of each swath, so any instrument of the GPM constellation can be read. For the instruments covered and the quality flag see "PPS L1C" below |
| `GMI-L1B` | `GPMCOR_GMI_*_1B*`, `1B.GPM.GMI.*` | The files from JAXA G-Portal are the NASA 1B files; brightness temperature Tb, 9 channels in S1 and 4 in S2 |
| `TMI-L1B` | `1B.TRMM.TMI.*.nc` / `.HDF5` | TRMM TMI 1B from NASA PPS (`.../gpmdata/YYYY/MM/DD/1B/`) or GES DISC; same structure as GMI 1B: S1 10.65 GHz, S2 19.35 / 21.3 / 37.0 GHz, S3 85.5 GHz |
| `AMSR2-L1B` | `GW1AM2_*_L1SGBTBR_*.h5` | 6.925\~89 GHz. The file only has the latitude / longitude of 89A / 89B; those of each low-frequency band are derived from 89A with the co-registration parameters A1, A2 in the file header (the formula is in the reader's notes) |
| `AMSRE-L1B` | `PM1AME_*_L1SGBTBR_*.h5` | L1B of Aqua AMSR-E, the version reprocessed by JAXA into the **AMSR2 format** (G-Portal). The structure and the derivation of the low-frequency latitude / longitude are the same as for AMSR2; there is no 7.3 GHz; the A scan of 89 GHz failed from 2004/11, and only the B scan is used automatically. The old HDF4 format is not supported |
| `AMSR3-L1B` | `GGWAM3_*1BTBB*.nc` | 6.925\~183 GHz (adds 10.25, 165.5, 183.31±3 / ±7 GHz); each band carries its own latitude / longitude |
| `MWRI-L1` | `FY3?_MWRI?_*_L1_*.HDF` | FY-3A\~D MWRI: 10.65\~89 GHz V/H, 10 channels, one set of latitude / longitude. FY-3F/H MWRI-II: also 50.3\~53.75 GHz V/H and four 118.75 GHz channels (22 in total); FY-3G MWRI-RM: plus 165.5 and 183.31 GHz (26 in total); window and sounding channels each have their own latitude / longitude (see "Fengyun-3 MWRI" below) |
| `MWRI-L1C` | `FY3?_MWRI?_ORBT_L2_MIP_*.L1c` | FY-3D MWRI orbital L1c (NSMC data-assimilation product): binary, one 152-byte record per sample, 266 samples per scan; same 10 channels as MWRI-L1 |
| `SSMIS-SDR` | `NPR.SDRN.S[A-D].*` | Big-endian binary. Imager scenes (91, 150, 183 GHz, 180 points), environmental scenes (19, 22, 37 GHz, 90 points), LAS (50\~60 GHz sounding, 60 points, one scan in every 3) |
| `SSMI-SDR` | `NPR.SDRR.S[4-9].*`, `NPR.TDRR.S[4-9].*` | SDR brightness temperature and TDR antenna temperature use the same reader. 19 / 22 / 37 GHz 64 points; the A and B scans of 85 GHz are interleaved into 128 points |
| `SSMIS-TDR` | `NPR.TDRN.S[A-D].*` | Antenna temperature. Channels and scan structure as in SDR, but 19 / 22 GHz and 37 GHz each have their own latitude / longitude (SDR shares the 37 GHz position) |
| `ATMS-SDR` | `SATMS_*.h5`, `TATMS_*.h5` | SDR brightness temperature and TDR antenna temperature of JPSS ATMS (S-NPP, NOAA-20, NOAA-21), 22 channels (including the 50\~57 GHz temperature sounding channels). A file is only 32 seconds long; files of the same satellite that are contiguous in time are joined into one entry automatically. **The latitude / longitude are in separate `GATMO_*.h5` files, which must be in the same folder** (see below) |
| `AMSUA-L1B` | `AMSA_xxx_1B_M0x_*` (`.nat` or `.zip`) | EUMETSAT EPS native format of MetOp-A/B/C AMSU-A; the zip downloaded from the Data Store can be read directly. 15 channels (23.8, 31.4, 50.3\~57.29 GHz sounding, 89 GHz), 30 fields of view, 48 km at nadir. Radiances are converted to brightness temperature with the Planck function |

Details of the file structure of each format are written in the notes at the top of each reader. TDR is antenna temperature (no antenna pattern correction; it differs from brightness temperature by a few K),
and such channels are marked `TA` in the menu; PCT and RGB can still be made, but the values are slightly off.

### Fengyun-3 MWRI

Formats follow the MWRI "L1 data instruction" documents of each satellite on the NSMC document page
(<https://satellite.nsmc.org.cn/DataPortal/cn/support/document.html>), the FY-3A/B L1 Product Format and the FY-3D MWRI orbital L1C format document.

| Satellite | Instrument | Brightness temperature | Lat / lon | Time |
|---|---|---|---|---|
| FY-3A/B | MWRI | top-level `EARTH_OBSERVE_BT_10_to_89GHz` (10, scans, samples) | top level | `Scan_Time_and_Period` columns 1, 2 (days, ms) |
| FY-3C | MWRI | `Data/EARTH_OBSERVE_BT_10_to_89GHz` | `Geolocation/` | `Data/Scan_daycnt`, `Scan_mscnt` (lower case) |
| FY-3D | MWRI | `Calibration/EARTH_OBSERVE_BT_10_to_89GHz` | `Geolocation/` | `Calibration/Scan_Daycnt`, `Scan_Mscnt` |
| FY-3F/H | MWRI-II | `Calibration/EARTH_OBSERVE_BT` in `Window Channel/` and `Sounding Channel/` (scans, 520, channels) | `Geolocation/` of each group | `Scan_Daycnt`, `Scan_Mscnt` (x 0.1 ms) of each group |
| FY-3G | MWRI-RM | `S1/Data/EARTH_OBSERVE_BT_10_to_89GHz`, `S2/Data/EARTH_OBSERVE_BT_50_to_183GHz` (scans, 492, channels) | `S1/`, `S2/Geolocation/` | `Scan_Daycnt`, `Scan_Mscnt` (x 0.1 ms) of each group |

* Time is always days since 2000-01-01 12:00 UTC plus milliseconds since 12:00 of that day; checked against the observation start / end times in the file headers of the samples.
* Brightness temperature = value x 0.01 + 327.68 K (from the Slope / Intercept of the data set). The V / H order of each group was confirmed with sea-surface medians of the samples (V is higher in every pair).
* Sounding channels (50\~183 GHz) only have single-channel brightness temperature (no PCT / RGB); 50\~54 GHz use the "other" colormap group, 165.5 / 183.31 GHz the 150 GHz colormaps.
* The previous `MWRI-L1` only knew FY-3D's `Scan_Daycnt`, so FY-3B (`Scan_Time_and_Period`) and FY-3C (`Scan_daycnt`) files raised `KeyError`; this is fixed.
* **Verification** (2026/11): samples of FY-3B, FY-3C, FY-3D L1c, FY-3F, FY-3G and FY-3H were all drawn; the eye of Typhoon Goni in the FY-3D L1c is at the same place as in MERSI-II at the same time.
  **FY-3A has no samples** (240 samples per scan, written from the format document). FY-3D GBAL HDF files are read as before (verified earlier).

### ATMS geolocation files

`SATMS_*` (SDR) and `TATMS_*` (TDR) from NOAA CLASS contain only brightness temperature; the latitude / longitude are in the matching `GATMO_*` (ATMS-SDR-GEO).
The reader pairs them by the `<satellite>_d<date>_t<start>_e<end>_b<orbit>` part of the file name (the creation time `c…` need not be the same); files whose geolocation file is not found are skipped with a message. How to get them:

* When ordering ATMS SDR / TDR from CLASS, also tick Geolocation (GATMO);
* or download them from NOAA's open data (AWS S3, no account needed), e.g. for NOAA-20:
  `https://noaa-nesdis-n20-pds.s3.amazonaws.com/ATMS-SDR-GEO/YYYY/MM/DD/GATMO_j01_d…_t…_e…_b…_c…_oeac_ops.h5`
  (`noaa-nesdis-snpp-pds` for S-NPP, `noaa-nesdis-n21-pds` for NOAA-21; `ATMS-SDR/` and `ATMS-TDR/` in the same bucket are SATMS and TATMS).

The CLASS files that pack geolocation and data into one file, `GATMO-SATMS_*` and `GATMO-TATMS_*`, are also supported (the datasets are the same), but there are no samples and they have not been tested.
SDR and TDR have the same output file name (they differ only in the format column of the data listing), so when both are plotted in the same run the latter overwrites the former.

### AMSU-A

* The sampling is very coarse (30 fields of view, 48 km at nadir, over 130 km at the edge), so only the rough structure of a typhoon can be made out; its main uses are the large-scale distribution at 23.8 / 31.4 / 89 GHz
  and the 50\~57 GHz temperature sounding channels (warm core). There is no dual polarization, only single-channel brightness temperature.
* Failed channels: a channel with fewer than half of its samples valid is treated as having no data at all (89 GHz of MetOp-B in the sample); a channel that is very noisy but still has values is shown as it is
  (54.4 GHz of MetOp-B in the sample).
* A file is one complete orbit, so its beginning and end overlap at high latitudes; in the overlap the nearest sample is taken, and there is a seam at the boundary.

### Instruments covered by PPS L1C

The 1C products on the NASA PPS archive (`arthurhouhttps.pps.eosdis.nasa.gov/gpmdata/YYYY/MM/DD/1C/`; user name and password are both the registered email)
are all in one and the same format, and the single reader `PPS-L1C` covers them all (the directories for 1987\~2026 were surveyed in 2026/11):

| Instrument | Satellite | Years | Channels | Tested |
|---|---|---|---|---|
| SSM/I | F08, F10, F11, F13, F14, F15 | 1987\~ | 19V/H, 22V, 37V/H, 85V/H | F14 1998/10/13 |
| TMI | TRMM | 1997\~2015 | 10V/H, 19V/H, 21V, 37V/H, 85V/H | 1998/10/13 |
| AMSR-E | Aqua | 2002\~2011 | 10\~89 GHz V/H | 2008/09/12 |
| SSMIS | F16, F17, F18, F19 | 2005\~ | 19\~183 GHz | F16 2017/09/06, F18 2026/10/20 |
| AMSR2 | GCOM-W1 | 2012\~ | 10\~89 GHz V/H | 2025/09/22 |
| GMI | GPM | 2014\~ | 10\~183 GHz (also 1C-R) | 2026/10/29 (1C and 1C-R) |
| MWI | WSF-M1 | 2024\~ | 10\~89 GHz V/H | 2026/11/17 |
| AMSR3 | GOSAT-GW | 2025\~ | 10\~183 GHz | 2026/10/20 |
| AMSU-B | NOAA-15, 16, 17 | 2000\~2010 | 89, 150, 183±1/3/7 | NOAA-15 2008/09/12 |
| MHS | NOAA-18, 19, MetOp-A, B, C | 2005\~ | 89, 157, 183±1/3, 190 | NOAA-18 2008/09/12, MetOp-B/C 2026/10/20 |
| ATMS | S-NPP, NOAA-20, NOAA-21 | 2011\~ | 23.8, 31.4, 88.2, 165.5, 183±1\~7 | S-NPP, NOAA-20 2026/10/20 |
| SAPHIR | Megha-Tropiques | 2011\~2021 | 183.31±0.2\~11, 6 in total | 2017/09/06 |

* The cross-track scanning AMSU-B, MHS, ATMS and SAPHIR have no dual polarization, only single-channel brightness temperature (no PCT or RGB).
* AMSU-B channels are written in the file as "89.0 +/- 0.9 GHz" and "150.0 +/- 0.9 GHz" (the double sidebands of the window channels); the menu uses only the center frequencies `89` and `150`;
  the ± is kept only for the channels near the absorption line (183.31±7…).
* **The A scan of AMSR-E 89 GHz failed in 2004/11**, and after that only the B scan has data: the merged `89V` / `89H` automatically use only the B scan
  (10 km along-track sampling, not 5 km).
* **GMI 1C-R**: the version in which the high frequencies (S2: 166, 183 GHz) are re-registered to the sampling positions of the low frequencies (S1); S1 is exactly the same as in 1C; the format column of the listing shows `L1C-R`.
  1C and 1C-R of the same orbit have the same output file name, so when they are plotted together the latter overwrites the former.
* Failed channels (150 GHz of F18 SSMIS, 183±3 of NOAA-15 AMSU-B in 2008…) still appear in the menu; selecting one shows
  "[skip] this channel has no valid data inside the map range".
* **Quality flag** (`[swath] mask_bad_quality`, default `false`): according to the L1C ATBD, when `Quality` is negative (serious error) the brightness temperature of the channel in error has already been set to
  missing, so it is enough to remove the missing values when reading. `Quality` is shared by a group of channels, and many negative values actually refer to only one of the channels
  (`-4` one channel missing data, SSMIS `-125` 150H failed, AMSU-B `-100`\~`-104` a certain channel unusable), while the other channels are still valid;
  setting it to `true` discards these samples entirely (the 183 GHz of F18 SSMIS and the 89 / 150 GHz of NOAA-15 AMSU-B disappear completely).
* Satellites that were not tested (F08\~F13, F15, F17, F19, NOAA-16/17/19/21, MetOp-A) have the same format as the tested instruments of the same type.

## Adding a format

Add `<instrument>-<format>.py` in `sensors/mw/` (the file name must contain `-`), write a class that inherits from `satima.mw.reader.SwathReader`,
with the last line `READER = <class name>`:

| Item | Contents |
|---|---|
| `NAME` `SENSOR` `DESCRIPTION` `PATTERNS` | Format name (= file name), instrument, description for `--list`, file search patterns (glob) |
| `scenes(folder)` | Returns a list of `Scene(satellite, instrument, start time, end time, file=, orbit=, format=)`, one per file |
| `channels(scene)` | Returns a list of `Channel(name, frequency GHz, polarization, geolocation group, key=for the reader's own use)` (in menu order); the name can be made with `band_name(frequency, polarization, offset)` |
| `geo(scene, group)` | Returns `[Geo(lat, lon, time of each scan, key for identification)]`; channels with merged A / B scans return two |
| `tb(scene, channel)` | Returns `[brightness temperature array]` (K, float32, NaN where invalid); the number and shapes correspond to `geo()` |

Everything else (listing, menus, interpolation, PCT / RGB, plotting, saving, multi-file merging) is handled by the main workflow.

## Verification (2026/11)

Samples: 9 files in `L1C-sample`, 26 files in `L1B-sample`.

**Reading**
* 7 of the 8 readers were checked by producing images from real files: AMSR2 2025-09-22 (Ragasa), SSMIS F16 2017-09-06 (Irma),
  SSMIS F17 2019-12-24 (Phanfone; the same orbit as the old `PHANFONE.png`, structure and time 1026Z agree),
  SSM/I F14 1998-10-13 (Zeb) and F15 2013-11-06 (Haiyan), the typhoon of 2026-10-20 with AMSR3 / SSMIS / ATMS / MHS,
  and global maps from FY-3D MWRI, GMI and WSF-M MWI.
* Binary formats: 15 SSMIS SDR files were parsed record by record to the end of the file, and the position matched the file size exactly; the size of the 6 SSM/I files = header + a whole number of records.
* L1B versus L1C for the same orbit:

  | Instrument | Latitude / longitude difference (maximum) | Brightness temperature L1B − L1C (mean) |
  |---|---|---|
  | AMSR2 | low frequencies (derived from 89A) ≤ 2 m; 89A / 89B 0 | +0.3 \~ +2.7 K (the XCAL calibration adjustment) |
  | AMSR3 | ≤ 4 m | −0.6 \~ +3.3 K |
  | GMI | 0 | 0.000 K (GMI is the XCAL reference) |

**Interpolation**
* Values taken at the sample positions: all three index-space methods reproduce the original values exactly.
* Fractional coordinates: round-trip error 0 at the samples and < 1e-8 at the cell centers; interpolation from the grid nodes compared with solving pixel by pixel differs by a median of 2e-6 and a maximum of 0.04 sampling intervals.
* The numba and numpy versions give pixel-identical results; the 8-bit images from strips + multiple cores and from the whole image on a single core are pixel-identical.
* The interpolated values of cubic do not go outside the range of the original sample values (89H: interpolated 110.32\~293.98 K, original 110.12\~294.19 K).
* gauss: the sample cases of WSF-M MWI (ratio 3.5), GMI, AMSR2 (A/B merged and 37 GHz), SSMIS and MHS show no jagged edges or stripes;
  the numba and numpy versions are pixel-identical, strips and the whole image are identical, the values do not go outside the range of the original sample values, and the pixels that have values are the same as with the other methods.
* Multi-file merging: one file split into two halves overlapping by 80 scans and merged again is exactly the same as the original file; two identical files merged equal one file.
  (The samples contain no two files that are really contiguous in time, so this was not tested with real adjacent files.)

**L1B readers added on 2026/11/19** (compared with the PPS L1C of the same orbit, taking the nearest point for each sample)
* `AMSRE-L1B` (2009/10/17, Hurricane Rick): the positions of 10.65\~36.5 GHz and 89B differ from L1C by 0.000 km (the derivation of the low-frequency latitude / longitude is correct);
  brightness temperature difference: mean −1.8\~+1.0 K, standard deviation 0.01\~1.0 K (XCAL calibration).
* `TMI-L1B` (1998/10/13): the positions of the 9 channels are the same; brightness temperature difference: mean −0.8\~+0.7 K, standard deviation 0.02\~0.5 K.
* `ATMS-SDR` (NOAA-20 2025/09/28, Hurricane Humberto): the positions of the 7 channels that SDR has in common with L1C are the same, and the brightness temperature difference is a constant (−1.5\~+1.1 K, standard deviation 0.00 K);
  TDR differs from L1C by −1.7\~+0.5 K, standard deviation 0.1\~0.8 K (antenna pattern correction). L1C does not have the 50\~57 GHz sounding channels, so they cannot be compared.
* `AMSUA-L1B` (MetOp-B, MetOp-C 2026/04/12): there is no product to compare against. The geolocation was confirmed with the coastline (the Iberian Peninsula, France, the UK and north-west Africa all line up);
  the brightness temperature was only checked for plausibility (sea surface at 23.8 GHz about 180\~220 K, land about 280 K, sounding channels 213\~254 K); **there is no numerical comparison with other products**.
  There is no sample from MetOp-A (M02).
* `SSMIS-TDR` (F16 2009/10/18 02:27Z, Hurricane Rick; against the `NPR.SDRN` of the same orbit, matching scan by scan on the scan time):
  the scan times are identical to SDR (the port of the old program took the scan counter for the time, which has been fixed; the zero padding at the end of the file up to a multiple of 512 is also handled);
  the correlation between TA and TB is above 0.996 for all 19 channels; TA − TB has a mean of −0.9\~−6.9 K (largest at 19V), a standard deviation of 0.4\~2.2 K at 91 / 150 / 183 GHz
  and 0.3\~0.7 K for the sounding channels. The positions of 37 / 91 / 150 / 183 GHz and of the sounding channels differ from SDR by 1.1 km (0.01 degree, i.e. rounding);
  the 19 / 22 GHz positions differ from SDR by about 16 km: TDR has its own 19 / 22 GHz latitude / longitude, and interpolating with them to the SDR positions gives a smaller standard deviation
  (19H 4.5 → 3.1 K, 19V 2.8 → 1.9 K, 22V 2.0 → 1.4 K), which shows that **the NOAA SDR resamples 19 / 22 GHz to the 37 GHz positions**
  (which is why the 19 / 22 GHz standard deviation of SDR against L1C below cannot be brought down). There is only this one sample (F16); F17\~F19 were not tested.

**PPS L1C versus NOAA SDR** (same orbit; the geolocation and calibration processing of the two products differ, so only the order of magnitude can be confirmed)
* SSM/I F14 1998/10/13: the standard deviation of Tc − Tb for each channel is 0.5\~0.9 K (mean −1.8\~+0.7 K, i.e. the XCAL calibration adjustment);
  the geolocation of the two differs by a median of 3.7 km (L1C uses re-geolocated data). The interleaving order of the A/B scans of 85 GHz agrees on both sides.
* SSMIS F16 2017/09/06: standard deviation 0.4\~0.7 K at 37 GHz and 1.4\~2.3 K at 91 / 150 / 183 GHz; the geolocation differs by 4\~11 km (in L1C each channel group has its own set of
  latitude / longitude, while the environmental channels of the NOAA SDR share one set). The standard deviation at 19 / 22 GHz is 2\~4 K and cannot be brought down however the data are aligned;
  presumably the NOAA SDR applies an additional resampling to these channels; this was not confirmed further.

**Parts that could not be verified**
* Colormap `rammb89`: the values were derived backwards from the color bar of a finished RAMMB image; they are not an officially published definition.
* TDR antenna temperature is not converted to brightness temperature.

## Speed (for reference)

2200 × 2333 pixels, AMSR2 89 GHz (A/B merged, cubic): sampling and coloring take about 2.4 seconds (multiple cores, numba);
on the first run numba needs about 5 seconds to compile (cached afterwards). A 2000 × 778 global map takes about 1 second (not counting reading the coastline).
gauss and cubic take about the same time (WSF-M MWI 89H at 2048 × 2057: sampling, coloring and saving take 1.3 seconds in total);
without numba the gauss sampling takes about 19 seconds on a single core (3.6 seconds on a single core with numba).

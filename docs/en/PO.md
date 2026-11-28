# Satima: Polar-Orbiting Imagers (PO)

> This is the detailed documentation for **PO (polar-orbiting imagers)**, rewritten from the original Satima-PO. For how to run the integrated version, the data listing, the layered config files
> (`Satima.toml` + `Satima-PO.toml`), the interface language, the image modes (native / dpi) and the continuity rules for polar orbiters, see [README.en.md](../../README.en.md).
> PO is the family that changed most of the three: the differences from the original Satima-PO are listed under "Differences from the original three programs" in the README.

Single channels (reflectance, brightness temperature with colormaps, DNB) and true color from VIIRS, MODIS, AVHRR, METimage, OLCI, SLSTR, SGLI and FY-1C/1D MVISR.
The 8 scripts, originally independent of each other, were first integrated into Satima-PO and then merged into Satima; the calibration and algorithms of each format remain in its own script.
FY-1C/1D MVISR (`MVISR-L1B`) was rewritten from `SatimaV3-FY1D-NSMC.py` and NOAA AVHRR in NSMC format (`AVHRR-NSMC`) from `SatimaV3-NOAA-NSMC_accnumpy.py`,
both added in 2026/11 (see "FY-1C/1D MVISR" and "NOAA AVHRR/3 (NSMC format)").

## Usage

```
python Satima.py -t po               list polar-orbiting imager data only
python Satima.py -f VIIRS-HDF5       search one format only
```

```
1   S-NPP     VIIRS    HDF5   2023/05/26 04:03-04:10Z   NOAA       22 bands   5 files
2   METOP-B   AVHRR3   NAT    2026/04/12 09:31-11:13Z   EUMETSAT   IR/VIS
Input the index (several: (3,5) or [1,3,5]): 1
Projection: 1.MERC  2.CYL  3.LCC
Input the projection (Enter = MERC):
Map center (lat lon): 18 121
Range (degrees of latitude / longitude from the center; one number = both): 5 8
(channel / product menu of each format)
Select colormap (several allowed, 0 = all): 0
Image width (…; Enter = native sampling interval …):
Image title: Haikui
```

* One row of the listing is **one continuous segment**: short granules such as VIIRS (about 85 seconds per file), METimage (about 1 minute) and MODIS (5 minutes), when of the same satellite and
  contiguous in time, are joined into one entry automatically, and the notes column shows the number of files and of missing files. Whole-orbit files such as AVHRR each get their own row. For the rules see "Continuity rules for polar orbiters" in the README.
* VIIRS is no longer one row per band: after the data is selected, the menu lists the channels that every granule of this segment has, plus true color when M3 / M4 / M5 (plus I1) are all present.
* Map range: `Satima-PO.toml` defaults to `[map] range_input = "center"` (the center is asked first, then the range, as in the original Satima-PO);
  4 numbers are taken as south north west east. The center longitude may be given in -180\~180 or 0\~360, and ranges across the date line can be plotted directly.
* **Center time**: the time on the image and in the file name is the actual observation time of the pixel closest to the center point. Center point = the center you input; when the four edges are input it is the midpoint of the range.
* For infrared channels the highest brightness temperature near the center point is printed (`[eyetemp]`), so when plotting a typhoon put the center on the eye.
* Colormaps, true color styles and DNB scalings can all be multi-selected; `0` = all.
* The band name in the image annotation and file name follows GEO: `BAND-04`, `BAND-03A` (`B04`, `B03A` before 2026/11); the native channel names of VIIRS (`I05`), SGLI (`VN11`) and METimage are unchanged.
* The title (formerly Picname) is asked by default (`[title] ask = true` in `Satima-PO.toml`) and may be left empty; it is put into `{title}` of the file name. PO defaults to `[title] draw = false`, so the title is only used in the file name and not drawn on the image; with `true` it is also drawn at the lower left corner (the common default, as in GEO and MW).

## Config file highlights

| Item | Setting | File |
|---|---|---|
| Data / output directories | `[general] data_dir`, `output_dir`; relative paths are relative to the folder containing `Satima.py`; the output directory may contain fields, e.g. `plot/{sensor}/{time:%Y%m%d}` | common |
| Coastline | `[coastline] draw / color / linewidth`; the resolution is in `[map] resolution` | common (PO overrides it to green 0.5) |
| Annotation | `[annotation] draw = false` turns it off; or individually `[[output.image]] annotation = false`; text template `[annotation] text` | common |
| File name | `[output] filename`; empty fields are omitted automatically; PO also has `{res}` (MODIS, SGLI) | common |
| Image mode | `[image] mode = "native"` (native resolution) or `"dpi"` | common |
| dpi | One image per `[[output.image]]`: `dpi`, `suffix` (used only with `mode = "dpi"`) | common |
| Background color | `[image] background` (outside the swath); OLCI and SGLI set it to white in `[readers.<format>.image]` | common / PO |
| Range warning | `[general] range_warning`: when the half-width exceeds this value, a reminder about memory use is given | PO |
| Merging of long files | `[swath] merge`, `merge_gap` | PO |
| Eye temperature | `[eyetemp] enable / radius_deg / legacy_counts` | PO |
| Colormap menu | `[colormap] ir / wv`; items can be added, removed or reordered | PO |
| Per-format settings | `[readers.<format>.<section>]`, e.g. `[readers.MODIS-HDF4.general] data_dir = "D:/modis"` | any |

Examples of default file names: `Haikui.JPSS-1.VIIRS.I5.BD.202608190028Z.png`; the second image in dpi mode `S-NPP.VIIRS.I1.202305260403Z.500DPI.png`

## Sampling and image size

Every sample of the swath carries its own latitude / longitude. The two image modes draw in the same way (`satima/po/canvas.py`); only the way the image size is determined differs:

| | `mode = "native"` | `mode = "dpi"` |
|---|---|---|
| Image size | One pixel at the map center = the native sampling interval (`[image] size` can enlarge it by an integer factor or set the width) | The dpi of `[[output.image]]` x the matplotlib default size in inches (a square range at 1000 dpi is about 3700 pixels high) |
| Number of images | One (`[image] suffix` fills `{suffix}` in the file name) | One per `[[output.image]]` |
| Values | Nearest neighbour: each pixel takes the nearest valid sample; values are not interpolated | Same as left (pcolormesh is no longer used) |

- The **native sampling interval** is measured from the latitude / longitude of the data being plotted (the finer of the along-scan and along-track intervals inside the range), so no per-format setting is needed; the 250 m / 375 m of sharpened true color
  is also handled automatically. When the map range does not include nadir, what is obtained is the finest interval inside the range (pixels at the swath edge are larger anyway).
- `[image] size`: `"native"`, `"x2"`\~`"x8"` (integer multiples) or a fixed width (pixels). When the width is below `min_size` (default 2000) it is automatically enlarged by an
  integer factor of at least `min_scale` (default 4); when the long side exceeds `max_size` the largest integer factor that fits is used instead; `min_size = 1` turns the automatic enlargement off.
- The output grid is the regular grid of the shared map (MERC / CYL / LCC, WGS84). The nearest sample is found with a KD-tree on the map plane: Mercator and LCC are
  conformal projections, so nearest on the plane is nearest in reality; for equidistant cylindrical the east-west direction is first multiplied by the cosine of the central latitude to make it approximately isotropic (with a very large latitude range there is a slight deviation at the edges).
- Pixels farther from the nearest sample than about half a sample cell at that location are treated as outside the swath and show `[image] background`. Samples with invalid data take no part in the sampling, so
  bow-tie overlaps and VIIRS bow-tie deletion areas are filled in by valid pixels of adjacent scans; the gaps between scans near nadir
  (Aqua MODIS after its orbit was lowered) are filled in too. The gap from missing files exceeds the distance limit and is left blank, not filled.
- When the output pixels are coarser than the sampling (a large range in dpi mode, or a smaller width was specified), it amounts to subsampling.
- When "0 = all colormaps" is selected the sampling is computed only once; when dpi mode outputs several images each one is sampled once, while the reader reads and calibrates the files only once.
- The LCC frame is the smallest rectangle covering the whole latitude / longitude range, so its four corners extend beyond the input range; when cropping the data the reader uses the latitude / longitude range actually covered by the frame.

## Formats at a glance

| Script | Satellite / sensor | Files | Products |
|---|---|---|---|
| VIIRS-HDF5 | S-NPP, NOAA-20/21 VIIRS | `SV*_ops.h5` + GMTCO / GITCO / GDNBO (or GMODO / GIMGO); CLASS packaged files can also be read | Single channels, DNB (4 scalings), TrueColor-M / TrueColor-I |
| MODIS-HDF4 | Terra, Aqua MODIS | `MOD02*/MYD02*` + MOD03/MYD03 | Thermal infrared, visible, true color (RGB / CORR / WV) |
| AVHRR-NAT | Metop-A/B/C AVHRR/3 | `AVHR_*_1B_*.zip` (EPS native) | CH1–CH5 |
| AVHRR-NC | Metop-A/B/C AVHRR/3 | `W_XX-EUMETSAT-...METOP*.nc` | CH1–CH5 |
| METimage-NC | Metop-SG METimage (VII) | `*VII-1B*.zip / .nc` | 20 channels, true color (General / CIRA / Worldview) |
| OLCI-SAFE | Sentinel-3 OLCI | `S3*_OL_1_*.SEN3` (zip / tar / folder) | Oa01–Oa21, true color (RGB / CORR / WV) |
| SLSTR-SAFE | Sentinel-3 SLSTR | `S3*_SL_1_RBT*.SEN3` (zip / tar / folder) | S1–S6 reflectance, S7–S9 brightness temperature |
| SGLI-HDF5 | GCOM-C SGLI | `GC1SG1_*_1BSG_*.h5` | VN01–VN11, SW01–SW04 reflectance, TI01–TI02 brightness temperature, true color (RGB / CORR / WV) |
| MVISR-L1B | FY-1C/1D MVISR | `FY1?_AVHRR_HRPT_L1_*.L1B` (NSMC HRPT 1B) | CH1–CH10: CH3–CH5 brightness temperature, the others reflectance; true color (RGB / CORR / WV) |
| AVHRR-NSMC | NOAA-15~19 AVHRR/3 | `NA??_AVHRR_HRPT_L1_*.L1B` (NSMC HRPT 1B) | CH1–CH5 (CH3A / CH3B switched by scan line) |
| AVHRR-L1B | TIROS-N~NOAA-19, MetOp AVHRR/1, /2, /3 | NOAA Level 1b (POD / KLM): CLASS `NSS.GHRR/LHRR/HRPT/FRAC.*`, Chiba University HRPT (including `.l1b.bz2`), also `.gz` / `.bz2` / `.xz`; **detected by content, not by file name** | CH1–CH5 (CH3A / CH3B switched by scan line on AVHRR/3; no CH5 on AVHRR/1) |
| MERSI-HDF5 | FY-3A/B/C MERSI-1, FY-3D MERSI-II, FY-3E MERSI-LL, FY-3F/H MERSI-III, FY-3G MERSI-RM | `FY3?_MERSI_*_L1_*_{1000M,0250M,0500M}_*.HDF` + GEO1K / GEOQK / GEOHK (NSMC) | single channels (250 m / 500 m / 1 km), low-light LL (3 scalings), true color (General / CIRA / Worldview) |
| VIRR-HDF5 | FY-3A/B/C VIRR | `FY3?_VIRRX_GBAL_L1_*_1000M_MS.HDF` (+ GEOXX for FY-3C) | CH1–CH10: CH3–CH5 brightness temperature, the rest reflectance; true color (General / CIRA / Worldview) |
| AVHRR-L0 | TIROS-N~NOAA-19 AVHRR/1, /2, /3 | HRPT raw frames (Level 0) recorded by a receiving station: `.raw`, `.raw.xz` etc., Chiba University `.l0.bz2` / `.gi.bz2`, also inside a `.tar`; **detected by content**; navigation needs TLEs | as AVHRR-L1B |

VN07 / VN10 of SGLI are high-gain (ocean color) channels that saturate over cloud at a reflectance of about 0.16 / 0.13, and VN08 saturates at about 0.49; to plot clouds use VN11 (868.5 nm).

Start and end times of each format: for VIIRS, METimage, AVHRR-NAT, OLCI and SLSTR they are taken from the file name; MODIS, MERSI and VIRR file names only have the start time, so the nominal 5 minutes is assumed;
AVHRR-NC and SGLI file names have no end time, so the listing shows only the start time (files within 10 minutes of each other are still joined into one segment). AVHRR-NAT and AVHRR-NC read only one file at a time.
For MVISR, AVHRR-NSMC and AVHRR-L1B the start and end times are taken from the header of the file (to the millisecond). For AVHRR-L0 they come from the time codes and the number of frames.

### VIIRS files

* File name `<product[-product…]>_<satellite>_d<date>_t<start>_e<end>_b<orbit>_c<creation time>_<origin>_ops.h5`. An ordinary file has only one product (`SVM15`, `GMTCO`…);
  NOAA CLASS packaged files list the contained products at the start of the file name (`GMTCO-SVM03-SVM04-SVM05_…`, `GDNBO-SVDNB_…`). The two kinds can be mixed:
  the program records which channels and geolocation each file contains, and each channel and each geolocation points to the file that contains it. The group names inside the HDF5 are the same in both kinds of file.
* Geolocation: M channels `GMTCO` → `GMODO`, I channels `GITCO` → `GIMGO`, DNB `GDNBO`. Left of the arrow is the terrain-corrected version; the right one is used only when the left one is not found,
  and this is noted on screen. The Rayleigh correction of true color uses sea level when the geolocation file has no height data.
* DNB: when `GDNBO` has `Latitude_TC` / `Longitude_TC` (terrain-corrected), they are preferred.
* When there are several files for the same product of the same granule, the one with the latest creation time (`c…`) is used.
* Once the range is known, only the granules that intersect the map are read (judged from one point every 64 pixels, with 2 degrees added outside the range), plus one extra before and after.
* Channels that only some of the granules in a segment have are also listed in the menu (with a notice); after one is selected, granules without that channel are skipped. When a channel needed for true color is missing in any granule, the program stops and points out which one.
* Channels that have data but no geolocation file are not offered in the menu; the channel column of the data listing says which kind is missing (e.g. `0 bands, no GITCO`).
* When only one product can be selected, the program does not ask.
* The DNB ERF + saturation correction (scaling 3) counts the saturated fraction over the loaded area, so the result depends on which granules were loaded (now the granules that intersect the map plus one before and after).
* **Verification status** (2026/11, real files of S-NPP, NOAA-20 and NOAA-21): for the I5 colormap, I1, M5, DNB (ERF), TrueColor-M and TrueColor-I, the center time and eye temperature
  are the same as the original Satima-PO and the images are visually consistent; granule continuity, reading only the intersecting granules, selecting several segments at once, the missing-geolocation notice and `Latitude_TC` were all run with real files.
  Packaged files were tested with files assembled from real data (`tests/pack_viirs.py`); the output is pixel-identical to that from the separate files.
  **No real files to test yet**: original packaged and aggregated files downloaded from CLASS, `GMODO` / `GIMGO`, M-band brightness temperatures (M12\~M16), removal of duplicate files and marking of missing files
  (the latter items were only run with the synthetic files of `tests/fake_viirs.py` to confirm that the workflow runs through).

### Fengyun-3 MERSI (MERSI-HDF5)

* One reader for five MERSI types: FY-3A/B/C MERSI-1, FY-3D MERSI-II, FY-3E MERSI-LL, FY-3F/H MERSI-III and FY-3G MERSI-RM; the instrument column
  of the listing shows which one. Format and calibration follow the "L1 data instruction" documents of each satellite on the NSMC document page
  (<https://satellite.nsmc.org.cn/DataPortal/cn/support/document.html>).
* Files: `FY3X_MERSI_{GBAL|GRAN}_L1_YYYYMMDD_HHmm_<kind>_{MS|Vn}.HDF`, 5 minutes per granule. Data `1000M` / `0250M` / `0500M` (FY-3G),
  geolocation `GEO1K` / `GEOQK` / `GEOHK` (FY-3G) / `GEOXX` (FY-3A/B 250 m), in the same folder. Contiguous granules are joined automatically;
  when a granule has several versions of the same kind (`V0`, `V1`) the higher version is used.
* The menu is grouped by resolution (250m / 500m / 1km) into reflectance, brightness temperature and low light; a channel available at two
  resolutions is listed at both. True color (channels 3, 2, 1 of MERSI-1/II/III) is offered at 250 m and 1 km. When segments of different MERSI
  types are selected at once, later segments are matched by kind and wavelength (e.g. MERSI-II CH24 to MERSI-LL CH06).
* Latitude / longitude: geolocation file of the same resolution → (250 m) GEO1K interpolated scan by scan → tie points in the data file (every
  5 / 20 pixels); the last two are noted on screen. Compared on FY-3H (which has both GEOQK and GEO1K): GEO1K interpolation differs from GEOQK by
  33 m (median); tie points differ by only 3 m over the sea but are **not terrain corrected**: 0.4 km (median) and 2.3 km (99%) above 1500 m.
* Angles for true color and low-light scaling come from GEO1K / GEOHK (interpolated from 1 km for 250 m); without them (e.g. only 0250M +
  GEOQK) solar angles are computed from time and lat/lon and satellite angles from the sub-satellite track (where the along-scan spacing of the
  tie points is smallest) and the satellite height. For FY-3D 1 km true color, computed instead of GEO1K angles change the CREFL result by at
  most 1/255.
* Only the needed part is read: once the range is known, the tie points tell which scans and columns are inside the map (plus 0.3 degrees), and
  only that block is read (a whole 250 m granule is 8000 x 8192).
* Reflectance: the NSMC calibrated reflectance Ref = (k0 + k1 dn + k2 dn²) / 100 (MERSI-RM: k0 + k1 dn), i.e. π L / E0 without dividing by
  cos(solar zenith); single channels are drawn with this value (gray 0\~1). True color multiplies by the squared Earth-Sun distance ratio,
  applies the sun-zenith correction and then CREFL (there are no official FY-3 coefficients: O3 / H2O are borrowed from the nearest MODIS bands,
  the Rayleigh optical depth is computed from the equivalent central wavelength).
* Brightness temperature: radiance → inverse Planck at the equivalent central wavenumber Te → Tbb = A Te + B (the NSMC formula). The equivalent
  central wavelength comes from the file (`Effect_Center_WaveLength`), A and B from `TBB_Trans_Coefficient(_A / _B)`; the document tables are used
  when they are missing or implausible. satpy uses the nominal wavelength and (Te − B) / A and differs from this by 0.1\~1.4 K (FY-3D CH24:
  1.4 K; equivalent central wavenumber 933.364 cm⁻¹ vs 925.9 cm⁻¹ for the nominal 10.8 µm). MERSI-1 channel 5 (11.25 µm) follows satpy: the Slope
  is mistakenly 1 (should be 0.01), A, B = 1.0047, −0.8549.
* Low light (MERSI-LL LL, 0.5\~0.9 µm, similar to the VIIRS DNB): 1. ERF Dynamic 2. ERF + saturation correction (as the VIIRS DNB; radiance in
  W m⁻² sr⁻¹ converted to W cm⁻² sr⁻¹) 3. Linear (π L / E0, for day and twilight). The lunar phase is computed with Meeus' low-precision formula
  (`satima/po/dnb.py`, no ephem needed), the lunar zenith angle comes from GEO1K (90 degrees when missing). FY-3E is in a dawn-dusk orbit, so
  night granules often cross the terminator; ERF adapts pixel by pixel to the solar zenith angle.
* MERSI-II 250 m infrared shows horizontal stripes over cold clouds; they are in the data (no destriping is done here).
* **Verification** (2026/11): samples of FY-3D (Typhoon Goni 2020), FY-3E (Milton, Melissa, Kong-rey), FY-3F, FY-3G (Melissa) and FY-3H were all
  drawn (250 m / 500 m / 1 km single channels, true color in three styles, low light in three scalings); coastlines match the clouds, and
  reflectances are pixel-identical to satpy. **MERSI-1 (FY-3A/B/C) has no samples**: written from the format documents and satpy, and only run
  with synthetic files assembled from FY-3D data (`tests/fake_fy3.py`).

### Fengyun-3 VIRR (VIRR-HDF5)

* The Visible and Infrared Radiometer of FY-3A/B/C: 10 channels at 1.1 km, 5 minutes per granule (1800 scans x 2048 pixels). **No samples**:
  written from the NSMC FY-3A/B L1 Product Format and FY-3C VIRR L1 format documents and satpy `virr_l1b`, and only run with synthetic files
  (`tests/fake_fy3.py`).
* Files: `FY3X_VIRRX_GBAL_L1_YYYYMMDD_HHmm_1000M_MS.HDF`; FY-3C also has a `…_GEOXX_MS.HDF` geolocation file, FY-3A/B keep lat/lon, angles and
  time in the data file.
* Reflectance % = DN x k + b (`RefSB_Cal_Coefficients`; the FY-3B values in the files are wrong, satpy's coefficients are used). Infrared: the
  per-scan `Emissive_Radiance_Scales` / `Offsets` give radiance, inverted with `Emissive_Centroid_Wave_Number` (spelled `Emmisive_…` on FY-3A/B);
  `Emissive_BT_Coefficients` are applied as Tbb = A Te + B when the A / B layout can be recognized (otherwise no correction).
* Products: single channels (CH3–CH5 brightness temperature, the rest reflectance), true color (CH1, CH9, CH7; styles as MERSI).

### FY-1C/1D MVISR

* Files: NSMC HRPT 1B, fixed-length records of 28400 bytes (big-endian); the format is in NSMC's "FY-1 (02) batch polar-orbiting meteorological satellite data format" §3.1.
  The file name says AVHRR, but the instrument is the 10-channel MVISR. One file is one received pass (about 9\~20 minutes); when several files contiguous in time are selected they are merged according to `[swath] merge`.
* Once the range is known, the scan lines that intersect the map are first picked out with the tie points (3 degrees added outside the range), and only those lines are decoded and interpolated.
* **Calibration**: every scan line carries calibration coefficients for the 10 channels. The format document says "(2,10), (1,10) holds the slopes and (2,10) the intercepts"; these are Fortran
  subscripts, and the actual storage order is slope1, intercept1, slope2, intercept2… The original script took the first 10 values as slopes and the last 10 as intercepts, so every pair was wrong; this was the main reason
  there was "no good calibration method". Correctly paired, the meaning is the same as NOAA pre-KLM 1b:
  - Visible / near infrared (CH1, CH2, CH6\~CH10): albedo (%) = slope / 2^30 × count + intercept / 2^22; reflectance = albedo / 100, not divided by cos(solar zenith angle) (same as AVHRR-NAT).
    The zero point falls on the space count (12\~29).
  - Infrared (CH3\~CH5): radiance (mW/(m² sr cm⁻¹)) = slope / 2^30 × count + intercept / 2^22; the slope is negative and is updated for every scan line from the onboard blackbody.
    Brightness temperature T\* = c2·ν / ln(1 + c1·ν³ / L), T = (T\* − A) / B; the central wavenumber ν and the band correction A, B are taken from FILTER_FUNCTIONS of the RTTOV coefficient files
    (filter functions provided by NSMC; `rtcoef_fy1_4_mvisr.dat`, `rtcoef_fy1_3_mvisr.dat`):

    | Channel | FY-1D ν (cm⁻¹) | A | B | FY-1C ν (cm⁻¹) | A | B |
    |---|---|---|---|---|---|---|
    | CH3 | 2699.586 | 3.3894 | 0.995275 | 2691.460 | 3.6385 | 0.994675 |
    | CH4 | 925.553 | 0.7552 | 0.997388 | 912.074 | 0.8458 | 0.996996 |
    | CH5 | 850.420 | 0.7725 | 0.997133 | 829.868 | 0.7262 | 0.997266 |

  - Check (5 FY-1D files, 2003\~2006): brightness temperatures about 190\~302 K; the split-window difference CH4 − CH5 has a median of 1.5\~2.3 K over warm ocean (CH4 > 290 K)
    and about 1 K over cold cloud tops (CH4 < 215 K), both in the normal range. With 1/wavelength as the central wavenumber (12.0 µm = 833 cm⁻¹) the cold-cloud difference is about 3.5 K, which is not plausible.
  - Limitations: no nonlinearity correction (neither the format document nor the files have coefficients for it). The visible coefficients are identical in all files from 2003\~2006 (presumably pre-launch),
    so sensor degradation is not reflected and absolute reflectances are only indicative (images are not affected).
* **Geolocation**: the 51 tie points are at pixels 25, 65, …, 2025 (1-based, as in NOAA LAC 1b), verified with scan angles derived from the satellite zenith angles
  (first point about −54.0°, spacing 2.16° = 40 pixels × 0.0541°). The original script spread the 51 points evenly over the whole scan line, putting the swath edges about 24 pixels off.
  Along the scan, interpolation is linear in Earth-centered Cartesian coordinates with the Earth central angle as the coordinate, extrapolating 24 pixels at both ends; leave-one-out estimates are about 0.4 km inside and about 1\~4 km at the edges
  (the tie points themselves are quantized to 1/128°, about 0.9 km). In a night infrared image with the Kanto region of Japan at the eastern swath edge, the coastline matches the land–sea temperature boundary.
* Latitude / longitude are signed integers; the original script read them as unsigned, so scan lines in the southern and western hemispheres were skipped.
* A few pixels at the edges of some channels are zeros or constant values filled in during co-registration (pixels 1\~2 of CH5/9/10, 1\~4 of CH7, the last 2 of CH3/CH8); they are set invalid.
* The last 4-byte word of each scan line holds CH9, CH10, 0 in that order (actual data; the format document says 0, CH9, CH10). The scan line count in the data header does not match the actual number of records, so it is computed from the file size.
* **True color** (0 in the channel menu): R = CH1 (0.58–0.68 µm), G = CH9 (0.53–0.58), B = CH7 (0.43–0.48); modes as in SGLI: 1.RGB, 2.CORR, 3.WV.
  - Sensor degradation: the visible coefficients are identical in all files from 2003\~2006, and have been updated in the 2011 files (CH1 slope 0.0893 → 0.1292). The TOA reflectance of thick ice cloud should be nearly neutral in the visible
    (estimated with the CREFL forward model: with 200\~400 hPa of Rayleigh above the cloud top, CH7/CH1 stays within 1.00 ± 1%), but the measured CH7/CH1 decreases year by year:

    | Date | 2003/09 | 2004/10 | 2005/07 | 2005/09\~10 (3 passes) | 2006/04 | 2006/10 | 2011/05 (new coefficients) |
    |---|---|---|---|---|---|---|---|
    | CH9/CH1 | 1.00 | 0.99 | 0.97 | 0.95 | 0.92 | 0.94 | 0.93 |
    | CH7/CH1 | 0.95 | 0.94 | 0.84 | 0.82\~0.83 | 0.77 | 0.72 | 0.85 |

    A direct composite turns yellow, more so in later years. True color therefore first applies **inter-channel gains**: the median TOA ratio of thick ice cloud over the whole pass (CH4 < −30°C, CH1/cos(θs) > 0.5, θs < 75°,
    satellite zenith < 50°) brings CH9 and CH7 in line with CH1; three passes from the same period agree within 2%. With fewer than 100 thick ice cloud samples (near the terminator),
    the measurements with the same calibration set are interpolated by date. The gains and their basis are printed at run time. Single-channel output is not affected; the absolute degradation of CH1 itself cannot be determined this way.
  - Rayleigh / ozone: there are no spectral response functions for the visible channels, so the band limits of the format document (boxcar) are averaged with solar-spectrum weighting:
    τR = 0.0577 / 0.0945 / 0.2133, aO3 = 0.0912 / 0.0911 / 0.0054 (CH1 / CH9 / CH7). Recomputing the CREFL coefficients of MODIS B1/B3/B4 and VIIRS M3/M4/M5/I1 the same way gives τR within 4% and ozone mostly within 6%.
    Ozone 0.319 cm-atm, sea level, no water vapor correction.
  - Geometry: solar zenith, satellite zenith and relative azimuth of the 51 tie points in the file, interpolated with the Earth central angle as for geolocation. Compared with astronomical calculation by pyorbital, the solar zenith differs by < 0.02°;
    relative azimuth = solar azimuth − satellite azimuth (0° = backscatter, the CREFL convention).
  - Known limitations: thick clouds are slightly cream after correction (B/R about 0.9; CREFL corrects high clouds for the full atmosphere, MODIS does the same); sun glint turns khaki at large solar zenith angles
    (CREFL assumes a Lambertian surface); clear ocean in 2003\~2006 is short of blue and looks olive (very sensitive to the CH7 calibration and τR;
    with τR = 0.18 the ocean looks normal in these three scenes, but with the new 2011 coefficients the current τR already gives a blue ocean, so this is not handled by adjusting τR).
* **Verification status** (2026/11): 13 FY-1D files (2003\~2011), all 10 channels readable; CH4 colormap, CH2 and true color (2003/09 Typhoon Dujuan, 2005/09 Typhoon Longwang, 2011/05 Yellow Sea) images visually normal.
  No FY-1C samples (coefficients included, not tested).

### NOAA AVHRR/3 (NSMC format)

* Files: NSMC's NOAA-15\~19 HRPT 1B (`NA16_AVHRR_HRPT_L1_ORB_MLT_NUL_…L1B`), 22016 bytes per record; the format is in NSMC's
  "NOAA series polar-orbiting meteorological satellite data format" §3.2. It is close to the 16-bit NOAA KLM 1b format but differs in a few places: the scaling of the visible coefficients (slope ×10¹⁰, intercept ×10⁷),
  the order of the infrared coefficients (quadratic, linear, constant) and the satellite codes (NOAA-16/17/18 = 3/11/13); the radiance conversion constants in the header are wrong (a NOAA-18 file holds NOAA-16 / 17 values) and are not used.
* Calibration: visible channels are dual-gain albedo (the intersection field is 0 before 2008 and is computed from the two lines; slope 2 of NOAA-17 CH3A overflows int32 and is derived from continuity),
  reflectance = albedo / 100. The linear infrared radiance already includes the space radiance Ns (checked with a space count of about 990), and CH4 / CH5 get the NOAA KLM nonlinearity correction;
  brightness temperature T = (T\* − A) / B with ν, A, B from the NOAA KLM User's Guide (identical to the RTTOV NOAA-18 coefficients).
* Channel 3A / 3B: the channel 3 flag in the scan line bit field has opposite meanings before 2008 and from 2010 on, and is inconsistent within a file too, so the sign of the linear term of the 3B calibration coefficients is used instead
  (negative = 3B; in 3A mode that field holds the 3A albedo coefficients), verified with the count distributions.
* Quality of the calibration coefficients: the infrared coefficients of single scan lines jump (about 1\~4 K); lines that differ from the 51-line running median by more than a threshold use the median.
  The infrared coefficients of NOAA-17 2010/10/28 are invalid for the whole file (period of AVHRR failure).
* Damaged data, not marked by the quality flags: misaligned lines (more than 10% of pixels with implausible CH4 − CH5, or discontinuous with the adjacent lines) are masked entirely,
  and runs of at most 2 such lines are interpolated from the lines above and below; single bit errors only mask the CH4 / CH5 pixels with an implausible split-window difference (so they are not taken as the maximum eye temperature);
  single-line streaks in the visible channels are replaced by the mean of the lines above and below. NOAA-16 2005/09/28 (whole pass) and 2005/07/16 north of 32°N are misaligned garbage and are masked (no data inside the range).
* Geolocation: the 51 tie points are at pixels 25\~2025 (NSMC format document), interpolated with the Earth central angle (as in MVISR-L1B). In CH2 images, the coastlines of Taiwan (scan angle about 34°) and Hainan (near the swath edge)
  match.
* **Verification status** (2026/11): 23 NOAA-16 / 17 / 18 files (2004\~2016) readable; CH4 − CH5 over warm ocean about 0.8\~3.5 K. Meranti (2016/09/13) CH4,
  Utor (2013/08/11) CH1 and NOAA-17 CH3A images are visually normal. No NOAA-15 or NOAA-19 samples.

### NOAA AVHRR Level 1b (CLASS, Chiba University)

* Files: the two NOAA Level 1b formats, documented in the POD Guide (TIROS-N\~NOAA-14) and the NOAA KLM User's Guide (NOAA-15 onwards) from NCEI;
  the PDFs are downloaded in `D:\github\HRPT\doc`.
  - POD: GAC records are 3220 bytes (409 pixels, data from the 3rd record), LAC/HRPT 14800 bytes (2048 pixels).
  - KLM: GAC 4608 and LAC/HRPT/FRAC 15872 bytes (10-bit).
  - Every file in the folder is identified by its content (file size and record length, header fields, a 0 / 122-byte TBM / 512-byte prefix),
    so renamed files are read too.
  - The satellite comes from the dataset name inside the file (`NSS.HRPT.NN.…`), otherwise from the spacecraft ID in the header
    (POD and KLM use different numbers).
  - Chiba University files (processing block `B0000000`) are listed with the source Chiba. Chiba stores NOAA-18 in the POD format,
    and its header ID 7 means NOAA-9 in POD.
  - Compressed files `.gz` / `.bz2` / `.xz` (e.g. Chiba's `n1297121720.l1b.bz2`) need not be unpacked: the listing decompresses only the
    start of the file, and plotting decompresses the whole file into memory (no decompressed file is left on disk; a 70 MB .bz2 takes
    about 2\~3 s). The decompressed size of a .bz2 is not known in advance, so the number of scan lines is first taken from the header
    and replaced by the actual length when the file is read. Files whose inner extension is on the skip list (e.g. HSD `.DAT.bz2`) are not
    checked; when both the compressed and the uncompressed version of a file are present, only the uncompressed one is listed.
* Channels: AVHRR/3 as in AVHRR-NSMC. CH3A / CH3B follows the scan-line bit field; when a KLM satellite is stored in the POD format,
  it follows the sign of the channel-3 calibration slope. AVHRR/2 (NOAA-7, 9, 11, 12, 13, 14) has CH1\~CH5;
  AVHRR/1 (TIROS-N, NOAA-6, 8, 10) has CH1\~CH4.
* Visible calibration: reflectance = albedo / 100, not divided by cos(solar zenith angle).
  - KLM: the operational coefficients in the records (dual gain, updated by NOAA).
  - POD: the single slope in the records (mostly pre-launch values; NOAA-14 monthly from 1996/11).
  - KLM satellite stored in the POD format (or a KLM file without usable coefficients): the NOAA STAR in-orbit coefficients
    (see the visible calibration of AVHRR-L0 below), otherwise the pre-launch dual-gain coefficients of Appendix D.
* IR calibration: brightness temperature T = (T\* − A) / B.
  - KLM: N = a0 + a1·C + a2·C². The coefficients in the records include the space radiance and the nonlinearity correction;
    a2 is scaled by 10⁶ in format versions <= 2 and by 10⁷ for CH4 / CH5 in versions >= 3.
  - POD: the 1b coefficients are only a two-point linear calibration (space and blackbody), and the space radiance depends on the satellite
    and the period (0 or negative; NOAA-9 changed several times in 1986\~87). They are therefore recomputed with the KLM model:
    the linear position from the space and blackbody counts in the telemetry, the blackbody radiance from the 1b coefficients,
    and then Walton's nonlinearity correction. For NOAA-14 this is identical to the guide's quadratic; for the other satellites
    it differs from the guide's Ns = 0 correction tables by about 0.1\~0.4 K.
  - Station products (both KLM and POD formats; the Chiba NOAA-12 1b coefficients are 0.3\~0.7 K below the PRTs), IR coefficients that are
    constant over the pass (Chiba NOAA-17), and KLM satellites stored in the POD format (the Chiba NOAA-18 IR slopes are off by about
    3.5 times): calibrated from the PRT, blackbody and space counts following KLM Guide 7.1.2.4.
    Verified on CLASS NOAA-15/18/19: within ±0.05 K of the operational coefficients.
    One scan line in five carries the PRT reference reading; the guide says 0, but in practice it is often 1\~5 (never 0 on POD
    satellites), so a line whose three readings have a median < 20 is taken as the reference. Station data often have bit errors in
    single readings, so each line uses the median of its three readings, and each of the four PRTs is cleaned of outliers and
    smoothed before averaging (the blackbody temperature of the Chiba NOAA-17 / 18 files used to jump by several K; it is now smooth).
* Navigation: 51 tie points (LAC/HRPT pixels 25\~2025, GAC pixels 5\~405) interpolated in earth-central angle.
  Scan lines with invalid navigation, or navigation discontinuous with the neighbouring lines, are masked.
  LAC/HRPT station products (Chiba etc.) are navigated from the TLE when one is available (as in AVHRR-L0; the navigation inside the files
  is 4\~11 km off the coastlines); `[tle] station_l1b = false` in `Satima-PO.toml` keeps the navigation in the file.
* Time: damaged time codes are replaced by the expected time. They are common at the start of Chiba files; CLASS files have
  occasional errors too, e.g. a day number of 19. Lines around data gaps are not flagged; scan-line numbers are unreliable at gaps
  and are not used.
* Merging sources: LAC, HRPT and FRAC of the same satellite have the same resolution and are merged into one image when
  selected together in the list (e.g. `[10,11]`); GAC forms its own group. Different sources of the same pass
  (CLASS LAC and station HRPT, CLASS and Chiba) are merged as follows:
  - Scan lines are aligned by time, then checked with the channel-4 counts: the same data should be almost bit-identical,
    and a time offset of a whole number of scan lines is corrected from the counts.
  - Each scan line uses the source with the lowest fraction of missing and spike pixels. The source that is best overall
    in the overlap is preferred; another source is used only when it is clearly better.
  - Missing or spike pixels in the chosen line are filled with good values from the other sources at the same position.
  - Example: the LAC (20:28–20:34) and HRPT (20:31–20:45) of Mitch (1998/10/26 NOAA-14) overlap by 903 lines with 99.94%
    identical counts; the merged image covers 1\~58N and the stripes at the start of the HRPT are replaced by the LAC.
    Synthetic tests: a 1-second time offset is aligned automatically, and 287 of 300 spikes in the primary source
    are filled from the other source.
* **Verification** (2026/11): all 30 sample files are read (NOAA-9\~19, GAC / LAC / HRPT, CLASS and Chiba).
  - Over warm ocean CH4 − CH5 is 1.3\~3.4 K, over cold cloud 0\~1.9 K.
  - Navigation and eye positions: the eyes of Haiyan (2013/11/07 GAC), Mitch (1998/10/26 NOAA-14) and Utor (2013/08/11 Chiba NOAA-18)
    match the best tracks, and the coastlines line up.
  - Bit-error speckle in station data (implausible split-window difference in about 10% of the pixels of the Chiba NOAA-17 file and
    about 2% of the NOAA-18 file; < 0.1% in CLASS files) and the stripes at the start and end of a pass are drawn as they are for now.
  - No samples: TIROS-N, NOAA-6/7/8/13, MetOp (NOAA 1b), FRAC.
  - Chiba `.l1b.bz2` (two NOAA-12 passes on 1997/12/17, NOAA-19 on 2014/10/07): identical to the decompressed files (navigation, time,
    calibration); the eye of Super Typhoon Paka (1997/12/17 21Z, 14.1N 140.7E) is clear. About 30% of the 2014 NOAA-19 scan lines are
    empty and some sections are very noisy; the L0 of the same pass is bit-for-bit the same, so this is a reception problem.
  - USGS EarthExplorer AVHRR Orbital Segments (M2M dataset `avhrr_orbseg`, 1992/04\~1999/09, not only North America):
    `ao14101797050456_012257.l1b` (NOAA-14 1997/10/17 05:04\~05:37, 33 continuous minutes) is POD with a 122-byte TBM and is read as is;
    Super Typhoon Ivan (05:21Z, 14.5N 130.6E) is complete (the Japanese station data of the same pass are broken near the typhoon,
    where the satellite was too low). It is listed with the source USGS: the file name is an EarthExplorer Entity ID (`AO` + 14 digits),
    or the file has a 122-byte TBM and the station code of its dataset name is `SF` (Sioux Falls; renamed files are recognized too).
    Chiba = station products such as processing block `B0000000` (a 122-byte TBM counts only when the file has no dataset name);
    everything else is CLASS.

### NOAA AVHRR HRPT raw frames (Level 0)

* Files: HRPT minor frames recorded directly by a receiving station, one frame per scan line (11090 10-bit words, each stored in 16 bits;
  little- or big-endian, detected automatically).
  - Uncompressed or `.xz` / `.gz` / `.bz2`; an uncompressed `.tar` can also be used as is (every file in it is listed, no need to unpack).
    Subfolders up to two levels below the data folder are searched too. A file present both in a `.tar` and in a folder is listed once.
  - Detected by the frame sync, not by the file name. The satellite comes from the spacecraft address in the ID word (AAPP's table,
    TIROS-N\~NOAA-19; addresses used by more than one satellite are resolved by year); `noaa18`, `n14` and the like in the path take precedence.
  - Raw frames carry no year: it is taken from a date in the path (`2003-02-04/`, `20130811.063140.noaa18.raw.xz`), otherwise from the
    modification time of the file (or tar member).
  - Chiba University CEReS L0 (`20141007163202.n19.l0.bz2`, `n1297121720.gi.bz2`): 16-bit big-endian words after a 240-byte ANSI tape label
    (`vol1` + `hdr1` + `hdr2`, 80 bytes each), which is skipped. The year comes from the file name: a date followed directly by the time
    (`yyyymmddhhmmss`), or `n` + satellite number + `yymmddhh` (used only when none of the original forms matches, so paths that were
    parsed before give the same result); `n12…` also identifies the satellite. Bytes 13\~16 of vol1 hold the size of the whole file and give
    the number of frames and the end time of a .bz2 (ignored when they do not fit the frame length).
  - Samples: the three tars in `D:\github\L0` (2003/02/04 NOAA-14\~17, 2013/08/11 and 2013/11/07 NOAA-15/16/18/19; 42 passes).
* TLEs: raw frames carry no orbit data, so the satellite position is computed from TLEs (SGP4; needs the `sgp4` or `pyorbital` package).
  - TLEs are read from `[tle] dir` in `Satima-PO.toml` (default `tle/`, including subfolders and .zip) and from the `.tle` / `.txt` files in
    the data folder; the set with the epoch closest to the observation is used, with a warning beyond `[tle] max_age` days (default 5).
    There is also a warning when the three nearest sets disagree by more than 20 km (historical archives occasionally contain
    contradictory elements, e.g. two interleaved families 160 km apart for NOAA-9 in 1988 at CelesTrak).
  - TLEs are not shipped with the program: run `python tools/fetch_spacetrack_tle.py` to download them from Space-Track (free account;
    put the account and password in `.spacetrack` in your home folder (`identity = …` and `password = …`, or just the account and the
    password on two lines), or in the environment variables `SPACETRACK_IDENTITY` / `SPACETRACK_PASSWORD`). The files go to
    `tle/spacetrack/` (kept out of version control per the Space-Track terms of use). By default it downloads every year of TIROS-N\~NOAA-19
    from launch to end of operation (about 195 queries, 10 minutes); `--from` / `--to` or satellite names restrict it, and years already
    downloaded are skipped. For 1980\~2004 the CelesTrak historical archives also work (no account; put the zips from
    https://celestrak.org/NORAD/archives/ in `tle/`).
* Navigation: the same geometry as AAPP (geodetic nadir, scan plane perpendicular to the inertial velocity, 25 µs per sample, intersection
  with the WGS84 ellipsoid), computed at the 51 tie points and interpolated as in AVHRR-L1B, plus a roll, yaw and scan step per satellite
  (`avhrr_hrpt.NAV_PARAMS`).
  - The parameters are fitted to coastlines: the IR image (land colder than sea on winter nights, the opposite on summer days) is matched
    with Satima's 30-arc-second land/sea mask in 128×128 blocks, and the roll, yaw and step of each satellite are solved by least squares
    (with one time offset per pass to absorb clock and TLE errors). Data: 14 ascending and descending passes near Japan on 2003/02/04 and
    2013/11/07; the RMS residual drops from 0.9\~3.3 km to 0.3\~0.5 km.
  - Results (sign convention of this program: roll > 0 towards pixel 2048): NOAA-18 roll −0.18°, NOAA-19 yaw +0.33°, NOAA-15 yaw +0.18°,
    NOAA-16 yaw −0.13° with a step of 0.05396° (AAPP: 0.053955°), the others 0.05406\~0.05409° (AAPP derives 0.0540865° from the sampling
    rate); all satellites share a roll of about −0.03°. Without them the scan edges are off by 5\~10 km. They are of the same size as the
    constant attitude corrections in NOAA 1b headers (NOAA-18 roll 0.15°, NOAA-19 yaw 0.25°) and AAPP's default biases. NOAA-14 has only two
    damaged descending passes (roll −0.076°, yaw −0.056°). NOAA-12 is fitted on two Chiba descending passes of 1997/12/17: roll +0.295°,
    yaw +0.12° (AAPP lists roll −4.9 mrad with the opposite sign), residual from 5.8 km down to 0.6 km; a constant along-track offset of about
    −3\~−3.8 km remains (both passes are descending, so a clock error cannot be told from a north-south offset).
    Satellites without a fit use AAPP's step and a roll of −0.03°.
  - Nadir: the KLM Guide and AAPP both use the ellipsoid normal (geodetic nadir). Against the coastlines the north-south error is about −1 km
    with the geodetic nadir and −3.5 km with the geocentric one (it changes sign between ascending and descending passes, so it is separable
    from clock errors). NOAA's operational navigation from 2013 on also matches the geodetic nadir, while the 1998 CLASS 1b is closer to
    the geocentric one.
  - Independent check (2013/08/11 passes, not used in the fit): most blocks improve from 1\~5 km to within 1 km.
  - No clock correction: the raw time codes lack NOAA's clock correction, leaving a constant along-track offset. Against the CLASS GAC of
    2013/11/07, NOAA-18 runs 0.225 s fast, about 1.5 km along track (the 2013/08/11 coastlines show about −1.8 km); NOAA-19 0.05 s,
    NOAA-15 0.003 s.
  - The navigation inside the Chiba L1B files (same passes) against the coastlines: NOAA-12 off by 4.5\~6 km across track, NOAA-19 by
    8\~11 km at the scan edges (no yaw correction); the TLE navigation (with the parameters above) is within 0.2 km across track. AVHRR-L1B
    therefore also navigates station 1b files from the TLE by default (`[tle] station_l1b`).
* Time: damaged time codes are repaired as in AVHRR-L1B. These frames are treated as missing: fewer than 4 matching sync words, flagged as
  test signal in the KLM ID word, or a repeat of the previous frame (same time code and minor frame number; as in AAPP).
  Reception gaps (e.g. about 90 frames missing between 00:59:59 and 01:00:14 in the 2003/02/04 NOAA-17 pass) show up as a band without data.
* Calibration:
  - IR always from the PRT blackbody temperature and the blackbody and space counts (as above). PRT coefficients: KLM Appendix D from
    NOAA-15 on; pygac values for older satellites (the same as AAPP's avhcal.txt), verified on CLASS 1b (NOAA-9/10/11/14): the blackbody
    temperature agrees with that of the 1b coefficients at the blackbody counts to within ±0.012 K (median). The NOAA-14 band corrections
    (ν, A, B) of pygac and AAPP differ, but the brightness temperatures agree to within 0.01 K.
  - End to end: the 2013/11/07 raw frames averaged to GAC pixels (GAC pixel n = LAC samples 5n\~5n+3; counts differ by 0.25) compared with
    the operational calibration of the CLASS GAC: CH4 median difference −0.020 K for NOAA-18, −0.033 K for NOAA-19, −0.022 K for NOAA-15
    (5th\~95th percentile about ±0.15 K).
  - The Chiba L0 and L1B of the same pass have 100% identical counts. The NOAA-12 PRT coefficients are the same as in AAPP's avhcal.txt
    (and the same set verified on CLASS for NOAA-11/14); the blackbody temperature implied by the Chiba 1b coefficients is 0.3\~0.7 K
    lower (the brightness temperatures likewise, e.g. Paka's eye 21.9°C against 23.5°C from the PRTs), so AVHRR-L1B now calibrates POD
    station products from the PRTs as well. The L0 and L1B now agree exactly in navigation, IR and visible.
  - Visible: NOAA-16\~19 use the NOAA STAR in-orbit coefficients (F. Wu's vicarious calibration over the Libyan desert, one set per month;
    the files shipped with AAPP, in `conversion_table/AVHRR/`), taking the set whose Center date is closest to the observation. Coverage:
    NOAA-16 2003/05\~2014, NOAA-17 2004/02\~2010, NOAA-18 2005/10\~2017, NOAA-19 2009\~2017. More than 120 days before the first set,
    the coefficients are interpolated between the launch date (pre-launch values) and the first set; more than 400 days from any set, the
    pre-launch values are used. The difference is large: in 2013 the NOAA-18 CH2 is 24.5% above pre-launch and NOAA-16 CH2 32%.
    NOAA-15 and older satellites have no STAR coefficients and use pre-launch values (dual gain from Appendix D from NOAA-15 on,
    POD Guide Table 3.3.2-1 for older satellites), so their reflectances are low.
  - 3A / 3B: as in AAPP, from the channel-3 space counts (at least 5 of 10 samples < 100 for 3A, > 900 for 3B) when they agree with bit 0
    of the ID word (1 = 3A); lines that disagree or are undecided are interpolated from their neighbours, then an 11-line median is taken.
* Comparison with AAPP 8.15 (source in `D:\github\L0\AAPP_8.15.tgz`): the geodetic nadir, scan step, address table, repeated frames,
  3A/3B and STAR coefficients above come from AAPP. AAPP does not handle the displaced and repeated segments caused by AVHRR scan-motor
  jitter (next item), and it times every scan line from the first one with the nominal 1/6 s.
* **Damage in the data itself** (2026/11 samples): in the 2003/02/04 passes about 22% of the NOAA-14 pixels and 2\~4% of NOAA-15/16 are
  displaced or repeated segments (identical runs of values appear in neighbouring frames at drifting positions, giving stair-like stripes);
  NOAA-17 and the 2013 files < 0.5% (NOAA-15 in 2013 about 3%). The pattern fits the MIRP re-sending old scan lines when the AVHRR scan
  motor jitters, not a decoding problem; it is drawn as is for now and can raise the reported eye temperature.

## True color correction check

| | Solar zenith angle correction | Rayleigh (CREFL) | O3 | H2O | Terrain height | Sharpening | Enhancement |
|---|---|---|---|---|---|---|---|
| VIIRS | already included in SDR | ✔ | climatology | ✔ | ✔ GMTCO Height | I1 ratio sharpening | Simple(sRGB) / WorldView / CIRA |
| MODIS | ✔ (tapering off from 88°) | ✔ | climatology | ✔ | ✔ MOD03 Height | B1 ratio sharpening | CIRA / WorldView |
| METimage | ✔ | ✔ (ISRF-weighted τ) | climatology | ✔ | — sea level | — | CIRA / WorldView |
| OLCI | ✔ | ✔ (including sea-level pressure) | ECMWF | — | ✔ altitude | unsharp mask | CIRA (white-point normalization) / WorldView |
| SGLI | ✔ | ✔ (τR from the JAXA ATBD) | climatology 300 DU (JAXA Koz) | — | — sea level | — (all 250 m) | CIRA / WorldView |
| MVISR | ✔ | ✔ (τR estimated from the band limits) | climatology 0.319 cm-atm | — | — sea level | — | CIRA / WorldView; plus inter-channel gains from thick ice cloud |

SGLI true color uses JAXA's L1B RGB combination R = VN08, G = VN05, B = VN03. VN08 / VN05 saturate over bright cloud (at about 0.49 / 0.67); saturated pixels are filled in from
VN11 (868.5 nm, low gain) to keep clouds from turning cyan; this fill is Satima's own processing, not an official JAXA method.

The CREFL `chand`, `csalbr`, the compositing formula, the solar zenith angle correction and the enhancement curves were originally copied separately in each of four scripts, mathematically exactly the same, and have been gathered into
`satima/po/truecolor.py`. The original scripts and the new version were compared with the same input and the results agree (differences within 1e-6). The coefficients specific to each satellite, the water vapor formulas,
terrain corrections, sharpening methods and enhancement fine-tuning remain in the individual scripts.

## Bugs of the original scripts that were fixed (Satima-PO period)

- VIIRS: the channel name `M10` was read as `M1`; M8–M11 and I3 were treated as brightness temperature, so no data could be read; I4, M12–M14 and M16 had no colormap menu (crash); for single channels, fill rows were removed separately from the data rows and from the latitude / longitude rows, so the two row counts could disagree
- AVHRR-NAT: visible channels raised a NameError; the radiance scaling of CH3A/3B should be 1e4; the scaling of the CH3B central wavenumber should be 1e2 (originally the CH3B brightness temperature was about -212°C and the whole image was blank); 3A/3B were not distinguished by scan line; the number of VEADR records and the dummy MDR were hard-coded; the radiance was printed once for every scan line
- AVHRR-NAT / AVHRR-NC: the band correction formula for brightness temperature was inverted. The official one (AVHRR L1B Product Guide §6.3.4) is T = A + B·T\*, but it was written as (T\* − A)/B, giving brightness temperatures that were too high: about 1.5–2.5 K for CH3B, about 0.1–0.5 K for CH4 and about 0–0.35 K for CH5, the difference growing as cloud tops get colder. Verified by integrating the AVHRR/3 spectral response functions (the A and B fitted from the SRF agree with the coefficients in the file; the error of the official formula is < 0.003 K)
- METimage: channel 11 (2.25 µm reflectance) showed the colormap menu and output 13 identical images; the last channel in the menu was misprinted as 18; the tar branch was missing an import
- SLSTR: visible channels changed from the approximation πL/1500 to πL/(F0·cosθs) (F0 taken from the product per detector, θs interpolated from the tie-point grid, verified by astronomical calculation); the ESA vicarious calibration coefficients are applied only to old-version products (new versions already include them); the processing that drew dark sea surface as white was removed; the channel format in the file names of the two dpi outputs was inconsistent (B8 / B08)
- SGLI: missing values (Error_DN) were treated as the maximum value (about 14% outside the scan range was drawn as bright gray); SW01–SW04 wrongly used the thermal infrared brightness temperature conversion table; the 1 km channels (SW01/02/04) were geolocated incorrectly; latitude / longitude now use spline interpolation in Cartesian coordinates (error about 1/14 of the original, about 10 times faster)
- Center time (the 7 scripts other than MODIS): it used to be estimated linearly from the start and end times of the file according to the middle row of the cropped range; now the actual observation time of the pixel closest to the center point is taken (what MODIS did originally), and the time scales of the formats (TAI / GPS / UTC) are all converted to UTC (`satima/po/obstime.py`). Measured errors of the old method: 24 s when the cropped range is cut off by the file boundary (VIIRS); 50 minutes when the beginning and end of a full-orbit NAT file overlap at high latitudes; 50 s at the edge of the SLSTR conical scan
- Shared: the colormap menu item `10.NHC` should be `10.AVN`; the range check for file selection used the length of the file name `len(h5f)`; the cropped range was one row and one column short; it crashed when there was no data inside the range; `eval(input())` was replaced by the safe `ast.literal_eval`
- MVISR (`SatimaV3-FY1D-NSMC.py`, when it was added in 2026/11): slopes and intercepts of the calibration coefficients were paired wrongly; for visible channels the albedo was treated as a radiance and divided by an assumed solar irradiance;
  brightness temperatures used 1/wavelength as the central wavenumber without band correction; the pixel positions of the tie points were wrong (edges about 24 pixels off); latitude / longitude were read as unsigned integers (scan lines in the southern and western hemispheres were skipped);
  the end-time field of the data header had the wrong type; the number of records to read was taken from the scan line count in the data header (which does not match the actual count); the whole data array was printed every time
- AVHRR-NSMC (`SatimaV3-NOAA-NSMC_accnumpy.py`, when it was added in 2026/11): the band correction of brightness temperature was written as T = A + B·T\* (NOAA KLM uses (T\* − A) / B; NOAA-18 CH4
  is about 0.3 K too warm at 200 K and CH3B about 1.6 K too warm at 300 K); for visible channels the albedo was multiplied by 100π / F and divided by 100, scaling reflectance up by about 2.4 (saturating above 41% albedo);
  the menu offered channel 3 but selecting it crashed; the tie points were spread evenly over pixels 0\~2047 (should be 25\~2025), using interp2d which scipy 1.14 removed; the QC dropped coefficients outside the 2nd / 98th percentiles
  (only for CH1 and CH4, always changing 4% of the scan lines of every file); all scan lines used the date of the first one; it crashed when there was no data in the range.
  The extra division by 1000 / 10 of the visible coefficients, the order of the infrared coefficients and the nonlinearity correction were correct (the NSMC format differs from NOAA KLM)

## Adding a satellite

Put a `<sensor>-<format>.py` in `sensors/po/`, structured like the other scripts: name the single-file class `<sensor>File` and the reader
`<sensor>Reader` (e.g. `AvhrrFile` / `AvhrrReader`).

```python
from satima.core.colormaps import ask_colormaps
from satima.core.i18n import tr
from satima.po import SatReader, ui

class XxxFile:                    # a single file
    def __init__(self, fullname):
        self.fullname, self.satellite = fullname, 'XXX'
        self.time = ...           # start time (UTC datetime)
        self.end = ...            # end time; leave it unset if the file name has none (or set GRANULE_SECONDS = nominal length in the reader)
    def extract(self, band): ...  # reflectance 0–1 / brightness temperature °C
    def geocoord(self): ...       # lats, lons (same size as the data)
    def row_times(self, nrows):   # observation time of each row (UTC datetime64; for the conversions see satima/po/obstime.py);
        ...                       # if the file has no per-row times, obstime.linear(start, end, nrows) can be used

class XxxReader(SatReader):       # the reader
    NAME = 'XXX-FMT'
    SENSOR = 'XXX'
    DESCRIPTION = '中文說明'
    DESCRIPTION_EN = 'English description'
    PATTERNS = ('XXX_*.h5',)
    SOURCE = 'XXX'                # source column of the data listing

    def make_file(self, name):
        return XxxFile(name)
    def select_band(self):
        self.band = ui.ask_int(tr('輸入頻道：', 'Input the band: '), (1, 2))
    def select_variants(self):    # infrared returns a list of colormaps, the others [None]
        return ask_colormaps(self.cfg, 'ir') if self.band == 2 else [None]
    def prepare(self):            # self.files = the files of this segment (sorted by time); when there are several, concatenate along track (self.concat)
        f = self.files[0]
        data = f.extract(self.band)
        lats, lons = f.geocoord()
        if self.slice_swath(lats, lons) is None:    # sets self.lats / self.lons / self.box
            return False
        self.data = self.crop(data)
        self.time = self.center_time(f.row_times(lats.shape[0]))   # observation time of the pixel closest to the center point
        return True
    def plot(self, canvas, variant):
        if variant is None:
            canvas.mesh(self.lons, self.lats, self.data, cmap='gray', vmin=0., vmax=1.)
            cm = ''
        else:
            cm = canvas.bt_mesh(self.lons, self.lats, self.data, variant)
        return dict(band='B%02d' % self.band, colorbar=cm)

READER = XxxReader
```

* The data listing, continuity rules, data selection, projection and range, colormap menu, coastline, annotation, file name and saving are all handled by the main workflow; nothing extra needs to be written.
* `canvas.mesh()` / `bt_mesh()` / `rgb()` only hand the data to be drawn this time over to the canvas; sampling is done only after the main workflow has decided the image size,
  so `plot()` is called only once for each colormap / style.
* The `band` and `colorbar` returned by `plot()` are filled into the file name and the annotation below the image.
* Formats with one file per band (like VIIRS) must override `segments()` to organize the files into granules and override `begin()` to pick the files for the selected channels; see `sensors/po/VIIRS-HDF5.py` for how it is written.
* Formats that can only read one file at a time set `MULTI_FILE = False`.
* When one reader reads several instruments (like MERSI-HDF5 with MERSI-II / MERSI-LL / …), put the instrument name in `Segment.extra['sensor']` (or override `sensor_label()`); the instrument column of the listing, the annotation and `{sensor}` in file names then use it.

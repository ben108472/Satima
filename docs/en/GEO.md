# Satima: Geostationary Satellites (GEO)

> This is the detailed manual for **GEO (geostationary satellites)**; its content comes from the original Satima-GEO. For how to run the integrated version, the data list, the layered config files
> (`Satima.toml` + `Satima-GEO.toml`), the interface language, the image modes (native / dpi) and the continuity rules for polar orbiters, see [README.en.md](../../README.en.md).
> Of the settings mentioned below, the common sections (`[general]`, `[map]`, `[coastline]`, `[image]`, `[title]`, `[annotation]`, `[output]`…)
> are written in `Satima.toml` and the GEO-only ones in `Satima-GEO.toml`; any section can be overridden in `Satima-GEO.toml`.
> The "native resolution" behavior described here refers to `[image] mode = "native"`; with `mode = "dpi"` the image size is determined by the dpi instead, and the resampling method is unchanged.


Supports Himawari-8/9 (AHI), GOES-16\~19 (ABI), GK-2A (AMI), MTG-I (FCI), MTSAT-1R/2 (Imager), COMS-1 (MI),
GMS-1\~5 / GOES-9 (VISSR), GOES-8\~15 (Imager, GVAR / AREA), SMS-1/2, GOES-1\~7 (VISSR / VAS, McIDAS AREA),
Meteosat-8\~11 (MSG SEVIRI, HRIT / native / NetCDF), Meteosat-2\~7 (MFG MVIRI, FCDR / NetCDF):
single bands, IR / water vapor colormap enhancement, CIRA true color (True Color), and GeoColor, which adds night-time clouds and city lights;
MTG's FDHSI and HRFI are listed separately, and the HRFI TrueColor is sharpened with the 0.5 km red band.

## Running

```
python Satima.py                  # search the data in data_dir of the config file
python Satima.py -d D:\data\h9    # specify the data directory
python Satima.py -f AHI-HSD       # search one format only
python Satima.py -o D:\png        # specify the output directory
python Satima.py --list           # list the available formats
```

Workflow: select data (several time slots allowed) → map projection (Enter = MERC) → map range (filled in automatically for regional scans; not needed for FD) →
products (several allowed) → colormaps → image width (Enter = native resolution; `x3` = 3 times the native resolution).

Settings are in `Satima.toml` (common) and `Satima-GEO.toml` (GEO-only items and overrides; UTF-8, Chinese text is allowed); items not written there use the defaults in `satima/core/config.py` and `satima/geo/config.py`.

Directories in the config files (`data_dir`, `output_dir`, `colormap_dir`, `conversion_table_dir`) may be absolute or relative paths;
relative paths are always relative to the directory containing Satima.py (by default `data/` holds the data and `plot/` the images), regardless of where the program is run from.
Relative paths given to `-d` / `-o` on the command line follow the usual convention and are relative to the current directory.

The map package is selected with `[map] backend`: `basemap` or `cartopy`. The resampling of the satellite image is exactly the same;
only the data source of the map lines differs: GSHHS for basemap, Natural Earth for cartopy.
`[map] resolution` accepts either `c/l/i/h/f` or `110m/50m/10m`; they are mapped automatically.

### Map projection

Chosen at run time (just pressing Enter = `[map] projection`, default merc; with `ask_projection = false` it is not asked):

| Menu | Projection | Description |
|---|---|---|
| MERC | Mercator | Default, same as the old version (no projection added to the file name) |
| CYL | Equidistant cylindrical | Equal spacing in longitude and latitude |
| LCC | Lambert conformal conic | The standard parallels `lcc_parallels` default to 1/6 and 5/6 of the north-south range, and the central longitude `lcc_central_longitude` defaults to the middle of the range; the frame is the smallest rectangle covering the whole lon/lat range |
| FD | Full disk (satellite view) | No lon/lat range needs to be entered; the sub-satellite point and satellite altitude are taken from the data |

The Earth model for MERC and LCC is always the WGS84 ellipsoid (same as PROJ `ellps=WGS84`; the lon/lat from satellite navigation are WGS84 geodetic latitudes anyway):
* basemap defaults to a sphere of radius 6370997 m; Satima-GEO changes it to the same ellipsoid as cartopy.
* The grids and images of the two map backends are therefore pixel-identical (Mercator, LCC and CYL each tested; lon/lat difference < 1e-12°); only the coastline data differ.
* Images from the old version (basemap sphere): the height of a Mercator image of the same range differs (e.g. 1692 -> 1685 rows for 5\~35°N), and positions differ by up to about half a pixel (0.9 km).

The FD frame is determined by the band of each product:
* For data in the standard geostationary projection (Himawari, GOES-R, GK-2A, MTG, MTSAT, MSG, MFG FCDR…), with a full-disk scan the frame is the original image itself
  (e.g. MTSAT VIS 11000x11000, IR 2750x2750); at native-resolution output, output pixels map one-to-one onto original pixels, which amounts to drawing the original image directly
  (verified pixel by pixel). In this case row/column numbers are converted linearly from scan angles without going through lon/lat, so even large images are fast.
* Regional scans (half disk, RSS, mesoscale etc.) still draw the whole disk; the frame is the Earth edge plus 1%.
* For old spin-scan satellites (GMS VISSR, GOES-1\~15 etc.) the original image is not an ideal geostationary projection (attitude, different east-west sampling interval);
  the frame is an ideal projection of the Earth edge plus 1%, and each output pixel is still sampled with that satellite's exact navigation.
* basemap does not support the sweep = x of GOES-R, so FD uses sweep = y instead (still sampled with exact navigation, just no longer one-to-one with original pixels);
  cartopy does support it.

File name field `{projection}`: empty for MERC (same file names as the old version), otherwise CYL, LCC or FD, e.g.
`Himawari-9.AHI.BAND-13.BD.FD.202611080400Z.png`.

## Directory layout

```
Satima.py                     entry point (shared by the three families)
Satima.toml, Satima-GEO.toml  common config file, GEO config file
satima/core/              core shared by the three families (satellite-independent)
  config.py                 reading of the layered config files, common defaults, per-format overrides [readers.<format>]
  session.py                main flow: searching for data, the sectioned data list, selecting data
  ui.py                     interactive input
  maps.py                   map projections (MERC / CYL / LCC / FD) and the map backend: basemap or cartopy ([map] backend)
  canvas.py                 image size (native / dpi), title and annotation, saving
  colormaps.py              colormaps (.cmap)
  naming.py                 file name templates
  parallel.py               multi-core (thread parallelism, [memory] workers)
  pngio.py                  multi-core PNG writing (parallel chunked compression)
satima/geo/               GEO code
  config.py                 GEO's own defaults and checks
  family.py                 the GEO flow: rows of the list, projection and range, product menu, drawing and saving each image
  reader.py                 reader base class GeoReader, Channel, Scene
  nav.py                    navigation: GeosNav (geostationary projection), ProjGridNav (regular projected grid), LonLatNav (per-point lon/lat)
  resample.py               output grid and resampling (reads only the window covering the map range)
  products.py               product menu and single-band / colormap drawing
  landsea.py                grayscale basemaps: land/sea mask (rasterized GSHHS) and compositing of clouds + the CWA land/sea / LIGHTS night-lights basemaps
  bluemarble.py             grayscale basemaps NASA / NICT: Blue Marble 2004 / 2002 (download, JPEG tiled GeoTIFF) and compositing with semi-transparent white clouds
  _resample_nb.py           numba kernels for resampling (block averaging, value lookup, full-disk row/column numbers and lon/lat)
  truecolor/                true color (Rayleigh correction, solar zenith angle correction, hybrid green, sharpening)
    kernel.py                 numba kernel for per-pixel compositing ([truecolor] engine)
    geocolor.py               GeoColor (night-time IR high clouds, low clouds, city lights and day/night blending)
    background.py             city lights (built-in tiled GeoTIFF reader, no GDAL needed; single-band DEFLATE, RGB JPEG) and land/sea mask (built into basemap)
sensors/geo/              readers for each satellite, file name <sensor>-<format>.py
  AHI-HSD.py                Himawari-8/9 Himawari Standard Data (.DAT / .DAT.bz2)
  ABI-NC.py                 GOES-16~19 L1b Radiances NetCDF
  AMI-NC.py                 GK-2A L1B NetCDF
  FCI-NC.py                 MTG-I FCI L1c FDHSI / HRFI (read directly from the zip)
  MTSAT-HRIT.py             MTSAT-1R / MTSAT-2 Imager JMA HRIT (.gz, or read directly from the whole .tar)
  MI-HDF.py                 COMS-1 MI L1B HDF-EOS5 (.he5) — rewritten from the old COMS_HDF.py, not yet validated with real data
  GMS-VISSR.py              GMS-1~5, GOES-9 JMA VISSR (VS/IR*.Z?? and VISSR_*.A.IMG, gzip, may be inside a tar)
  GMS-TAPE.py               GMS-1 VISSR tape images (FGGE 1978~1979; 36-bit words, file name is the volume number, e.g. VS2303)
  GMS-AREA.py               GMS McIDAS AREA (GMSX navigation, e.g. GMS_IR_1_1979_007_1533.area)
  vissr_common.py           (helper) field definitions of the JMA VISSR parameter blocks, calibration tables, navigation parameters and Earth-edge correction, shared by GMS-VISSR and GMS-TAPE
  vissr_nav.py              (helper) JMA VISSR navigation MGIVSR and Earth-edge correction, shared by GMS-VISSR, GMS-TAPE and GMS-AREA
  vissr_nav_cy.pyx          (helper) Cython version of the above; build_vissr_nav_cy.py compiles it into a .pyd
  GOES-GVAR.py              GOES-8~15 Imager raw GVAR (NOAA CLASS goes??.yyyy.ddd.hhmmss, NSMC GOS*_GVR_*.dat)
  GOES-AREA.py              GOES-8~15 Imager McIDAS AREA (CLASS GVAR_IMG *.BAND_0?)
  goes_imager.py            (helper) GOES I-M O&A navigation (same as McIDAS nvxgvar), Gould / BCD numbers, calibration
  goes_nav_cy.pyx           (helper) Cython version of the GOES navigation
  gvar_format.py, gvar_area_format.py  (helper) format definitions copied verbatim from the old GVAR_format.py and GVAR_AREA_format.py
  GOES-PREGVAR.py           SMS-1/2, GOES-1~7 VISSR / VAS McIDAS AREA (goesNN.yyyy.ddd.hhmm.AREA_IR / AREA_VIS)
  mcidas_goes_nav.py        (helper) McIDAS 'GOES' navigation (spin scan, same as nvxgoes / VisAD GOESnav) and exact inverse
  mcgoes_nav_cy.pyx         (helper) Cython version of the above
  pregvar_area_format.py    (helper) format definitions copied verbatim from the old AREA_format.py
  SEVIRI-HRIT.py            MSG SEVIRI EUMETSAT HRIT (H-000-MSG*, including -C_ compressed and .bz2)
  SEVIRI-NAT.py             MSG SEVIRI Level 1.5 native (MSG?-SEVI-MSG15-*.nat, or the Data Store .zip without unpacking)
  SEVIRI-NC.py              MSG SEVIRI Level 1.5 NetCDF (W_XX-EUMETSAT-Darmstadt,...SEVIRI_C_EUMG_*.nc) — no sample yet to validate with
  seviri_common.py          (helper) SEVIRI channels, calibration (nominal / GSICS / Meirink), navigation (including the HRV upper / lower windows)
  seviri_native_hdr.py      (helper) copy of the old prologue.py (satpy's native header structures), with np.bool changed to np.bool_
  MVIRI-FCDR.py             MFG MVIRI FCDR (MVIRI_FCDR-EASY/FULL_*.nc or .zip, FIDUCEO_FCDR_*.nc)
  MVIRI-NC.py               MFG MVIRI EUMETSAT NetCDF (...MVIRI_C_EUMS_*.nc, the format of the old MFG.py) — no sample yet to validate with
colormap/                 colormap files
conversion_table/         calibration tables (GK-2A etc.)
nightlights/              city lights for GeoColor: VIIRS 2024 annual mean night lights, 500 m (Cloud-Optimized GeoTIFF, 61 MB)
landsea/                  land/sea mask for the CWA basemap (GSHHS full resolution, 30 arc-seconds; generated automatically if missing)
bluemarble/               Blue Marble 2004 (about 1 km) / 2002 (about 1.85 km) tiled GeoTIFFs for the NASA / NICT basemaps (downloaded automatically if missing)
```

## How it works

Longitude/latitude of each output pixel → that band's navigation `nav.to_rowcol()` → row/column number in the original image → take the value.

* Only the window covering the map range is read (with a 2-pixel margin at the edges); neither the lon/lat of the whole image nor a KD-tree is needed.
  For a small area of GOES the amount read is tens to hundreds of times less than in the old version.
* The original pixel corresponding to each output pixel is stored as an integer index (the nearest pixel center), and lon/lat are computed in pieces, so large images (full-disk VIS etc.) need no whole-image lon/lat arrays.
  This is independent of the projection: MERC / CYL / LCC / FD all use the same resampling, and the coastlines are drawn with the same projection, so the two always match.
* Coloring is computed by Satima-GEO in float64 from the color-level definitions (colormaps.colorize); the 8-bit image is then handed to matplotlib to overlay coastlines and annotation.
  When matplotlib does the coloring itself it first scales in float32, so with very fine colormaps (e.g. BD, about 50,000 levels) a few pixels occasionally differ by one level; grayscale and true color are exactly the same.
* When the output is coarser than the native resolution, block averaging is done first; block origins are aligned to multiples of the factor, so the 2x2 blocks of the 0.5 km red band line up naturally with the 1 km blue pixels (no color fringing in true color).
* Row/column numbers always start from 1, following the CGMS convention (row 1, column 1 = upper-left corner).

### Taking values when enlarging (`[image] interpolation`, default bilinear)

When output pixels are finer than or similar to the original pixels (block-averaging factor = 1), there are two ways to take the value:

* **Nearest neighbour (nearest)**: takes the nearest original pixel directly. The magnification is usually not an integer and also changes with latitude (Mercator) and the satellite's oblique viewing angle,
  so some original pixels are drawn repeatedly as 2 output pixels and some as only 1. The picture therefore shows irregular stair-steps and a patchwork look, most visible at the color-level edges of colormap enhancements (BD etc.).
  * Example: GOES-19 C13 (nominal 2 km), Mercator map of 5\~35°N, 2612 pixels wide at native resolution.
  * Near the eye (17°N) about 8% of the original pixels are duplicated.
  * Near 30°N output pixels are about 1.84 km and original pixels about 2.1\~2.4 km; about 40% of the original pixels are repeated 2\~3 times.
* **Bilinear (bilinear)**: interpolates between the 4 neighbouring original pixels (offsets are quantized to 1/254 pixel and stored as int8; the nearest-pixel index is unchanged); edges are smooth.
  * Where a neighbour is missing (Earth edge, space) the nearest pixel is used instead, so boundaries are the same as with nearest neighbour.
  * Single-pixel extremes are diluted a little by their neighbours. The maximum / minimum cloud-top temperature is still computed from the original pixels, the same as with nearest neighbour.
  * When output pixels map one-to-one onto original pixels (Himawari, MTSAT, SEVIRI full disk at native resolution), the result is pixel-identical to nearest neighbour.

This applies to single bands, colormap enhancement and the night-time IR of GeoColor. The daytime bands of true color stay nearest neighbour: the 1 km bands are ratio-sharpened with the 0.5 km red band,
and the detail is restored by the red band. Results of strip processing, multi-core, and the numba and numpy versions are all pixel-identical.

#### Showing the original pixels: multiples of the native resolution (`[image] size = "x4"`; `x3` can also be entered at run time)

The old version's `m.pcolormesh` is equivalent to nearest neighbour: at the same image width, 91\~99% of output pixels take the same original pixel, and the rest differ by 1 pixel only at the boundary between two original pixels.
The old version does not look patchy because its width was fixed at 2000 pixels (reduced). Taking GOES-19 5\~35°N as an example, about 36% of the original pixels near the eye are not drawn at all.

Once projected onto a Mercator or similar map, satellite pixels vary in size, so a uniform output grid cannot have exactly one cell per original pixel: pixels are either duplicated or dropped, or else interpolated.
To see every original pixel with unchanged values, enlarge the image width to 3\~4 times the native resolution and use nearest neighbour (`interpolation = "nearest"`).
Each original pixel is then drawn as a block of about 3x3\~4x4 (i.e. the effect of a high-resolution pcolormesh):

| Image width (GOES-19 C13, 5\~35°N, 100\~50°W) | Original pixels not drawn near the eye | Block-size unevenness (eye / 30°N) |
|---|---|---|
| 2000 (old version) | about 36% | 0% / 22% |
| 2612 (native resolution) | about 0 | 24% / 44% |
| 10448 (x4) | about 0 | 15% / 13% |

Block-size unevenness = standard deviation / mean of the number of output pixels occupied by each original pixel.

* `size` can be `"native"`, `"x2"`\~`"x8"` (the width used when you just press Enter at run time) or a fixed width. Entering `x3` (or `3x`) at the image-width prompt applies to that run only.
* Native resolution is always enlarged by an **integer factor**, so that every original pixel has the same size:
  * When the width is below `min_size`, the smallest factor that reaches the lower limit is used, e.g. 357 pixels -> x6 = 2142. Previously it was enlarged straight to 2000 (5.6 times), and original pixels were drawn irregularly as 5 or 6 output pixels.
  * When the width is below `min_size` it is enlarged by at least `min_scale` (default 4), e.g. 1500 pixels -> x4 = 6000 rather than x2 = 3000;
    with `min_scale = 1` only the smallest factor that reaches the lower limit is used. Nothing changes when the native width already reaches `min_size`.
  * When the long side exceeds `max_size`, the largest factor that fits is used instead. Only if even the native resolution does not fit is the image shrunk to the limit.
* The resampling of basemap and cartopy is exactly the same (the satellite image is resampled by Satima-GEO; the map package only draws coastlines), so the duplication on Mercator and similar maps is the same for both:
  near 30°N both have 40\~50% of the original pixels duplicated. The only difference is the GOES-R full disk: basemap does not support sweep = x and is not aligned with the original pixels
  (5489 wide at native resolution), whereas cartopy maps one-to-one onto the original 5424x5424 image.
  * Previously basemap's Mercator was a sphere and cartopy's the WGS84 ellipsoid; the grids were offset by about half a pixel and the duplicated rows fell at different places, so the same spot looked slightly different.
    Both now use WGS84 (see "Map projection") and are pixel-identical.
* When the full disk is aligned with the original pixels (Himawari, GK-2A, MTG, MTSAT, SEVIRI; GOES-R needs cartopy), native resolution already has one output pixel per original pixel and no factor is applied.
  Enlarging would only copy each original pixel into N x N; the 0.5 km and 1 km bands would also hit the long-side limit and become a non-integer factor (e.g. 22000 x4 -> 25000, 1.14 times),
  bringing back irregular duplication. GOES-R with basemap, and the full disks of old satellites such as GMS and GOES-1\~15, are not aligned and are still enlarged.
* Satellites whose original pixels are not square (scanning instruments such as GMS, SMS / GOES-1\~15): the east-west sampling interval is finer than the north-south scan-line spacing,
  e.g. GOES-8\~15 IR is 2.29 km east-west × 4.01 km north-south (visible 0.57 × 1.00 km), and GMS-5 IR is 3.43 × 5.02 km.
  * Output pixels are square. "Native resolution" is based on the **finer direction** (2.29 km for GOES IR), so all original pixels are kept; in the coarser direction each original pixel is drawn over several output pixels
    (about 1.75 for GOES, about 1.47 for GMS-5). The sizes in both directions and the resolution adopted are printed at run time.
  * Previously the nominal north-south resolution (4 km for GOES IR, 5 km for GMS) was used and the east-west direction was compressed: a full disk used only 57% of the original pixels (GOES); now it is 100%,
    and the image width changes from 2744 to 4793 (GOES IR), from 2195 to 3204 (GMS-5 IR) and from 10977 to 19172 (GOES visible).
  * The pixel size is measured from the navigation (the scan-angle spacing of adjacent pixels near the image center, see `maps.source_pixel_km`), so nothing has to be set separately in each reader;
    if a file such as a McIDAS AREA has been subsampled, what is measured is the actual spacing in the file. Pixels of the standard geostationary projection (Himawari, GOES-R, MTSAT, SEVIRI etc.) are square and are not affected.
  * The native resolution of every projection (Mercator, LCC etc.) is based on this. The number of duplicated rows in the north-south direction is not an integer; to make every original pixel the same size, additionally use `x4`
    (GOES: each original pixel is about 4 × 7 output pixels).
* The number of pixels grows with the square of the factor. The x4 example above is 10448x6770, takes 10.6 s, and the PNG is about 11 MB.
* Coastlines and text are scaled in proportion to the image width (the figure is fixed at 5 inches wide), so their thickness relative to the satellite pixels does not change. If they look too thick when zoomed in, reduce `[coastline] linewidth`.

### True color (`[truecolor]`)

Follows the daytime true color of CIRA GeoColor (SHAC True Color; Miller et al. 2016 BAMS, 2020 JTECH):

1. Solar zenith angle correction.
2. Rayleigh scattering correction (`rayleigh_correction`). Blue, green, red and near-IR are all corrected, and the path reflectance above cloud tops is reduced according to the 10.3 µm brightness temperature (no reduction at 283 K, 70% reduction at 233 K).
   * `"full"` (default): after pyspectral's path reflectance is subtracted, the result is divided by the downward and upward transmittances and the spherical albedo of the atmosphere is taken into account. This is the SeaDAS / CREFL approach, and CIRA's is also adapted from SeaDAS.
   * `"subtract"`: the old version's approach, which only subtracts the path reflectance. Clouds turn yellowish, clear sky darker, and areas near the terminator orange.
   * `cloud_rayleigh_transmittance = true` (default): the transmittance and spherical albedo over high clouds also use the same reduction factor (only that fraction of the atmosphere remains above the cloud top).
     If only the path reflectance is reduced while the transmittance still uses the whole atmosphere, the blue band of high clouds is divided by too small a transmittance (at a solar zenith angle of 80°, about 0.60 for blue and about 0.85 for red), and cloud tops turn bluish.
   * `ozone = 0.30` (default, atm-cm; 0 or false turns it off): absorption correction for the ozone Chappuis band; reflectance is divided by the ozone transmittance.
     * Ozone absorbs red the most and blue the least. The absorption coefficients are Leckner's (1978) table averaged over each band's spectral response; for AHI 0.47 / 0.51 / 0.64 µm they are
       0.011 / 0.040 / 0.082 per atm-cm; the slant-path factor uses spherical geometry with an ozone-layer height of 22 km (SPCTRAL2).
     * The ozone amount is a fixed value. The real amount is about 0.24\~0.45, lower in the tropics and higher at mid and high latitudes in spring, so there are residuals at high latitudes.
   * These two items are not in the CIRA papers; they were added because cloud tops turned bluish at sunrise and sunset. Blue/red ratio of deep convective clouds (10.4 µm < −55°C) after correction (Himawari-9, gain not applied):

     | Solar zenith angle | < 50° | 50\~70° | 70\~78° | 78\~82° | 82\~86° |
     |---|---|---|---|---|---|
     | Both off | 1.16 | 1.17\~1.26 | 1.35 | 1.42 | 1.53 |
     | Both on | 1.06 | 1.06\~1.07 | 1.09 | 1.10 | 1.12 |

     With both off, the output is pixel-identical to that before they were added. The comparison table with CIRA SLIDER below is the result with both off.
3. Green
   * AHI / AMI / FCI: 0.91 × 0.51 µm + 0.09 × 0.86 µm (hybrid green, `hybrid_green_fraction`).
     The CIRA paper says about 7%. Compared with CIRA's actual images, at 7% the vegetation of all four satellites is redder than CIRA's; at 9% the hue is closest, and deserts become only 1\~2 levels greener.
   * ABI has no green band; `abi_green = "ahi"` handles it the CIRA way: first synthesize the 0.51 µm green band based on AHI, then blend in 0.86 µm in the same way.
     * CIRA's lookup table is not public; here G0.51 = 0.976 × blue^0.668 × red^0.347 was fitted from two Himawari-9 time slots (0\~50 degrees N). The median error against AHI's actual green band is less than 1 level.
     * It can also be set to `[blue, red, near-IR]` linear weights, e.g. the old version's CIMSS `[0.45, 0.45, 0.10]`. With these weights land is noticeably too green.
4. Ratio sharpening with the 0.5 km red band, terminator fade-out, log10 stretch (reflectance 0.025\~1.20 → [-1.6, 0.176], the same as the CIRA paper).
5. Tone curve (`tone`) and saturation (`saturation`); with the defaults `tone = "cira"`, `saturation = 1.0` the output is bit-identical to that before these two items were added.
   * The original CIRA formula has two inherent limitations:
     * Reflectance is clipped at 1.2 but the white point is at 10^0.176 = 1.5, so the brightest output only reaches 241; thick clouds (reflectance around 0.75) are about 212 and look grayish,
       and reflectances 0.6\~1.0 span only 32 levels, leaving very little gradation in typhoon cloud tops.
     * Reflectance below 0.025 is always 0. After Rayleigh correction the red band over clear-sky ocean is about 0.02, so red is 0 for 50\~60% of ocean pixels and the sea appears as blotchy blue-black.
   * `tone = "enhanced"`: the same log curve; the midtones are almost unchanged and only the two ends are modified.
     * The white point is lowered to `white_point` (default 1.25), and display values above `shoulder` (0.95) are compressed smoothly towards 1, so the brightest cloud tops do not become a flat blown-out white.
     * In the shadows a smooth floor (`shadow_toe`, default 0.006) replaces the hard clip, so reflectances below 0.025 still have gradation.
   * `saturation`: after the stretch, color differences are amplified about the luminance. Below a display luminance of 0.50 the configured factor applies; it decreases above that, and above 0.85 there is no amplification, so clouds are not tinted.
   * `tone = "enhanced"`, `saturation = 1.25` compared with the original CIRA formula (median 8-bit levels):

     | Satellite (time UTC) | Curve | Clear-sky ocean | Ocean with red = 0 | Vegetation | Bright cloud | Image maximum |
     |---|---|---|---|---|---|---|
     | Himawari-9 (2026-10-21 01:00) | cira | 0, 28, 57 | 52% | 66, 85, 64 | 212, 213, 218 | 241 |
     | | enhanced | 9, 33, 68 | 0% | 66, 91, 65 | 222, 224, 228 | 254 |
     | GOES-16 (2019-09-01 16:00) | cira | 0, 37, 67 | 57% | 70, 94, 76 | 209, 214, 218 | 241 |
     | | enhanced | 6, 43, 80 | 5% | 69, 101, 77 | 219, 224, 228 | 251 |

   * The pixel-by-pixel comparisons with CIRA SLIDER below are all results with `tone = "cira"`; with enhanced, clouds and ocean are brighter than CIRA's.

**Relative band gain for Himawari-9 (`band_gain`)**: Himawari-9 clouds are cyan-greenish in true color (before the fix, the cloud levels were about +5 for green−red and about +9 for blue−red).
Each band was checked with deep convective clouds (10.4 µm < −50°C, red reflectance > 0.8); such clouds are almost white in the visible:

| Satellite (time UTC) | Blue | Green | Red | Near-IR | Blue/Red | Green/Red |
|---|---|---|---|---|---|---|
| Himawari-8 (2020-09-04 03:00) | 0.927 | 0.902 | 0.897 | 0.919 | 1.035 | 1.008 |
| Himawari-9 (2026-10-21 01:00) | 0.992 | 0.934 | 0.886 | 0.910 | 1.124 | 1.057 |
| GOES-19 (2026-10-28 13:00) | 0.885 | — | 0.882 | 0.926 | 1.006 | — |

(Reflectances are medians where the solar and satellite zenith angles are < 40°.) Red and near-IR agree across the three satellites; only B01 and B02 of Himawari-9 are too bright.
Two other Himawari-9 time slots (2023-07-01 03:00, 2026-11-09 05:00) also give blue/red of 1.11\~1.12 and green/red of 1.04\~1.05, so this is not a one-off.
Using the updated calibration coefficients in the file header does not remove the difference either.
(The table above is TOA reflectance, without ozone or Rayleigh correction.)

After adding the ozone and high-cloud transmittance corrections, the corrected ratios for the same deep convective clouds (solar zenith angle < 50°, gain not applied):

| Satellite (time UTC) | Blue/Red | Green/Red |
|---|---|---|
| Himawari-8 (2020-09-04 03:00) | 0.98 | 0.98 |
| Himawari-9 (2026-10-20 22:00, 10-21 01:00) | 1.06, 1.07 | 1.02, 1.02 |
| GOES-16 (2019-09-01 16:00) | 1.06 | — |
| GOES-19 (2025-10-28, 2026-10-28 13:00) | 0.95, 0.95 | — |

The calibration differences between satellites remain, so `band_gain` is kept. By default true color / GeoColor multiply Himawari-9 B01 by 0.94 and B02 by 0.98;
after this the two ratios are both 1.00 for solar zenith angles < 70°, and 1.03 and 1.01 at 80° (the old values, with the two corrections off, were 0.92 and 0.95).
`band_gain` is written as `{ satellite = { band = factor } }`; other satellites can be added too; single-band images are not affected; set `band_gain = false` to turn it off.

Pixel-by-pixel comparison with CIRA SLIDER GeoColor; the figures are the median absolute level difference (0\~255), given as "new version / old version":

| Satellite (time UTC) | Red | Green | Blue | Vegetation hue (green−red) |
|---|---|---|---|---|
| GOES-19 CONUS (2026-11-08 18:01) | 4 / 3 | 4 / 8 | 4 / 11 | 0 / +8 |
| Himawari-9 full disk (2026-11-09 05:00) | 2 / 4 | 3 / 5 | 5 / 7 | 0 / −6 |
| GK-2A full disk (2026-11-08 04:00) | 3 / 4 | 4 / 5 | 8 / 5 | +2 / −5 |
| MTG-I1 full disk (2026-11-09 13:20) | 2 / 4 | 3 / 10 | 5 / 20 | −4 / −10 |

* The vegetation hue column is "our (green−red) − CIRA's (green−red)": negative means redder and not green enough compared with CIRA, positive means greener.
* Cloud and land hues agree: clouds are white, and the arid surfaces of the western United States, the Sahara and Australia are tan-brown just as in CIRA.
* Remaining differences:
  * CIRA's products for the different satellites are not fully consistent with each other. For example, the clear-sky blue of GK-2A is about 12 levels darker than this program's, while that of MTG is about 10 levels brighter.
  * Over shallow water CIRA is more cyan-green. CIRA's documentation also notes that its lookup table may be noisy or incorrect over shallow water.
  * These differences probably come from CIRA's unpublished SeaDAS parameters and lookup tables.

### GeoColor (`[geocolor]`, the second true color style)

Follows CIRA GeoColor V2.0 (Miller et al. 2020 JTECH, Section 5). When the data have the 3.9 µm and 10.4 µm bands (AHI, ABI, GK-2A, MTG),
the RGB group of the product menu gets a GeoColor entry after TrueColor:

* **Daytime**: the same as TrueColor, including the terminator fade-out at SZA 80°\~90° (CIRA's daytime layer also darkens in the twilight zone).
* **Night-time**: three layers stacked from top to bottom.
  1. High clouds: the colder the 10.4 µm brightness temperature, the whiter.
  2. Low clouds: the 10.4 − 3.9 µm brightness temperature difference is drawn in light blue; land and water use different thresholds, and the land/sea mask comes from basemap's built-in data.
  3. Surface: city lights are drawn in the orange-yellow of sodium lamps; other land is dark purple and water is black.
* **Terminator**: day and night are blended according to cos(solar zenith angle).
  * At sunrise and sunset 3.9 µm contains reflected sunlight, so low clouds fade briefly; during this time the daytime layer dominates, and CIRA shows the same thing.
* **City lights**: the VIIRS 2024 annual mean night lights in `nightlights/`, at the full 500 m resolution; the overview level inside the file is chosen automatically according to the output pixel size.
  * Data source: EOG VNL v2, taken from the Zenodo mirror (Hengl, OpenGeoHub, doi:10.5281/zenodo.17294744, CC BY 4.0).
  * Read with the built-in tiled reader; GDAL is not needed. The GEOS bundled with basemap conflicts with GDAL's version, and GDAL cannot be imported when basemap is loaded first.
  * Set `nightlights = ""` to draw no lights.
* **Blue floor for the ABI synthetic green (0.3 × red)**: the synthetic green is fitted from AHI. Near the terminator the Rayleigh correction removes almost all of the blue; without a floor green would also approach 0 and the twilight zone would turn red.
  This also affects TrueColor.

The default parameters are taken from the images currently on CIRA's SLIDER (after the 2025 update) and differ from the paper in three places:
* IR range `ir_range = [178, 301]` K (paper: [200, 280]). In both, the cold end increases by 20 K with latitude.
* Day/night blending `twilight_cos = [0.0, 0.18]`, `twilight_power = 1.75`, a transition over SZA of about 80°\~90° (paper: [0.1, 0.3], 1.5, about 72°\~84°).
  CIRA's description is "narrowing the twilight zone, extending daytime by up to about one hour".
* The water background is black. The gray of CIRA's clear-sky ocean at night comes from the IR layer.

These parameters were inferred pixel by pixel from GOES-19 (2026-11-12 23:01 UTC, sunset crossing the United States) and MTG-I1 (2026-11-09 13:20 UTC); the two satellites give consistent results.
Level comparison of each GOES-19 night-time layer with CIRA (medians):

| Region | This program RGB | CIRA RGB |
|---|---|---|
| Clear-sky land | 54, 47, 64 | 60, 51, 65 |
| Clear-sky ocean | 14, 14, 14 | 13, 13, 13 |
| Low clouds (land) | 176, 211, 251 | 177, 212, 252 |
| High clouds (< 230 K) | 165, 165, 166 | 166, 166, 166 |
| City lights | 201, 176, 144 | 191, 167, 135 |
| Twilight SZA 70°\~84° | 108, 99, 98 | 104, 104, 105 |

City lights are about 10 levels brighter than CIRA's, possibly because CIRA's night-lights data are from an older year (2015 in the paper).
For the look of the original paper, just change the three parameters above back to the values in parentheses.

### Grayscale basemaps (second level under GRAY)

After `1.GRAY` is chosen as the colormap, the basemap is asked next (several allowed, 0 = all):

```
BAND-13 colormaps: 0.ALL 1.GRAY 2.BD 3.COLOR ...
  Select colormaps (several allowed, 0 = all): 1
BAND-13 GRAY basemaps: 0.ALL 1.NONE 2.CWA 3.LIGHTS 4.NASA 5.NICT
  Select basemaps (several allowed, 0 = all):
```

* `NONE`: ordinary grayscale image (no basemap).
* `CWA`: land/sea basemap (`[landsea]`, see the next section). Available for both IR and visible / near-IR (visible has grayscale only, so the basemap is asked directly).
* `LIGHTS`: night-lights basemap (`[nightmap]`), IR only.
* `NASA`, `NICT`: Blue Marble basemaps (`[bluemarble]`, `[nictmap]`), available for both IR and visible / near-IR
  (the menu for visible is `1.NONE 2.CWA 3.NASA 4.NICT`).
* Water vapor does not see the surface, so it has no basemap and nothing is asked.
* The name is the `{colorbar}` field of the file name, e.g. `BAND-13-CWA`, `BAND-13-LIGHTS`, `BAND-13-NASA`, `BAND-13-NICT`; it can be changed with `name` in each section.

### Land/sea basemap CWA (`[landsea]`)

Grayscale clouds are overlaid on single-color land (olive green) and ocean (dark blue), imitating the "color infrared cloud image" on the CWA official website.

* Compositing: output = a × gray + (1 − a) × base color.
  * IR: gray is a piecewise-linear mapping of brightness temperature (`ir_gray`, +80°C black \~ −95°C white); a = clip((26.1 − T) / 70.8, 0, 1)^0.8,
    fully transparent above about +26°C and opaque below −45°C.
  * Visible: gray = reflectance; a ramps from reflectance 0.10 (fully transparent) to 0.50 (opaque).
* There are only two base colors, blended at coasts according to the land fraction. The terrain texture over land on the official website is actually the IR itself (mountain surfaces are colder at night, and the alpha and gray level change accordingly),
  not a separate hillshade image.
* The defaults were inferred pixel by pixel from the official website's Himawari-9 full-disk "color infrared" and "grayscale infrared" images of 2026-11-13 11:30 UTC against the B13 brightness temperature of the same time slot
  (alpha fitting error 0.009). An image made by this program for the same time slot differs from the official one, pixel by pixel, by a median of R 3, G 3, B 5 levels
  (the official image is a JPEG and has yellow coastlines and lat/lon lines).
* Land/sea mask (`landsea/gshhs_f_30s.tif`): basemap's built-in GSHHS full-resolution coastlines (lakes count as water, islands in lakes as land)
  are rasterized with a numba scan-line algorithm into land fraction at 30 arc-seconds (about 0.93 km) (4× vertical sampling, exact horizontal coverage fraction; coasts are smooth without jaggies)
  and stored as a tiled GeoTIFF with overview levels at 2, 4, 8… times. When sampling, the coarsest level that is not coarser than the output pixel is chosen and bilinearly interpolated; across the date line only the tiles actually used are read.
  If the file does not exist it is generated automatically on first use (about 30 s with numba; without numba matplotlib is used instead, about 12 minutes).
* By default coastlines, borders and lat/lon lines are not drawn (`lines = false`); set it to true to draw lines according to `[coastline]` and the related settings.

### Night-lights basemap LIGHTS (`[nightmap]`)

Clouds are overlaid on the VIIRS annual mean night lights (city lights), in the same way as GeoColor's night-time layers but without the low-cloud layer that needs 3.9 µm:

* Base color: `land_color` (dark purple) / `water_color` (black) of `[geocolor]`, blended at coasts according to the land fraction of `[landsea] mask`;
  where there are city lights it becomes the orange-yellow of sodium lamps (same as GeoColor). The lights file is `[geocolor] nightlights`; if it is not found this product is skipped.
* Clouds: output = N_IR × white + (1 − N_IR) × base color, N_IR = 1 − N(brightness temperature)[IRmin, IRmax], i.e. the alpha of GeoColor's high-cloud layer
  (`[geocolor] ir_range`, default 178\~301 K, the cold end increasing by 20 K between latitudes 30°\~60°).
* Changing the `[geocolor]` parameters affects both GeoColor and LIGHTS. `[nightmap]` has only two items, `name` and `lines`.

### Blue Marble basemaps NASA / NICT (`[bluemarble]`, `[nictmap]`)

White clouds, semi-transparent according to brightness temperature / reflectance, are overlaid on NASA Blue Marble, imitating the view of each band over Blue Marble on NICT's "Himawari-8 Real-time Web"
(https://himawari8.nict.go.jp/). The two differ only in the basemap; the cloud compositing is the same (its settings are all in `[bluemarble]`):

| Menu | Config section | Basemap | Resolution | Appearance |
|---|---|---|---|---|
| `NASA` | `[bluemarble]` | Blue Marble Next Generation 2004 (NASA GIBS) | about 0.98 km | Darker: ocean nearly black, snow cover on the Tibetan Plateau |
| `NICT` | `[nictmap]` | Blue Marble 2002 `land_shallow_topo` (the version used by the NICT website) | about 1.85 km | Brighter: ocean a uniform blue-violet (RGB 10, 10, 51), shallow water cyan-green, no snow on the plateau |

* This view on the NICT website has two layers: the cloud image of each time slot is a **pure white + alpha** PNG overlaid on a fixed Blue Marble.
  So the compositing is: output = a × white + (1 − a) × Blue Marble. Pixels without data show only the basemap, and outside the Earth it is transparent.
* Alpha: inferred pixel by pixel from the NICT website's cloud images (10d level, 5500-pixel full disk) against Himawari-9 HSD of the same time slot:
  * IR (`ir_alpha`): B13, the 2026-10-20 21:00 UTC full disk and the northern half at 2026-11-09 05:00 UTC.
    Opacity is linear in **10.4 µm radiance** (not linear in brightness temperature): fully transparent above 10.5°C, 230.5 (about 90%) at −80°C, about 235 at −84°C.
    Both time slots differ from NICT by only 0.4 level (8-bit quantization error); if it were linear in brightness temperature the difference would be up to 28 levels.
    All IR bands use the same brightness temperature → opacity curve.
  * Visible (`vis_alpha`): B03, 2026-11-09 05:00 UTC. Opacity = reflectance × 255 (NICT's 2 km is a 4×4 average of 0.5 km,
    with no solar zenith angle correction), a difference of 0.3 level. This is exactly the default of `[visible]` (no correction, factor 1, no gamma);
    changes to `[visible]` are applied here too.
* The basemaps (`image`) are all NASA public-domain images (credit NASA Earth Observatory / Goddard Space Flight Center when citing),
  stored as tiled JPEG GeoTIFFs (with overview levels; open-ocean tiles with identical content share one copy of the data; readable by GDAL).
  When sampling, the coarsest level that is not coarser than the output pixel is chosen and bilinearly interpolated (same as the land/sea mask).
  If the file does not exist it is downloaded automatically on first use; if the download fails this product is skipped and the reason printed, with no retry within the same run.
  * **NASA**: Blue Marble Next Generation (MODIS 2004 composite), taken from the NASA GIBS WMTS layer `BlueMarble_NextGeneration`
    (EPSG:4326). Level 6 is about 0.98 km, 40960×20480, 4273 tiles, 67.5 MB on disk, about 6 minutes as measured.
    The 512×512 JPEG tiles from GIBS are used unchanged as the tiles of the GeoTIFF (no recompression), and the overview levels (about 2, 4, 8… km) are also taken directly from the coarser GIBS levels.
    For a finer basemap (level 7, about 0.5 km, an estimated 250 MB and about half an hour to download) run
    `python -c "from satima import bluemarble; bluemarble.build('bluemarble/BlueMarble_NextGeneration_GIBS_500m.tif', 7)"`,
    then point `image` to that file.
  * **NICT**: Blue Marble 2002 `land_shallow_topo` (land, shallow water and terrain shading), taken from Wikimedia Commons'
    [Land shallow topo 2048.jpg](https://commons.wikimedia.org/wiki/File:Land_shallow_topo_2048.jpg)
    (despite the 2048 in the file name it is actually 21600×10800, 1/60°, 57.6 MB). It is cut into 512×512 tiles and recompressed as JPEG (quality 92), 44.9 MB on disk;
    the conversion temporarily needs about 1.5 GB of memory.
    Wikimedia sometimes rate-limits (HTTP 429 Too many requests): try again later, or download the file with a browser first and then run
    `python -c "from satima import bluemarble; bluemarble.build_bm2002('bluemarble/BlueMarble_2002_land_shallow_topo_2km.tif', r'downloaded_file.jpg')"`.
* **Comparison with the NICT website's basemap** (NICT's `FULL_24h/BlueMarble`, full-disk 4d level and East Asia 10d level):
  * NICT uses Blue Marble 2002: the mean color differences over land and ocean are both within 1 level (land −0.2, ocean R +0.9 / B −1.0).
    NASA (the 2004 version), by contrast, is 11\~26 levels darker over land and 33 levels darker in the blue of the ocean.
  * The pixel-by-pixel texture still differs (land rms about 15 levels): NICT's basemap is offset from the true position by about 2\~3 pixels (4\~6 km in East Asia, towards the northwest),
    and after shifting, the rms drops to 9. Satima-GEO's Blue Marble 2002 is aligned with the GSHHS coastlines to within half a pixel (about 1 km), so it does not follow NICT's offset;
    the clouds in both are placed at the correct position by satellite navigation.
  * eoimages.gsfc.nasa.gov, where NASA hosts the original Blue Marble files, could not be reached during testing, so the files are obtained from GIBS and Wikimedia Commons instead.
* By default coastlines, borders and lat/lon lines are not drawn (`lines = false`, set separately in `[bluemarble]` and `[nictmap]`).

### Low-memory mode (`[memory] low_memory`, on by default)

The output image is processed in strips: each strip computes its own row/column numbers, block averaging, sampling, true color compositing and coloring, and is released as soon as it is done;
only the final 8-bit image is kept. The source data of each band are still read only once (the read window is estimated by sparse sampling, and an extra read is made when a strip falls outside it),
so the readers need no changes and the result is pixel-identical to whole-image processing (compared for MERC / CYL / LCC / FD, true color, colormaps and each satellite).

Himawari full disk (peak of the Python process, measured on a 16 GB computer):

| Product | Whole image | In strips |
|---|---|---|
| True color 8000x8000 | 9.4 GB | 4.3 GB |
| B03 8000x8000 | 6.2 GB | 2.9 GB |
| B01 native resolution 11000x11000 | 6.3 GB | 2.5 GB |

A full disk still requires the whole source image to be read in (B03 22000x22000 in float32 is about 1.9 GB); this part is not affected by strip processing.

### Multi-core (`[memory] workers`, default auto) and numba acceleration (`[truecolor] engine`)

* **Strips in parallel**: after the output image is split into strips, several threads process them at the same time (satima/core/parallel.py). numpy, numba, zlib and bz2 release the GIL while computing,
  so threads are enough to use multiple cores and large arrays need not be copied. The result is pixel-identical to single-core.
* **numba kernels**: the per-pixel true color compositing is fused into a single loop (truecolor/kernel.py), including solar / satellite angles, Rayleigh lookup-table interpolation and correction,
  green, sharpening, fade-out and stretch. Block averaging, value lookup, and full-disk row/column numbers and lon/lat also use numba kernels (_resample_nb.py).
  * The numpy version is limited by memory bandwidth and the GIL and is only about 2 times faster on 8 cores; the numba version speeds up almost linearly with the number of cores.
  * Differences from the numpy version come only from the order of floating-point operations: at most 1\~2 levels, in only 0.001%\~0.08% of the pixels.
  * The first run has to compile (about 10\~20 s; the result is cached in `__pycache__`). Without numba, numpy is used automatically.
* **Parallel decompression of Himawari segments**: the 10 HSD segments of one band are decompressed (bz2) and calibrated at the same time.
* **Saving**:
  * matplotlib draws the coastlines and annotation only once, on a transparent canvas, which is then composited with the image (for both basemap and cartopy). The original `savefig(bbox_inches='tight')`
    drew twice, resampling the whole image each time.
  * PNGs are compressed on multiple cores (satima/core/pngio.py, chunked compression like pigz).
  * Compared with the original output there are only occasional differences of 1\~3 levels where lines overlap.
  * **With cartopy the image was originally drawn by matplotlib, which shifted pixels**: the coordinate transform of GeoAxes is non-affine, so matplotlib resamples the image once more when drawing it.
    Even with the canvas the same size as the image, 54% of the pixels of a GOES-19 Mercator image were replaced by the adjacent pixel above, below, left or right (an offset of about half a pixel, 0 or 1 pixel after rounding).
    * Color-level edges therefore looked broken up and less like stair-steps, but pixel positions were off by up to 1 pixel (about 2 km).
    * Now cartopy also goes through compositing: the image is unchanged pixel for pixel, the same as basemap. Lines differ from those drawn by matplotlib by at most 1 level (0.03\~0.06% of the pixels).
    * On a full disk the image outside the Earth is transparent anyway (background color); the Earth edge is decided by Satima-GEO's per-pixel visibility and is no longer clipped a few extra pixels by cartopy's disk polygon.
* **The image height is now the integer number of grid rows**: previously the image height was a floating-point number, so matplotlib scaled the image by about 0.4 pixel and one row was lost in nearest-neighbour sampling.
  Now the image, the coastlines and the grid are aligned one-to-one.
  * The floating-point error of inches x dpi (e.g. 1691.9999999999998) would be truncated by Agg to one row less, so 0.001 pixel is added to both the width and the height of the figure.

Himawari-9 true color, 500 m native-resolution full disk (22000x22000; only 4 of the 10 segments were on hand, about 40% of the full disk),
measured on an i7-9700 (8 cores) / 16 GB:

| | Before | Now |
|---|---|---|
| Reading and computing | 467.5 s | 41.3 s |
| Drawing and saving | 103.5 s | 9.2 s |
| Total | 571 s (9.5 min) | 50.6 s |
| Peak memory | 9.8 GB | 8.3 GB |

A complete full disk (10 segments) is about 2.5 times that: about 24 minutes before, about 2 minutes now.
Each additional thread uses several hundred MB more memory; if memory runs short, reduce `workers`.

## Adding a satellite

Add `<sensor>-<format>.py` in `sensors/geo/` (the file name must contain `-`; a .py without `-` is treated as a helper module and is not loaded),
write a class that inherits `GeoReader`, and end with the line `READER = ClassName`. Only the satellite-specific parts need to be written:

| Item | Content |
|---|---|
| `NAME` `SENSOR` `DESCRIPTION` | Format name (= file name), sensor name, description for `--list`. The data list shows, in order, satellite, `SENSOR`, format (the part of `NAME` after '-', e.g. HSD, GVAR; a scene can customize it with `Scene.format`, e.g. NSMC GVAR), time, area |
| `PATTERNS` | File search patterns (glob) |
| `CHANNELS` | `{id in the data: Channel(id, center wavelength µm, sub-satellite point resolution km, 'vis'/'nir'/'swir'/'wv'/'ir', display name, functional nickname)}`; the display name is used in menus, annotation and file names (e.g. `BAND-13`); for the functional nickname see the table below |
| `TRUECOLOR` | (optional) `{'TrueColor': dict(red=, green=, blue=, nir=, ir=, lo_blur=)}`; `green=None` when there is no green band (ABI synthetic green) |
| `scenes(folder)` | Group the files into `Scene`s (all bands of the same time slot and the same area) |
| `nav(scene, channel)` | Return the navigation object (see below) |
| `read(scene, channel, (r0, r1, c0, c1))` | Read a window and calibrate it: reflectance 0\~1 / brightness temperature °C, invalid values NaN |
| `satellite_position(scene)` | (sub-satellite point longitude, satellite altitude km), used by true color |
| `platform(scene)` | (optional) pyspectral's satellite name; defaults to `scene.satellite` |

### Band names

* Menus, image annotation and file names always use `BAND-01`\~`BAND-16` (each agency's channel numbers; MTG in order of wavelength).
  MTSAT has no official channel numbers, so JMA's `VIS`, `IR1`\~`IR4` are kept (IR1 10.8, IR2 12.0, IR3 6.75 water vapor, IR4 3.75 µm);
  COMS keeps KMA's `VIS`, `IR1`, `IR2`, `WV`, `SWIR`;
  GMS-1\~4 have `VIS`, `IR`, and GMS-5 / GOES-9 have `VIS`, `IR1`, `IR2`, `IR3` (water vapor) (JMA names).
  MTG's FDHSI and HRFI are two separate data entries, and the annotation and file name get a dataset tag appended: `BAND-14 HRFI` / file name `BAND-14_HRFI`,
  `TrueColor FDHSI` / file name `TrueColor_FDHSI`. HRFI bands keep the number of the same band (BAND-03, BAND-08, BAND-09, BAND-14).
* Note: MTG has an extra 0.9 µm band (BAND-05), so the numbers after it are shifted relative to the other three satellites; e.g. MTG's BAND-13 is 9.7 µm ozone and
  BAND-14 is the 10.5 µm window (BAND-13 is the window on the other three). The menus all show the wavelength; compare with the table below.
* Internally (file reading, pyspectral) the ids in the data are still used: B13, C13, IR105, ir_105_hr.
* Functional nicknames (used in menus, in English) follow NOAA's official ABI nicknames; other satellites are mapped by function. When adding a satellite please keep to the same set of names:

| Nickname | GOES ABI | Himawari AHI | GK-2A AMI | MTG FCI |
|---|---|---|---|---|
| Blue | BAND-01 0.47 | BAND-01 0.47 | BAND-01 0.47 | BAND-01 0.444 |
| Green | — | BAND-02 0.51 | BAND-02 0.51 | BAND-02 0.51 |
| Red | BAND-02 0.64 | BAND-03 0.64 | BAND-03 0.64 | BAND-03 0.64 (also in HRFI) |
| Veggie | BAND-03 0.86 | BAND-04 0.86 | BAND-04 0.86 | BAND-04 0.865 |
| Near-IR WV | — | — | — | BAND-05 0.914 |
| Cirrus | BAND-04 1.37 | — | BAND-05 1.38 | BAND-06 1.38 |
| Snow/Ice | BAND-05 1.6 | BAND-05 1.6 | BAND-06 1.61 | BAND-07 1.61 |
| Cloud Particle Size | BAND-06 2.2 | BAND-06 2.3 | — | BAND-08 2.25 (also in HRFI) |
| Shortwave Window | BAND-07 3.9 | BAND-07 3.9 | BAND-07 3.83 | BAND-09 3.8 (also in HRFI) |
| Upper-level WV | BAND-08 6.2 | BAND-08 6.2 | BAND-08 6.24 | BAND-10 6.3 |
| Mid-level WV | BAND-09 6.9 | BAND-09 6.9 | BAND-09 6.95 | — |
| Lower-level WV | BAND-10 7.3 | BAND-10 7.3 | BAND-10 7.34 | BAND-11 7.35 |
| Cloud-top Phase | BAND-11 8.4 | BAND-11 8.6 | BAND-11 8.59 | BAND-12 8.7 |
| Ozone | BAND-12 9.6 | BAND-12 9.6 | BAND-12 9.63 | BAND-13 9.66 |
| Clean LW Window | BAND-13 10.3 | BAND-13 10.4 | BAND-13 10.4 | BAND-14 10.5 (also in HRFI) |
| LW Window | BAND-14 11.2 | BAND-14 11.2 | BAND-14 11.21 | — |
| Dirty LW Window | BAND-15 12.3 | BAND-15 12.4 | BAND-15 12.36 | BAND-15 12.3 |
| CO2 LW | BAND-16 13.3 | BAND-16 13.3 | BAND-16 13.31 | BAND-16 13.3 |

Wavelengths (µm) are from each agency's official table: NOAA, JMA, KMA GK2A Fact Sheet, EUMETSAT.

### Navigation

The navigation is chosen according to the projection; there is no need to squeeze every projection into one piece of code:

* Geostationary projection: `GeosNav(...)`, or from CGMS parameters `GeosNav.from_cgms(shape, cfac, lfac, coff, loff, ...)`
* Other regular-grid projections (Lambert, Mercator, equal lat/lon…): `ProjGridNav(shape, proj parameters, x0, dx, y0, dy)`
* Old satellites that only have per-point lon/lat arrays: `LonLatNav(lon, lat)`
* If none of these applies, inherit `Nav` and implement `to_rowcol()`, `to_lonlat()` and `key()`
  (e.g. `VissrNav` in `sensors/geo/vissr_nav.py`, the spin-scan navigation of GMS)

### GMS (VISSR) navigation

* Follows MGIVSR of JMA MSC Technical Note No.23; each point is computed from the attitude and orbit predictions at its own observation time (no interpolation);
  the inverse solves for the scan angle directly (the original program's inverse uses a small-angle approximation with VMIS), and the observation time is corrected iteratively until convergence.
  The forward computation differs from the old Cython version (gms_nav.pyx) by 0.0000 km, and the line/pixel number of the sub-satellite point differs by 0.01 from the simplified coordinate table in the header.
* `[navigation] engine`: `numba` (default, no compilation needed), `cython` (compile once with `python sensors/geo/build_vissr_nav_cy.py`),
  `numpy`; the three give the same result (difference < 1e-10 pixel). A 2000x2000 output takes about 0.4\~0.5 s.
* Line numbers use the LCW line number (starting from 1, same as JMA); satpy adds 1 line too many. For GOES-9 IR the center pixel number must have the header's
  "Pixel difference of VISSR center" (1.5) subtracted, as verified with the Earth edge (without the subtraction the east-west offset is -1.54 pixels).
* `[navigation] earth_edge_correction` (on by default): the Earth-edge correction of MSC Technical Note No.26, which compares the actual east and west Earth edges in the LCW of each scan line
  with the predicted edges and solves by least squares for the north-south and east-west offsets and the atmosphere width. Comparison with clear-sky coastlines of Western Australia:
  off by 2\~5 km before correction, about 1.5\~3 km after. GMS-AREA has no per-line data, so this correction is not applied.
* Full header: `reader.header(scene, 'IR1')` returns all records (names and fields as in the old _VISSR_*.py), including the Line_Control_Word of every line.

### GMS-1 tape images (GMS-TAPE)

GMS-1 VISSR tape images from the FGGE period (1978~1979), one file per tape (the file name is the volume number: 6 characters starting with `V`, no extension,
e.g. `VS2303`, `V10585`); visible and IR are both in the same file.

* The content is simply the JMA VISSR format (same as GMS-VISSR), except that each 4-byte word is stored as one 36-bit word:
  integers are 36 bits; reals are a 27-bit mantissa + 9-bit exponent (double precision 63 + 9); text is 9-bit EBCDIC;
  IR is stored at 9 bits per pixel and visible at 6 bits, packed continuously. Each record is preceded by a 4-byte length, and the beginning and end are EBCDIC tape labels.
  The decoding is taken from the vissr.c / gms_area.c supplied with the tapes, and the field positions correspond one-to-one to JMA's "FORMAT OF VISSR ARCHIVE DATA".
* Calibration and navigation both use the tape's own parameters: the IR brightness temperature table, the albedo table of each visible sensor, the coordinate conversion parameters, and the attitude and orbit predictions;
  navigation and Earth-edge correction are the same code as GMS-VISSR (`sensors/geo/vissr_common.py`), and line numbers use the LCW line number.
* Only image scan lines (Line name = 1) are taken; blackbody, test, annotation and gray-scale records are not used.
* Missing visible lines: these tapes used the redundant sensor set, in which visible sensor 2 had failed (`VIS redundant 2ch` = 0 in the Mode Block),
  so 1 of every 4 scan lines is missing. `[visible] fill_missing_lines = true` (default) fills it with the mean albedo of the lines above and below;
  `false` leaves it as no data (horizontal lines in the background color on the image). What is filled in is the albedo after each line's own calibration, not an observation of sensor 2.
* Full header: `reader.header(scene, 'IR')`, names and fields as in GMS-VISSR; values keep the precision of the tape (integers int64, reals float64),
  plus `Labels` (the tape labels).
* Third parameter block: its position corresponds to ⑧ (Data segment 8 of the status table) and ⑨ β-angle samples in the JMA document, but the document says "Omitted" for the content of both,
  with no field definitions. On the tapes only ⑧ has content (S/DB frame settings such as operation mode, line-number range, sensor failure table and frame constants, mostly duplicating the Mode Block);
  neither calibration nor navigation needs it, so it is kept as raw words (`Data_segment_8`); ⑨ is all 0.

### GOES-8\~15 (GVAR) reading and navigation

* Band names `BAND-01`\~`BAND-06` (= GVAR channels 1\~6; GOES-8\~11 have 1\~5, GOES-12\~15 have 1, 2, 3, 4, 6).
* Error tolerance for raw GVAR (cases where the CRC or the Data Valid flag is wrong but the image is fine are very common, so neither is used as grounds for discarding data):
  the 3 copies of the block header are majority-voted byte by byte (recognizable even if the first copy is bad); when the word count is wrong, the next block header is searched for again starting from the previous block;
  NSMC files are read in fixed 26624-byte slots with the starting point detected automatically; fields of the Line Documentation such as channel and detector
  are corrected with the whole-file mode for "the n-th record within the block"; each record locates its line number from its own scan count RISCT, cross-checked against Block 0;
  the navigation parameters are the mode of the O&A parameters of all Block 0s.
* Navigation is the same as McIDAS nvxgvar (0.0000 km difference from the old version's le2ll); when IMC is off (e.g. GOES-9 at 155°E in 2003\~2005)
  it is computed with the scan time of each line (the old version and McIDAS use a single time, giving differences of up to about 18 km within one image).
  The center of each IR line is in the middle of its 4 visible lines (+1.5 lines; the old version and McIDAS do not add this, and a 10.7 µm coastline comparison confirms that +1.5 is correct).
* Calibration: visible = bias + gain x count (+ quadratic term), then multiplied by the albedo coefficient (the old version's formula was (count + bias) x gain, which made space too bright);
  IR follows the tables in conversion_table/GOES.
  The GVAR of GOES-14/15 has no visible coefficients, so the per-detector pre-launch coefficients published by NOAA OSPO are used (Weinreb & Han 2009; GOES-14 no longer borrows those of GOES-15).
  `[visible] post_launch_calibration` (on by default) additionally multiplies by NOAA's post-launch calibration (sensor degradation): GOES-11/12/13/15 use STAR's monthly correction values
  (conversion_table/GOES/GOES_Imager_Vis_OpCal.txt, effective from the second Tuesday of the following month), GOES-8/10 use OSPO's linear degradation formula,
  and GOES-14 a fixed x 1.0695 (from comparison with MODIS during the post-launch test in 2009/12, found only in a PLT weekly report summary); GOES-9 has no official value.
  The visible coefficients in GVAR Block 0 have always been the pre-launch values (unchanged for GOES-13 over 2013\~2017), so the correction is not applied twice.
  The IR tables (NOAA OSPO) differ by < 0.02 K from the latest coefficients published by OSPO (the STAR-corrected version for GOES-13 CH6, RevH + STAR for GOES-14/15);
  the tables only list 180\~340 K, so the warm end (3.9 µm can reach 343 K) is extended with T = a + b x Teff(n, R) fitted to the table itself, while the cold end (space) is not extended.
  Detectors follow the GVAR specification: odd LIDET = north (DET1 in the table), even = south (DET2); LSIDE 1023 = side 2.
  The NOAA coefficients for GOES-12 CH6 (13.3 µm) have no spectral response correction; the GSICS comparison with IASI shows a bias of about 3 K, which is currently not corrected.
* Coastline comparison (clear-sky desert coasts): visible within 1 km; the east-west registration between IR and visible differs by about 0.4 IR pixel (varies by satellite).

### SMS, GOES-1\~7 (pre-GVAR) reading and navigation

* One data entry = the AREA_VIS and AREA_IR of the same satellite and time; the instrument column of the menu shows VISSR (SMS, GOES-1\~3) or VAS (GOES-4\~7).
  Bands: `VIS`, `IR` (the IR of VISSR mode), and `BAND-01`\~`BAND-12` of VAS (according to the file's filter map, e.g. 7, 8, 10).
* Image sources: VISR (McIDAS brightness values, 1 byte per pixel), VAS (mode AA, before 1987/03) and AAA (GOES-7), 2 bytes with up to 3 bands interleaved;
  for each line the position of the band is found from the band table.
* Calibration: visible is a square-root-encoded 6-bit value x 4, reflectance = (B/252)^2 (the albedo table of the same-type GMS VISSR is exactly this relation, rms 0.002;
  SMS / GOES have no absolute calibration. The old version displayed B/255 directly, which is equivalent to gamma already applied; for the same look turn on `[visible] gamma`);
  for IR brightness values, 418 - B when B >= 176, otherwise 330 - B/2 (K); VAS mode AA radiance = (P - YSUBZ) x 2^(F + DELTAF - 15);
  AAA radiance = (IAB2 x P/32 - IAB1) / 2^(15 - IFAB); radiance -> brightness temperature T = (FK2/ln(FK1/R + 1) - TC1)/TC2
  (the old version wrote it as ((FK2 - TC1)/TC2)/ln(...), about 0.3 K too warm). FK1\~TC2 keep the old version's single set of VAS coefficients (shared by GOES-4\~7; there are no per-satellite values).
* Navigation is the same as McIDAS nvxgoes (VisAD GOESnav): the 8 visible lines of one revolution (spin) share the same time; Kepler orbit, spin axis, attitude;
  the inverse solves the rotation directly (McIDAS's inverse interpolates the time from the beta records and is not the exact inverse of its forward computation) and is the exact inverse of the forward computation (1e-7 line).
  It differs from the old version's compiled nav.pyd by 0\~0.5 km (McIDAS rounds the seconds of the epoch time, the old version truncated them); the old version's Python source nav.py
  used floating-point division for ietimh, NUMSEN and the scan count and differs from McIDAS by about 11 km on average (what the old reader actually loaded was nav.pyd).
* IR is sampled once per revolution; the AREA IR line number = the first visible line of that revolution, and the pixel center is taken at +3.5 lines (the middle of the 8 visible lines).
  Comparing visible / IR cloud edges of the same time slot in 27 regions gives +3.25\~+8 lines in every case; east-west it is +7 pixels in 1988 and -3 pixels in 1994, which is the registration error of the time and is not corrected.
* The navigation of this era itself has errors of several km (visible coastlines: north-south +2\~+6 km, east-west -7\~+2 km, varying from image to image); the McIDAS formulas are followed with no further correction.
* Files whose navigation block has no orbit parameters (e.g. GOES-1 and GOES-3 in 1979, which McIDAS cannot navigate either) are skipped after the reason is listed.
* Full header: `reader.header(scene, 'B08')` returns the AREA directory, navigation block, AAA calibration block, per-line prefix data (_Data_Block) and the band table.

### MSG SEVIRI (Meteosat-8\~11)

* Bands `BAND-01`\~`BAND-12` (= EUMETSAT channels 1\~12: VIS006 … IR_134, HRV; in the menu BAND-12 is labelled "HRV (High Resolution
  Visible)", while the annotation and file name remain BAND-12); the data list shows FD / RSS / ROI and IODC.
* The three formats give the same result (values differ from satpy by < 1e-4, lon/lat by 0.0000 km; HRIT was validated with satpy's public Zenodo sample
  MSG4 2018/02/28 15:00, native with a Meteosat-9 IODC full disk and Meteosat-11 RSS).
* Image row 0 is the southernmost and column 0 the easternmost; data before 2017/12 (Earth model 1) have the 1.5 km offset corrected as described by EUMETSAT.
* The HRV (1 km) full disk is split into a lower and an upper window, each located by its own east boundary (the data are not padded to 11136 columns wide).
  The upper window actually starts at the line LowerNorthLineActual of the epilogue (one line earlier than UpperSouthLineActual;
  verified by image continuity, same as satpy). The IODC HRV has a slanted no-data area where the two windows meet; the original data are like that.
* Calibration: `[seviri] calibration` = nominal (default) / gsics (IR) / meirink (visible); reflectance includes the Sun-Earth distance correction,
  and brightness temperature uses EUMETSAT's α, β or the quadratic formula depending on the radiance type (effective / spectral). Pixels with zero radiance are NaN (satpy would give -276 °C).
* Compressed HRIT (-C_) needs pyPublicDecompWT, or xRITDecompress specified with `[seviri] xrit_decompress`; decompression is done only in memory / a temporary directory.
* zip of native: listing only decompresses the header; the .nat is decompressed into memory when the image is first read (about 270 MB; only one is kept).
* The np.bool in the old version's prologue.py was removed in numpy 1.24 and later, so the old MSG_HRIT fails when reading the prologue; fixed in the new version.

### MFG MVIRI (Meteosat-2\~7)

* Bands `VIS` (5000x5000), `IR` 11.5 µm, `WV` 6.4 µm (2500x2500).
* FCDR: navigation is the standard geostationary projection with an 18-degree field of view (i.e. the old version's regfeo), with the projection longitude taken from the file name (E5700 = 57.0);
  it is the same as satpy (lon/lat 0.0000 km, brightness temperature difference < 1e-4), and the image orientation was additionally verified with the satellite zenith angle in the file (0.2 degree for south-to-north, 8 degrees for the reverse).
  The VIS FCDR value is the bidirectional reflectance factor (already divided by cos of the solar zenith angle); it is multiplied back by cos(solar zenith angle) to give a TOA reflectance consistent with the other satellites.
  In some files the time variable is one day later than the file name; the file name is taken as correct.
* MVIRI-NC (the format of the old version's MFG.py) is ported from the old version; geolocation uses the lon/lat in the file and brightness temperature uses EUMETSAT's A, B coefficients; no sample yet to validate with.

## Validation (2026/11)

* New readers compared with the old version's readers (with the 1-based fix applied): AHI, AMI, FCI navigation 0.0000 pixel, values identical;
  ABI navigation differs from satpy by < 2 m (the old GOES16_NC had 0.01\~0.04 pixel error because it computed in float32).
* For true color, the block correspondence of each band on the output grid agrees 99.98\~100%, and the jointly fitted shift is 0.
* GMS-VISSR / GMS-AREA: for the same time slot (GMS-1 1979/10/11 02:33) the brightness temperatures of the two formats differ by < 0.26 K;
  redrawing IR197901071533.png with GMS-AREA gives the same result.
* GMS-TAPE (four tapes: VS2283, VS2288, VS2303, V10585):
  * Compared with the same tape (VS2303) previously converted to JMA format and then read with GMS-VISSR: visible pixel-identical, IR difference < 0.00003 K,
    lon/lat difference < 0.00007° (rounding from reals being reduced to 32 bits in the conversion).
  * Compared with the McIDAS AREA of the same time slot (GMS_VIS_1_1979_282_2333.area, VS2288): on observed lines the AREA brightness = count x 4 (all identical),
    and on missing lines the AREA is also the mean of the lines above and below; the AREA's GMSX navigation parameters (center line 5064.5, center pixel 6688.5, starting line number) are the same as the tape's,
    and the lon/lat at the same line/pixel number differ by 0.55 km on average (GMSX stores the stepping angle with only 4 significant digits).
  * Earth-edge correction amounts: north-south -0.2~+1.3 lines, east-west +0.2~+1.6 pixels (IR); the comparison with Australian coastlines agrees.
  * After the shared code was factored out of GMS-VISSR, the headers, navigation and image values of 11 bands in total from GMS-1, GMS-3, GMS-5 and GOES-9 are bit-identical to those before the change.
* Not yet validated for lack of samples: GK-2A regional data in Lambert projection (lc).

# Satima: Polar-Orbiting Imagers (PO)

> This is the detailed documentation for **PO (polar-orbiting imagers)**, rewritten from the original Satima-PO. For how to run the integrated version, the data listing, the layered config files
> (`Satima.toml` + `Satima-PO.toml`), the interface language, the image modes (native / dpi) and the continuity rules for polar orbiters, see [README.en.md](../../README.en.md).
> PO is the family that changed most of the three: the differences from the original Satima-PO are listed under "Differences from the original three programs" in the README.

Single channels (reflectance, brightness temperature with colormaps, DNB) and true color from VIIRS, MODIS, AVHRR, METimage, OLCI, SLSTR and SGLI.
The 8 scripts, originally independent of each other, were first integrated into Satima-PO and then merged into Satima; the calibration and algorithms of each format remain in its own script.

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
* The title (formerly Picname) is asked by default (`[title] ask = true` in `Satima-PO.toml`) and may be left empty; it is put in the file name and drawn at the lower left corner of the image.

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

VN07 / VN10 of SGLI are high-gain (ocean color) channels that saturate over cloud at a reflectance of about 0.16 / 0.13, and VN08 saturates at about 0.49; to plot clouds use VN11 (868.5 nm).

Start and end times of each format: for VIIRS, METimage, AVHRR-NAT, OLCI and SLSTR they are taken from the file name; MODIS file names only have the start time, so the nominal 5 minutes is assumed;
AVHRR-NC and SGLI file names have no end time, so the listing shows only the start time (files within 10 minutes of each other are still joined into one segment). The two AVHRR formats read only one file at a time.

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

## True color correction check

| | Solar zenith angle correction | Rayleigh (CREFL) | O3 | H2O | Terrain height | Sharpening | Enhancement |
|---|---|---|---|---|---|---|---|
| VIIRS | already included in SDR | ✔ | climatology | ✔ | ✔ GMTCO Height | I1 ratio sharpening | Simple(sRGB) / WorldView / CIRA |
| MODIS | ✔ (tapering off from 88°) | ✔ | climatology | ✔ | ✔ MOD03 Height | B1 ratio sharpening | CIRA / WorldView |
| METimage | ✔ | ✔ (ISRF-weighted τ) | climatology | ✔ | — sea level | — | CIRA / WorldView |
| OLCI | ✔ | ✔ (including sea-level pressure) | ECMWF | — | ✔ altitude | unsharp mask | CIRA (white-point normalization) / WorldView |
| SGLI | ✔ | ✔ (τR from the JAXA ATBD) | climatology 300 DU (JAXA Koz) | — | — sea level | — (all 250 m) | CIRA / WorldView |

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

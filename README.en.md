# Satima — Weather Satellite Imaging

The integrated version of Satima-GEO (geostationary satellites), Satima-MW (microwave) and Satima-PO (polar-orbiting imagers): one main program, one user interface,
layered config files. 中文版: [README.md](README.md)

| Family | Data | Readers | Details |
|---|---|---|---|
| **GEO** geostationary satellites | Himawari, GOES, GK-2A, MTG, MTSAT, COMS, GMS, Meteosat… | `sensors/geo/` (18 formats) | [docs/en/GEO.md](docs/en/GEO.md) |
| **MW** microwave radiometers | AMSR2/3, GMI, SSMIS, SSM/I, ATMS, AMSU-A, MWRI (FY-3A\~H), NASA PPS L1C… | `sensors/mw/` (13 formats) | [docs/en/MW.md](docs/en/MW.md) |
| **PO** polar-orbiting imagers | VIIRS, MODIS, AVHRR, METimage, OLCI, SLSTR, SGLI, FY-1 MVISR, FY-3 MERSI / VIRR | `sensors/po/` (12 formats) | [docs/en/PO.md](docs/en/PO.md) |

## Running

```
python Satima.py                  search for data in the config file's data_dir (default: data next to Satima.py)
python Satima.py -d D:\data       set the data directory
python Satima.py -o D:\png        set the image output directory
python Satima.py -t geo           list one family only: geo / mw / po
python Satima.py -f AHI-HSD       search one format only
python Satima.py -c my.toml       set the common config file (family files are looked up in the same folder)
python Satima.py --lang en        interface language zh / en
python Satima.py --mode dpi       image mode native (native resolution) / dpi
python Satima.py --list           list the available formats
```

The workflow is the same for all three families: **select data → projection → map range → products → colormap / style → image width → title → plot and save**.
Only the menu contents differ by family (the full disk projection exists only for GEO, PCT and RGB only for MW…).

### Data listing

Data of all three families can be kept in the same folder. The listing is split into sections by family with continuous numbering; when only one family has data, the section headers are not printed:

```
[GEO]
  1   MTSAT-1R     IMAGER   HRIT   2013/11/07 11:30Z   FD DK01   5 bands
  2   Himawari-9   AHI      HSD    2026/10/20 21:00Z   FLDK      1 bands
[MW]
  3   F16       SSMIS   SDR   2017/09/06 07:25-09:20Z   CLASS   19 ch
  4   GCOM-W1   AMSR2   L1B   2025/09/22 04:34-05:24Z   JAXA    14 ch
  5   NOAA-20   ATMS    SDR   2025/09/28 05:50-06:01Z   CLASS   22 ch   21 files
[PO]
  6   S-NPP     VIIRS    HDF5   2023/05/26 04:03-04:10Z   NOAA       2 bands   4 files (1 missing)
  7   METOP-B   AVHRR3   NAT    2026/04/12 09:31-11:13Z   EUMETSAT   IR/VIS
Input the index (several: (3,5) or [1,3,5]):
```

* Columns: satellite, instrument, format, time, observation area (GEO) or data source (MW, PO), number of bands / channels, **notes**.
* Notes column: when several files are joined into one entry it shows the number of files, and adds `(N missing)` only when files are missing in between; for GEO regional observations it shows the latitude / longitude range of the image.
* `3` = a single entry; `(3,5)` = entries 3 through 5; `[1,3,5]` = the listed entries. Only one format can be selected at a time; if you select different formats you are asked to select again.
* GEO: selecting several entries = one image for each time; for MW and PO, selected entries that are contiguous in time can be joined into one image (see "Continuity rules for polar orbiters").

### Map range

`[map] range_input` decides which one is asked first, but both input methods can be used at any time:

* Input **4 numbers** = south north west east, e.g. `20 27 117 124` (across the date line use 0\~360, e.g. `10 30 170 200`)
* Input **2 numbers** = center latitude longitude, then the range is asked (degrees of latitude and longitude outward from the center; one number = both the same)

For polar orbiters (PO) the center time and the eye temperature are based on the center point: when the four edges are input it is the midpoint of the range; to specify the typhoon center, input 2 numbers.
GEO regional observations (target area, mesoscale) are not asked and use the whole image; full disk (FD) needs no range.

Projections: MERC (Mercator), CYL (equidistant cylindrical), LCC (Lambert conformal conic); GEO also has FD (full disk).
Mercator and LCC always use the WGS84 ellipsoid as the Earth model, so **GEO, MW and PO images of the same range and the same size are aligned pixel for pixel** and can be overlaid directly.

## Config files

| File | Contents |
|---|---|
| `Satima.toml` | Common settings: language, directories, map, coastline, image mode and size, title, annotation, file name |
| `Satima-GEO.toml` | Geostationary satellites: true color, GeoColor, base map, infrared colormaps… |
| `Satima-MW.toml` | Microwave: interpolation, PCT, RGB recipes, colormaps… |
| `Satima-PO.toml` | Polar-orbiting imagers: eye temperature, colormaps… |

* Priority: built-in defaults < `Satima.toml` < `Satima-<family>.toml` < `[readers.<format>.<section>]` < command line.
* **Each setting only needs to be written in one place.** A family file can override any section of the common file (write only the keys you want to change); when both have the same key, the family file wins.
  For example, the coastline in the common file is white 0.2; writing `[coastline] color = "#00BB00"` in `Satima-PO.toml` makes only the polar orbiters use green.
* A line such as `[readers.AHI-HSD.general] data_dir = "D:/data/himawari"` overrides further for a single format; it can be written in the common file or in a family file.
* Directories may be absolute or relative paths; **relative paths are always relative to the folder containing `Satima.py`**, regardless of where the program is run from
  (relative paths given to `-d` / `-o` on the command line follow the usual convention and are relative to the current directory).
* Every key is explained in the comments of the config files. Versions with English comments are in `config-en/` (same keys and values); **the program only reads the four files next to `Satima.py`**,
  and `config-en/` is never read; to use the versions with English comments, copy them out over the originals. `python tests/check_config_sync.py` checks whether the keys and values of the two sets agree.

### Language

`[general] language = "zh"` (Traditional Chinese) or `"en"` (English); `--lang` on the command line overrides it temporarily; only one language is shown at a time.
It only affects the messages on screen; product and band names, colormap names, file names and the annotations on the image are the same in both languages, so the files produced in the two languages are exactly the same.

### Image mode (`[image] mode`)

Applies to all three families; `--mode` on the command line overrides it temporarily, and it can also be set individually in a family file or in `[readers.<format>.image]`.

| | `mode = "native"` (native resolution) | `mode = "dpi"` |
|---|---|---|
| Image size | One pixel at the map center = the native resolution of this product (for polar orbiters and microwave, the native sampling interval) | Number of pixels = inches x dpi |
| Number of images | One | One per `[[output.image]]` (default 1000 dpi and 500 dpi) |
| At run time | Asks for the image width: Enter = native resolution, `x3` = 3 times, or input the number of pixels directly | Not asked |
| Related settings | `[image] size / max_size / min_size / min_scale / suffix` | `[[output.image]] dpi / suffix / annotation` |

* The size in inches in dpi mode follows the matplotlib defaults (a 6.4 x 4.8 inch figure, axes about 4.96 x 3.70 inches): the map is fitted into the axes area keeping its aspect ratio, and the white margins are trimmed.
  A wide range at 1000 dpi is about 4960 pixels wide, a square range about 3700 pixels high. The defaults are read from the matplotlib settings, so they change if you edit `matplotlibrc`.
* Font sizes and line widths are in points: native-resolution mode treats the image as 5 inches wide, dpi mode uses the actual size in inches, so whatever the number of output pixels,
  text and lines take up the same proportion of the image.
* **In both modes the image is first sampled onto the grid of the output map and is not resampled by matplotlib**; pcolormesh is no longer used either
  (the original dpi mode of PO was pcolormesh). GEO samples by analytic navigation, MW by swath interpolation (`[swath] interpolation`),
  and PO takes the nearest valid sample. To see "one block per sample", set `[image] interpolation = "nearest"` for GEO and
  `[swath] interpolation = "nearest"` for MW.
* When dpi mode outputs several images they are told apart by `{suffix}` or `{dpi}` in the file name; when the template has neither field, images that have a suffix automatically get `.suffix` appended to the file name.
* GEO and MW sample once for each dpi, so outputting two images takes about twice as long as one.

## Continuity rules for polar orbiters

Polar-orbiter files vary greatly in length: an ATMS file is 32 seconds, METimage about 1 minute, VIIRS about 85 seconds, MODIS 5 minutes, JAXA L1B half an orbit,
SSMIS a whole orbit. The program treats "adjacent" and "same segment" as two separate things (`satima/core/segments.py`):

* **Adjacent** (used to count missing files): start of the next file − end of the previous file ≤ the adjacency threshold. The threshold follows the granule length: half the median length of the files in the group,
  limited to 2\~120 seconds (ATMS 16 seconds, VIIRS about 43 seconds, files of 5 minutes or longer 120 seconds). The gap left by one missing file equals the length of one file, so it always exceeds the threshold.
* **Same segment** (one row of the data listing): files of the same satellite, same instrument and same format whose gaps do not exceed 10 minutes are joined into one segment; they are joined even when files are missing in between,
  and the notes column of the listing shows `(N missing)`. Two overpasses of the same place are at least one orbit apart (about 100 minutes), so within 10 minutes it cannot be a different overpass;
  the gap is left blank when sampling and is not smoothed over by interpolation.
* The orbit number is not used to split segments: it increments at the ascending-node equator crossing, which would cut the overpass of a low-latitude typhoon exactly in two.
* **Short granules** (a file shorter than 15 minutes: VIIRS, ATMS, METimage, MODIS…) are joined into one entry automatically, without asking.
  **Long files of a whole / half orbit** each get their own row; when several rows that are contiguous in time are selected together, `[swath] merge` decides whether you are asked to join them into one image (`"ask"` / `true` / `false`).
* When a granule has duplicate files (downloaded again, different creation times), only the newest is used.
* VIIRS: once the range is known, only the granules that intersect the map are read, plus one extra before and after.
* The time on the image and in the file name is the time the swath **passes over the map range**, not the start time of the file.

Current coverage: ATMS (MW) and all PO readers use these rules. The other MW readers (one file = a whole or half orbit) keep the original
`[swath] merge` / `merge_gap`: when several files that are contiguous in time are selected, you are asked whether to merge them.

### VIIRS and ATMS files

* **NOAA CLASS packaged files** (geolocation and data packed in one file, with the contained products listed at the start of the file name, e.g. `GMTCO-SVM03-SVM04-SVM05_j01_…`,
  `GATMO-SATMS-TATMS_j01_…`) can be read, and separate files and packaged files can be mixed; an aggregated file that combines several granules into one file is treated as one longer granule according to the start and end times in its file name.
* **VIIRS geolocation files**: the terrain-corrected `GMTCO` / `GITCO` are preferred; when they are not found, `GMODO` / `GIMGO` without terrain correction are used instead, and this is noted on screen.
  Over the sea the two are the same; high mountains are displaced at the swath edge (by about terrain height × tan(view angle); a 3 km mountain at a 60° view angle is off by about 5 km).
* **VIIRS DNB**: when `GDNBO` has terrain-corrected latitude / longitude (`Latitude_TC` / `Longitude_TC`), they are preferred.
* **ATMS geolocation files**: when `GATMO` has the per-band beam latitude / longitude (`BeamLatitude`), they are used; old-format files without them fall back to the single set `Latitude` / `Longitude`.

> **VIIRS: verified with real files** (I1, I5, DNB, M3 / M4 / M5 of S-NPP, NOAA-20 and NOAA-21, 2026/10): for the I5 colormap, I1 reflectance,
> DNB (ERF), M5, TrueColor-M and TrueColor-I, the center time and eye temperature are numerically the same as the original Satima-PO and the images are visually consistent; automatic joining of granules, reading only the intersecting granules,
> selecting one segment each of two satellites at once, the notice when a geolocation file is missing, and the DNB `Latitude_TC` were all run with real files. Packaged files were tested with files assembled from real data
> (`tests/pack_viirs.py`); the output is pixel-identical to that from the separate files.
>
> **Fengyun-3: verified with real files** (2026/11): single channels and true color of MERSI-II (FY-3D), MERSI-LL (FY-3E), MERSI-III (FY-3F, FY-3H)
> and MERSI-RM (FY-3G), the MERSI-LL low-light band; MWRI (FY-3B, FY-3C), MWRI-II (FY-3F, FY-3H), MWRI-RM (FY-3G) and the FY-3D MWRI orbital L1c.
> Reflectances are pixel-identical to satpy; brightness temperatures follow the NSMC documents (equivalent central wavenumber and Tbb = A Te + B)
> and differ from satpy (which uses the nominal wavelength) by 0.1\~1.4 K. FY-3A/B/C MERSI-1 and VIRR have no samples: they were written from the
> format documents and only run with synthetic files assembled from FY-3D data (`tests/fake_fy3.py`); FY-3A MWRI has no samples either.
>
> **Parts not yet verified with real files**: original packaged and aggregated files downloaded from CLASS, VIIRS `GMODO` / `GIMGO`, M-band brightness temperatures (M12\~M16),
> ATMS packaged files and old-format geolocation files (the VIIRS items among these were only run with the synthetic test files of `tests/fake_viirs.py` to confirm that the workflow runs through);
> MODIS, METimage, OLCI, SLSTR, SGLI and AVHRR-NC have no samples, and it was only confirmed that the code loads.
> Take care the first time you run these with real data; if there is a problem, keep the error message and the file names.

## Differences from the original three programs

**Common to all three families**
* One main program `Satima.py`; config files split into a common file + three family files; messages switchable between Chinese / English.
* `[image] mode` (native / dpi) is available for all three families; the size rules of dpi mode are the same for all three.
* The data listing is shown in sections, with an added notes column.
* The default file name template gains `{suffix}` (omitted when empty; in native-resolution mode the file names are the same as before).

**GEO, MW**: in native-resolution mode the output is pixel-identical to the original programs (see "Verification"). GEO gains the center + range input method.

**PO** (the most changes)
* pcolormesh removed: both modes sample onto the map grid (the approach of the original "native-resolution mode").
* Mercator changed from a sphere to the WGS84 ellipsoid: the image height for the same range differs slightly (e.g. 6918 rows become 6896 rows), so images are not pixel-identical to earlier ones.
* CYL and LCC projections added; the range can be input as four edges; messages are available in Chinese.
* Relative paths are now relative to the folder containing `Satima.py` (formerly the current directory at run time), and images are saved to `plot/` by default (formerly next to the data).
* The title (Picname), besides being put in the file name, can also be drawn at the lower left corner of the image (`[title] draw`; the common default is true, the PO config file defaults to false, i.e. file name only); the annotation below the image now uses the `[annotation] text` template.
* The listing changed from "one row per file" to "one row per continuous segment"; VIIRS is no longer one row per band: the channels or true color are chosen from a menu after the data is selected.
* Config keys: `[output] facecolor` → `[image] background`, `[output] ask_title` / `title` → `[title] ask` / `default`,
  `[coastline] resolution` → `[map] resolution` (old keys give a message asking you to rewrite them).

## Directory layout

```
Satima.py                 main program
Satima.toml               common config file
Satima-GEO.toml  Satima-MW.toml  Satima-PO.toml     config file of each family
config-en/                config file templates with English comments (not read by the program)
satima/
  core/                   shared by the three families
    session.py            command line, reader loading, sectioned data listing, data selection
    config.py             reading of the layered config files, common defaults and checks
    i18n.py               interface language tr(Chinese, English)
    ui.py                 interactive input, menus, listings
    maps.py               map projections (MERC / CYL / LCC / FD) and map backends (basemap / cartopy)
    canvas.py             image size (native / dpi), coastline, annotation, saving
    naming.py             output file name template
    colormaps.py          colormaps
    segments.py           continuity rules for polar orbiters
    pngio.py  parallel.py multi-core PNG output
  geo/                    geostationary satellites: reader, nav, resample, products, truecolor/, landsea, bluemarble
  mw/                     microwave: reader, swath (interpolation), grid, products
  po/                     polar-orbiting imagers: reader, canvas (nearest-neighbour sampling), swath, truecolor, obstime, tiepoint, archive
  <family>/family.py      workflow of the family (data search, rows of the listing, range, products, plotting)
  <family>/config.py      the family's own defaults and checks
sensors/geo/  sensors/mw/  sensors/po/    readers of each format (.py files with '-' in the name are loaded automatically)
colormap/                 colormaps (.cmap)
conversion_table/         calibration tables
bluemarble/ landsea/ nightlights/   base map data for GEO
data/  plot/              default data and output directories
docs/                     detailed documentation of each family
tests/                    regression tests and checking tools
```

Required packages are listed in `requirements.txt`. Python 3.11 or later (the built-in tomllib is used). When a package is missing, only the formats that use it are skipped.

### Adding a reader

Put a `<instrument>-<format>.py` in `sensors/<family>/`, with the last line `READER = <reader class>`. Base class and how to write one, for each family:

| Family | Base class | Details |
|---|---|---|
| GEO | `satima.geo.reader.GeoReader` | [docs/en/GEO.md](docs/en/GEO.md) "Adding a satellite" |
| MW | `satima.mw.reader.SwathReader` | [docs/en/MW.md](docs/en/MW.md) "Adding a format" |
| PO | `satima.po.SatReader` | [docs/en/PO.md](docs/en/PO.md) "Adding a satellite" |

Write messages as `tr('中文', 'English')` (`from satima.core.i18n import tr`); for the reader description use `DESCRIPTION` (Chinese) and `DESCRIPTION_EN`.

## Verification (2026/11)

`tests/regress.py` runs the original three programs and the integrated version with the same input and compares the output pixel by pixel:

```
python tests/regress.py run orig      run the originals (needs the original three projects and the sample data; the paths are at the top of the script)
python tests/regress.py run new       run the integrated version
python tests/regress.py cmp orig new  compare
```

| Case | Contents | Result |
|---|---|---|
| geo_ahi_merc | Himawari-9 AHI Mercator: BAND-13 colormap, TrueColor | 2 images pixel-identical to the original |
| geo_ahi_fd | AHI full disk: BAND-13 grayscale, CWA land-sea base map | 2 images identical |
| geo_mtsat_lcc | MTSAT-1R LCC: VIS, IR1 colormap | 2 images identical |
| mw_amsr2 | AMSR2 L1B: 89V colormap, PCT89, NRL-37 | 3 images identical |
| mw_atms | ATMS SDR, 21 granules joined into one entry automatically: 88.2 GHz | 1 image identical |
| mw_ssmis_l1c | PPS L1C (TMI): 10.65V | 1 image identical |
| po_avhrr | MetOp-B AVHRR EPS native: CH4 colormap | Image produced; because of the switch to WGS84 and the shared canvas it is not pixel-identical to the original (visually consistent) |
| po_viirs_i5 | NOAA-20 VIIRS I5 colormap, 3 granules joined into one entry automatically | Center time and eye temperature the same as the original; image 2897 x 2985 (original 2894 x 3000), visually consistent |
| po_viirs_dnb | NOAA-21 VIIRS DNB (ERF), only the 5 intersecting granules of 7 are read | Center time the same as the original, visually consistent (the new version uses `Latitude_TC`) |
| po_viirs_tc | NOAA-21 VIIRS TrueColor-I (M3 / M4 / M5 + I1) | Center time the same as the original, visually consistent |

dpi mode: for each of the three families one sample was output as two images at 1000 / 500 dpi. The 1000 dpi AVHRR image is 2830 x 3696 pixels; the pcolormesh of the original PO
gave 2818 x 3696 (same height; the 12-pixel difference in width comes from the change from sphere to ellipsoid).

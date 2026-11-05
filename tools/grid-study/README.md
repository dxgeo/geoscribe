# grid-study

The measurements behind
[`plans/transformation-data.org`](../../plans/transformation-data.org)
§ MEASURED (2026-11-05): how far PROJ's datum-shift grids shrink, what
spherical-harmonic geoid models cost, and whether NADCON5 can be rebuilt from
NGS's own control points. These are study scripts, not part of the crate's
build; the numbers they printed are recorded in the plan, and
`inventory.csv` is the per-file table behind its inventory totals.

## Needs

- PROJ-data: the grid GeoTIFFs and `proj.db`, with `PROJ_DATA` pointing at
  them (`projinfo --searchpaths` lists where PROJ looks).
- GDAL (`gdal_translate`, `gdalinfo`), `sqlite3`, `zstd`, `xz`.
- GMT, for the NADCON5 rebuild only (`gmt surface`).
- Python 3 with numpy.
- Network, for the two `fetch_*` scripts only: ICGEM (GFZ) for the EGM
  coefficients, NGS for the NADCON5 build.

Everything works in a *work directory* you name, outside this repository:
exported grids, downloads and intermediate files go there, the scripts stay
here. Downloaded files are data from third parties, under their own terms
(plan D-T7); keep them out of the repository. Run Python with `-I`.

```sh
export PROJ_DATA=/path/to/proj
W=/some/scratch/dir
sh setup.sh $W            # grids.txt and the sample grids, band-sequential
sh fetch_egm.sh $W        # EGM96 and EGM2008 coefficients (~74 MB), parsed
sh fetch_nadcon5.sh $W    # NADCON5's 2007 -> 2011 CONUS build files (~10 MB)
sh nadcon5_surface.sh $W  # NGS's gridding command, rerun with GMT
```

`setup.sh` forces band-sequential output on purpose. GDAL keeps a source's
interleave, and NADCON5's GeoTIFFs are pixel-interleaved: read as
band-sequential they scramble, and look 3-40x less compressible than they are.
The first measurements made that mistake.

## Which script gave which result

| result in § MEASURED | script | run as |
|---|---|---|
| the whole PROJ dataset (312 files), `inventory.csv` | `inventory.py` | `python3 -I inventory.py $W` |
| per-file ratios at two steps | `gridpack.py` | `QH=1e-5 QV=1e-4 DROP_ACC=1 python3 -I gridpack.py $W/grids GRID...` |
| coarse-to-fine (SZ3-style) predictor | `interp.py` | `python3 -I interp.py $W/grids GRID...` |
| lossless at native step; entropy-coder bounds | `entropy.py` | `python3 -I entropy.py $W/grids 1e-6 1e-5 GRID...` |
| native steps, roughness, upsampling | `probe.py` | `python3 -I probe.py $W/grids GRID...` |
| metric lattice, filler, where the bits go | `structure.py` | `python3 -I structure.py $W/grids GRID...` |
| national geoids as residuals from EGM2008's grid | `xgrid.py` | `python3 -I xgrid.py $W/grids GRID...` |
| EGM96 as coefficients, against PROJ's grid | `egm96_test.py` | `python3 -I egm96_test.py $W` |
| EGM2008 against its grid; uniform step | `egm2008_test.py` | `python3 -I egm2008_test.py $W` (minutes) |
| why EGM2008's step must shrink with degree | `egm_errcheck.py` | `python3 -I egm_errcheck.py $W egm2008_test.py` |
| EGM2008 with a degree-dependent step | `egm2008_step.py` | `python3 -I egm2008_step.py $W` |
| NADCON5: grid against its source points; their size | `nadcon5_points.py` | `python3 -I nadcon5_points.py $W` |
| NADCON5: GMT's rebuild against PROJ's grid | `cmp_surface.py` | `python3 -I cmp_surface.py $W` |
| NADCON5: the same, by distance to the nearest point | `nadcon5_distance.py` | `python3 -I nadcon5_distance.py $W` |
| compression, predictors, interpolation; licences | `survey_tiffs.sh` | `sh survey_tiffs.sh` |

`egm_parse.py` is called by `fetch_egm.sh`. Sizes depend a little on the
`zstd` version; the plan's figures came from zstd at `-19` on 2026-11-05,
with PROJ 9.9.0's data and GMT 6.7.0.

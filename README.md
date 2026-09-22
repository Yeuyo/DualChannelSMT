# DualChannelSMT

MATLAB code for two-channel single-molecule tracking (SMT) analysis. Single molecules are localised and tracked in one channel, and each localisation is then scored against a second, co-registered channel that is segmented into intensity tiers. The output describes how single-molecule binding relates to local concentration in the second channel.

The entry point is **`concentrationCompare_SlidingScan.m`**. Everything else in the repository is a supporting library.

> **Status:** research code. Paths, intensity thresholds and several constants are hard-coded for one dataset and must be changed before use on new data. See [Known limitations](#known-limitations).

---

## What the pipeline does

For each matched pair of channel-1 and channel-2 image stacks:

1. **Localisation and tracking (channel 1).** The stack is read with `tiffread`, localised with `localizeParticles_ASH`, and linked into trajectories with `buildTracks2_ASH` (MTT algorithm). Trajectories are flattened into a table with columns `[time, x, y, frame, trajectoryID, intensityTier, isLongBinding]`.

2. **Spatial clustering of localisations.** A greedy `pdist2` loop repeatedly picks the localisation with the most neighbours within 4 px, assigns that neighbourhood a cluster ID, and removes it from the pool. For each cluster, occupancy is computed as the number of distinct frames in which the cluster is occupied, divided by `totalFrames`. Written to `*_Occupancy.xls`.

3. **Intensity tiering (channel 2).** Each channel-2 frame is quantised into 21 levels with `multithresh(frame, 20)` / `imquantize`. Levels are then binned into three tiers using the indices `lI`, `mI`, `hI` (default 8 / 13 / 16), interpreted as low / medium / high local concentration. Every channel-1 localisation is assigned the tier of the channel-2 pixel it sits on, in the same frame.

4. **Long-binding filter.** A trajectory is called "long binding" if it persists for at least `tLong` seconds (i.e. `tLong / (ExposureTime/1000)` frames). The fraction of long-binding localisations falling in each tier is computed and normalised to sum to 1 (`realL`, `realM`, `realH`).

5. **Per-site temporal analysis.** Pixels that are in the high tier for at least `fL` frames define candidate sites. Sites at least 4 px inside these regions (via `bwdist`) are taken as seeds. For each seed, the script counts, frame by frame, how many channel-1 localisations fall within `distTol` pixels, producing an occupancy trace, a per-site trajectory count and a per-site frame count. Traces are saved as PNGs; counts go to `*_Summary.xls`.

---

## Repository layout

```
concentrationCompare_SlidingScan.m   Main analysis script — run this
tiffread.m                           TIFF stack reader (v2.4, F. Nédélec, EMBL)
Batch_MTT_code/
  localizeParticles_ASH.m            Localisation on an in-memory 3-D stack
  buildTracks2_ASH.m                 Trajectory building on an in-memory stack
  detect_et_estime_part_1vue*.m      MTT detection / Gaussian estimation core
  carte_H0H1_1vue.m                  GLRT detection map
  estim_param_part_GN.m              Gauss–Newton parameter estimation
  main.m                             Original standalone MTT batch example
  localizeParticles*.m, buildTracks2.m   File-based variants (not used by the main script)
  bfmatlab/                          Bio-Formats MATLAB toolbox (bundled)
  saveastiff_4.0/                    TIFF writer (Y. Tak, BSD — see its license.txt)
```

The `*_ASH` variants take a 3-D image matrix directly instead of reading a TIFF from disk; these are the ones the main script calls.

---

## Requirements

- **MATLAB R2019a or later** (`writecell` was introduced in R2019a).
- **Image Processing Toolbox** — `multithresh`, `imquantize`, `bwdist`, `imregionalmax`, `adapthisteq`, `fspecial`, `imfilter`, `bwareaopen`.
- **Statistics and Machine Learning Toolbox** — `pdist`, `pdist2`, `squareform`.
- **Parallel Computing Toolbox** — only if you use `localizeParticlesPar.m`; the main script does not.

---

## Input data

Two folders of multi-page `.tif` stacks:

| Folder | Contents |
| --- | --- |
| `c1Path` | Single-molecule channel (e.g. 640 nm), one stack per cell/field |
| `c2Path` | Second-channel masks / concentration images, one stack per cell/field |

**File naming matters.** Files are paired by sorted order, and the pairing is only accepted if the first two underscore-separated tokens of the two filenames are identical:

```
640/   Cell01_Pos1_640.tif
Masks/ Cell01_Pos1_mask.tif     →  paired on "Cell01_Pos1"
```

Both folders must contain the same number of files, and the two channels must already be spatially registered (see the note on the hard-coded offset below).

---

## Quick start

1. Clone the repository and open MATLAB **with the repository root as the working directory** — the script calls `addpath(genpath(['.' filesep 'Batch_MTT_code' filesep]))`, which is relative to `pwd`.
2. Open `concentrationCompare_SlidingScan.m` and edit the user block at the top (lines 4–17).
3. Run the script.

### User-editable parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `c1Path` | `G:\...\640` | Folder of single-molecule stacks |
| `c2Path` | `G:\...\Masks` | Folder of second-channel stacks |
| `resultPath` | `G:\...\Result` | Output folder (must exist) |
| `tLong` | `3` | Seconds a trajectory must persist to count as long binding |
| `distTol` | `4` | Radius (px) for associating a localisation with a site |
| `lI, mI, hI` | `8, 13, 16` | Quantisation level indices defining the low / medium / high tiers |
| `fL` | `1` | Frames a pixel must be in the high tier to seed a site |
| `totalFrames` | `300` | Denominator for cluster occupancy |
| `pxSize` | `0.11` | µm per pixel |

### Acquisition and tracking parameters

Set below the `%% Parameters` divider. The defaults correspond to a 500 ms exposure, 580 nm emission, NA 1.49 objective, 0.11 µm pixels, `Dmax = 0.3 µm²/s`, one gap allowed, no deflation loops. `impars.FrameRate` is the frame *interval* in seconds despite the name.

---

## Outputs

Written to `resultPath`:

| File | Contents |
| --- | --- |
| `<name>_Occupancy.xls` | Occupancy (fraction of frames occupied) for each spatial cluster |
| `<name>_Summary.xls` | Per site: total trajectories detected, number of frames with a detection |
| `<name>_<site>.png` | Binary occupancy trace vs frame for each site |

Tier fractions (`realL`, `realM`, `realH`) are printed to the console but **not** saved to file.

---

## Known limitations

These are properties of the current code, not of the method. They are worth resolving before the repository is used by anyone else.

**Hard-coded and dataset-specific:**

- Absolute Windows paths must be replaced.
- `coord = [row + 30, col - 30]` (line 246) applies a fixed **+30/−30 pixel offset** between the two channels before distance calculations. This is a field-of-view correction to account for acquisition from two different cameras setup; it will be wrong for any other setup. Ideally, this should be derived from a bead-based registration.
- The tier indices `lI/mI/hI` are positions in a per-frame 21-level quantisation, so they correspond to **different absolute intensities in every frame and every cell**. This makes "high concentration" a relative, image-dependent call rather than a calibrated one. If the comparison of interest is between conditions or concentrations, consider fixing thresholds.
- `totalFrames = 300` is fixed, while the actual stack length (`nbImages`) is read from the file. Occupancy will be wrong for stacks of any other length.

**Interpretation:**

- Association of a localisation with a site is a proximity call within `distTol` pixels, comparable to the localisation precision and PSF width. Co-localisation at this scale does not establish molecular interaction. A randomised-position or rotated-mask control gives the expected rate under a no-association null and is worth reporting alongside the measured rate.
- The "long binding" cut-off at `tLong` seconds is not corrected for photobleaching, so the measured residence times are a lower bound. Bleaching-rate measurements under matched illumination, or a bleaching-corrected fit, would be needed for absolute dwell times.
- Per-site and per-localisation counts are not independent observations. Statistics should be computed at the level of the biological replicate (independent cells and independent experimental days), not pooled localisations.

---

## Attribution

This repository bundles third-party code:

- **MTT / SLIMFAST localisation and tracking** (`Batch_MTT_code/`) — core detection and estimation routines by the MTT authors, batch wrapper by M. Mir, in-memory variants (`*_ASH`) by Anders Sejr Hansen (2016–2017). If you publish with this code, cite the original MTT method paper.
- **`tiffread` v2.4** — François Nédélec, EMBL, 1999–2006.
- **`saveastiff` v4.0** — YoonOh Tak, 2016, BSD licence (see `Batch_MTT_code/saveastiff_4.0/license.txt`).

---

## Citation

* McSwiggen, D. T. et al. (2019). **Evidence for DNA-mediated nuclear compartmentalization distinct from phase separation.** *eLife*, 8, e47098. https://doi.org/10.7554/eLife.47098

* Sergé, A., Bertaux, N., Rigneault, H., & Marguet, D. (2008). Dynamic multiple-target tracing to probe spatiotemporal cartography of cell membranes. *Nature Methods*, 5(8), 687–694. https://doi.org/10.1038/nmeth.1233
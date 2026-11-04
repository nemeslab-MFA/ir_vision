# Imaging protocol

How to take VIS and IR images for `residuals.py`, step by step, in one session, and how to run the tool on them with the notebook `demo.ipynb` (section 9). The short version for ticking off at the microscope is [`MEASUREMENT_CHECKLIST.md`](MEASUREMENT_CHECKLIST.md).

**Why the care.** Rhombohedral graphite is only 0.4–2 % darker in IR than hexagonal graphite of the same thickness. Uneven illumination, camera offsets and VIS colour errors are just as large. The tool can remove them only if they are measured in the same session and every image is taken with the same settings.

---

## 1. Requirements

The calibration shipped with the tool (`calib/calib.npz`) is valid for:

- **90 nm SiO₂/Si** substrate;
- the same objective, IR camera and IR filter as the calibration measurements;
- flakes from **dark blue to light purple** in the VIS image (VIS red contrast about −0.70 to −0.15). Very thick flakes (light blue, yellow, white) fall outside it.

A different substrate, objective or camera needs a new calibration.

## 2. Start of the session

1. Switch on the lamp and let it **warm up for at least 20 minutes**. The lamp brightness drifts while it warms up, and the flat field (step 4) would no longer match the flake images.
2. Clean the sample surface if needed. Dust shows up as dark dots in IR.

## 3. Set the cameras (once per session)

Set everything on a piece of bare substrate next to the flakes, then **don't change any of it** until the session ends. If a setting must change, take new dark and flat frames (step 4) for the new setting.

**IR camera**

| Setting | Target | Why |
|---|---|---|
| file format | raw **16-bit TIFF**, no scale bar or text burned in | the tool needs the raw camera counts |
| black level / brightness | one fixed value | it adds an offset to every pixel; the dark frame measures it |
| exposure, gain | substrate at **50–80 % of full scale**; no flake pixel at 0; substrate not saturated | linear signal on both flake and substrate |
| frame averaging | on, if the camera has it (e.g. 16 frames) | lower noise, most important on thin flakes |
| field diaphragm | open enough for flake plus substrate margin | fixed, because the flat field depends on it |

**VIS camera**

| Setting | Target | Why |
|---|---|---|
| file format | **TIFF** (8 or 16 bit, RGB), no JPEG, no overlays | lossless colours |
| exposure | manual; **every substrate channel ≤ ~220 of 255** | saturated channels compress the colour contrast |
| white balance | manual, set on the **same white reference** every session; bare substrate ≈ **R:G:B 227:145:175** | the hexagonal prediction is calibrated on these colours |

**Both:** auto-exposure, auto white balance, HDR, auto-contrast, gamma and image enhancement **off**. Keep the objective and lamp setting the same as in the calibration.

Write all settings into `notes.txt` (step 7).

## 4. IR calibration frames (start of the session)

Take these with exactly the settings from step 3.

1. **Dark frames (required):** block the lamp or close the shutter, then take **3 IR frames**. They measure the camera offset and hot pixels.
2. **Flat frames (strongly recommended):** move to **bare substrate** with no flakes and no dust in view, then take **5–8 IR frames**, moving the stage slightly between frames. They record the vignetting, the bright spot in the centre of the field and fixed dust on the optics. Moving the stage lets the tool separate dust on the sample (which moves) from the illumination pattern (which doesn't).

## 5. Each flake

1. **Position:** centre the flake. Leave **bare substrate around it on at least 25 % of the field**, and keep the flake **away from the field edge**. The substrate is the reference for every contrast value and for removing the background.
2. **Focus** on the flake surface. The VIS image must be sharp; a blurred VIS image breaks the registration and the colour-to-IR prediction.
3. **VIS image.**
4. **IR image 1**, without moving the stage.
5. **Shift:** move the stage by about **¼ of the field** in one direction. Don't refocus unless the image is clearly soft.
6. **IR image 2** at the shifted position. This is optional but recommended: the tool averages the two, which halves the noise, and any flat-field error shows up as a mismatch between them.
7. **Check the histograms** before moving on:
   - IR: no spike at 0, substrate not saturated;
   - VIS: substrate ≤ ~220 in every channel.

## 6. End of the session

Repeat step 4: **3 dark frames and 5–8 flat frames**. If they differ clearly from the start frames, the lamp drifted during the session; use the frames closer in time to each flake.

## 7. Files and notes

Store the session like this:

```
<session>/
  calibration/   dark_start_1..3.tif   flat_start_1..8.tif
                 dark_end_1..3.tif     flat_end_1..8.tif
  <flake A>/     vis.tif   ir1.tif   ir2.tif
  <flake B>/     vis.tif   ir1.tif   ir2.tif
  notes.txt
```

`notes.txt` should contain:

- date, sample, substrate;
- objective and lamp setting;
- IR exposure, gain, black level, frame averaging, filter;
- VIS exposure, white balance and which white reference;
- times of the start and end calibration frames;
- anything unusual (refocused, setting changed, lamp flicker).

## 8. Before leaving the microscope

- [ ] Dark and flat frames taken at the start **and** at the end.
- [ ] Every flake has a VIS image and IR image 1 (and preferably IR image 2).
- [ ] Histograms checked for every flake.
- [ ] `notes.txt` complete.

## 9. Running the tool

The easiest way is the notebook **`demo.ipynb`** in the tool's folder (setup on Windows: see the README, *Quick start in Jupyter*):

1. Open it in JupyterLab and run the first cells (**Shift + Enter**): they check the installation and run the demo flakes.
2. In section 4 of the notebook, set:
   - `FLAKE_FOLDER` to a flake folder, e.g. `r"C:\data\<session>\<flake A>"`;
   - `CALIBRATION_FOLDER` to the session's `calibration` folder, so the dark and flat frames don't have to be copied into every flake folder. Leave it `None` if they are in the flake folder (names starting with `dark` and `flat`);
   - `PART = "thick"` (or `"thin"`) for extra figures of one part of the flake.
3. Run the cell. The figures appear in the notebook; the files are saved in `<flake A>/residuals/`.

For all flakes of the session at once, set `PARENT_FOLDER` to the session folder in section 5 of the notebook.

From the command line (on Windows: Anaconda Prompt), the same is:

```bash
python residuals.py "<session>/<flake A>" --dark "<session>/calibration/dark*.tif" --flats "<session>/calibration/flat*.tif"
python residuals.py "<session>" --batch --dark "<session>/calibration/dark*.tif" --flats "<session>/calibration/flat*.tif"
```

`dark*` and `flat*` take both the start and end frames; to use only one set, write `dark_start_*` and `flat_start_*`.

## 10. After running the tool

Check the warnings. They are printed in the notebook (above the figures) and saved in `<name>_summary.txt`. Each one points to something to fix in the next session:

| Warning | Fix next time |
|---|---|
| VIS substrate channel near saturation | lower the VIS exposure |
| VIS white balance differs from the calibration | set the white balance on the calibration's white reference |
| no flat field / only 2 flat frames | take 5–8 flat frames, stage moved between them |
| only one IR image | take the shifted second IR image |
| bare substrate covers only … % of the VIS image / almost no substrate in the IR field | choose a position with more substrate around the flake |
| … % of the flake outside the calibrated colour range | flake too thick, or VIS settings differ from the calibration |

# Measurement checklist for `residuals.py`

A one-page checklist to use at the microscope. Each item is something the tool needs or checks.

The calibration in `calib/calib.npz` only holds for **90 nm SiO₂/Si**, the same objective and the same IR camera and filter. Other setups need a new calibration.

---

## 1. Settings: set once, never change during the session

- [ ] Lamp switched on and warmed up for **≥ 20 min**.
- [ ] Objective, lamp setting, IR filter: the same as for the calibration. Write them down.
- [ ] **IR:** exposure, gain and **black level / brightness** fixed. Substrate at **50–80 % of full scale**; no flake pixel at 0; substrate not saturated.
- [ ] **VIS:** manual exposure. Every substrate channel **≤ ~220 of 255** (the tool warns above ~245).
- [ ] **VIS:** manual white balance on the **same white reference** as the calibration. Bare substrate should read about **R:G:B = 227:145:175** (the tool warns when off by > 25 %).
- [ ] HDR, auto-contrast, auto-exposure, gamma and image enhancement **off** in both cameras.
- [ ] Field diaphragm set and **not touched again**.

## 2. Calibration frames: same settings as the flakes

- [ ] **IR dark:** 1–3 frames, lamp blocked or shutter closed (**required**; the tool stops without it).
- [ ] **IR flat:** 3–8 frames of **bare substrate**, stage moved a little between frames. These remove vignetting, the central bright spot and dust. With only 2 frames the tool still works but warns; with none it only fits a plane.
- [ ] Take the dark and flat frames again **at the end** of the session, to check for lamp drift.

## 3. Each flake

- [ ] Flake centred, **not touching the field edge**, with **bare substrate on ≥ 25 % of the field** in both VIS and IR. With less, the background can only be corrected by a plane and the residual shows gradients.
- [ ] Flake colour between **dark blue and light purple**. Very thick (light blue, yellow) parts fall outside the calibration and are greyed out.
- [ ] Focus on the flake surface; the VIS image must be **sharp**.
- [ ] **VIS image:** RGB, raw TIFF (8 or 16 bit), no overlays (scale bar, text).
- [ ] **IR image 1:** raw 16-bit TIFF, no overlays. If the camera can average frames, use it (e.g. 16 frames); thin flakes need the low noise most.
- [ ] Move the stage by **about ¼ of the field** in one direction, without refocusing.
- [ ] **IR image 2** at the shifted position, same settings (recommended: halves the noise and checks the flat field).
- [ ] Histograms checked: no IR spike at 0, IR substrate not saturated, VIS substrate ≤ ~220.

## 4. Files for the tool

One folder per flake. The tool recognises files by name: `dark*`/`black*` are dark frames, `flat*` are flat frames, the RGB image is the VIS image, and the rest are IR images (at most two).

```
my_flake/
  vis.tif
  ir1.tif
  ir2.tif                 (shifted, optional)
  dark1.tif ... dark3.tif
  flat1.tif ... flat8.tif
```

- [ ] The session's dark and flat frames are **copied into every flake folder**, or passed with `--dark "…/dark*.tif" --flats "…/flat*.tif"`.
- [ ] No other images (screenshots, exports) in the flake folder.
- [ ] A `notes.txt` with: date, objective, IR exposure / gain / black level / filter, VIS exposure and white balance (which white reference), lamp setting, substrate.

## 5. After running the tool

- [ ] Read `<name>_summary.txt`. Each warning points to a setting from sections 1–3 to fix next time.
- [ ] Compare blue areas only with neighbouring areas of the same VIS colour (thickness).

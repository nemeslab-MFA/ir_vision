# Imaging protocol: VIS + IR images for rhombohedral graphite detection

Version 2026-10-21. Follow the steps in order during one session. The background is in `measurement_checklist.md`.

**Why it matters.** The RG signal in IR is small: about 0.004–0.02 in contrast for purple to dark-blue flakes, and less on thin flakes. Uneven illumination, camera offsets and VIS colour errors are as large as that, unless they are measured in the same session and the images are taken the same way every time.

---

## 1. Before the session

- [ ] Choose the flakes (section 7). Keep **one flake blind**: measure its Raman map, but don't open the map.
- [ ] Use the same substrate as before (90 nm SiO₂/Si). Write it down if it is different.
- [ ] Switch on the lamp and let it **warm up for at least 20 minutes**.

## 2. Fix the settings (once, for the whole session)

Set everything once and write it into `notes.txt`. **Don't change any of these during the session.** If something has to change, repeat steps 3 and 6 for the new settings.

| | Setting | Target |
|---|---|---|
| IR | exposure, gain, **brightness / black level** | substrate at about 50–80 % of full scale (≈ 11 000 of 16 383 is good); no flake pixel at 0; thin parts not saturated |
| IR | field diaphragm | open enough for the flake plus a substrate margin; **don't touch it again** |
| IR | filter / band | always the same; write it down |
| VIS | exposure | manual; **red channel of the substrate at most ~200 of 255** (now ≈ 225, too bright: thick parts saturate) |
| VIS | white balance | manual, set on the **same white reference sample** every session |
| both | objective, lamp setting | always the same; write them down |
| both | HDR, auto-contrast, gamma, enhancement | **off** |
| both | file format | **16-bit TIFF**, raw camera image, no overlays (scale bar, text) |

## 3. Calibration frames at the START

Take all of these with the settings from step 2.

- [ ] **IR dark:** 3 frames with the lamp blocked or the shutter closed.
- [ ] **IR flat:** 5–8 frames of **bare substrate** (no flakes, no dust in view). Move the stage a little between frames.
- [ ] **VIS flat (new, required):** 3–5 frames of bare substrate, stage moved between frames, same VIS exposure and white balance.

## 4. Each flake

1. Centre the flake. Leave **substrate all around it, at least ~25 % of the field**. The flake must not touch the field edge.
2. Focus on the flake surface. VIS has to be **sharp** (375 was out of focus and failed).
3. **VIS image**, position A.
4. **IR, position A: 16 frames** (or as many as practical). Averaging them cuts the noise 4×, which is most important for **thin flakes**.
5. Move the stage by about **¼ of the field** (60–90 IR px) in one direction. Don't refocus unless the image is clearly soft.
6. **IR, position B: 16 frames.** Also take a **VIS image at position B** (it checks the VIS flat field).
7. Check the histograms before moving on:
   - IR: no spike at 0, no saturated substrate.
   - VIS: substrate not saturated, red channel ≤ ~200.

## 5. Raman maps (calibration flakes only)

- [ ] Map the whole flake plus a little substrate.
- [ ] Export as PNG, ideally with a **transparent background**; the colour-bar layout of IR v4 also works. Raw data values are a bonus.
- [ ] Write down the **scale-bar length** and the step size (µm per pixel) for every map.
- [ ] **Check that the map's file name matches the flake** (`66-03` vs `66_01` caused doubts).
- [ ] Blind flake: save its map, but don't open it.

## 6. Calibration frames at the END

- [ ] Repeat step 3: IR dark (3), IR flat (5–8), VIS flat (3–5).

Comparing the start and end frames shows whether the lamp or the illumination drifted during the session.

## 7. Which flakes to image

- [ ] **Thin flakes and thin parts (light purple, a few to ~10 layers)** with hexagonal and RG areas, plus Raman maps. These matter most scientifically, and the calibration has almost no data there. Image them with the most IR frames (16 or more).
- [ ] Purple to dark-blue flakes with **RG and hexagonal areas at the same thickness**, side by side.
- [ ] A spread of thicknesses, including some thick **hexagonal** areas.
- [ ] Flakes at different positions in the field (the flake itself centred, but different flakes at different stage positions).
- [ ] Few wrinkles and folds on the areas of interest.
- [ ] Note in `notes.txt` which flakes have light-blue (very thick) parts.

## 8. Files and notes

```
<session folder>/
  calibration/  ir_dark_start_1..3.tif   ir_flat_start_1..8.tif   vis_flat_start_1..5.tif
                ir_dark_end_1..3.tif     ir_flat_end_1..8.tif     vis_flat_end_1..5.tif
  VIS/          <flake>_VIS_A.tif   <flake>_VIS_B.tif
  IR/           <flake>_IR_A_01..16.tif   <flake>_IR_B_01..16.tif
  raman/        <flake>.png   (later: <flake> copy.png with RG painted pure green, unsure pure blue)
  notes.txt
```

`notes.txt` should contain:
- date and objective
- IR exposure, gain, brightness / black level, filter
- VIS exposure and white balance (and which white sample)
- lamp setting and substrate
- the times of the start and end calibration frames
- which flakes are thick, and which flake is blind
- for each Raman map: scale bar and step size

## 9. Before leaving the microscope

- [ ] All calibration frames are taken, at the start **and** at the end.
- [ ] Every flake has VIS A, VIS B, IR A (×16) and IR B (×16).
- [ ] The histograms were checked (step 4.7).
- [ ] `notes.txt` is complete.
- [ ] One flake is kept blind.

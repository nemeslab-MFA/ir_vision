"""IR/VIS residual maps for finding rhombohedral (ABC) graphite in few-layer graphite flakes.

Residual = IR contrast - IR contrast that HEXAGONAL graphite of the same VIS colour (thickness) would have.
Rhombohedral graphite is darker in IR, so blue (negative) areas are RG candidates. See README.md.

Quick start (one folder per flake, containing vis*.tif, ir*.tif, dark*.tif, flat*.tif):
    python residuals.py path/to/flake_folder
    python residuals.py path/to/parent_folder --batch        # every sub-folder that contains images
Explicit files:
    python residuals.py --vis V.tif --ir IR.tif [--ir2 IR_shifted.tif] --dark "dark*.tif" --flats "flat*.tif"
In Python / Jupyter (see demo.ipynb):
    import residuals
    r = residuals.run("path/to/flake_folder", part="thick")
    residuals.show(r)
"""
import argparse, glob, os, time
import numpy as np
import tifffile
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from skimage import transform, morphology
from skimage.registration import phase_cross_correlation

HERE = os.path.dirname(os.path.abspath(__file__))
WARN = []                 # warnings collected during a run (printed and written to the summary)
FLAT = {"median": None, "used_for_ir": False, "sub_share": 1.0, "vis_sub_share": 1.0}   # substrate shares: IR field / whole VIS image   # median flat level (dark subtracted): absolute substrate reference for flake-filled fields
VIS_MAX = 1600            # VIS long side used by the calibration geometry
K = 50                    # kNN neighbours


class InputError(Exception):
    """Missing or unusable input files (a clear message instead of a traceback)."""


def reset_state():
    WARN.clear()
    FLAT.update(median=None, used_for_ir=False, sub_share=1.0, vis_sub_share=1.0)


# ---------------- loading / corrections ----------------
def load_vis(path):
    v = tifffile.imread(path).astype(np.float32)[..., :3]
    sat = np.any(v >= (255 if v.max() <= 255 else 65535) - 1, -1)
    f = VIS_MAX / max(v.shape[:2])
    if f < 1:
        v = transform.resize(v, (round(v.shape[0] * f), round(v.shape[1] * f), 3), order=1, anti_aliasing=True, preserve_range=True)
        sat = transform.resize(sat.astype(float), v.shape[:2], order=1) > 0.5    # majority saturated (noise specks ignored)
    return v, sat


def combine_flats(frames):
    """Combine bare-substrate frames. Dust on the substrate only DARKENS a frame (small dark specks) and moves between
    frames when the stage is moved. >= 3 frames: the median rejects it. 2 frames: the median is just the mean, so
    wherever one frame has a local dark speck that the other does not, the other frame's value is used."""
    if len(frames) != 2:
        return np.median(frames, 0)
    a, b = frames
    ra, rb = a / ndi.median_filter(a, 9) - 1, b / ndi.median_filter(b, 9) - 1       # local relative deviation
    speck_a = ndi.binary_dilation((ra < -0.02) & (ra < rb - 0.02), iterations=2)
    speck_b = ndi.binary_dilation((rb < -0.02) & (rb < ra - 0.02), iterations=2)
    out = 0.5 * (a + b)
    out = np.where(speck_a & ~speck_b, b, out)
    return np.where(speck_b & ~speck_a, a, out)


def ir_session(dark_glob, flat_glob, ir_example=None):
    """Dark + flat field (one or more frames each). Missing flats: a plain field mask from the IR image (warning)."""
    pick = lambda x: list(x) if isinstance(x, (list, tuple)) else (sorted(glob.glob(x)) if x else [])
    darks, flats = pick(dark_glob), pick(flat_glob)            # a glob pattern or a list of files
    if not darks:
        raise InputError("no dark frame found (dark*.tif, lamp blocked, same IR settings) - it is required.")
    dark = np.median([tifffile.imread(f).astype(np.float32) for f in darks], 0)
    if len(flats) == 2:
        WARN.append("only 2 flat-field frames: dust specks were removed by comparing the two frames; 3-8 frames "
                    "(stage moved between them) remove dust more reliably.")
    if not flats:
        WARN.append("no flat field (flat*.tif): illumination is only corrected by a plane -> expect slow gradients "
                    "of ~0.01-0.04 under the flake. Take 3-8 bare-substrate IR frames next time.")
        I = tifffile.imread(ir_example).astype(np.float32) - dark
        field = I > 0.5 * np.median(I)
        return dark, np.where(field, 1.0, np.nan).astype(np.float32), morphology.binary_erosion(field, morphology.disk(12))
    hot = dark > np.median(dark) + 10 * 1.4826 * np.median(np.abs(dark - np.median(dark)))
    dark = np.where(hot, ndi.median_filter(dark, 5), dark)
    F = combine_flats([tifffile.imread(f).astype(np.float32) - dark for f in flats])
    field = F > 0.5 * np.median(F)
    FLAT["median"] = float(np.median(F[field]))
    Fn = np.where(field, F / np.median(F[field]), np.nan)
    field = morphology.binary_erosion(field, morphology.disk(12))
    return dark, Fn, field


def load_ir(path, dark, Fn):
    I = tifffile.imread(path).astype(np.float32)
    sat = I >= 16383
    return (I - dark) / Fn, sat


# ---------------- geometry ----------------
def vis_to_ir(img, tf, shape, order=1):
    """Resample a VIS-frame image (H, W[, C]) into the IR frame; tf maps VIS (x, y) -> IR (x, y)."""
    inv = transform.SimilarityTransform(matrix=tf).inverse
    sig = 0.5 / abs(tf[0, 0]) if order == 1 else 0
    if img.ndim == 3:
        return np.stack([vis_to_ir(img[..., k], tf, shape, order) for k in range(img.shape[2])], -1)
    src = ndi.gaussian_filter(img, sig) if sig else img
    return transform.warp(src, inv, output_shape=shape, order=order, cval=np.nan, preserve_range=True)


def register(vis_gray, ir, field, tf0, search=120):
    """Translation only (scale/rotation fixed by tf0): phase correlation of gradient magnitudes inside the field."""
    v = vis_to_ir(vis_gray, tf0, ir.shape)
    ok = field & np.isfinite(v) & np.isfinite(ir)
    win = ndi.gaussian_filter(morphology.binary_erosion(ok, morphology.disk(15)).astype(float), 8)   # soft edge: the
    def g(a):                                                                                          # field border must
        gm = ndi.gaussian_gradient_magnitude(np.where(ok, a, np.nanmedian(a[ok])), 1.5)                # not dominate
        gm = gm / np.nanpercentile(gm[ok], 99)
        return (gm - np.mean(gm[ok])) * win
    shift, _, _ = phase_cross_correlation(g(ir), g(v), upsample_factor=10, normalization=None)
    dy, dx = shift
    if abs(dy) > search or abs(dx) > search:
        raise RuntimeError(f"registration shift too large ({dx:.0f}, {dy:.0f}) px - check the images")
    tf = np.array([[1, 0, dx], [0, 1, dy], [0, 0, 1.0]]) @ tf0
    return tf, (dx, dy)


# ---------------- contrasts ----------------
def poly_bg(img, mask, order):
    H, W = img.shape
    yy, xx = np.mgrid[:H, :W]; x = xx / W - 0.5; y = yy / H - 0.5
    A = np.stack([x ** i * y ** j for i in range(order + 1) for j in range(order + 1 - i)], -1)
    m = mask & np.isfinite(img)
    cf = np.linalg.lstsq(A[m], img[m], rcond=None)[0]
    return A @ cf


def substrate_level(g, min_share=0.05):
    """Brightness of the bare substrate: the BRIGHTEST histogram peak holding >= min_share of the pixels (within
    +-4 %). Graphite on SiO2 is darker than the substrate, so this works even when a flake fills most of the image."""
    h, e = np.histogram(g, bins=128); c = 0.5 * (e[:-1] + e[1:])
    hs = ndi.gaussian_filter1d(h.astype(float), 1.5)
    peaks = [c[i] for i in range(1, len(hs) - 1) if hs[i] >= hs[i - 1] and hs[i] >= hs[i + 1]]
    big = [p for p in peaks if np.mean(np.abs(g / p - 1) < 0.04) >= min_share]
    h2, e2 = np.histogram(g, bins=256); common = e2[np.argmax(h2)]     # the level the calibration was built with
    level = max(big) if big else common
    return common if abs(level / common - 1) < 0.04 else level      # normal images: unchanged; flake-filled: brightest peak


def vis_contrast(vis, vsat):
    """Per-channel VIS contrast in the VIS frame. The substrate is found on the whole VIS image (much larger than the
    IR field, so plenty of substrate even under a big flake): dominant brightness peak, then iterated with a
    3rd-order background per channel. Returns rgb contrast (H, W, 3) and the substrate mask."""
    g = vis.mean(-1); ok = ~vsat & (g > 0)
    mode = substrate_level(g[ok])
    sub = ok & (np.abs(g / mode - 1) < 0.06)
    FLAT["vis_sub_share"] = float(sub.mean())
    order = 3 if sub.mean() >= 0.25 else 1        # substrate only at the borders (flake fills the image): a plane,
    if order == 1:                                # not a 3rd-order surface that would bend under the flake
        WARN.append(f"bare substrate covers only {sub.mean():.0%} of the VIS image: VIS background reduced to a plane.")
    for _ in range(3):
        rgb = np.stack([vis[..., k] / poly_bg(vis[..., k], sub, order) - 1 for k in range(3)], -1)
        new = morphology.binary_erosion(ok & (np.abs(rgb).max(-1) < 0.03), morphology.disk(3))
        if new.sum() < 3000:                      # never let the substrate mask collapse
            break
        sub = new
    return rgb, sub


def contrasts(rgb, sub, ir, field, sat):
    """rgb, sub already in the IR frame. IR: plane fitted on the substrate. Returns rgb, ir_c, substrate, flake."""
    ok = field & np.all(np.isfinite(rgb), -1) & ~sat
    sub = sub & ok
    FLAT["used_for_ir"] = False
    FLAT["sub_share"] = float(sub.sum() / max(field.sum(), 1))
    if sub.sum() >= 500:
        ir_c = ir / poly_bg(ir, sub & np.isfinite(ir), 1) - 1
    elif FLAT["median"]:                           # flake fills the IR field: the flat field IS the substrate reference
        ir_c = ir / FLAT["median"] - 1
        FLAT["used_for_ir"] = True
        WARN.append(f"almost no bare substrate inside the IR field ({sub.sum()} px): IR contrast taken directly from the "
                    "flat field (assumes identical IR settings and a stable lamp between flats and flake image).")
    else:
        ir_c = ir / poly_bg(ir, sub & np.isfinite(ir), 0) - 1
        WARN.append("almost no bare substrate inside the IR field and no flat field: IR contrast is unreliable.")
    flake = ok & (rgb.mean(-1) < -0.04)
    flake = morphology.remove_small_objects(morphology.binary_opening(flake, morphology.disk(1)), 60)
    return rgb, ir_c, sub, flake


# ---------------- hexagonal IR model (hybrid) ----------------
class Hybrid:
    def __init__(self, cal):
        self.X, self.y = cal["X"], cal["y"]
        self.sd = self.X.std(0) + 1e-9
        self.tree = cKDTree(self.X / self.sd)
        self.p = cal["poly"]; self.d50, self.d95 = float(cal["d50"]), float(cal["d95"])

    def predict(self, X):
        dd, nn = self.tree.query(X / self.sd, k=K)
        pk = np.median(self.y[nn], 1); dist = dd[:, -1]
        pp = np.polyval(self.p, X[:, 0])
        w = np.clip((dist - self.d50) / (self.d95 - self.d50), 0, 1)
        return (1 - w) * pk + w * pp


def branch_dist(rgb, pg, pb):
    R = rgb[..., 0]
    return np.hypot(rgb[..., 1] - np.polyval(pg, R), rgb[..., 2] - np.polyval(pb, R))


# ---------------- small holes ----------------
def fill_small_holes(img, inside, max_px=30):
    """Fill small NaN holes (<= max_px pixels, e.g. isolated saturated or defect pixels) inside `inside` by a
    normalised Gaussian average of the valid neighbours. Larger holes stay NaN. Returns (filled image, filled mask)."""
    holes = inside & ~np.isfinite(img)
    lab, n = ndi.label(holes)
    if n == 0:
        return img, np.zeros(img.shape, bool)
    sizes = np.bincount(lab.ravel()); small = (sizes <= max_px); small[0] = False
    fill = small[lab]
    ok = np.isfinite(img)
    out = img.copy()
    for sig in (1.0, 2.0, 4.0):                   # grow the neighbourhood until every small hole has neighbours
        num = ndi.gaussian_filter(np.where(ok, img, 0.0), sig); den = ndi.gaussian_filter(ok.astype(float), sig)
        todo = fill & ~np.isfinite(out) & (den > 1e-3)
        out[todo] = (num / np.maximum(den, 1e-9))[todo]
    return out, fill & np.isfinite(out)


# ---------------- one flake ----------------
def process(vis_path, ir_path, dark, Fn, field, cal, ir2_path=None, tf=None, no_thick=False):
    """Returns dict with residual, ir_c, rgb, flake, excluded, tf (all in the IR frame of ir_path)."""
    vis, vsat = load_vis(vis_path)
    ir, isat = load_ir(ir_path, dark, Fn)
    tf0 = cal["tf0"] if tf is None else tf
    if tf is None:
        tf, shift = register(vis.mean(-1), ir, field, tf0)
    rgb_v, sub_v = vis_contrast(vis, vsat)
    rgb_ir = vis_to_ir(rgb_v, tf, ir.shape)
    sub_ir = vis_to_ir(sub_v.astype(np.float32), tf, ir.shape, order=0) > 0.5
    sat = (vis_to_ir(vsat.astype(np.float32), tf, ir.shape, order=0) > 0.5) | isat
    if ir2_path:                                    # shifted IR: register it to the same VIS, map into this frame
        ir2, isat2 = load_ir(ir2_path, dark, Fn)
        tf2, _ = register(vis.mean(-1), ir2, field, tf)
        M = transform.ProjectiveTransform(matrix=tf2 @ np.linalg.inv(tf))
        ir2w = transform.warp(np.where(field, ir2, np.nan), M, output_shape=ir.shape, order=1, cval=np.nan, preserve_range=True)
        s2 = transform.warp(isat2.astype(float), M, output_shape=ir.shape, order=0) > 0.5
        a = np.where(field, ir, np.nan)
        ir = np.where(np.isfinite(a) & np.isfinite(ir2w), 0.5 * (a + ir2w), a)
        sat |= s2
    rgb, ir_c, sub, flake = contrasts(rgb_ir, sub_ir, ir, field, sat)
    thick = np.zeros(flake.shape, bool)
    if np.any(cal["pg"]) and not no_thick:         # calibration mode (make_calib.py) fits the branch afterwards
        thick = branch_dist(rgb, cal["pg"], cal["pb"]) > 0.10
        thick = morphology.binary_opening(ndi.median_filter(thick, 3) & flake, morphology.disk(1))
    excluded = (thick | morphology.binary_dilation(sat, morphology.disk(1))) & flake
    use = flake & ~excluded & np.all(np.isfinite(rgb), -1) & np.isfinite(ir_c)
    res = np.full(ir.shape, np.nan, np.float32)
    if cal is not None and "X" in cal:
        res[use] = ir_c[use] - Hybrid(cal).predict(rgb[use])
    res, filled = fill_small_holes(res, flake)    # isolated blank pixels -> neighbours' average
    excluded = excluded & ~filled
    ir_c = np.where(filled & ~np.isfinite(ir_c), fill_small_holes(ir_c, flake)[0], ir_c)
    return dict(residual=res, filled=filled, ir_c=ir_c.astype(np.float32), rgb=rgb.astype(np.float32), flake=flake,
                excluded=excluded, substrate=sub, tf=tf, thick_share=float(thick[flake].mean()) if flake.any() else 0.0,
                vis_ir=vis_to_ir(vis / np.percentile(vis, 99.5), tf, ir.shape))



# ---------------- input finding and checks ----------------
def find_inputs(folder):
    """Classify the TIFFs of one flake folder by name and shape."""
    if not os.path.isdir(folder):
        raise InputError(f"folder not found: {folder}")
    files = sorted({os.path.normcase(f): f for e in ("*.tif", "*.tiff", "*.TIF", "*.TIFF")
                    for f in glob.glob(os.path.join(folder, e))}.values())   # Windows globs ignore case: no duplicates
    low = lambda f: os.path.basename(f).lower()
    dark = [f for f in files if low(f).startswith(("dark", "black"))]
    flat = [f for f in files if low(f).startswith("flat")]
    rest = [f for f in files if f not in dark + flat]
    vis, ir = [], []
    for f in rest:
        with tifffile.TiffFile(f) as t:
            shp = t.pages[0].shape
        (vis if len(shp) == 3 and shp[-1] in (3, 4) else ir).append(f)
    if len(vis) != 1:
        raise InputError(f"{folder}: expected exactly one VIS (RGB) image, found {len(vis)}: {[os.path.basename(v) for v in vis]}")
    if not ir:
        raise InputError(f"{folder}: no IR image found")
    if len(ir) > 2:
        WARN.append(f"more than two IR images; using {os.path.basename(ir[0])} and {os.path.basename(ir[1])}, "
                    f"ignoring {[os.path.basename(i) for i in ir[2:]]} (use --ir/--ir2 to choose)")
    return dict(vis=vis[0], ir=ir[0], ir2=ir[1] if len(ir) > 1 else None, dark=dark, flat=flat)


def check_vis(path):
    """Warn about VIS settings that break the colour calibration (saturation, white balance)."""
    v = tifffile.imread(path)[..., :3]
    full = 255 if v.max() <= 255 else 65535
    ok = np.all(v < 0.98 * full, -1)                          # saturated (very thick) flakes are not substrate
    if 1 - ok.mean() > 0.2:                                   # the substrate itself is saturated: no colour reference
        WARN.append(f"VIS image overexposed: {1 - ok.mean():.0%} of the pixels are saturated, probably including the bare "
                    "substrate. The substrate reference is then wrong and the residual is NOT usable. Retake the VIS image "
                    "with a lower exposure (substrate channels <= ~220).")
    g = v.astype(np.float32).mean(-1); h, e = np.histogram(g[ok], bins=256); mode = e[np.argmax(h)]
    sub = ok & (np.abs(g / mode - 1) < 0.05)
    lev = [float(np.median(v[..., k][sub])) / full * 255 for k in range(3)]
    if max(lev) > 245:
        WARN.append(f"VIS substrate channel near saturation (R,G,B = {lev[0]:.0f}, {lev[1]:.0f}, {lev[2]:.0f} of 255): "
                    "contrasts are compressed. Use a lower VIS exposure (all channels <= ~220).")
    ref = np.array([227, 145, 175.])                      # substrate levels of the calibration sessions
    if np.max(np.abs(np.array(lev) / lev[1] - ref / ref[1])) > 0.25:
        WARN.append(f"VIS white balance differs from the calibration (substrate R:G:B = {lev[0]:.0f}:{lev[1]:.0f}:{lev[2]:.0f}, "
                    "calibration ~227:145:175). Set the white balance on the same white reference as the calibration.")
    return lev


# ---------------- output ----------------
def save(out, name, r, vis_orientation=True):
    """Outputs: the three panels (visible, corrected IR, residual) as single SVGs, drawn exactly as in the combined
    SVG, plus the combined figure; raw float data in data/. vis_orientation: rotate by 180 deg (the IR camera is
    rotated ~180 deg relative to VIS)."""
    data = os.path.join(out, "data")
    os.makedirs(data, exist_ok=True)
    rot = (lambda a: np.rot90(a, 2)) if vis_orientation else (lambda a: a)
    res, ir, flake, exc, vis = (rot(r[k]) for k in ("residual", "ir_c", "flake", "excluded", "vis_ir"))
    outline = rot(r["outline"]) if "outline" in r else None
    if outline is not None:                       # part figures: show only the part in the IR panel too
        ir = np.where(outline, ir, np.nan)
    tifffile.imwrite(os.path.join(data, f"{name}_vis.tif"), np.clip(np.nan_to_num(vis, nan=0) * 255, 0, 255).astype(np.uint8))
    tifffile.imwrite(os.path.join(data, f"{name}_ir_corrected.tif"), ir.astype(np.float32))
    tifffile.imwrite(os.path.join(data, f"{name}_residual.tif"), res.astype(np.float32))
    if "filled" in r:
        tifffile.imwrite(os.path.join(data, f"{name}_filled_pixels.tif"), (rot(r["filled"]) * 255).astype(np.uint8))
    from matplotlib.figure import Figure          # no pyplot: works the same in scripts and Jupyter, any backend
    ys, xs = np.nonzero(flake); sl = (slice(max(ys.min() - 15, 0), ys.max() + 15), slice(max(xs.min() - 15, 0), xs.max() + 15))
    med = float(np.nanmedian(res))
    grey = np.zeros(flake.shape + (4,)); grey[exc] = (0.5, 0.5, 0.5, 1)

    def panel(ax, kind):                          # one drawing routine for the combined AND the single images
        if kind == "visible":
            ax.imshow(np.clip(np.nan_to_num(vis[sl], nan=1.0), 0, 1)); ax.set_title("visible")
            if outline is not None:
                ax.contour(outline[sl], [0.5], colors="yellow", linewidths=1.2); ax.set_title("visible (part outlined)")
        elif kind == "ir":
            im = ax.imshow(ir[sl], cmap="gray", vmin=np.nanpercentile(ir[flake], 1), vmax=np.nanpercentile(ir[flake], 99.5))
            ax.figure.colorbar(im, ax=ax, fraction=0.046)
            ax.set_title("IR, corrected (dark, flat field, background)\ncontrast to the substrate")
        else:
            im = ax.imshow(res[sl], cmap="RdBu_r", vmin=med - 0.01, vmax=med + 0.01); ax.imshow(grey[sl])
            ax.figure.colorbar(im, ax=ax, fraction=0.046)
            ax.set_title(f"residual = IR − predicted hexagonal IR\nblue = RG candidate (scale centred on {med:+.3f})")
        ax.axis("off")

    kinds = (("visible", "1_visible"), ("ir", "2_ir_corrected"), ("residual", "3_residual"))
    fig = Figure(figsize=(20, 6.5)); ax = fig.subplots(1, 3)
    for a_, (k, _) in zip(ax, kinds):
        panel(a_, k)
    fig.suptitle(name, fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(os.path.join(out, f"{name}_residual.svg"))
    for k, fn in kinds:                           # the same panels, one per file, same size as in the combined figure
        fig = Figure(figsize=(20 / 3, 6.5)); a_ = fig.subplots()
        panel(a_, k)
        fig.suptitle(name, fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(os.path.join(out, f"{name}_{fn}.svg"))


def select_part(r, part, by="auto"):
    """Mask of the thick (or thin) part of the flake: the dark (or light) group of a 2-class split (Otsu).
    by = "ir": corrected IR contrast (thicker = darker); "vis": VIS brightness (thicker = darker); "auto": IR, but VIS
    when bare substrate covers < 50 % of the VIS image (flake fills the frame -> IR gradient). Small holes
    (excluded pixels) are bridged so the part comes out as whole regions."""
    from skimage.filters import threshold_otsu
    if by == "auto":                               # flake fills the frame -> IR gradient not removable -> use VIS
        by = "vis" if (FLAT["used_for_ir"] or FLAT["vis_sub_share"] < 0.5) else "ir"
    region = ndi.binary_fill_holes(morphology.binary_closing(r["flake"], morphology.disk(3)))
    img = r["ir_c"] if by == "ir" else np.nanmean(r["vis_ir"], -1)
    img = np.where(region & np.isfinite(img), img, np.nan)
    fill = ndi.median_filter(np.nan_to_num(img, nan=float(np.nanmedian(img))), 7)
    sm = ndi.gaussian_filter(np.where(np.isfinite(img), img, fill), 3)
    thr = threshold_otsu(sm[region])
    m = region & ((sm < thr) if part == "thick" else (sm >= thr))
    m = ndi.binary_fill_holes(morphology.remove_small_objects(morphology.binary_opening(m, morphology.disk(3)), 400))
    WARN.append(f"--part {part}: selected by {'IR contrast' if by == 'ir' else 'VIS brightness'} (2-class split at "
                f"{thr:+.3f}); {int(m.sum())} px. The {part}-part figures centre the residual scale on this part only.")
    return m


def run_one(name, inp, out, cal, no_thick=False, vis_orientation=True, part=None, part_by="auto"):
    t0 = time.time()
    check_vis(inp["vis"])
    dark, Fn, field = ir_session(inp["dark"], inp["flat"], inp["ir"])
    r = process(inp["vis"], inp["ir"], dark, Fn, field, cal, inp.get("ir2"), no_thick=no_thick)
    if not no_thick and r["thick_share"] > 0.5:
        WARN.append(f"{r['thick_share']:.0%} of the flake is outside the calibrated colour range (too thick, or VIS settings "
                    "differ) -> recomputed with --no-thick. The absolute residual is then offset: compare areas of the "
                    "same VIS colour (last panel).")
        r = process(inp["vis"], inp["ir"], dark, Fn, field, cal, inp.get("ir2"), no_thick=True)
        no_thick = True
    if not inp.get("ir2"):
        WARN.append("only one IR image: a second IR image with the flake shifted by ~1/4 field halves the noise and "
                    "checks the flat field.")
    save(out, name, r, vis_orientation)
    figures = [os.path.join(out, f"{name}_residual.svg")]
    if part:                                       # the same figures again, for the thick (or thin) part only
        m = select_part(r, part, part_by)
        rp = dict(r, flake=m & r["flake"], excluded=r["excluded"] & m,
                  residual=np.where(m, r["residual"], np.nan).astype(np.float32), outline=m)
        save(out, f"{name}_{part}", rp, vis_orientation)
        figures.append(os.path.join(out, f"{name}_{part}_residual.svg"))
        tifffile.imwrite(os.path.join(out, "data", f"{name}_{part}_mask.tif"),
                         ((np.rot90(m, 2) if vis_orientation else m) * 255).astype(np.uint8))
    f = r["flake"] & ~r["excluded"] & np.isfinite(r["residual"])
    res = r["residual"][f]; med = float(np.median(res)) if f.any() else float("nan")
    lines = [f"flake {name}"] + [f"  {k}: {v}" for k, v in inp.items()] + [ f"time {time.time() - t0:.1f} s", f"flake pixels {int(f.sum())}",
             f"residual median {med:+.4f}, robust std {1.4826 * float(np.median(np.abs(res - med))):.4f}",
             f"colour cut {'off (--no-thick)' if no_thick else 'on'}",
             f"small holes filled from neighbours: {int(r['filled'].sum())} px (see data/{name}_filled_pixels.tif)", "warnings:"] + [f"  - {w}" for w in WARN]
    with open(os.path.join(out, f"{name}_summary.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"{name}: done in {time.time() - t0:.1f} s -> {figures[0]}  (residual median {med:+.4f})")
    for w in WARN:
        print(f"  WARNING: {w}")
    return dict(name=name, out=out, figures=figures, residual_median=med, warnings=list(WARN), result=r)


# ---------------- entry points for Python / Jupyter ----------------
def run(folder=None, *, vis=None, ir=None, ir2=None, dark=None, flats=None, out=None, name=None, part=None,
        part_by="auto", no_thick=False, ir_orientation=False, calib=None):
    """Process one flake. Either `folder` (files found by name, see README) or explicit files: vis, ir, optional ir2,
    dark and flats (a file, a list of files, or a pattern like "C:/data/dark*.tif"). With a folder, dark and flats
    (if given) replace the folder's own dark*/flat* files, e.g. one calibration folder for the whole session.
    part: None, "thick" or "thin". Outputs go to <folder>/residuals unless `out` is given.
    Returns a dict: name, out, figures (SVG paths), residual_median, warnings, result (the arrays)."""
    reset_state()
    cal = dict(np.load(calib or os.path.join(HERE, "calib", "calib.npz")))
    if vis:
        if not (ir and dark):
            raise InputError("vis= needs ir= and dark= (and preferably flats=)")
        inp = dict(vis=vis, ir=ir, ir2=ir2, dark=dark, flat=flats)
        name = name or os.path.splitext(os.path.basename(ir))[0]
        out = out or "residuals"
    elif folder:
        folder = os.path.normpath(folder)
        inp = find_inputs(folder)
        if dark:                                    # shared session calibration frames instead of the folder's own
            inp["dark"] = dark
        if flats:
            inp["flat"] = flats
        name = name or os.path.basename(folder)
        out = out or os.path.join(folder, "residuals")
    else:
        raise InputError("give a flake folder, or vis=, ir= and dark=")
    return run_one(name, inp, out, cal, no_thick, not ir_orientation, part, part_by)


def run_batch(parent, out=None, **kw):
    """Process every sub-folder of `parent` that contains TIFFs. A flake with missing files is skipped with a message.
    Keyword arguments as in run(). Returns the list of run() results."""
    subs = sorted(d for d in glob.glob(os.path.join(parent, "*")) if os.path.isdir(d) and find_tifs(d))
    print(f"{len(subs)} flake folders")
    done = []
    for d in subs:
        try:
            done.append(run(d, out=out and os.path.join(out, os.path.basename(d)), **kw))
        except InputError as e:
            print(f"{os.path.basename(d)}: skipped - {e}")
    return done


def find_tifs(folder):
    return [f for e in ("*.tif", "*.tiff", "*.TIF", "*.TIFF") for f in glob.glob(os.path.join(folder, e))]


def show(results, width=1100):
    """Show the figures of run() / run_batch() results in a Jupyter notebook."""
    import re
    from IPython.display import HTML, display
    for r in results if isinstance(results, list) else [results]:
        for f in r["figures"]:
            with open(f, encoding="utf-8") as fh:
                svg = fh.read()
            svg = svg[svg.index("<svg"):]                                      # drop the XML header
            svg = re.sub(r'(<svg[^>]*?)\swidth="[^"]*"\s+height="[^"]*"', r'\1 width="100%"', svg, count=1)  # fit the page
            display(HTML(f"<div style='max-width:{width}px'>{svg}</div>"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", nargs="?", help="flake folder (or parent folder with --batch)")
    ap.add_argument("--batch", action="store_true", help="process every sub-folder of FOLDER that contains TIFFs")
    ap.add_argument("--vis"); ap.add_argument("--ir"); ap.add_argument("--ir2", help="same flake, shifted (optional)")
    ap.add_argument("--dark", help='file or glob, e.g. "dark*.tif" (with a folder: replaces its own dark files)'); ap.add_argument("--flats", help='glob, e.g. "flat*.tif" (with a folder: replaces its own flat files)')
    ap.add_argument("--out", help="output folder (default: <flake folder>/residuals)")
    ap.add_argument("--name", help="flake name used in the output files (default: folder name)")
    ap.add_argument("--calib", default=f"{HERE}/calib/calib.npz", help="calibration file (default: calib/calib.npz)")
    ap.add_argument("--no-thick", action="store_true",
                    help="keep pixels outside the calibrated colour range (switched on automatically if > 50 %% of the flake)")
    ap.add_argument("--part", choices=["thick", "thin"], help="also write the figures for the thick (or thin) part only")
    ap.add_argument("--part-by", choices=["auto", "ir", "vis"], default="auto",
                    help="how the part is found: IR contrast, VIS brightness, or auto (IR; VIS when the flake fills the field)")
    ap.add_argument("--ir-orientation", action="store_true", help="keep the IR camera orientation (default: rotate to VIS orientation)")
    a = ap.parse_args()
    kw = dict(part=a.part, part_by=a.part_by, no_thick=a.no_thick, ir_orientation=a.ir_orientation, calib=a.calib)
    try:
        if a.vis:                                          # explicit files
            run(vis=a.vis, ir=a.ir, ir2=a.ir2, dark=a.dark, flats=a.flats, out=a.out, name=a.name, **kw)
        elif a.folder and a.batch:
            run_batch(a.folder, out=a.out, dark=a.dark, flats=a.flats, **kw)
        elif a.folder:
            run(a.folder, out=a.out, name=a.name, dark=a.dark, flats=a.flats, **kw)
        else:
            ap.print_help()
    except InputError as e:
        raise SystemExit(f"ERROR: {e}")

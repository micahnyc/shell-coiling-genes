"""
Measure Raup parameters (W, D, T, S) from a photo of a shell by clicking landmarks.

Two kinds of photo are supported:

  spire        A snail photographed in "apertural view": the opening faces the
               camera and the tip of the spire points up. Gives W, D, T and S.

  planispiral  A flat coil (nautilus, ammonite, ram's-horn snail) photographed
               looking straight down the coiling axis. Gives W and D. T is 0
               by definition; S cannot be seen from this view, so pass --S if
               you measured it from another photo (default 1).

Usage
-----
    python measure.py photo.jpg --mode spire
    python measure.py photo.jpg --mode planispiral --S 1.3
    python measure.py --landmarks photo.landmarks.json   # recompute, no clicking
    python measure.py --selftest                         # check the maths

Clicking: left-click to place a point, right-click (or Backspace) to undo the
last point, Enter to finish a step that takes several points.

Each run saves:
  <photo>.landmarks.json   the clicked points (so you can recompute later)
  <photo>.fit.png          the photo with landmarks next to the fitted model
  measurements.csv         one row per specimen (appended)

How the numbers are computed (apex/centre = origin of the logarithmic spiral)
---------------------------------------------------------------------------
  W  successive whorl edges along one ray shrink by a factor W per turn, so
     W = exp(-slope) of log(distance from apex) against whorl number.
  D  inner / outer distance of the opening from the coiling axis.
  S  height / width of the opening.
  T  (distance of the opening's centre along the axis from the apex) /
     (distance of the opening's centre from the axis).

Caveats: a photo is a 2D projection, so the opening is measured as it appears
in the photo's plane. Real openings are often tilted, and in tightly coiled
(involute) shells the inner edge of the opening is hidden, which makes D come
out too large. Treat values as estimates and be consistent between specimens.
"""

import argparse
import csv
import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

from raup import draw

STEPS = {
    "spire": [
        ("apex", 1, "Click the APEX (very tip of the spire)"),
        ("base", 1, "Click the BASE of the coiling axis (bottom end of the columella)"),
        ("ap_outer", 1, "Click the OUTER edge of the opening (furthest from the axis)"),
        ("ap_inner", 1, "Click the INNER edge of the opening (nearest the axis)"),
        ("ap_top", 1, "Click the TOP of the opening (nearest the apex)"),
        ("ap_bottom", 1, "Click the BOTTOM of the opening"),
        ("whorls", None, "Click the outer edge of each earlier whorl on the SAME side as "
                         "the opening, going up the spire. Enter when done (at least 1)"),
    ],
    "planispiral": [
        ("centre", 1, "Click the CENTRE of the coil"),
        ("ap_outer", 1, "Click the OUTER edge of the opening (furthest from the centre)"),
        ("ap_inner", 1, "Click the INNER edge of the opening (nearest the centre)"),
        ("whorls", None, "Click the outer edge of each earlier whorl along the SAME line "
                         "from the opening to the centre, going inward. Enter when done"),
    ],
}


# ---------------------------------------------------------------- maths

def fit_W(distances):
    """distances[k] = distance from the apex/centre of the same point k whorls ago."""
    d = np.asarray(distances, float)
    if len(d) < 2 or np.any(d <= 0):
        raise ValueError("need the opening edge plus at least one earlier whorl")
    k = np.arange(len(d))
    slope = np.polyfit(k, np.log(d), 1)[0]
    return float(np.exp(-slope))


def compute_spire(lm):
    apex = np.array(lm["apex"][0], float)
    u = np.array(lm["base"][0], float) - apex
    u /= np.linalg.norm(u)

    def rz(p):
        v = np.asarray(p, float) - apex
        z = v @ u
        r = abs(v[0] * u[1] - v[1] * u[0])
        return r, z

    r_out, _ = rz(lm["ap_outer"][0])
    r_in, _ = rz(lm["ap_inner"][0])
    _, z_top = rz(lm["ap_top"][0])
    _, z_bot = rz(lm["ap_bottom"][0])

    D = r_in / r_out
    S = (z_bot - z_top) / (r_out - r_in)
    r_c, z_c = (r_out + r_in) / 2, (z_top + z_bot) / 2
    T = z_c / r_c
    pts = [lm["ap_outer"][0]] + list(lm["whorls"])
    W = fit_W([np.linalg.norm(np.asarray(p, float) - apex) for p in pts])
    return dict(W=W, D=D, T=T, S=S)


def compute_planispiral(lm, S=1.0):
    c = np.array(lm["centre"][0], float)
    dist = lambda p: float(np.linalg.norm(np.asarray(p, float) - c))
    r_out, r_in = dist(lm["ap_outer"][0]), dist(lm["ap_inner"][0])
    W = fit_W([r_out] + [dist(p) for p in lm["whorls"]])
    return dict(W=W, D=r_in / r_out, T=0.0, S=S)


def compute(lm, mode, S=None):
    if mode == "spire":
        return compute_spire(lm)
    return compute_planispiral(lm, S if S is not None else 1.0)


def estimate_turns(lm, W):
    """Number of whorls to draw: enough to shrink the opening to ~1% of its size."""
    n_clicked = len(lm["whorls"])
    return float(max(n_clicked + 1, min(15, np.log(100) / np.log(max(W, 1.01)))))


# ---------------------------------------------------------------- clicking

def collect(image_path, mode):
    img = plt.imread(image_path)
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.imshow(img)
    ax.set_axis_off()
    lm = {}
    for key, n, prompt in STEPS[mode]:
        ax.set_title(prompt, fontsize=11, wrap=True)
        fig.canvas.draw_idle()
        print(prompt)
        pts = plt.ginput(n=n if n else -1, timeout=0, show_clicks=True)
        if n and len(pts) < n:
            sys.exit("cancelled")
        lm[key] = [list(map(float, p)) for p in pts]
        xs, ys = zip(*pts) if pts else ([], [])
        ax.plot(xs, ys, "o", ms=6, mec="k", mfc="cyan" if key != "whorls" else "yellow")
    plt.close(fig)
    return lm


# ---------------------------------------------------------------- output

def save_figure(image_path, lm, mode, params, out_png):
    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(1, 2, 1)
    if image_path and os.path.exists(image_path):
        ax1.imshow(plt.imread(image_path))
    else:
        ax1.invert_yaxis()
        ax1.set_aspect("equal")
    for key, pts in lm.items():
        if pts:
            xs, ys = zip(*pts)
            ax1.plot(xs, ys, "o", ms=6, mec="k",
                     mfc="yellow" if key == "whorls" else "cyan")
    ax1.set_axis_off()
    ax1.set_title(os.path.basename(image_path) if image_path else "landmarks")

    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    p = dict(params, turns=estimate_turns(lm, params["W"]),
             view=(90, -90) if mode == "planispiral" else (0, 90))
    ax2.set_proj_type("ortho")  # photos are (nearly) flat projections
    draw(ax2, p, "fitted model")
    if mode == "spire":
        ax2.invert_zaxis()  # photos have the apex at the top
        ax2.invert_xaxis()  # face the opening, coiled right-handed (dextral)
    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)


def append_csv(path, row):
    new = not os.path.exists(path)
    with open(path, "a", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(row))
        if new:
            w.writeheader()
        w.writerow(row)


# ---------------------------------------------------------------- self-test

def synthetic_landmarks(mode, W, D, T, S, n_whorls=3, rot_deg=7, scale=120, shift=(400, 80)):
    """Landmarks of a perfect model shell, placed into 'pixel' coordinates."""
    a = 1.0
    Rc = a * (1 + D) / (1 - D)
    if mode == "spire":
        # (r, z) in the plane of the photo, apex at origin, z pointing to the base
        pts = dict(apex=[(0, 0)], base=[(0, 2 * T * Rc + 3)],
                   ap_outer=[(Rc + a, T * Rc)], ap_inner=[(Rc - a, T * Rc)],
                   ap_top=[(Rc, T * Rc - a * S)], ap_bottom=[(Rc, T * Rc + a * S)],
                   whorls=[((Rc + a) * W ** -k, T * Rc * W ** -k) for k in range(1, n_whorls + 1)])
    else:
        pts = dict(centre=[(0, 0)], ap_outer=[(Rc + a, 0)], ap_inner=[(Rc - a, 0)],
                   whorls=[((Rc + a) * W ** -k, 0) for k in range(1, n_whorls + 1)])
    th = np.radians(rot_deg)
    R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    return {k: [list(R @ np.array(p) * scale + shift) for p in v] for k, v in pts.items()}


def selftest():
    cases = [("spire", dict(W=2.2, D=0.05, T=1.8, S=1.1)),
             ("spire", dict(W=1.25, D=0.1, T=6.0, S=1.0)),
             ("planispiral", dict(W=3.0, D=0.0, T=0.0, S=1.3)),
             ("planispiral", dict(W=1.9, D=0.35, T=0.0, S=1.0))]
    ok = True
    for mode, true in cases:
        lm = synthetic_landmarks(mode, **true)
        got = compute(lm, mode, S=true["S"])
        err = max(abs(got[k] - true[k]) for k in "WDTS")
        ok &= err < 1e-6
        print(f"{mode:12s} true {true}  ->  "
              + "  ".join(f"{k}={got[k]:.4g}" for k in "WDTS")
              + ("   OK" if err < 1e-6 else f"   FAIL (err {err:.2g})"))
    lm = synthetic_landmarks("spire", W=2.2, D=0.05, T=1.8, S=1.1)
    params = compute(lm, "spire")
    save_figure(None, lm, "spire", params, "selftest.fit.png")
    print("saved selftest.fit.png")
    return ok


# ---------------------------------------------------------------- main

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("image", nargs="?", help="photo of the shell")
    p.add_argument("--mode", choices=STEPS, default="spire")
    p.add_argument("--landmarks", help="recompute from a saved .landmarks.json")
    p.add_argument("--S", type=float, help="opening shape (planispiral mode only)")
    p.add_argument("--specimen", help="name for the CSV row (default: file name)")
    p.add_argument("--csv", default="measurements.csv")
    p.add_argument("--selftest", action="store_true")
    args = p.parse_args()

    if args.selftest:
        sys.exit(0 if selftest() else 1)

    if args.landmarks:
        saved = json.load(open(args.landmarks, encoding="utf-8"))
        lm, mode, image = saved["landmarks"], saved["mode"], saved.get("image")
        S = args.S if args.S is not None else saved.get("S")
        stem = args.landmarks.replace(".landmarks.json", "")
    elif args.image:
        mode, image, S = args.mode, args.image, args.S
        lm = collect(image, mode)
        stem = os.path.splitext(image)[0]
        json.dump(dict(image=image, mode=mode, S=S, landmarks=lm),
                  open(stem + ".landmarks.json", "w", encoding="utf-8"), indent=1)
    else:
        p.error("give a photo, --landmarks, or --selftest")

    params = compute(lm, mode, S)
    print("  ".join(f"{k}={params[k]:.3f}" for k in "WDTS"))
    save_figure(image, lm, mode, params, stem + ".fit.png")
    append_csv(args.csv, dict(specimen=args.specimen or os.path.basename(stem), mode=mode,
                              n_whorls_clicked=len(lm["whorls"]),
                              **{k: round(params[k], 4) for k in "WDTS"}))
    print(f"saved {stem}.fit.png and a row in {args.csv}")


if __name__ == "__main__":
    main()

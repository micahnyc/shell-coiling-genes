"""
Raup shell model: generate coiled tubes (shells, horns) from four parameters.

Based on Raup, D. M. (1966) "Geometric analysis of shell coiling: general
problems", Journal of Paleontology 40(5): 1178-1190.

Parameters
----------
W : whorl expansion rate. How much the tube's size multiplies per full turn.
    Snails ~2, nautilus ~3, clams 100s-10,000s (the valve is less than one turn).
D : distance of the tube from the coiling axis, as (inner edge / outer edge)
    of the opening. 0 = tube touches the axis; near 1 = loose, open coil.
T : translation rate. How far the tube moves along the axis relative to
    how far it moves outward. 0 = flat (planispiral); large = tall spire.
S : shape of the opening, as height/width of an ellipse. 1 = circle.

The model: take the opening (the "generating curve"), then rotate it around
the axis while scaling it by W ** (turns). Every point of the tube is the
same opening at a different size and angle, which is why one small set of
numbers can describe nautilus, clams, snails and horns.

Usage
-----
    python raup.py                       # grid of all presets -> shells.png
    python raup.py --preset snail        # one preset
    python raup.py --W 2 --D 0.2 --T 1.5 --S 1.2 --turns 5   # custom shell
    python raup.py --list                # show the presets
"""

import argparse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Approximate parameter values, chosen by eye to give recognisable shapes.
# They are starting points, not measurements of real specimens.
PRESETS = {
    "nautilus": dict(view=(80, -90), W=3.0, D=0.0, T=0.0, S=1.3, turns=3.0,
                     note="flat, tightly coiled, each whorl covers the last"),
    "ammonite": dict(view=(80, -90), W=1.9, D=0.35, T=0.0, S=1.0, turns=5.0,
                     note="flat but open coil; whorls only just touch"),
    "rams-horn snail": dict(view=(80, -90), W=2.0, D=0.25, T=0.0, S=1.0, turns=4.5,
                            note="Planorbidae: flat coil like a ram's horn"),
    "snail": dict(view=(15, -60), W=2.2, D=0.05, T=1.8, S=1.1, turns=5.0,
                  note="garden-snail-like spire"),
    "tower snail": dict(view=(10, -60), W=1.25, D=0.1, T=6.0, S=1.0, turns=12.0,
                        note="Turritella-like tall spire"),
    "clam": dict(view=(25, -40), W=10000.0, D=0.0, T=0.3, S=1.0, turns=0.5,
                 note="one valve: huge expansion, less than one turn"),
    "ram's horn": dict(view=(55, -70), W=6.0, D=0.45, T=0.35, S=0.9, turns=1.0,
                       note="bighorn-like curl: open coil drifting sideways"),
    "markhor horn": dict(view=(15, -60), W=1.6, D=0.6, T=4.0, S=1.0, turns=3.0,
                         note="corkscrew: open coil with strong translation"),
}


def raup_surface(W, D, T, S, turns, n_theta=400, n_phi=60):
    """Return X, Y, Z arrays of the tube surface and the growth fraction.

    The tube ends at theta = 0 with an opening of radius 1, and is traced
    backwards in time (towards the tip/apex) for `turns` full turns.
    """
    theta = np.linspace(-2 * np.pi * turns, 0, n_theta)  # coiling angle
    phi = np.linspace(0, 2 * np.pi, n_phi)  # angle around the opening
    theta, phi = np.meshgrid(theta, phi)

    a = 1.0  # radius of the final opening (sets overall scale)
    D = min(D, 0.999)
    # Centre of the opening: D = (Rc - a) / (Rc + a) solved for Rc
    Rc = a * (1 + D) / (1 - D)

    growth = W ** (theta / (2 * np.pi))  # size relative to the final opening

    radial = (Rc + a * np.cos(phi)) * growth
    axial = (T * Rc + a * S * np.sin(phi)) * growth

    X = radial * np.cos(theta)
    Y = radial * np.sin(theta)
    Z = axial
    age = (theta - theta.min()) / (theta.max() - theta.min())  # 0 old, 1 new
    return X, Y, Z, age


def draw(ax, params, title=None):
    X, Y, Z, age = raup_surface(**{k: params[k] for k in "WDTS"},
                                turns=params["turns"])
    colors = cm.YlOrBr(0.25 + 0.6 * age)
    ax.plot_surface(X, Y, Z, facecolors=colors, rstride=1, cstride=1,
                    linewidth=0, antialiased=True, shade=True)
    spans = [np.ptp(X), np.ptp(Y), np.ptp(Z)]
    ax.set_box_aspect([max(s, 1e-3 * max(spans)) for s in spans])
    ax.set_axis_off()
    elev, azim = params.get("view", (25, -60))
    ax.view_init(elev=elev, azim=azim)
    label = (f"W={params['W']:g}  D={params['D']:g}  "
             f"T={params['T']:g}  S={params['S']:g}")
    ax.set_title(f"{title}\n{label}" if title else label, fontsize=9)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--preset", choices=PRESETS.keys())
    p.add_argument("--W", type=float)
    p.add_argument("--D", type=float)
    p.add_argument("--T", type=float)
    p.add_argument("--S", type=float)
    p.add_argument("--turns", type=float)
    p.add_argument("--out", default=None, help="output image (PNG)")
    p.add_argument("--show", action="store_true", help="open an interactive window")
    p.add_argument("--list", action="store_true", help="list presets and exit")
    args = p.parse_args()

    if args.list:
        for name, v in PRESETS.items():
            print(f"{name:16s} W={v['W']:<8g} D={v['D']:<5g} T={v['T']:<5g} "
                  f"S={v['S']:<5g} turns={v['turns']:<5g} {v['note']}")
        return

    custom = any(getattr(args, k) is not None for k in ("W", "D", "T", "S", "turns"))

    if args.preset or custom:
        params = dict(PRESETS[args.preset]) if args.preset else \
            dict(W=2.0, D=0.1, T=1.5, S=1.0, turns=5.0)
        for k in ("W", "D", "T", "S", "turns"):
            if getattr(args, k) is not None:
                params[k] = getattr(args, k)
        fig = plt.figure(figsize=(6, 6))
        draw(fig.add_subplot(projection="3d"), params, args.preset or "custom")
        out = args.out or f"{(args.preset or 'custom').replace(' ', '_')}.png"
    else:
        fig = plt.figure(figsize=(16, 8))
        for i, (name, params) in enumerate(PRESETS.items()):
            draw(fig.add_subplot(2, 4, i + 1, projection="3d"), params, name)
        fig.suptitle("Raup shell model: one tube, four parameters", fontsize=14)
        out = args.out or "shells.png"

    fig.tight_layout()
    fig.savefig(out, dpi=130)
    print(f"saved {out}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()

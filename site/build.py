"""
Build the results page from analysis/results.

Writes docs/index.html (GitHub Pages) with the chart data inlined.
The page text and tables live in site/page_src.html; the chart data
(genome scan and LG6 positions) are read from the result files here.

Usage:  python site/build.py
"""

import json
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
RES = os.path.join(ROOT, "analysis", "results")

HEAD = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<style>:root{padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}
body{margin:0;font:14px system-ui,sans-serif}img{max-width:100%}[hidden]{display:none!important}</style>
</head>
<body>
"""


def plot_data():
    lods = pd.read_csv(os.path.join(RES, "lod_scores.csv.gz"))
    lods["bin"] = lods.cM.round(0)
    out = {}
    for t in ["W", "D_like", "S", "shellLength"]:
        g = lods.groupby(["LG", "bin"])[t].max().reset_index()
        out[t] = [[int(a), float(b), round(float(c), 2)] for a, b, c in g.values]
    out["lens"] = {int(k): round(v, 1) for k, v in lods.groupby("LG").cM.max().items()}
    pos = pd.read_csv(os.path.join(RES, "marker_v2_positions.csv")).merge(lods[["marker", "S"]], on="marker")
    lg6 = pos[pos.LG == 6]
    out["lg6"] = [[round(b / 1e6, 2), round(c, 1), round(s, 2)] for b, c, s in lg6[["bp", "cM", "S"]].values]
    return out


def main():
    src = open(os.path.join(HERE, "page_src.html"), encoding="utf-8").read()
    page = src.replace("/*DATA*/", json.dumps(plot_data()))
    os.makedirs(os.path.join(ROOT, "docs"), exist_ok=True)
    with open(os.path.join(ROOT, "docs", "index.html"), "w", encoding="utf-8") as f:
        f.write(HEAD + page + "\n</body>\n</html>\n")
    print("wrote docs/index.html")


if __name__ == "__main__":
    main()

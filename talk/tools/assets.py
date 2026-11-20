import json, sys
g = json.load(open(sys.argv[1])); out = sys.argv[2]
INK, RULE, VIR = "#110F0E", "#D5D4D1", "#0B6E58"
for s in (39, 93, 100, 102):
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 400 260" width="400" height="260">'
           f'<path d="M0 {g["y0"]}H400" stroke="{RULE}" stroke-width="1" fill="none"/>'
           f'<path d="{g[f"faint{s}"]}" stroke="{INK}" stroke-opacity="0.22" stroke-width="5" stroke-linecap="round" fill="none"/>'
           f'<path d="{g["model"]}" stroke="{VIR}" stroke-width="2.4" stroke-linejoin="round" fill="none"/>'
           f'<path d="{g[f"bins{s}"]}" stroke="{INK}" stroke-width="7.5" stroke-linecap="round" fill="none"/>'
           f'</svg>')
    open(f"{out}/transit-{s}.svg", "w").write(svg)
svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1664 300" width="1664" height="300">'
       f'<path d="{g["cover"]}" stroke="{INK}" stroke-width="11" stroke-linecap="round" fill="none"/></svg>')
open(f"{out}/folded.svg", "w").write(svg)
svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 624 624" width="624" height="624" shape-rendering="crispEdges">'
       f'{g["pixels"]}</svg>')
open(f"{out}/star-pixels.svg", "w").write(svg)

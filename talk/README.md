# Talk: Discovering an exoplanet

A 5–6 minute talk for a general audience on how TIC 402710280.01 was found. The deck is a Claude Slides artifact at <https://claude.ai/artifact/BCwWDprCzFEZiNFWLzbwzq> (private to its owner). This folder holds its source as of 2026-11-20, so the deck can be rebuilt if the artifact is lost.

## The story, one sentence per slide

| # | Slide file | Sentence on the slide | Backdrop | Time |
|---|---|---|---|---|
| 1 | `cover` | Discovering an exoplanet | The four transits overlaid, as dots | 0:15 |
| 2 | `light` | The light in this picture left its star in 1803. | The real TESS image, the star marked as on a finder chart | 0:30 |
| 3 | `dimmed` | In June 2021, the star dimmed. | One light curve, circled by hand; 0.2%, 5 hours | 0:30 |
| 4 | `transit` | A planet crossing its star does exactly this. | Star, planet and path to scale; the planet moves | 0:35 |
| 5 | `archive` | The dip sat in NASA's public archive. No one had flagged it. | 10,324,827 | 0:25 |
| 6 | `claude` | I searched the archive with Claude. | Who did what | 0:40 |
| 7 | `timeline` | Four years later, it dimmed again. | Time line, three click steps | 0:50 |
| 8 | `planet` | It fits a planet nearly four times as wide as Earth. | Earth, the candidate and Neptune to scale | 0:30 |
| 9 | `impostors` | Most candidates turn out to be something else. | Three impostors and why each fails | 0:25 |
| 10 | `triceratops` | So we asked TRICERATOPS for the odds. | The acronym spelled out, 0.3%, a doodle | 0:25 |
| 11 | `next` | If the planet is real, the star dims again on 19 February 2027. | The sky again; blog and paper addresses | 0:30 |

Things to say that are not on the slides:

- Slide 6: we found time stamps running about two minutes early in an MIT data product (QLP, Sectors 100 and 102) and reported it.
- Slide 7: the 2026 dips come from our own photometry of the raw full-frame images. The three clicks add the 2025 dip, the two 2026 dips, and the 45 crossings nobody measured.
- Slide 11: the TESS Science Office approved two-minute cadence for the 2027 transit as DDT 106 (`proposals/tess-ddt/`).

Wordings that are deliberate:

- "TESS caught four cleanly": a fifth crossing, in Sector 66, was recorded but swamped by scattered light.
- "No one had flagged it": as of 2026-10-18 the star was not a TOI or community TOI and had no entry with planet parameters on ExoFOP, the Exoplanet Archive or ADS.
- "0.3%": the TRICERATOPS false-positive probability, 0.003 ± 0.004 over ten runs. Without high-resolution imaging the object stays a candidate.
- "1803": the June 2021 image, minus 218 light-years (66.83 pc).

## Files

- `deck/deck.json`, `deck/deck.css`, `deck/slides/*.html`: the deck in the Slides artifact's own format (one `<section>` per file, 1920×1080 canvas). In the artifact these live under `project/`.
- `assets/*.svg`: the images the slides reference. In the artifact each is an uploaded asset, addressed as `/_blob/<id>`:

| File | Asset id |
|---|---|
| `transit-39.svg` | `9c60989b2b496bed04b0998d6c13db0e` |
| `transit-93.svg` | `66f445a845c45439fe22db7d75044656` |
| `transit-100.svg` | `fca25a41645452f7227dfa5ae44c32b9` |
| `transit-102.svg` | `4df9f9db08545d31f8cf3e1bb9215086` |
| `folded.svg` | `579841f83000deff968100d47cc61726` |
| `sky.svg` | `ac6cca3efc84afed6fc98f1ae37d7162` |

A new artifact assigns new ids, so the `src` attributes in the slides must be updated after re-uploading.

## Rebuilding the figures

Every chart is drawn from the published data record (doi:10.5281/zenodo.23390273) and the star image from a TESScut cutout of Sector 39. Requires `numpy` and `astropy`.

```sh
# 1. Cutouts of the Sector 39 full-frame images: 13×13 pixels, and 43×19 for the full-bleed sky
curl -L -o cut13.zip "https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra=248.5236&dec=-63.5497&y=13&x=13&sector=39"
curl -L -o cut43.zip "https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra=248.5236&dec=-63.5497&y=19&x=43&sector=39"

# 2. Median out-of-transit image of each cutout
python tools/pix.py <13x13 cutout>.fits img13.npy
python tools/pix.py <43x19 cutout>.fits img43.npy

# 3. Chart geometry from the transit windows of the data record, then the SVG files
python tools/geom.py <record>/data/transit_windows img13.npy geom.json
python tools/assets.py geom.json assets
python tools/sky.py img43.npy assets/sky.svg
```

`assets.py` also writes `star-pixels.svg`, a 13×13 view the deck no longer uses.

## Presenting and exporting

- Motion plays only in the artifact's Present mode: the hand circle on slide 3, the moving planet on slide 4 and the click steps on slide 7. Exports to Google Slides, PowerPoint and PDF are static.
- The deck is set in Google Sans at its four named weights (400, 500, 600, 700) so that an export to Google Slides finds the same face. An uploaded font falls back to Arial there.
- Palette: paper `#FAF9F7`, ink `#110F0E`, secondary ink `#575552`; viridian `#0B6E58` for anything computed or fitted, ochre `#CE871B` for marks made by hand, ultramarine `#2256BB` for the current position, cinnabar `#D02E1E` (`#EE6951` on dark) for something to follow. Dark ground `#131110`.

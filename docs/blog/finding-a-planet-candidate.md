---
title: "Finding a Planet Candidate in Public TESS Data"
date: 2026-10-15
description: "A week-long search of ten million TESS light curves turned up a possible sub-Neptune 218 light-years away. How it was found, what an audit broke, and the February test that could rule it out."
tags: [astronomy, exoplanets, tess, ai]
---

Last week a search of public NASA data turned up a planet candidate. It orbits a star called TIC 402710280, about 218 light-years away in the southern constellation Triangulum Australe, once every 36.34 days. If it is real, it is a little smaller than Neptune and sits at around 500 K. It may not be real. TESS will look at the star again on 19 February 2027, and I have already published the exact test that will decide whether the planet survives that observation.

This post explains how it was found, in plain terms. The honest version includes several things we got wrong along the way, which turned out to be the most useful part.

First, credit. Claude, Anthropic's AI model running in Claude Code, did nearly all of the technical work. It wrote the search pipeline, ran it in parallel cloud sessions, vetted the candidates, fitted the transits and drafted the paper. I set the goal, asked questions, pushed back, arranged independent checks, and did the parts that need a human, like corresponding with archives and filling in Zenodo forms. For the first three days, two of my handful of messages were "Continue?" and "Try again". When I write "we" below, read it generously.

## How it started

On 9 October I asked Claude whether we could search the TESS archive for something nobody had reported yet, the way amateur astronomers sometimes do. I pointed it at [Paras Chopra's tess-planet-candidates](https://github.com/paraschopra/tess-planet-candidates) as an example, and said I'd be fine if we found nothing. Being an astronomer, even an amateur one, has been a childhood wish of mine.

Six days later there was a preprint with a DOI. Most of this post is about how the candidate was checked in between.

## Why a planet like this gets missed

TESS, the Transiting Exoplanet Survey Satellite, looks for transits: the small, regular dips in a star's brightness when a planet passes in front of it. It observes the sky in strips called sectors, about 27 days each, and then moves on. A typical star gets a sector, then nothing for a year or two, then another sector.

A planet with a 36-day orbit therefore shows at most one transit per sector. Planet searches work by finding dips that repeat, so a single dip in 2021 and another in 2025 tend to stay unconnected. TIC 402710280 had a second disadvantage. It has only ever been recorded in TESS's full-frame images, the wide pictures of each camera's whole field, and never at the two-minute cadence that the mission's main search pipeline processes.

## The search

Claude searched 10,324,827 public light curves (a light curve is a star's brightness over time) for individual transit-shaped dips. It slid box-shaped templates of ten different durations, from 45 minutes to 18 hours, along each light curve and flagged anything more than six times the local noise.

A search that large finds a great deal of rubbish, so every search was also run on the same light curves turned upside down. An inverted light curve contains no transits, because dips become bumps, so anything the pipeline "finds" there is a false alarm. That gives a measured false-alarm rate to compare against, instead of a guessed one.

After removing glitches shared by many stars at once, dips at the edges of data gaps, and shapes that don't look like transits, 2,263 stars still had at least one convincing single dip. The 142 of them brighter than TESS magnitude 10.5 were checked individually, at the pixel level and against catalogues of known planets and binaries. Twenty-seven survived with no known explanation. For those, Claude made its own light curves from the 2026 full-frame images, which no public pipeline had processed yet, and looked for more transits.

## Finding the period

In the public data, TIC 402710280 had two dips: one on 19 June 2021 (Sector 39) and one on 11 June 2025 (Sector 93). They had the same depth and the same duration, and were 1,453.58 days apart.

Two transits don't give you a period. They give you a list of candidates: 1,453.58 days, or half that, or a third, or 1,453.58 divided by any whole number. A transit in 2023 (Sector 66) might have settled it, but the data at the predicted times were swamped by sunlight scattered off the Earth and Moon, which is common in TESS images.

Claude's photometry of the 2026 images found two more transits, on 21 February and 29 March 2026 (Sectors 100 and 102), 36.36 days apart. Because those two are on consecutive orbits, with data in between showing no transit, the period can't be a multiple of 36.36 days. It could still be a fraction of it, such as 18.2 or 12.1 days, so Claude checked every fraction from a half to a fifth. Each one predicts transits at times when TESS was watching and saw nothing. The weakest of those non-detections is at 18.6σ; for half the period it is 28.5σ. Only one period fits all the data: 36.33986 ± 0.00051 days, a precision of about 44 seconds. The gap between the 2021 and 2025 dips is exactly 40 orbits.

![The four transits](fig1.png)
*The four transits used in the fit, with the Sector 66 window on the right. Grey points are individual measurements, blue points are 20-minute averages, and the red curve is the transit model.*

## What it would be

The transit removes 0.21% of the star's light (2,092 ± 47 parts per million). The depth tells you the planet's size relative to the star: here, 4.2% of the star's radius. The star is a K dwarf, a bit smaller and cooler than the Sun (0.82 solar radii, about 5,180 K), so the planet would be 3.75 ± 0.22 times the radius of Earth. Neptune is 3.88.

It would orbit at 0.21 AU, about half Mercury's distance from the Sun, and receive about ten times the sunlight Earth does. Its equilibrium temperature, assuming it reflects nothing and spreads the heat evenly, is about 500 K. That number is a rough guide and not a surface temperature, since a planet this size probably has a thick atmosphere.

Planets between Earth and Neptune in size are the most common kind found so far, and the Solar System has none, so nobody knows quite what they're like inside. A typical planet of this radius would have about 14 Earth masses, and would pull its star back and forth at about 2.8 metres per second. Current spectrographs can measure that for a star this bright (V = 10). Measuring the mass would give the planet's density, which is the first real clue to its composition.

None of this makes it special. Plenty of sub-Neptunes are known. Ones with orbits longer than a month around stars this bright and nearby are less common, and they're the ones where follow-up measurements are practical.

## Reasons it might not be a planet

Most signals that look like planets in TESS data turn out to be something else. The usual impostor is a blended eclipsing binary: two stars eclipsing each other somewhere near the target, whose deep eclipses get diluted by the target's light into a shallow, planet-sized dip. TESS pixels are 21 arcseconds across, so many stars can share one.

The first check is where the missing light comes from. If you subtract an image taken during a transit from one taken just outside it, what's left is a picture of the dimming. Fitting its position for all four transits puts the source 1.5 ± 1.3 arcseconds from the target. Before trusting that, Claude tested the method on 40 confirmed TESS planets. It put 30 of them on their host star, couldn't decide for 5, failed or gave an odd answer for 5, and never put one on the wrong star.

That rules out every catalogued star bright enough to cause the signal except one. A faint background star sits 3.5 arcseconds away, about 350 times fainter than the target and roughly 20 times farther away. It's too close for the image test to separate. To produce a 0.21% dip in the combined light, though, it would have to lose about 74% of its own light at each eclipse, and an eclipse that deep between two ordinary stars has a different shape from what we see. Claude fitted the transits with every eclipse that kind of pair could produce. All of them fit worse than a planet, by a wide margin (Δχ² ≥ 142, roughly 12σ). The only blends that come close would need an unusual pairing, a small hot star eclipsed by a larger, much dimmer one. A ground-based telescope in ordinary seeing can separate two stars 3.5 arcseconds apart, so one light curve of a future transit would settle this.

The target itself could be an eclipsing binary. If so, there would usually be a second eclipse half an orbit later. There isn't one, down to 77 parts per million. The Gaia satellite also sees the star as a single, well-behaved point of light.

The standard tool for weighing all of this is TRICERATOPS, which estimates how likely each false-positive scenario is given the data and the stars nearby. It gives a false-positive probability of 0.003. The usual threshold for calling a candidate a validated planet is 0.015, but that threshold assumes someone has imaged the star at high resolution to rule out a companion within about an arcsecond. Nobody has, so this stays a candidate.

## The audit

By 13 October I was excited, and asked Claude why it kept describing the find so cautiously when we seemed to have discovered a planet. It ran the false-positive numbers and I got more excited. I wrote to ExoFOP, the archive where TESS candidates are posted, to ask for an account. The reply politely asked for my institutional affiliation, which I don't have.

Then I had two independent audits run on Claude's report. I read the first one feeling fairly deflated, because it was mostly right:

- The analysis had counted a fifth transit, in the Sector 66 data swamped by scattered light. It shouldn't have counted, and including it had pulled the measured depth down. Without it the depth went from 1,655 to about 2,100 ppm, and the planet grew from 3.4 to 3.75 Earth radii.
- The transit shape had been fitted by stacking all the transits on one strict period. Because the transits don't arrive exactly on schedule (more on that below), the stack was smeared, and the fitted shape was wrong. Refitting with each transit centred on its own measured time fixed it.
- The period's uncertainty was quoted as ±6 seconds. Once the timing scatter was accounted for, it is ±44 seconds.
- Predicted transit times labelled as UTC were in fact in a barycentric time scale, about 6 minutes off.
- The neighbouring star's separation had been taken from catalogue coordinates for the year 2000. The target moves across the sky by 61 milliarcseconds a year, so during the transits the stars were 3.5 arcseconds apart, not the 2.2 we had been using.
- ExoFOP changed its rules on 19 August 2026: a new candidate now needs a peer-reviewed paper before it can be posted. Claude had told me to upload soon without checking that page.

Every one of these has been fixed, and the candidate came out of it with stronger evidence than before. The corrected transit shape now matches the density expected for this kind of star, and the false-positive probability dropped from 0.015 to 0.003. The second audit had measured the transit times with its own photometry and got the same pattern of early and late transits we did, which was reassuring in an odd way.

The lesson I'd take from this is about cost. Claude produced a careful-looking analysis faster than I could read it. The mistakes were the kind a good referee catches, and catching them took a second, independent pass that cost very little. I'd now treat that second pass as part of the analysis.

## The part we can't explain

A planet on a simple orbit is a good clock, and this one isn't keeping time. Against a constant period, the four transits arrive 4 minutes late, 23 minutes early, 7 minutes late and 36 minutes late, while each time is measured to within 2 or 3 minutes. The two 2026 transits, on consecutive orbits, are 29 minutes further apart than the period says they should be.

Claude looked hard for an instrumental cause. Changing how the star's baseline brightness is modelled moves no transit time by more than 2.3 minutes. Every version of the public Sector 93 light curve gives the same time to within 2 minutes. One real problem did turn up: the time stamps on the raw full-frame images are corrected for the light-travel time to a reference point on the detector, rather than to this star. Fixing that moved two transits by 48 and 20 seconds, which is far too little to explain the scatter.

That leaves two possibilities. Either something in the data shifts the times in a way we haven't found, or another planet in the system is tugging this one ahead and behind on its orbit (called transit-timing variations). Four transits can't tell those apart, so the paper reports the scatter and doesn't pick one. The next few transits should.

There's also a second, fainter signal in the light curve: a dip of 0.04% every 3.78 days, at the same depth in all eight sectors. On the target it would be a planet about 1.75 times Earth's size. But its transits are too short for a planet crossing the middle of this star, and the same faint neighbour could produce it with a 14% eclipse. I don't know which it is, and the paper doesn't claim it.

## What happens next

Since ExoFOP now wants a refereed paper first, the route to an official candidate goes through a journal, realistically with a professional co-author. In the meantime the work is public in three records on Zenodo: the [preprint](https://doi.org/10.5281/zenodo.23390317), the [data and code](https://doi.org/10.5281/zenodo.23390274), and a [pre-registered test](https://doi.org/10.5281/zenodo.23390325). I borrowed that last idea from someone on r/ClaudeAI who [published a candidate](https://doi.org/10.5281/zenodo.22967456) found the same way a few weeks earlier, along with predictions fixed before the new data arrived.

The pre-registration exists because every choice in the analysis (the period, the transit shape, how the data were cleaned) was made after looking at the transits. However careful you are, it's easy to fool yourself with data you've already seen. Predicting data that doesn't exist yet avoids that. The record contains the predicted times, the code that will analyse the new data (frozen, with checksums), and four possible outcomes:

- **Recovered:** a dip at the predicted time, detected at Δχ² ≥ 25, with a depth between 1,046 and 3,138 ppm.
- **Inconsistent:** a clear dip with the wrong depth.
- **Excluded:** good data covering the window, with no dip anywhere near the expected depth. This would mean the period is wrong.
- **Insufficient data:** anything else.

Two transits are registered. The first is on 19 February 2027, during TESS Sector 114, predicted at BJD 2461455.98 (about 11:30 UTC) ± 30 minutes. The second is on 3 May 2027 at 03:41 UTC ± 31 minutes, when the whole transit is in darkness from La Silla in Chile. I'll publish the result as a new version of the same record whichever way it goes.

## If you have a telescope in the southern hemisphere

The star is at declination −63.5°, so it's out of reach from most of the northern hemisphere; from 20°N it never gets more than about 6° above the horizon. From the south it's easy to find: V = 10.0, about 6° north of Atria (Alpha Trianguli Australis).

The dip is about 2.3 millimagnitudes deep and lasts 5.3 hours. That is at the hard end of ground-based photometry. You'd want a good telescope, a steady sky, and the whole transit plus an hour or two either side. Most amateur setups won't separate the 3.5″ neighbour, so a light curve wouldn't settle the blend question, but it would show whether the transit is there and roughly when it happened. Given the timing puzzle, that's useful.

These are the transits in the next year that happen at night, with the star more than 30° up, from at least one well-known southern site:

| Date (UTC) | Predicted mid-transit (UTC) | Where it is fully or mostly observable |
|---|---|---|
| 27 March 2027 | about 19:35 | Western Australia |
| 3 May 2027 | 03:41 | Chile, southern Brazil |
| 8 June 2027 | about 11:55 | Eastern Australia, New Zealand |
| 14 July 2027 | about 20:05 | South Africa |
| 20 August 2027 | about 04:15 | Chile (about 60% of the transit) |

Each prediction is uncertain by about ±30 minutes, so start early. The next three transits (2 November, 8 December and 14 January) are badly placed for every southern site I checked. The data record lists the next 30 predicted transits. If you try one, I'd very much like to hear about it; my email is in the paper.

## A closing thought

A few years ago this project would have needed someone who already knew TESS's data products, its time systems and its failure modes. This week it needed curiosity, an AI model, and a willingness to have the result audited. I think the result is as valid as one produced the old way, provided it's checked as thoroughly. Producing a confident-looking analysis has become much cheaper, and checking one hasn't kept pace. Pre-registration and independent audits are cheap ways to close that gap, and I'd recommend both to anyone trying this.

I'll write again after 19 February.

---

*Preprint: [doi:10.5281/zenodo.23390317](https://doi.org/10.5281/zenodo.23390317). Data and code: [doi:10.5281/zenodo.23390274](https://doi.org/10.5281/zenodo.23390274). Pre-registered test: [doi:10.5281/zenodo.23390325](https://doi.org/10.5281/zenodo.23390325). Code: [github.com/specbug/ad-astra](https://github.com/specbug/ad-astra).*

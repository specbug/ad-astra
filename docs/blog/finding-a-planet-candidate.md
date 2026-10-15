---
title: "My First Planet (Candidate)"
date: 2026-10-15
description: "Every 36.34 days, a star 218 light-years away dims by 0.21%. I think it's a planet nobody had noticed. How I found it, how I tried to break it, and the February test that could still rule it out."
tags: [astronomy, exoplanets, tess, ai]
---

I started out wanting to be an astronaut. Unusually, the ambition didn't fade as I grew up so much as narrow: first to astrophysics, then to physics, and finally to computer science, each step a little further from the sky and a little closer to what I turned out to be good at. I assumed that was where the astronomy ended.

Last week, computer science brought me back to it. There is a star in the southern constellation Triangulum Australe, 218 light-years away, that dims by 0.21% for a little over five hours every 36.34 days. As far as I can tell, nobody had noticed until I found it in public NASA data. If the dip is what I think it is, it's a planet slightly smaller than Neptune.

It is still only a candidate. On 19 February 2027 the TESS space telescope will look at the star again, and I have already published the test that will decide whether the planet survives that observation. It may not.

## Where it was hiding

TESS, the Transiting Exoplanet Survey Satellite, looks for transits: the small, regular dips in a star's brightness when a planet crosses in front of it. It photographs the sky in strips called sectors, about 27 days each, and then moves on. A typical star gets one sector, then nothing for a year or two, then another.

So a planet with a 36-day orbit shows at most one transit per sector. Planet searches work by finding dips that repeat, and a single dip in 2021 and another in 2025 tend never to be connected. My star, TIC 402710280, had a second disadvantage. TESS has only ever recorded it in its full-frame images, the wide pictures of each camera's whole field, and never at the two-minute cadence that the mission's main search pipeline processes. Long-period planets around stars like this are the ones most likely to have slipped through, so that is where I looked.

## Ten million light curves

A light curve is a star's brightness over time, and the public TESS archive holds a great many of them. I searched 10,324,827, sliding box-shaped templates of ten durations (45 minutes to 18 hours) along each one and flagging any dip more than six times the local noise.

A search that size finds a lot of rubbish, and the hard part is knowing how much. My favourite trick from the whole project is the answer to that: run every search a second time on the same light curves turned upside down. An inverted light curve contains no transits, since its dips have become bumps, so everything the search "finds" there is a false alarm. That gives you a measured false-alarm rate instead of a hopeful guess.

After removing glitches shared by many stars at once, dips at the edges of data gaps, and shapes that don't look like transits, 2,263 stars still had at least one convincing single dip. I checked the 142 brightest of them one by one, down to the pixels and against catalogues of known planets and binary stars. Twenty-seven had no known explanation. For those I made my own light curves from TESS's 2026 images, which no public pipeline had processed yet, and looked for more transits.

## Two more dips

In the public data, TIC 402710280 had two dips, one on 19 June 2021 and one on 11 June 2025: four years apart, almost to the week, with the same depth and the same duration.

Two transits don't give you a period. They give you a list of possible periods: 1,453.58 days, or half that, or a third, or 1,453.58 divided by any whole number. TESS also watched the star in 2023, but the transit it should have caught then landed in data swamped by sunlight scattered off the Earth and Moon.

Then my photometry of the 2026 images showed two more dips, on 21 February and 29 March 2026, 36.36 days apart. This is where I stopped being calm. Those two transits are on consecutive orbits, with data in between showing nothing, so the period can't be a multiple of 36.36 days. It could still be a fraction of it, like 18.2 or 12.1 days, so I checked every fraction down to a fifth. Each one predicts transits at times when TESS was watching and saw nothing; the weakest of those non-detections is at 18.6σ. One period fits everything: 36.33986 ± 0.00051 days, pinned down to about 44 seconds. The four years between the first two dips turned out to be exactly 40 orbits.

![The four transits](fig1.png)
*The four transits, with the unusable 2023 window on the right. Grey points are individual measurements, blue points are 20-minute averages, and the red curve is the transit model.*

## What it might be like

The depth of a transit tells you how big the planet is compared with its star. This one blocks 2,092 ± 47 parts per million of the light, so the planet's radius is 4.2% of the star's. The star is a K dwarf, a little smaller, cooler and more orange than the Sun, which makes the planet 3.75 ± 0.22 times the size of Earth. Neptune is 3.88.

It would orbit at 0.21 AU, about half of Mercury's distance from the Sun, with a year five weeks long. It gets about ten times the sunlight Earth does, enough for an equilibrium temperature of around 500 K (roughly 225 °C). A planet this size probably has a thick atmosphere, so that number is only a rough guide to what it's like there.

Planets between Earth and Neptune in size are the most common kind found so far, and the Solar System has none, so nobody really knows what they're like inside. A typical planet of this radius would weigh about 14 Earths and tug its star back and forth at about 2.8 metres per second, which current spectrographs can measure for a star this bright. That would give the planet's density, the first real clue to whether it's mostly gas, water or rock. I'd love for someone to do it.

There are plenty of known planets like this. What makes this one useful is the combination: a month-long orbit around a bright, nearby star, where follow-up measurements are practical.

## Trying to kill it

Most signals that look like planets in TESS data turn out to be something else, so I spent more time trying to break this one than I did finding it.

The usual impostor is a blended eclipsing binary: two stars eclipsing each other somewhere near your target, whose deep eclipses are diluted by your star's light into a shallow, planet-sized dip. TESS pixels are 21 arcseconds across, and plenty of stars can hide in one.

The first test is where the missing light comes from. Subtract an image taken during a transit from one taken just outside it, and what's left is a picture of the dimming. Across all four transits, it sits 1.5 ± 1.3 arcseconds from my star. Before trusting that, I ran the same method on 40 confirmed TESS planets. It put 30 on their host star, couldn't decide or failed on the other 10, and never put a planet on the wrong star.

That rules out every catalogued star bright enough to fake the signal except one: a faint background star 3.5 arcseconds away, about 350 times fainter than mine and roughly 20 times farther off. It's too close for the image test to separate. But to cause a 0.21% dip in the combined light it would have to lose about 74% of its own light in each eclipse, and an eclipse that deep between two ordinary stars has a different shape from the one I see. I fitted the transits with every eclipse such a pair could produce, and all of them fit worse than a planet by a wide margin (Δχ² ≥ 142, roughly 12σ). A ground-based telescope can easily separate two stars 3.5 arcseconds apart, so one light curve of a future transit would settle it.

My star could itself be an eclipsing binary. If it were, there would usually be a second, shallower eclipse half an orbit later. There isn't one, down to 77 parts per million, and the Gaia satellite sees the star as a single, well-behaved point of light.

The standard tool for weighing all this is TRICERATOPS, which estimates the probability of every false-positive scenario given the data and the stars nearby. It gives 0.003. The usual threshold for calling a candidate a validated planet is 0.015, but validation also needs a high-resolution image ruling out a companion within an arcsecond of the star. Nobody has taken one yet, so for now it stays a candidate.

## Getting audited

Once I was confident, I had two independent audits run on the analysis. The first was deflating to read, mostly because it was right. I had counted a fifth transit in the scattered-light data from 2023, which I shouldn't have, and it had dragged the measured depth down; without it the planet grew from 3.4 to 3.75 Earth radii. I had stacked the transits on a strictly regular period, which smeared their shape, because, as it turns out, they don't arrive on schedule. I had understated the period's uncertainty by a factor of seven. And my predicted transit times, labelled UTC, were in a different time scale and about six minutes off.

Fixing all of it took less than a day, and the candidate came out stronger. The corrected transit shape matches the density expected for this kind of star, and the false-positive probability fell from 0.015 to 0.003. The second audit had measured the transit times with its own photometry and found the same pattern of early and late arrivals I had. In the end I was glad to have been wrong where only two auditors could see it.

## The part I can't explain

A planet on a simple orbit keeps very good time, and this one doesn't. Against a constant period, the four transits arrive 4 minutes late, 23 minutes early, 7 minutes late and 36 minutes late, while each time is measured to within 2 or 3 minutes. The two 2026 transits, on consecutive orbits, are 29 minutes further apart than they should be.

I looked hard for an instrumental cause. Changing how the star's baseline brightness is modelled moves no transit by more than 2.3 minutes, and every version of the public 2025 light curve gives the same time to within 2 minutes. I did find one real problem: the raw TESS images are time-stamped for a reference point on the detector, not for this particular star. Correcting it moved two transits by 48 and 20 seconds, nowhere near enough.

So either something in the data shifts the times in a way I haven't found, or another planet in the system is tugging this one ahead and behind along its orbit. Astronomers call the second possibility transit-timing variations, and it would be by far the more exciting answer. Four transits can't tell the two apart. The next few will.

There's also a second, much fainter signal in the light curve: a 0.04% dip every 3.78 days, at the same depth in all eight sectors. If it comes from my star, it's a second planet about 1.75 times the size of Earth, on a much tighter orbit. But its transits are suspiciously short for a planet crossing the middle of this star, and the same faint neighbour could produce it with a 14% eclipse. I don't know which it is, so I haven't claimed it.

## The bet

Every choice in my analysis, from the period to the transit shape to how the data were cleaned, was made after looking at the transits. However careful you are, it's easy to fool yourself with data you've already seen. The remedy is to predict data that doesn't exist yet and write down in advance what would count as being wrong.

So before publishing, I registered a test. It's a [public, dated record](https://doi.org/10.5281/zenodo.23390325) holding the predicted transit times, the frozen code that will analyse the new data, and four possible outcomes:

- **Recovered:** a dip at the predicted time, detected at Δχ² ≥ 25, with a depth between 1,046 and 3,138 ppm.
- **Inconsistent:** a clear dip with the wrong depth.
- **Excluded:** good data covering the whole window and no dip anywhere near the expected depth. My period would be wrong.
- **Insufficient data:** anything else.

Two transits are registered. The first is on 19 February 2027 at about 11:30 UTC (±30 minutes), while TESS is looking at this part of the sky. The second is on 3 May 2027 at 03:41 UTC (±31 minutes), when the whole transit happens in darkness at La Silla in Chile. I'll publish the result in the same record either way.

After that comes a proper journal paper. ExoFOP, the archive where TESS candidates are listed, now requires a peer-reviewed paper before it adds a new candidate, so that is the route to making this official. Until then, the [preprint](https://doi.org/10.5281/zenodo.23390317) and all the [data and code](https://doi.org/10.5281/zenodo.23390274) are public.

## Go and look at it

If you have a telescope in the southern hemisphere, you can watch this yourself. The star is at declination −63.5°, about 6° north of Atria (Alpha Trianguli Australis), at magnitude 10. From the northern hemisphere you're out of luck: at 20°N it never rises more than about 6° above the horizon.

The dip is about 2.3 millimagnitudes deep and lasts 5.3 hours, which is at the hard end of what ground-based photometry can do. You'd want a good telescope, a steady sky, and the whole transit plus an hour or two either side. Most setups won't separate the 3.5″ neighbour, so a light curve won't settle the blend question, but it would show whether the transit is there and roughly when it happened. Given the timing puzzle, that alone would be valuable.

These are the transits over the next year that happen at night, with the star more than 30° up, from at least one well-known southern site:

| Date (UTC) | Predicted mid-transit (UTC) | Where it's fully or mostly visible |
|---|---|---|
| 27 March 2027 | about 19:35 | Western Australia |
| 3 May 2027 | 03:41 | Chile, southern Brazil |
| 8 June 2027 | about 11:55 | Eastern Australia, New Zealand |
| 14 July 2027 | about 20:05 | South Africa |
| 20 August 2027 | about 04:15 | Chile (about 60% of the transit) |

Each time is uncertain by about ±30 minutes, so start early. The next three transits, on 2 November, 8 December and 14 January, are badly placed for every southern site I checked. The data record lists the next 30. If you try one, please write to me; my email is in the paper.

I'll write again after 19 February, whichever way it goes. Until then the star is six degrees from Atria, dimming once every five weeks.

---

*Thanks to Claude, which wrote and ran most of the code behind this and kept me honest about the word "candidate", and to [Paras Chopra](https://github.com/paraschopra/tess-planet-candidates), whose TESS project made me think an amateur could try this. The idea of pre-registering predictions came from [another amateur's candidate](https://doi.org/10.5281/zenodo.22967456).*

*Preprint: [doi:10.5281/zenodo.23390317](https://doi.org/10.5281/zenodo.23390317). Data and code: [doi:10.5281/zenodo.23390274](https://doi.org/10.5281/zenodo.23390274). Pre-registered test: [doi:10.5281/zenodo.23390325](https://doi.org/10.5281/zenodo.23390325). Code: [github.com/specbug/ad-astra](https://github.com/specbug/ad-astra).*

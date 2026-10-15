---
title: "A Shadow Every Thirty-Six Days"
date: 2026-10-15
description: "How a star 218 light-years away gave up a planet nobody had noticed: four shadows, one dead decoy, a period proven by elimination, and a prediction I cannot take back."
tags: [astronomy, exoplanets, tess]
---

A star in the southern constellation Triangulum Australe loses two-tenths of one percent of its light for five hours and sixteen minutes, and then gets it back. It does this every 36.34 days. It has almost certainly been doing it for billions of years. Until this month, as far as any public record shows, no one had noticed.

The star has no name. It has a catalogue number, TIC 402710280, which is what a star gets when it is one of the hundreds of millions that the TESS space telescope keeps an eye on. The shadow has a number too: TIC 402710280.01, the first candidate planet on that star. It is probably a world a shade smaller than Neptune, wrapped in gas, on an orbit that would fit comfortably inside Mercury's. I have written a [preprint](https://doi.org/10.5281/zenodo.23390317) about it, which you are welcome to read. This post is the story the preprint is not allowed to tell: how you find something in data that ten million other searches have already been through, and how you convince yourself it is real when every instinct you own is lying to you.

## A sky photographed in strips

TESS, the Transiting Exoplanet Survey Satellite, does not stare at stars. It photographs the sky in strips, called sectors, each about 27 days long, then turns to the next strip. A star near the ecliptic might be watched for 27 days and then not again for two years.

Now do the arithmetic that the survey's design hides from you. A planet that orbits its star every 36 days cannot, by construction, be seen to transit twice in one 27-day sector. At most you catch one shadow. The pipelines that find planets in TESS data look for *repetition*. A single dip is just a dip, and the sky is full of dips: cosmic rays, spacecraft wobbles, sunlight scattered off the Earth and Moon, stars that flicker for reasons of their own.

So the planets with long years fall through. Not because anyone is careless, but because a search tuned to find repeating signals inside a sector is mathematically blind to a signal that repeats more slowly than the sector is long. If you want those planets, you have to stitch together sectors observed years apart, and you have to be willing to take single dips seriously.

That was the bet.

## Ten million light curves, and the same ones upside down

A light curve is a star's brightness over time. TESS has produced tens of millions of them. Mine was a search of 10,324,827 public light curves for individual transit-shaped dips, using box-shaped filters of ten durations, from 45 minutes to 18 hours, with a detection threshold set by each light curve's own measured noise.

Here is the problem with any search of that size. It *will* find things. With ten million light curves, a one-in-a-million coincidence happens ten times before breakfast. The question is never "did the detector fire?" but "how often does this detector fire on nothing?"

So every search was also run on light curves flipped upside down. A flipped light curve has no transits in it; a planet's dip becomes a bump. Whatever the pipeline "discovers" in flipped data is a false alarm by construction. That gives you a *measured* false-alarm rate instead of an assumed one. In two-minute TESS data the flipped light curves produced no planet-shaped dips at all. In the full-frame images, the coarser 10- and 30-minute data that cover most of the sky, the flipped rate was low but not zero, mostly because scattered light makes bumps that turn into dips when inverted. That told me exactly how much to distrust a single full-frame dip. The answer was: a lot.

The funnel narrowed quickly. After removing systematics that hit many stars at once, events at the edges of data gaps, and shapes that do not look like transits, 2,263 stars still had single transit-like events. Of those, 142 were bright enough (TESS magnitude below 10.5) to be worth checking by hand, at the pixel level and against every catalogue of known objects. Twenty-seven survived without a known explanation and without a known period.

TIC 402710280 was one of the twenty-seven.

## Two shadows, four years apart

In the public data, the star showed exactly two dips. One in Sector 39, in June 2021. One in Sector 93, in June 2025. Same depth, same duration. The gap between them: 1,453.58 days.

Two transits do not give you a period. They give you a menu. The true period could be 1,453.58 days, or half that, or a third, or any 1,453.58 / *n* for any whole number *n*. There were dozens of candidates on the menu, and the data could not choose among them.

There was a third chance, in Sector 66 in 2023, when a transit might have fallen. For most of the possible periods, it fell straight onto a ramp of scattered sunlight: the background brightened by a factor of five across the window, and the spacecraft's own pipeline flagged the data. That transit, if it exists, is lost in glare. I count it as neither for nor against.

Then I did the thing the public products could not. TESS had observed the star again in 2026, in Sectors 100 through 103, but no public light curves existed yet for those sectors. The raw, calibrated images did. So I made my own light curves from those images, pixel by pixel. In Sectors 100 and 102 there were two more transits, separated by 36.36 days.

## A period by elimination

Two transits on *consecutive* orbits change everything. If you watched the gap between them and saw no transit in it, the period cannot be 72.7 days or 109 days or any multiple of 36.36; those would have needed a missing transit in the middle. The only remaining possibilities are 36.36 days and its fractions: half of it, a third, a quarter, a fifth.

Every fraction makes a prediction. A period of 18.17 days, half the true one, predicts transits at three other times when TESS happened to be watching. They are not there. Each of those three predicted shadows would have been detected at more than 28 sigma, and none appears. A third of the period predicts eight transits in covered data; none appears. A quarter and a fifth predict thirty more between them; the least significant of all those absences is 18.6 sigma.

That is how you prove a period: not by finding it, but by killing every rival. The only period consistent with every hour of data TESS has ever taken of this star is 36.33986 ± 0.00051 days.

## What lives there

Numbers first, then what they mean.

The star is a K dwarf: a little smaller, cooler and redder than the Sun (0.82 solar radii, 5,180 K). It is 66.8 parsecs away, about 218 light-years. Gaia measures its distance to a tenth of a percent. At magnitude 10, it is invisible to the naked eye but easy in a small telescope.

The transit is 2,092 ± 47 parts per million deep. That depth is the ratio of two areas, the planet's silhouette and the star's disc, so it gives the planet's size relative to the star: 4.2% of the star's radius. With the star's size from catalogues, the candidate is 3.75 ± 0.22 times the radius of Earth. Neptune is 3.88. This is a world of Neptune's girth, a shade smaller.

It orbits at 0.206 astronomical units, about half Mercury's distance from the Sun. It moves at about 62 kilometres per second, twice Earth's orbital speed, and its year is five weeks long. It receives ten times the sunlight Earth does. If it reflected nothing and spread that heat evenly, it would sit at about 500 kelvin, around 225 °C.

Planets of this size are among the most common in the galaxy and among the least understood, because nothing like them orbits the Sun. Between about 1.5 and 2 Earth radii there is a known gap in the population, the *radius valley*. Planets above the gap seem to keep thick envelopes of hydrogen and helium; planets below it are bare rock. At 3.75, this one sits comfortably above. The simplest guess for its mass, from the statistical relation between size and mass, is about 14 Earth masses. That mass would make the star wobble towards and away from us at about 2.8 metres per second. Current spectrographs can measure that. It is a measurement someone could make next year.

And it orbits a star that has only ever been watched in full-frame images, which is why no automated pipeline ever stitched its transits together.

## The decoys

Here is the uncomfortable truth about transit astronomy: most things that look exactly like planets are not planets.

The commonest impostor is a *blended eclipsing binary*. Somewhere behind or beside your star, inside the same blurry patch of pixels, two other stars eclipse each other. Their combined light, diluted by your bright star, produces a shallow, planet-like dip. The camera cannot tell you which star dimmed, because TESS pixels are 21 arcseconds across, and at that scale dozens of stars blur into one.

So you go hunting for the decoy.

First, *where* does the light go missing? Subtract an image taken during transit from one taken just outside it, and the difference image is a picture of the missing light. Fit its centre carefully, transit by transit, and the four transits place the source 1.5 ± 1.3 arcseconds from our star. That rules out every catalogued star bright enough to fake the signal, each of which lies 21 sigma or more away. Every one but one.

The survivor is a faint background star 3.5 arcseconds away: magnitude 16.15, roughly 350 times fainter than our star, and, by its Gaia parallax, about 1,450 parsecs away, twenty times farther. At 3.5 arcseconds it sits well inside one TESS pixel, and no centroid can tell the two stars apart.

But arithmetic can. For that faint star to produce a 0.21% dip in the combined light, it would have to lose about 74% of its *own* light during each eclipse. That requires an eclipsing companion nearly as large as the star itself, and eclipses like that have a particular shape: long, sloping ingresses and egresses, nothing like the flat-bottomed U of a small planet crossing a large star. So I fitted the actual transit shape with every eclipse a pair of normal stars could produce at that dilution, across a range of sizes and three different models of how the star's brightness falls off toward its edge. Every such model fits worse than a planet, by a margin statisticians would call overwhelming (Δχ² ≥ 142). The only eclipsing systems that come closer need a hot, compact star eclipsed by a larger, darker one, a pairing that ordinary stars do not make.

The decoy is not dead. It is very badly wounded. A telescope on the ground can separate two stars 3.5 arcseconds apart in ordinary seeing, and a single light curve of a future transit, with an aperture that leaves out the neighbour, would finish it.

Second, is our star itself a binary? An eclipsing binary usually shows a second eclipse half an orbit later, when the brighter star passes in front of the dimmer one. There is none, down to 77 ppm at three sigma, and nothing anywhere between 0.3 and 0.7 of the orbit. Gaia sees a single, well-behaved point of light: its measure of wobble (RUWE) is a textbook 1.00.

Third, put it all together. TRICERATOPS, a tool built for this, weighs every false-positive scenario it knows against the planet hypothesis, given the stars around the target. It returns a false-positive probability of 0.003. The usual bar for calling something a "validated planet" is 0.015.

And yet I call it a candidate, not a planet, because the rule is that you need high-resolution imaging to rule out a companion hiding closer than about one arcsecond, too close for Gaia to resolve. Nobody has taken that image yet. Statistics computed without it are an upper bound on the evidence, not a verdict.

## The part that does not fit

If this were a clean story it would end there. It does not.

A planet on a perfect orbit is a clock. Once you know the period, every transit should arrive on schedule. Ours do not. Against the best straight-line ephemeris, the four transits arrive 4 minutes late, 23 minutes early, 7 minutes late, and 36 minutes late. The scatter is 22 minutes, while each time is measured to within two or three minutes. The two transits of 2026, on consecutive orbits, are 29 minutes further apart than the period says they should be.

There are two kinds of explanation, and they predict very different futures. One is instrumental: something in the data or the processing shifts the apparent times. I tried hard to find it. Changing how the brightness baseline around each transit is modelled moves no time by more than 2.3 minutes. Every flux version in the public Sector 93 light curve gives the same time within 1.9 minutes. While doing this I found and fixed a genuine time-stamp subtlety: the raw TESS images are time-stamped for the middle of the detector, not for the star, which shifted two transits by 48 and 20 seconds. Real, worth fixing, and nowhere near 22 minutes.

The other explanation is gravitational: another planet in the system tugging this one forward and back along its orbit, the effect called transit-timing variations. That would be the more exciting answer. Four transits cannot distinguish the two, so the preprint reports the scatter and refuses to choose. The next transits will.

There is one more loose thread. With the 36-day transits masked out, the light curve contains a second, much fainter signal: a 0.04% dip every 3.78 days, seen at the same depth in all eight sectors and independently in the old and new data. If it is on our star, it is a second planet about 1.75 times Earth's size, on a tighter orbit. But its transits are too short for a central crossing of this star, unless it is grazing the limb. And the faint neighbour, the same one, could produce it with a 14% eclipse of its own. I do not know which, and the preprint says so.

## Errors, mine

Science is supposed to be self-correcting. In practice the correcting is done by people, and it is done to your own work before you show it to anyone else, or it is done to you afterwards in public.

I am an amateur. I built this with Claude, an AI system from Anthropic, as a collaborator: it wrote most of the pipeline and ran most of the analysis, while I set the direction and checked the results. Then I had the work audited independently. The audits found real mistakes. An early version printed transit times in barycentric time, the clock of the solar system's centre of mass, while labelling them as local time at La Silla; they were about six minutes off. Another put the Sector 66 window, the one lost to scattered light, into the fit, which quietly distorted the transit's shape. A date conversion was wrong by three months. Every one is fixed, and every fix moved a number you can check.

I mention this because it is the most transferable lesson here. The danger in an amateur discovery is not a lack of talent. It is that nobody is paid to tell you that you are wrong. You have to arrange for that yourself, deliberately, before the result has your name on it.

## The bet

Every analysis like this has a quiet weakness. I chose the period, the transit shape and the way the data were cleaned *after* looking at the transits. However careful you are, you can always fool yourself with data you have already seen.

There is exactly one cure: predict data that does not exist yet, and say in advance what would count as being wrong.

So, before publishing, I registered a test. It is a [public, timestamped record](https://doi.org/10.5281/zenodo.23390325), and its rules cannot be changed. It fixes the predicted times of the next transits and the exact code that will analyse them. It also fixes the four possible verdicts: recovered, inconsistent, excluded, or insufficient data. And it sets the thresholds that separate them: a depth between 1,046 and 3,138 ppm and a detection at Δχ² ≥ 25, or no dip deeper than 1,046 ppm anywhere in the window.

Two transits are registered:

- **19 February 2027**, while TESS is watching this patch of sky again in Sector 114. Predicted mid-time BJD 2461455.9785, give or take 30 minutes.
- **3 May 2027**, at 03:41 UTC as seen from La Silla in Chile, where the whole transit falls in darkness. Give or take 31 minutes.

The uncertainty includes the unexplained 22-minute scatter, so even a late or early transit lands inside the window. But if TESS stares at the right patch of sky on 19 February and the star does not dim, the period is wrong, and I will say so in a new version of the same record. I committed to publish the result either way.

## What is left to do

The candidate needs three things to become a planet. It needs a high-resolution image, to exclude a companion hiding within an arcsecond. It needs a ground-based light curve that separates our star from the faint neighbour, which would settle both the 36-day and the 3.78-day signals at once. And it needs a radial-velocity measurement of that 2.8-metre-per-second wobble, which would give its mass and therefore its density. That density is the number that says whether this is a puffy gas world or something stranger.

None of that is in my reach. All of it is in reach of people with telescopes, and the preprint, the data and the code are public so that they can start. Once a refereed paper exists, the candidate goes to the community list of TESS planet candidates, which is where follow-up observers look.

Until then, there is a small, warm world going round an orange star every five weeks, and its shadow has crossed the star's face tens of billions of times without anyone writing it down. The light that carried the 2021 shadow left that star when Napoleon was still First Consul. The light of the next one is already on its way. It will reach Earth on the nineteenth of February.

Either the star dims or it does not. I have written down which.

---

*Preprint: [doi:10.5281/zenodo.23390317](https://doi.org/10.5281/zenodo.23390317). Data and code: [doi:10.5281/zenodo.23390274](https://doi.org/10.5281/zenodo.23390274). Pre-registered test: [doi:10.5281/zenodo.23390325](https://doi.org/10.5281/zenodo.23390325).*

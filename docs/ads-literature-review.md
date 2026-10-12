# Scoped ADS literature review

ADS access succeeded on 2026-10-12. This is an additional literature check of the first pilot's two variability leads and their known neighbors. It does not change the photometric findings or establish novelty. Earlier reports saying ADS was omitted describe the earlier review, before this retrieval.

The five initial requests all returned HTTP 200, with all reported hits returned in the first page. Their declaration, exact query list, timestamps, public responses, and response hashes are preserved under `data/ads-first-pilot/` and `data/ads-first-pilot-queries.json`. The query-list SHA-256 is `c7632fd5dcba6d9bad077f9971442e6469bc99ba4dc00d46282d88e68a8e778e`. The retrieval script reports that no credential was saved. The response hashes below were independently verified against the stored bytes.

| Query | Reported hits | Public response SHA-256 |
| --- | ---: | --- |
| `full:"149392671" OR full:"8883-946-1" OR full:"4757558780144090368"` | 0 | `2af5ba886b083125917ef0f07bf6a83501f942be7adc51af7482ccf81ec16786` |
| `full:"150104580" OR full:"8901-833-1" OR full:"5476797645662534144"` | 5 | `a75f333fd704c638d50ecd282b40ba8e4a52308b61deef978a3aaac40bf8a448` |
| `full:"LMC-RRLYR-34351" OR full:"4757558810204649856"` | 0 | `6453a596b98f77d9c0d21572913a136a8de2090286049c07870abc51d97f0ca5` |
| `full:"LMC-RRLYR-37553" OR full:"5476797641361944064"` | 0 | `8fcc2834ef36e862f158cc6aac42b2d002aee53f38e20cf43e0853dc2a5909ef` |
| `bibcode:2022A&A...658A..22R` | 1 | `f88acddd1922ea504e4951b49b1fa9436a5348fd4516e2e49092cdde78dbc5f7` |

## Assessment of the five TIC 150104580 search hits

The returned records do not identify TIC 150104580 as an astronomical object in their titles. They comprise two high-energy theory papers and three Earth-weather satellite records; two of the latter have the same title and DOI, representing duplicate-year records of one paper. The stored responses do not include abstracts for these five records. Their titles and subject matter supply no evidence of a previously reported stellar period or planetary candidate for this target.

| ADS bibcode | Returned title |
| --- | --- |
| `2024JHEP...08..200D` | Chaotic and thermal aspects in the highly excited string S-matrix |
| `2020JMetR..33.1113Y` | Capability of Fengyun-3D Satellite in Earth System Observation |
| `2019JMetR..33.1113Y` | Capability of Fengyun-3D Satellite in Earth System Observation |
| `2019AdAtS..36.1027Z` | Latest Progress of the Chinese Meteorological Satellite Program and Core Data Processing Technologies |
| `2017JHEP...12..073H` | Dissimilarities of reduced density matrices and eigenstate thermalization hypothesis |

The OR query does not identify which branch produced each hit. Bare numeric identifiers and a hyphenated alias without the `TYC` catalog prefix are insufficiently specific for interpreting these results. Punctuation normalization or unrelated occurrences of the numbers could explain the hits, but the saved responses do not establish the indexing mechanism. These records are therefore classified as unrelated apparent search matches, rather than evidence about the target.

## The white-dwarf system is already published

The direct bibcode query recovered [Raddi et al. (2022), *Kinematic properties of white dwarfs. Galactic orbital parameters and age-velocity dispersion relation*](https://ui.adsabs.harvard.edu/abs/2022A%26A...658A..22R/abstract), DOI `10.1051/0004-6361/202141837`. Its abstract describes 3,133 white dwarfs, including common-proper-motion pairs with nondegenerate companions, and explicitly points to the full CDS tables.

The target-level identity comes from the separately retrieved VizieR table, not from an identity match in this ADS abstract. As documented in [the stellar metadata report](qlp-stellar-metadata.md), table `J/A+A/658/A22/table2`, system 2655, pairs white dwarf Gaia DR3 4757558814501499392 (`SourceID1`) with the bright TIC 149392671 primary, Gaia DR3 4757558780144090368 (`SourceID2`). This confirms an existing published system. The known RR Lyrae at 29.25 arcseconds is a different source from the white dwarf at 6.20 arcseconds.

## Stronger follow-up queries

Follow-up should separate aliases so that the source of a hit is clear. Suitable starting queries are `collection:astronomy AND full:"TIC 150104580"` and `collection:astronomy AND full:"TIC 149392671"`, followed by separate full-prefixed Tycho phrases `collection:astronomy AND full:"TYC 8901-833-1"` and `collection:astronomy AND full:"TYC 8883-946-1"`. Prefix variants such as `TIC150104580` and `Gaia DR3 5476797645662534144` warrant separate searches too. Restricting to astronomy reduces incidental matches but can also omit relevant interdisciplinary records, so it is a precision check, not a completeness guarantee.

Zero identifier hits do not mean an object is new. Catalog rows and electronic supplements are not necessarily indexed as article full text, and identifiers may change between catalog releases. In particular, both RR Lyrae neighbors are already present in Gaia/OGLE/VSX metadata despite zero hits in these ADS queries. Coordinate-based bibliography, parent survey papers, supplements, alternate aliases, and the actual scientific claim must be reviewed before any novelty statement.

## Completed first-pilot refinement

All six requests in `data/ads-first-pilot-refined/` returned HTTP 200 and zero reported matches on 2026-10-12. The query-list SHA-256 is `1ab36c4d088dd2aa4e6007f70ee79a27ab6a95ea23fb763451fbffe9878ec123`. The stored response hashes were independently verified. These are scoped identifier nonmatches, not novelty clearance.

| Exact query | Public response SHA-256 |
| --- | --- |
| `collection:astronomy AND full:"TIC 149392671"` | `4b30bf43442a26c3413b5658fad246c02219a7d4f5b1120bc0ba21b37e5fc513` |
| `collection:astronomy AND full:"TIC 150104580"` | `a3b6ff76158b084b75f8cf42db85b5555eef40f1d1d6da8375edf8a9694f81f3` |
| `collection:astronomy AND full:"TYC 8883-946-1"` | `a946db98c19b3d80e5c001caa4e5bf5329bb2a5c5fa423bf9006a27d0eb80228` |
| `collection:astronomy AND full:"TYC 8901-833-1"` | `66e13d6db71aa26ee506ac05ed85634f8b8b0376954624303f7268af34c9147b` |
| `collection:astronomy AND full:"OGLE LMC-RRLYR-34351"` | `b2500158cc3fae2bc337432fe4e3d02ab76483a02af8c4d24908d476165f17bb` |
| `collection:astronomy AND full:"OGLE LMC-RRLYR-37553"` | `8db6bd57f325b5fa58ccd4b44c48effaa52ab19b0b8dc5051d4f67de3656fd0e` |

The full-prefixed Tycho searches remove the ambiguity of the initial bare alias query. The two full OGLE-identifier queries also return zero even though both RR Lyrae are independently established in the catalog evidence. This is a direct example of the limits of these exact-phrase ADS queries: known objects can lack searchable identifier matches. It does not determine whether punctuation, phrase conventions, supplement indexing, or another coverage limitation caused these particular nonmatches.

## Second-pilot lead queries

All six requests in `data/ads-second-pilot/` returned HTTP 200 and zero reported matches on 2026-10-12. The query-list SHA-256 is `b0c5523e8947aeb93f23514826c64422ad1ba5483e796bb9a8ce0a073f3bd625`. The stored response hashes were independently verified. These are scoped identifier nonmatches, not novelty clearance.

| Exact query | Public response SHA-256 |
| --- | --- |
| `collection:astronomy AND full:"TIC 355767420"` | `67a7f36a7cf66f12870f4cbcb851e2ffde1a7b48cdafc8c06d67d5f5f1f4f0bf` |
| `collection:astronomy AND full:"2MASS J00173070-6005040"` | `c82d4b0cac64a77878dd4f019eb32e19ea4d3a1faeee2c64048b6ad1f647fe74` |
| `collection:astronomy AND full:"TIC 201308024"` | `021c89ceb1615455f9dc1020445284ea7b6f038e4a1f320b3d8ef6a38be6a800` |
| `collection:astronomy AND full:"2MASS J00151270-5942190"` | `84847e7b01ee1b156fcfb93d195609f5a376d59f6a5b5db64d4e997e5ed3f7ff` |
| `collection:astronomy AND full:"TIC 281593180"` | `680dee852df4ddcb1854bc8b6ec9be95347b5f7c0d596faee751bde86df5e64d` |
| `collection:astronomy AND full:"2MASS J00321607-5732498"` | `06b1cb918b74d0d526be93fa44433778db34838288a10655c7aead7ce841ce8c` |

These six queries cover the TIC and 2MASS aliases of second-pilot inspection leads TIC 355767420, TIC 201308024, and TIC 281593180. A nonmatch neither establishes an intrinsic signal nor excludes prior reporting under other aliases or in electronic supplements. These searches are an early bibliography check; scientific validation and review of the specific proposed claim remain necessary.

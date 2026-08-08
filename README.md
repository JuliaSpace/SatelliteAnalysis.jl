SatelliteAnalysis
=================

[![CI](https://img.shields.io/github/actions/workflow/status/JuliaSpace/SatelliteAnalysis.jl/ci.yml?style=flat-square&logo=githubactions&logoColor=white&labelColor=475569&label=CI)](https://github.com/JuliaSpace/SatelliteAnalysis.jl/actions/workflows/ci.yml)
[![Codecov](https://img.shields.io/codecov/c/github/JuliaSpace/SatelliteAnalysis.jl?token=62H85L1AHF&style=flat-square&logo=codecov&logoColor=white&labelColor=475569)](https://codecov.io/gh/JuliaSpace/SatelliteAnalysis.jl)
[![docs-stable](https://img.shields.io/badge/docs-stable-16A34A?style=flat-square&logo=gitbook&logoColor=white&labelColor=475569)][docs-stable-url]
[![docs-dev](https://img.shields.io/badge/docs-dev-D97706?style=flat-square&logo=gitbook&logoColor=white&labelColor=475569)][docs-dev-url]
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495D1?style=flat-square&logo=julia&logoColor=white&labelColor=475569)](https://github.com/invenia/BlueStyle)
[![License](https://img.shields.io/github/license/JuliaSpace/SatelliteAnalysis.jl?style=flat-square&logo=readme&logoColor=white&labelColor=475569&color=0284C7)](https://github.com/JuliaSpace/SatelliteAnalysis.jl/blob/main/LICENSE.txt)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.10501188-DB2777?style=flat-square&logo=doi&logoColor=white&labelColor=475569)](https://zenodo.org/doi/10.5281/zenodo.10501188)

This package contains several functions to perform analysis related to satellites. Those
functions were split from the package
[SatelliteToolbox.jl](https://github.com/JuliaSpace/SatelliteToolbox.jl).

## Installation

This package can be installed using:

```julia-repl
julia> using Pkg
julia> Pkg.add("SatelliteAnalysis")
```

## Makie Theme

**SatelliteAnalysis.jl** ships a [Makie](https://makie.org) theme with coordinated dark and
light variants, designed for presentation slides and reports. It is provided by a package
extension that is loaded automatically when a Makie backend is available:

```julia
using CairoMakie   # Or GLMakie, WGLMakie, etc.
using SatelliteAnalysis

set_theme!(makie_theme(:dark))   # Use makie_theme() for the light variant.
scatter(rand(100))
```

The function `makie_palette` returns the categorical palette used by the theme, and the
keywords `fontscale` and `mono_ticklabels` adjust the font sizes and the tick label font.
See the [documentation][docs-makie-url] for a gallery and the complete API reference.

The theme bundles the IBM Plex Sans and IBM Plex Mono fonts, copyright © IBM Corp. and
distributed under the [SIL Open Font License 1.1](https://openfontlicense.org). The license
texts are available in `assets/fonts/`.

## Documentation

For more information, see the [documentation][docs-stable-url].

[docs-dev-url]: https://juliaspace.github.io/SatelliteAnalysis.jl/dev
[docs-makie-url]: https://juliaspace.github.io/SatelliteAnalysis.jl/dev/man/makie_theme/
[docs-stable-url]: https://juliaspace.github.io/SatelliteAnalysis.jl/stable

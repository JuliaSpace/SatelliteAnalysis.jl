SatelliteAnalysis.jl Changelog
==============================

Version 0.4.0
-------------

- ![Feature][badge-feature] We added the function `decay_analysis` to estimate the orbital
  decay lifetime of a satellite considering averaged perturbations (zonal harmonics with a
  J₂² correction, third bodies, atmospheric drag, and solar radiation pressure gated by the
  Earth shadow). The analysis is available when OrdinaryDiffEq.jl is loaded and returns a
  `DataFrame` with the mean orbital element evolution. (PR [#16][gh-pr-16])
- ![Feature][badge-feature] We added a Makie theme extension, moved from
  **SatelliteToolbox.jl**, where it was never released. Loading a Makie backend together
  with **SatelliteAnalysis.jl** provides the non-exported functions
  `SatelliteAnalysis.makie_theme` and `SatelliteAnalysis.makie_palette`, which style plots
  with coordinated dark and light variants using bundled IBM Plex fonts.
- ![Feature][badge-feature] The plotting functions `plot_world_map`, `plot_ground_track`,
  and `plot_ground_facility_visibility_circles` now apply the SatelliteAnalysis.jl Makie
  theme automatically and gained the keyword `theme` to select the variant (`:light` or
  `:dark`).
- ![Feature][badge-feature] We added the function `plot_decay_analysis`, available when
  **Makie.jl** is loaded, that plots the mean apogee and perigee altitude evolution
  computed by `decay_analysis` together with an information panel showing the satellite
  mass, the satellite mean area, and the estimated time to reenter. To support it,
  `decay_analysis` now records the metadata `Satellite Mass`, `Satellite Mean Area`, and
  `Terminate Altitude` in the output `DataFrame`.
- ![Bugfix][badge-bugfix] The functions `ground_repeating_orbit_adjacent_track_angle` and
  `ground_repeating_orbit_adjacent_track_distance` were ignoring the keyword `we`, and the
  functions `design_sun_sync_ground_repeating_orbit` and `sun_sync_orbit_inclination` were
  ignoring the keyword `R0` in part of the algorithm.
- ![Enhancement][badge-enhancement] We removed an allocation proportional to the analysis
  duration in `ground_track`.
- ![Info][badge-info] We fixed several typos and documentation errors, including wrong
  documented keyword defaults and signatures.
- ![Enhancement][badge-enhancement] We highly improved the `decay_analysis` performance
  (about 20x) by tuning the default integrator configuration for lifetime estimation and
  reducing the cost of the right-hand side. The integrator can now be configured through
  the keywords `solver`, `reltol`, and `abstol`.
- ![Enhancement][badge-enhancement] The keywords `Ap` and `F107` in `decay_analysis` now
  accept either a constant value or a function of time with the signature
  `(jd_utc::Number) -> Number`, allowing time-varying space index profiles. By default,
  the indices are obtained from **SpaceIndices.jl**, which must be initialized with
  `SpaceIndices.init()`.
- ![Info][badge-info] We added a test that validates the averaged decay dynamics against a
  full osculating (Cowell) reference propagation, bounding the neglected couplings.
- ![Info][badge-info] The plotting extension `SatelliteAnalysisPlottingExt` was renamed to
  `SatelliteAnalysisGeoMakieExt` since it only requires **GeoMakie.jl**. Extension names
  are not public API, so this change is internal.
- ![Info][badge-info] The decay analysis extension is now triggered by the package
  **OrdinaryDiffEqAdamsBashforthMoulton.jl**, which provides the default solver `VCABM`,
  adding support for **OrdinaryDiffEq.jl** v7, where this solver is no longer bundled.
  Loading **OrdinaryDiffEq.jl** v6 still activates the extension because it depends on
  that package, but users of **OrdinaryDiffEq.jl** v7 or newer must load
  **OrdinaryDiffEqAdamsBashforthMoulton.jl** explicitly.

Version 0.3.10
--------------

- ![Enhancement][badge-enhancement] Optimize analysis hot paths while preserving numerical
  precision.
- ![Bugfix][badge-bugfix] Fix world-map ticks to display the Equator.
- ![Bugfix][badge-bugfix] Harden frozen- and Sun-synchronous-orbit solvers, including degree
  handling, safeguards, and keyword compatibility.
- ![Info][badge-info] Expand orbit and analysis regression coverage.

Version 0.3.9
-------------

- ![Enhancement][badge-enhancement] Improve the default reduction function for ground
  facilities. ([#10][gh-pr-10])
- ![Enhancement][badge-enhancement] The GeoMakie.jl extension can now precompile. (PR
  [#9][gh-pr-9])

Version 0.3.8
-------------

- ![Info][badge-info] We updated the compat versions of the dependencies.
- ![Info][badge-info] The package is now being tested against Julia 1.10.

Version 0.3.7
-------------

- ![Bugfix][badge-bugfix] The interpretation of `duration` in `ground_track` was not
  correct. It must be the number of seconds the ground track will be computed **after** the
  initial time.

Version 0.3.6
-------------

- ![Feature][badge-feature] We added an in-place version of all plotting functions. Thus,
  the user can plot the analysis on top of existing figures, leading to better analysis
  options.
- ![Feature][badge-feature] We added the function `plot_world_map` to create a figure with
  only the world map.

Version 0.3.5
-------------

- ![Enhancement][badge-enhancement] We added the compatibility to GeoMakie 0.7.

Version 0.3.4
-------------

- ![Enhancement][badge-enhancement] The algorithm to compute the access to ground facilities
  now uses the local reference frame (NED) to compute the elevation angle, leading to a
  better precision. (PR [#5][gh-pr-5])

Version 0.3.3
-------------

- ![Bugfix][badge-bugfix] The function `ground_facility_gaps` was not taking into account
  the `step` parameter. (Issue [#3][gh-issue-3])

Version 0.3.2
-------------

- ![Enhancement][badge-enhancement] Minor source-code updates.
- ![Enhancement][badge-enhancement] We updated the documentation.

Version 0.3.1
-------------

- ![Feature][badge-feature] The package now contains extensions that are loaded when
  [GeoMakie.jl](https://github.com/MakieOrg/GeoMakie.jl) is loaded. In this case, we added
  functions to plot the ground tracks and the ground facility visibility circles.
- ![Feature][badge-feature] We added support to download the GeoJSON file with the
  countries' polygons that can be used, for example, to plot information in the world map.

Version 0.3.0
-------------

- ![BREAKING][badge-breaking] The algorithm to compute the beta angle was improved, but its
  behavior changed. Thus, this modification is breaking.
- ![BREAKING][badge-breaking] The algorithm to compute the eclipse time was simplified,
  changing its API.
- ![BREAKING][badge-breaking] The algorithm related to ground facility accesses and gaps was
  simplified, changing its API.
- ![Feature][badge-feature] We added algorithms to design Sun-synchronous orbits.
- ![Feature][badge-feature] We added a function to compute the ground facility visibility
  circle.
- ![Feature][badge-feature] We added a function to design frozen orbits.
- ![Feature][badge-feature] We added a function to compute the ground track inclination.
- ![Feature][badge-feature] We highly improved the package documentation.
- ![Enhancement][badge-enhancement] We improved many aspects of the package including
  comments, source code organization, allocations, and others.

Version 0.2.0
-------------

- ![Feature][badge-feature] The lightning analysis (eclipse) was added.
- ![Feature][badge-feature] Initial support for ground facility analysis.

Version 0.1.0
-------------

- Initial version.

[badge-breaking]: https://img.shields.io/badge/Breaking-DC2626?style=flat-square
[badge-deprecation]: https://img.shields.io/badge/Deprecation-D97706?style=flat-square
[badge-feature]: https://img.shields.io/badge/Feature-16A34A?style=flat-square
[badge-enhancement]: https://img.shields.io/badge/Enhancement-0284C7?style=flat-square
[badge-bugfix]: https://img.shields.io/badge/Bugfix-DB2777?style=flat-square
[badge-info]: https://img.shields.io/badge/Info-475569?style=flat-square

[gh-issue-3]: https://github.com/JuliaSpace/SatelliteAnalysis.jl/issues/3

[gh-pr-5]: https://github.com/JuliaSpace/SatelliteAnalysis.jl/pull/5
[gh-pr-9]: https://github.com/JuliaSpace/SatelliteAnalysis.jl/pull/9
[gh-pr-10]: https://github.com/JuliaSpace/SatelliteAnalysis.jl/pull/10
[gh-pr-16]: https://github.com/JuliaSpace/SatelliteAnalysis.jl/pull/16

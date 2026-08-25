SatelliteAnalysis.jl Changelog
==============================

Version 0.4.0
-------------

- ![Feature][badge-feature] We added the function `decay_analysis` to estimate the orbital
  decay lifetime of a satellite considering averaged perturbations (zonal harmonics with a
  J₂² correction, third bodies, atmospheric drag, and solar radiation pressure gated by the
  Earth shadow). By default, the input elements are treated as mean elements with respect
  to the averaged dynamics, following the same convention of semi-analytical tools such as
  STELA; the keyword `input_type` selects between `:mean` and `:osculating` input
  elements, converting the latter to mean elements before the propagation. The analysis
  is available when OrdinaryDiffEq.jl is loaded and returns a `DataFrame` with the mean
  orbital element evolution. (PR [#16][gh-pr-16])
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
  mass, the satellite mean area, and the estimated time to reenter. Keywords allow adding
  a mission name above the title, showing the absolute dates, and plotting the daily and
  the 81-day average F10.7 indices using a twin y-axis, extracted from the column
  `space_indices` through the keywords `f107_getter` and `f107_avg_getter`. To support it,
  `decay_analysis` now records the metadata `Satellite Mass`, `Satellite Mean Area`, and
  `Terminate Altitude` in the output `DataFrame`.
- ![Feature][badge-feature] We added the keyword `atmospheric_model` to `decay_analysis`,
  allowing the user to select the atmospheric density model used by the drag computation.
  It accepts any callable object, including callable structures carrying their own state,
  receiving the named tuple with the space indices provided by the keyword
  `space_indices`. By default, the analysis uses the NRLMSISE-00 model, which consumes the
  daily F10.7, the 81-day average F10.7, and the daily Ap. The keyword
  `atmospheric_model_name` overrides the model name recorded in the output metadata.
- ![Feature][badge-feature] We added the macros `@decay_analysis__jacchia77` and
  `@decay_analysis__jr1971` that provide keyword sets for `decay_analysis` selecting the
  Jacchia 1977 and the Jacchia-Roberts 1971 atmospheric models instead of the default
  NRLMSISE-00. Each macro expands to the keywords `atmospheric_model`,
  `atmospheric_model_name`, and `space_indices` in the keyword section of the call, and
  keywords passed after the macro override the ones it provides. The models consume the
  space indices `f107`, `f107_avg`, and `kp`, and the default source selected by the
  macros provides the observed values (space indices `F10obs`, `F10obs_avg_last81`, and
  `Kp_daily`), falling back to the predicted F10.7 and to Kp = 7 / 3 (equivalent to
  Ap = 9, as in STELA) outside the observed timespans.
- ![Feature][badge-feature] We added the keyword `verbose` to `decay_analysis`. When
  enabled, a progress interface is shown in `stderr` during the numerical integration: in
  interactive terminals, a live panel with a progress bar, the current perigee and apogee
  altitudes, the elapsed model time, and the elapsed wall time; otherwise, plain progress
  lines at every 10%. A summary line with the outcome and the wall time is printed at the
  end. Enabling the interface does not change the analysis result.
- ![Feature][badge-feature] We highly improved the `plot_decay_analysis` output for
  reports: the information panel now shows a card with the analysis assumptions
  (atmospheric model and drag and SRP coefficients), a
  dashed line marks the terminate altitude, the reentry is annotated next to the reentry
  marker, the keyword `show_dates` adds the absolute dates (the analysis timespan in the
  subtitle and the estimated reentry date in the panel and in the annotation), and the
  axis labels use human-readable unit names. The new keywords `title`, `subtitle`,
  `show_assumptions`, `show_dates`, `show_reentry_callout`, `fontscale`,
  `mono_ticklabels`, `panel_width`, `xlims`, and
  `ylims` control the figure. To support it, `decay_analysis` now records the metadata
  `Atmospheric Model`, `Drag Coefficient`, `Space Indices Source`, and `SRP Coefficient`
  in the output `DataFrame`.
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
- ![Enhancement][badge-enhancement] The keyword `space_indices` in `decay_analysis`
  provides the space indices required by the atmospheric model as a named tuple. It
  accepts either a constant `NamedTuple` or a function of time with the signature
  `(jd_utc::Number) -> NamedTuple`, allowing time-varying space index profiles, and the
  named tuple is passed to the atmospheric model and recorded in the column
  `space_indices` of the output `DataFrame`. By default, it provides the observed daily
  F10.7 (space index `F10obs`), the observed last-81-day average F10.7 (space index
  `F10obs_avg_last81`), and the observed daily geomagnetic index (space index `Ap_daily`)
  from **SpaceIndices.jl**, falling back to the predicted F10.7 (space index
  `F10predicted`) and to Ap = 9, as in STELA, outside the observed timespans. The required
  space index sets are initialized automatically on first use.
- ![Enhancement][badge-enhancement] We highly reduced the allocations of `decay_analysis`
  (about 95%) by reusing the Legendre buffers of the atmospheric and gravity models, using
  a static state vector in the numerical integration, caching the default gravity model,
  removing redundant conversions and dead code from the averaging loops, and forcing the
  specialization on the model callbacks.
- ![Enhancement][badge-enhancement] We highly reduced the compilation time at the first
  execution of `decay_analysis`: the extension now runs a precompilation workload
  covering the whole analysis pipeline, and the progress callback type no longer depends
  on the keyword `verbose`, sharing the solver specialization between the verbose and
  silent paths. The remaining first-call cost is dominated by loading the default EGM96
  gravity model.
- ![Bugfix][badge-bugfix] The decay analysis passed the Julian date to the gravity model
  where it expects elapsed seconds from the J2000.0 epoch. The error was harmless for the
  default EGM96 model, whose coefficients are static, but it would produce wrong results
  for ICGEM models with time-variable coefficients.
- ![Bugfix][badge-bugfix] The metadata `Description` of the `decay_analysis` output now
  propagates through DataFrame transformations, and the reported mean eccentricity is no
  longer clamped to 1e-6, removing a fake apogee and perigee split of about 13 m for
  circular orbits.
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

SatelliteAnalysis.jl Changelog
==============================

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

[badge-breaking]: https://img.shields.io/badge/BREAKING-red.svg
[badge-deprecation]: https://img.shields.io/badge/Deprecation-orange.svg
[badge-feature]: https://img.shields.io/badge/Feature-green.svg
[badge-enhancement]: https://img.shields.io/badge/Enhancement-blue.svg
[badge-bugfix]: https://img.shields.io/badge/Bugfix-purple.svg
[badge-info]: https://img.shields.io/badge/Info-gray.svg

[gh-issue-3]: https://github.com/JuliaSpace/SatelliteAnalysis.jl/issues/3

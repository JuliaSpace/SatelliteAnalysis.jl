SatelliteAnalysis.jl Changelog
==============================

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

# AGENTS.md

## Package Structure

- SatelliteAnalysis is a Julia package requiring Julia 1.10 or newer (`Project.toml`).
- The module entrypoint is `src/SatelliteAnalysis.jl`; add includes there and preserve its dependency order when new code relies on symbols from another source file.
- Source is organized by analysis domain, with shared helpers in `src/misc/`, ground-facility code in `src/ground_facilities/`, and plotting entrypoints in `src/plotting/`.
- Plotting is an optional `SatelliteAnalysisPlottingExt` extension under `ext/`, loaded only when `GeoMakie` is available. Exercise plotting changes with `GeoMakie` loaded.
- `test/runtests.jl` uses verbose domain-level testsets and includes each test file. Add or update the corresponding focused test file, then wire any new file into `test/runtests.jl`.
- Tests use `[extras]` and the `test` target in `Project.toml`; use `Pkg.test()` rather than a plain project session when test-only dependencies are needed.

## Commands

- Instantiate dependencies: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
- Run the complete test suite: `julia --project=. -e 'using Pkg; Pkg.test()'`
- Build the documentation locally: `julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()' && julia --project=docs docs/make.jl`
- Format all Julia code after installing JuliaFormatter in a non-package environment: `julia -e 'using JuliaFormatter; format(".")'`
- The initial dependency installation and test run can spend several minutes precompiling; use generous timeouts.

## Code Style

- Format with JuliaFormatter's BlueStyle configuration in `.JuliaFormatter.toml`; retain its alignment, whitespace, and long-function-definition settings.
- There is no formatter check in CI. Run formatting deliberately and inspect `git diff` before retaining formatter-only changes.
- Follow the existing docstring convention: document argument units, return values, keyword defaults, and any physical or reference-frame assumptions.
- Preserve the package's explicit units: angles are generally radians, distances/meters and durations/seconds are stated in APIs, and `angle_unit` options must remain explicit.
- Preserve DataFrame schemas and metadata, including documented unit and description metadata; test both values and schema/metadata when changing tabular output.
- Validate invalid physical inputs with clear `ArgumentError`s, matching the nearby API's behavior.

## Behavioral Constraints

- Ground facilities use WGS84 tuples `(latitude [rad], longitude [rad], altitude [m])`.
- Ground-facility access calculations can run concurrently; callbacks such as ECI-to-ECEF conversions must remain thread-safe.
- Access and gap results use UTC `DateTime` boundaries and duration metadata; do not silently change column names, units, or timestamp semantics.
- Numerical orbital calculations need tolerance-based regression tests and explicit angle-unit coverage where applicable.

## CI and Documentation

- CI builds and tests on Julia 1.10 and the latest stable Julia 1.x across supported Ubuntu, macOS, and Windows architecture combinations; standard CI also reports coverage.
- A nightly workflow builds and tests the package on the same platform matrix.
- Documentation uses the separate `docs/` environment with Documenter, CairoMakie, and GeoMakie; the docs workflow deploys through `julia-docdeploy`.

## Not Configured

- No root `Manifest.toml`, pre-commit configuration, or `deps/build.jl` exists.
- `Pkg.build()` is not required to reproduce the package test workflow.

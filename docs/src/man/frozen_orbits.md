Frozen Orbits
=============

```@meta
CurrentModule = SatelliteAnalysis
DocTestSetup = quote
    using SatelliteAnalysis
end
```

Due to the Earth's gravitational perturbation, the orbit of a salite will experience
secular changes in the argument of perigee. Hence, the satellite mean altitude per latitude
will differ during the mission. This effect can be problematic, especially if we must
compare images by a camera onboard the satellite in different periods. The altitude
variation will change the resolution, leading to some problems when comparing the data.

We can avoid this problem if we compute an eccentricity ``e``  and the argument of perigee
``\omega`` that yields theoretically:

```math
\begin{equation*}
  \frac{de}{dt} = 0,\ \frac{d\omega}{dt} = 0\ .
\end{equation*}
```

This orbit is called **frozen**. Refer to **[1]** for more information.

We can use the function:

```julia
frozen_orbit(a::Number, i::Number; kwargs...) -> Float64, Float64
```

to compute the eccentricity [ ] and argument of perigee [rad] that yield a frozen orbit when
the orbit has semi-major axis `a` [m] and inclination `i` [rad]. This function uses the
theory in **[1]**.

!!! note
    This function uses `BigFloat` internally to perform all computations, allowing very high
    degrees. However, the user must ensure that the default precision is enough for the
    required degree. Refer to the function `setprecision` for more information.

The following keywords are available:

- `gravity_model::Union{Nothing, AbstractGravityModel}`: Gravity model used to compute the
    frozen eccentricity. Refer to the object `AbstractGravityModel` of the package
    `SatelliteToolboxGravityModels.jl` for more information. If it is `nothing`, the system
    will automatically fetch and load the EGM96 gravity model. However, loading a gravity
    model can significantly decrease the performance. Thus, it is advisable to pass a
    gravity model here.
    (**Default** = `nothing`)
- `max_degree`: Maximum gravity model degree used to compute the frozen eccentricity. If it
    is equal to or lower than 0, the maximum degree in `grav_model` will be used. Otherwise,
    if it is lower than 3 or higher than the `grav_model` maximum degree, it will be clamped
    accordingly.
    (**Default** = 53)

## Examples

We will compute the eccentricity and argument of perigee that yields a frozen orbit using the
data from Amazonia-1 mission. First, we will use only 5 degrees, and the default gravity
model (EGM96):

```jldoctest
julia> frozen_orbit(7130.982e3, 98.410 |> deg2rad; max_degree = 5)
[ Info: Downloading the ICGEM file 'EGM96.gfc' from 'http://icgem.gfz-potsdam.de/getmodel/gfc/971b0a3b49a497910aad23cd85e066d4cd9af0aeafe7ce6301a696bed8570be3/EGM96.gfc'...
(0.0011108978494835141, 1.5707963267948966)

julia> e, ω = frozen_orbit(7130.982e3, 98.410 |> deg2rad; max_degree = 5)
(0.0011108978494835141, 1.5707963267948966)

julia> e
0.0011108978494835141

julia> ω |> rad2deg
90.0
```

If we want to use all the 360 terms in EGM96, we need to increase the precision of
`BigFloat` to keep the accuracy:

```jldoctest
julia> setprecision(1024)
1024

julia> e, ω = frozen_orbit(7130.982e3, 98.410 |> deg2rad; max_degree = 5)
(0.0011108978494835141, 1.5707963267948966)

julia> e
0.0011108978494835141

julia> ω |> rad2deg
90.0
```

We can use a different gravity model as follows:

```jldoctest
julia> jgm3 = GravityModels.load(IcgemFile, fetch_icgem_file(:JGM3))
[ Info: Downloading the ICGEM file 'JGM3.gfc' from 'http://icgem.gfz-potsdam.de/getmodel/gfc/a3375e01a717ac162962138a5e94f10466b71aa4a130d7f7d5b18ab3d5f90c3d/JGM3.gfc'...
IcgemFile{Float64}:
      Product type : gravity_field
       Model name  : JGM3
  Gravity constant : 3.986004415e14
            Radius : 6.3781363e6
    Maximum degree : 70
            Errors : formal
       Tide system : unknown
              Norm : fully_normalized
         Data type : Float64

julia> e, ω = frozen_orbit(7130.982e3, 98.410 |> deg2rad; max_degree = 70, gravity_model = jgm3)
(0.001163504566870769, 1.5707963267948966)

julia> e
0.001163504566870769

julia> ω |> rad2deg
90.0
```

## References

- **[1]** **Rosborough, G. W.; Ocampo, C. A (1991)**. _Influence of higher degree zonals on
  the frozen orbit geometry_. Proceedings of the **AAS/AIAA Astrodynamics Conference**,
  Durango, CO.

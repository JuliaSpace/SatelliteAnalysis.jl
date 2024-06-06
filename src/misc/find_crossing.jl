## Description #############################################################################
#
# Find crossing of a function.
#
############################################################################################

"""
    find_crossing(f::Function, t₀::Number, t₁::Number, s₀, s₁, vargs...; Δ = 1e-3, max = 100, kwargs...) -> T

Return the crossing time `tc` in which the function `f(t)` goes from the state `s₀` to the
state `s₁`. It is assumed that `f(t₀) = s₀` and `f(t₁) = s₁`.

The parameters in `vargs...` are passed to the function `f` after `t`, and the keywords
`kwargs...` are also passed to `f`. Hence, it will always be called as
`f(t, vargs...; kwargs...)`.

If the computed interval is smaller than `Δ`, or if the number of iterations is higher than
`max`, the algorithm stops.

!!! note

    The output type `T` is obtained by the type of `(t₁ + t₂) / 2`.

# Examples

```julia-repl
julia> SatelliteAnalysis.find_crossing(
    t -> (sin(t) > 0),
    -0.3,
    0.3,
    false,
    true;
    Δ = 1e-10
)
6.984919309616089e-11
```
"""
function find_crossing(
    f::Function,
    t₀::Number,
    t₁::Number,
    s₀,
    s₁,
    vargs::Vararg{Any, N};
    Δ = 1e-3,
    max = 100,
    kwargs...
) where N
    it = 0

    T = typeof((t₁ + t₀) / 2)

    ti₀ = T(t₀)
    ti₁ = T(t₁)

    while it <= max
        # Call the function at the middle of the interval.
        ti = (ti₁ + ti₀) / 2
        si = f(ti, vargs...; kwargs...)

        # Compute the new interval.
        if si == s₀
            ti₀ = ti
        elseif si == s₁
            ti₁ = ti
        else
            error("The function `f` returned an unexpected state.")
        end

        # If the interval is small enough, return.
        (ti₁ - ti₀) < Δ && break

        it += 1
    end

    return T(ti₁)
end

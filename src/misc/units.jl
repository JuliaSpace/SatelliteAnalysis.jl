## Description #############################################################################
#
# Functions to validate and convert the units selected by symbols in the API.
#
############################################################################################

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _angle_unit_factor(angle_unit::Symbol) -> Float64

Return the factor that converts an angle in radians to the unit `angle_unit`, which can be
`:rad` for radians or `:deg` for degrees.

# Extended help

## Throws

- `ArgumentError`: If `angle_unit` is not `:rad` or `:deg`.
"""
function _angle_unit_factor(angle_unit::Symbol)
    angle_unit == :rad && return 1.0
    angle_unit == :deg && return 180 / π

    return throw(
        ArgumentError("The angle unit `:$angle_unit` is not valid. Use `:rad` or `:deg`.")
    )
end

"""
    _distance_unit_factor(distance_unit::Symbol) -> Float64

Return the factor that converts a distance in meters to the unit `distance_unit`, which can
be `:m` for meters or `:km` for kilometers.

# Extended help

## Throws

- `ArgumentError`: If `distance_unit` is not `:m` or `:km`.
"""
function _distance_unit_factor(distance_unit::Symbol)
    distance_unit == :m && return 1.0
    distance_unit == :km && return 1 / 1000

    return throw(
        ArgumentError(
            "The distance unit `:$distance_unit` is not valid. Use `:m` or `:km`."
        ),
    )
end

"""
    _time_unit_factor(time_unit::Symbol[, valid_units::Tuple]) -> Float64

Return the factor that converts a time in seconds to the unit `time_unit`. The supported
units are `:s` for seconds, `:min` for minutes, `:h` for hours, `:d` for days, and `:y` for
Julian years. The tuple `valid_units` contains the units the caller accepts, and, if it is
omitted, only `:s`, `:min`, and `:h` are valid.

!!! note

    The symbol `:m` is not a valid time unit because it selects meters in the distance
    units.

# Extended help

## Throws

- `ArgumentError`: If `time_unit` is not in `valid_units`.
"""
function _time_unit_factor(time_unit::Symbol, valid_units::Tuple = (:s, :min, :h))
    if time_unit ∉ valid_units
        str_valid_units = join(("`:$u`" for u in valid_units), ", ", ", or ")

        # The minutes were selected by `:m` in previous versions.
        hint = time_unit == :m ? " Notice that the minutes are selected by `:min`." : ""

        throw(
            ArgumentError(
                "The time unit `:$time_unit` is not valid. Use $str_valid_units.$hint"
            ),
        )
    end

    time_unit == :s && return 1.0
    time_unit == :min && return 1 / 60
    time_unit == :h && return 1 / 3600
    time_unit == :d && return 1 / 86400

    # Julian year.
    return 1 / (365.25 * 86400)
end

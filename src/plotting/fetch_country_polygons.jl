## Description #############################################################################
#
# Fetch the GeoJSON data with the country polygons.
#
############################################################################################

export fetch_country_polygons

"""
    fetch_country_polygons(url = "https://datahub.io/core/geo-countries/r/countries.geojson"; kwargs...) -> String

Fetch the GeoJSON file with the country polygons in `url`. The algorithm stores the file in
a scratch space. The function returns a `String` with the file path.

!!! note
    If the file has already been downloaded, this function only returns its path. However,
    if the keyword `force_download` is `true`, the file is downloaded again from `url`.

# Keywords

- `force_download::Bool`: Download the file from `url` even if it already exists.
    (**Default** = `false`)
"""
function fetch_country_polygons(
    url = "https://datahub.io/core/geo-countries/r/countries.geojson";
    force_download::Bool = false
)

    url = "https://datahub.io/core/geo-countries/r/countries.geojson"
    filename = "countries.geojson"

    # Get the scratch space where the files are located.
    cache_dir          = @get_scratch!("geojson")
    filepath           = joinpath(cache_dir, filename)
    filepath_timestamp = joinpath(cache_dir, filename * "_timestamp")

    # If the file exists, only re-download if `force_download` is true.
    download_file = false

    if force_download ||
        isempty(readdir(cache_dir)) ||
        !isfile(filepath) ||
        !isfile(filepath_timestamp)

        download_file = true
    end

    if download_file
        @info "Downloading the file '$filename' from '$url'..."
        Downloads.download(url, filepath)
        open(filepath_timestamp, "w") do f
            write(f, string(now()))
        end
    end

    # Return the file path.
    return filepath
end

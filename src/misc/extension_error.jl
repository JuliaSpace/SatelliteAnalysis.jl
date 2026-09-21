## Description #############################################################################
#
# Function to throw the error related to functions provided by package extensions.
#
############################################################################################

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _extension_error(function_name::String, package::String, using_example::String, valid_call::String) -> Union{}

Throw the error of the fallback method of the function `function_name`, which is provided by
the package extension loaded together with `package`. This method is called if the
extension is not loaded or if the arguments are not valid. Hence, the message must help the
user in both cases: `using_example` is an example of a package that loads the extension,
and `valid_call` is the signature of the valid call.

# Extended help

## Throws

- `ErrorException`: Always.
"""
function _extension_error(
    function_name::String, package::String, using_example::String, valid_call::String
)
    return error(
        "The function `$function_name` is provided by a package extension. Load $package " *
        "(e.g., `using $using_example`) to use it. If $package is already loaded, check " *
        "the arguments: the valid call is `$valid_call`.",
    )
end

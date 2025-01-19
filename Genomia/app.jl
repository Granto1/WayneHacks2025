module App
# == Packages ==
# set up Genie development environment. Use the Package Manager to install new packages
using GenieFramework
@genietools

# == Code import ==
# add your data analysis code here or in the lib folder. Code in lib/ will be
# automatically loaded
function mean_value(x)
    sum(x) / length(x)
end

# == Reactive code ==
# add reactive code to make the UI interactive
@app begin
    # == Reactive variables ==
    # reactive variables exist in both the Julia backend and the browser with two-way synchronization
    # @out variables can only be modified by the backend
    # @in variables can be modified by both the backend and the browser
    # variables must be initialized with constant values, or variables defined outside of the @app block
    @in N = 10
    @out x = collect(1:10)
    @out y = randn(10) # plot data must be an array or a DataFrame
    @out msg = "The average is 0."
    @in shift = false

    # == Reactive handlers ==
    # reactive handlers watch a variable and execute a block of code when
    # its value changes
    @onchange N begin
        # the values of x, result and msg in the UI will
        # be automatically updated
        x = collect(1:N)
        y = rand(N)
        result = mean_value(rand(N))
        msg = "The average is $result."
    end
    # the onbutton handler will set the variable to false after the block is executed
    @onbutton shift begin
        y = circshift(y, 1)
    end
end

# == Pages ==
# register a new route and the page that will be loaded on access
@page("/", "app.jl.html")
end

# == Advanced features ==
#=
- The @private macro defines a reactive variable that is not sent to the browser.
This is useful for storing data that is unique to each user session but is not needed
in the UI.
    @private table = DataFrame(a = 1:10, b = 10:19, c = 20:29)

=#

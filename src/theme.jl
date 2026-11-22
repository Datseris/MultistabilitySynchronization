using CairoMakie

# decide theme:
ENV["COLORSCHEME"] = "JuliaDynamics"
ENV["BGCOLOR"] = :white
ENV["AXISCOLOR"] = :black

using MakieForProjects

"""
    record_interaction(file, figure; framerate = 30, total_time = 10)

As soon as you call this function, anything occurring on a live GLMakie window gets
recorded into the given file (`.mp4` or other video format).

## Keyword arguments

* `framerate = 30`.
* `total_time = 10`: Time to record for, in seconds.
* `sleep_time = 1`: Time to call `sleep()` before starting to record,
  so that you can go to the GLMakie window in time.
"""
function record_interaction(file, figure::Figure;
        framerate = 30, total_time = 10, sleep_time = 1,
    )
    ispath(dirname(file)) || mkpath(dirname(file))
    sleep(sleep_time)
    framen = framerate*total_time
    record(figure, file; framerate) do io
        for _ in 1:framen
            sleep(1/framerate)
            recordframe!(io)
        end
    end
    return
end
record_interaction(figure::Figure, file; kwargs...) =
record_interaction(file, figure; kwargs...)

export record_interaction

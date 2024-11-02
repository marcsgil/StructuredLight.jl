module StructuredLightMakieExt

import StructuredLight: visualize, save_animation
using Makie, LinearAlgebra

"""
    visualize(img;
    colormap=:jet, colorrange=Makie.automatic,
    share_colorrange=false, scaling=1)

Vizualize the image(s) `img`.

`img` can be a 2D, 3D or 4D array.

When using the 3D Array signature, the third dimension is interpreted as defining different images, which are displayed in a row.

When using the 4D Array signature, the third and fourth dimensions are interpreted as defining different images, which are displayed in a matrix.

## keyword arguments

`colormap` defines the colormap of the image(s).

`scaling` defines the maximum side length of the image, which is passed to `Makie.Figure`.

`share_colorrange` defines if the color range should be shared between all images. Has no effect for the 2D Array signature.

`colorrange` defines the color range of the image(s).
"""
visualize(img::AbstractArray{T}; kwargs...) where {T<:Real} = visualize(Array(img); kwargs...)

function visualize(img::Matrix{T};
    colormap=:jet, colorrange=Makie.automatic,
    share_colorrange=false, scaling=1) where {T<:Real}
    visualize(reshape(img, size(img)..., 1); colormap, colorrange, share_colorrange, scaling)
end

function visualize(img::Array{T,3}; colormap=:jet, colorrange=Makie.automatic,
    share_colorrange=false, scaling=1) where {T<:Real}
    visualize(reshape(img, size(img)..., 1); colormap, colorrange, share_colorrange, scaling)
end

function visualize(img::Array{T,4}; colormap=:jet, colorrange=Makie.automatic,
    share_colorrange=false, scaling=1) where {T<:Real}
    colorrange = share_colorrange ? extrema(img) : colorrange

    width, height = size(img)
    width = round(Int, width * scaling)
    height = round(Int, height * scaling)

    fig = Figure(figure_padding=0)
    for (m, column) in enumerate(eachslice(img, dims=4))
        for (n, _img) in enumerate(eachslice(column, dims=3))
            ax = Axis(fig[m, n]; aspect=DataAspect(), width, height)
            hidedecorations!(ax)
            hidespines!(ax)
            heatmap!(ax, _img; colormap, colorrange)
        end
    end
    colgap!(fig.layout, 0)
    rowgap!(fig.layout, 0)
    resize_to_layout!(fig)
    fig
end

"""
    save_animation(img::AbstractArray{T,3}, path;
    colormap=:jet, share_colorrange=false, scaling=1, framerate=16) where {T<:Real}

Save an animation of `img` at `path` with a `framerate`.

The colors are given by `colormap`.

`scaling` defines the maximum side length of the image, which is passed to `Makie.Figure`.

`share_colorrange` defines if the color range should be shared between all images.
"""
function save_animation(img::AbstractArray{T,3}, path;
    colormap=:jet, share_colorrange=false, scaling=1, framerate=16) where {T<:Real}
    colorrange = share_colorrange ? extrema(img) : Makie.automatic

    width, height = size(img)
    width = round(Int, width * scaling)
    height = round(Int, height * scaling)

    fig = Figure(figure_padding=0)
    ax = Axis(fig[1, 1]; aspect=DataAspect(), width, height)
    hidedecorations!(ax)
    hidespines!(ax)

    slices = eachslice(img, dims=3)
    hm = heatmap!(ax, first(slices); colormap, colorrange)
    resize_to_layout!(fig)

    record(fig, path, slices; framerate) do _img
        hm[3][] = _img
    end
end

end
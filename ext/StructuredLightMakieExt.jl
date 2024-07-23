module StructuredLightMakieExt

import StructuredLight: visualize, save_animation
using Makie, LinearAlgebra

"""
    visualize(img;
    colormap=:jet, colorrange=Makie.automatic,
    share_colorrange=false, max_size=1080)

Vizualize the image(s) `img`.

`img` can be a 2D, 3D or 4D array.

When using the 3D Array signature, the third dimension is interpreted as defining different images, which are displayed in a row.

When using the 4D Array signature, the third and fourth dimensions are interpreted as defining different images, which are displayed in a matrix.

## keyword arguments

`colormap` defines the colormap of the image(s).

`max_size` defines the maximum side length of the image, which is passed to `Makie.Figure`.

`share_colorrange` defines if the color range should be shared between all images. Has no effect for the 2D Array signature.

`colorrange` defines the color range of the image(s).
"""
visualize(img; kwargs...) = visualize(Array(img); kwargs...)

function visualize(img::AbstractMatrix{T};
    colormap=:jet, colorrange=Makie.automatic,
    share_colorrange=false, max_size=1080) where {T<:Real}
    visualize(reshape(img, size(img)..., 1); colormap, colorrange, share_colorrange, max_size)
end

function visualize(img::AbstractArray{T,3}; colormap=:jet, colorrange=Makie.automatic,
    share_colorrange=false, max_size=1080) where {T<:Real}
    visualize(reshape(img, size(img)..., 1); colormap, colorrange, share_colorrange, max_size)
end

function visualize(img::AbstractArray{T,4}; colormap=:jet, colorrange=Makie.automatic,
    share_colorrange=false, max_size=1080) where {T<:Real}
    colorrange = share_colorrange ? extrema(img) : colorrange

    width_factor = size(img, 1) * size(img, 4)
    height_factor = size(img, 2) * size(img, 3)

    if width_factor > height_factor
        width = max_size
        height = max_size * height_factor ÷ width_factor
    else
        height = max_size
        width = max_size * width_factor ÷ height_factor
    end

    fig = Figure(size=(width, height), figure_padding=0)
    for (m, column) in enumerate(eachslice(img, dims=4))
        for (n, _img) in enumerate(eachslice(column, dims=3))
            ax = Axis(fig[m, n], aspect=DataAspect())
            hidedecorations!(ax)
            heatmap!(ax, _img; colormap, colorrange)
        end
    end
    fig
end

"""
    save_animation(img::AbstractArray{T,3}, path;
    colormap=:jet, share_colorrange=false, max_size=1080, framerate=16) where {T<:Real}

Save an animation of `img` at `path` with a `framerate`.

The colors are given by `colormap`.

`max_size` defines the maximum side length of the image, which is passed to `Makie.Figure`.

`share_colorrange` defines if the color range should be shared between all images.
"""
function save_animation(img::AbstractArray{T,3}, path;
    colormap=:jet, share_colorrange=false, max_size=1080, framerate=16) where {T<:Real}
    colorrange = share_colorrange ? extrema(img) : Makie.automatic

    width_factor = size(img, 1)
    height_factor = size(img, 2)

    if width_factor > height_factor
        width = max_size
        height = max_size * height_factor ÷ width_factor
    else
        height = max_size
        width = max_size * width_factor ÷ height_factor
    end
    width, height

    fig = Figure(size=(height, width), figure_padding=0)
    ax = Axis(fig[1, 1], aspect=DataAspect())
    hidedecorations!(ax)
    hm = heatmap!(ax, rotr90(view(img, :, :, 3)); colormap, colorrange)

    record(fig, path, eachslice(img, dims=3); framerate) do _img
        hm[3][] = rotr90(_img)
    end
end

end
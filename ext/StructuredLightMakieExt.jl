module StructuredLightMakieExt

import StructuredLight: visualize, save_animation
using Makie, LinearAlgebra

"""
    visualize(img::AbstractMatrix{T}; colormap=:jet) where {T<:Real}
    visualize(img::AbstractArray{T,3}; colormap=:jet, max_size=1080, share_colorrange=false) where {T<:Real}
    visualize(img::AbstractArray{T,4}; colormap=:jet, max_size=1080, share_colorrange=false) where {T<:Real}

Vizualize the image(s) `img` with the `colormap`.

When using the 3D Array signature, the third dimension is interpreted as defining different images, which are displayed in a row.

When using the 4D Array signature, the third and fourth dimensions are interpreted as defining different images, which are displayed in a matrix.

`max_size` defines the maximum side length of the image, which is passed to `Makie.Figure`.

`share_colorrange` defines if the color range should be shared between all images.
"""
visualize(img; kwargs...) = visualize(Array(img); kwargs...)

function visualize(img::AbstractMatrix{T}; colormap=:jet) where {T<:Real}
    visualize(reshape(img, size(img)..., 1); colormap, max_size=1080)
end

function visualize(img::AbstractArray{T,3}; colormap=:jet, max_size=1080, share_colorrange=false) where {T<:Real}
    visualize(reshape(img, size(img)..., 1); colormap, max_size, share_colorrange)
end

function visualize(img::AbstractArray{T,4}; colormap=:jet, max_size=1080, share_colorrange=false) where {T<:Real}
    colorrange = share_colorrange ? extrema(img) : Makie.automatic

    width_factor = size(img, 1) * size(img, 4)
    height_factor = size(img, 2) * size(img, 3)

    if width_factor > height_factor
        width = max_size
        height = max_size * height_factor ÷ width_factor
    else
        height = max_size
        width = max_size * width_factor ÷ height_factor
    end

    fig = Figure(size=(height, width), figure_padding=0)
    for (m, column) in enumerate(eachslice(img, dims=4))
        for (n, _img) in enumerate(eachslice(column, dims=3))
            ax = Axis(fig[m, n], aspect=DataAspect())
            hidedecorations!(ax)
            heatmap!(ax, rotr90(_img); colormap, colorrange)
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
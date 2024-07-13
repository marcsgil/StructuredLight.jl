module StructuredLightMakieExt

import StructuredLight: visualize
using Makie, LinearAlgebra

function visualize(img::AbstractMatrix; colormap=:jet)
    visualize(reshape(img, size(img)..., 1); colormap)
end

function visualize(img::AbstractArray{T,3}; colormap=:jet, share_colorrange=false, max_size=1080) where {T}
    visualize(reshape(img, size(img)..., 1); colormap, share_colorrange, max_size)
end

function visualize(img::AbstractArray{T,4}; colormap=:jet, share_colorrange=false, max_size=1080) where {T}
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
            image!(ax, rotr90(_img); colormap, colorrange)
        end
    end
    fig
end

end
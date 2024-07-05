convert2color(x::Real, colormap, normalization) = convert(RGB{N0f8}, get(colorschemes[colormap], normalization(x)))

function convert2image!(dest, src, colormap, normalization)
    @tullio dest[n] = convert2color(src[n], colormap, normalization)
    reverse!(dest, dims=1)
end

Base.extrema(itr::AbstractArray{Complex{T}}) where {T<:Real} = extrema(abs2, itr)

"""
    visualize(ψ::AbstractArray{T,2}; colormap=:jet,ratio=1,range=extrema(ψ)) where T<: Real
    visualize(ψ::AbstractArray{T,3}; colormap=:jet,ratio=1,normalize_by_first=false) where T<: Real
    visualize(ψ::AbstractArray{T,4}; colormap=:jet,ratio=1,normalize_by_first=false) where T<: Real

Vizualize the beam described by array `ψ`. 

The colors are given by `colormap`. The full catalogue can be found at the [ColorSchemes.jl documentation](https://juliagraphics.github.io/ColorSchemes.jl/stable/catalogue/).

The image is rescaled by `ratio`.

When using the 3D Array signature, the third dimension is interpreted as defining different images, which are displayed in a row.

When using the 4D Array signature, the third and fourth dimensions are interpreted as defining different images, which are displayed in a matrix.
"""
function visualize(ψ::AbstractArray{T,2}; colormap=:jet, ratio=1, range=extrema(ψ)) where {T}
    img = similar(ψ, RGB{N0f8})
    convert2image!(img, ψ, colormap, scaleminmax(range...))
    imresize(img; ratio)
end

function _visualize(ψ::AbstractArray{T,3}; colormap=:jet, normalize_by_first=false) where {T}
    img = similar(ψ, RGB{N0f8})

    if normalize_by_first
        range = extrema(view(ψ, :, :, 1))
    end

    for n in axes(ψ, 3)
        if !normalize_by_first
            range = extrema(view(ψ, :, :, n))
        end
        convert2image!(view(img, :, :, n), view(ψ, :, :, n), colormap, scaleminmax(range...))
    end

    img
end

function visualize(ψ::AbstractArray{T,3}; colormap=:jet, ratio=1, normalize_by_first=false) where {T}
    img = _visualize(ψ; colormap, normalize_by_first)
    imresize(hcat((slice for slice in eachslice(img, dims=3))...); ratio)
end

function visualize(ψ::AbstractArray{T,4}; colormap=:jet, ratio=1, normalize_by_first=false) where {T}
    vcat((visualize(slice; colormap, ratio, normalize_by_first) for slice in eachslice(ψ, dims=4))...)
end

"""
    show_animation(ψs::AbstractArray{T,3}; colormap=:jet,ratio=1,fps=16,normalize_by_first=false) where T 

Open a gif showing a visualization of `ψs`with a framerate of `fps`. It displays properly on VSCode, but not on Jupyter.

The colors are given by `colormap`. The full catalogue can be found at the [ColorSchemes.jl documentation](https://juliagraphics.github.io/ColorSchemes.jl/stable/catalogue/).

The image is rescaled by `ratio`.

`normalize_by_first` defines if the intensities should be normalized by the first image.
"""
function show_animation(ψs::AbstractArray{T,3}; colormap=:jet, ratio=1, fps=16, normalize_by_first=false) where {T}
    Images.gif(imresize(_visualize(ψs; colormap, normalize_by_first), ratio=(ratio, ratio, 1)); fps)
end

"""
    save_animation(ψs::AbstractArray{T,3}, path; 
    colormap=:jet,
    ratio=1,
    fps=16, 
    encoder_options = nothing) where T

Save an animation of `ψs` at `path` with a framerate of `fps`.

The colors are given by `colormap`. The full catalogue can be found at the [ColorSchemes.jl documentation](https://juliagraphics.github.io/ColorSchemes.jl/stable/catalogue/).

The image is rescaled by `ratio`.

Follow the [VideoIO.jl documentation](https://juliaio.github.io/VideoIO.jl/stable/writing/) for more information on `encoder_options`.
"""
function save_animation(ψs::AbstractArray{T,3}, path; colormap=:jet, ratio=1, fps=16, encoder_options=nothing) where {T}
    imgstack = visualize.(eachslice(ψs, dims=3); colormap, ratio)
    if isnothing(encoder_options)
        VideoIO.save(path, imgstack, framerate=fps)
    else
        VideoIO.save(path, imgstack, framerate=fps, encoder_options=encoder_options)
    end
end
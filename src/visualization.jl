normalize(ψ::AbstractArray{T,2}) where T  = ψ/maximum(abs.(ψ))

function normalize(ψ::AbstractArray{T,3};normalize_by_first=false) where T
    result = similar(ψ)
    if normalize_by_first
        result = ψ/maximum(abs.(view(ψ,:,:,1)))
    else
        for n in axes(ψ,3)
            result[:,:,n] = view(ψ,:,:,n)/maximum(abs.(view(ψ,:,:,n)))
        end
    end
    result
end

function convert2image(ψ::AbstractArray{T,N}; colormap=:hot,ratio=1) where {T <: Real ,N}
    imresize(map( pixel -> convert(RGB{N0f8}, get(colorschemes[colormap], pixel)), Array(ψ) ), ratio=ratio)
end

function convert2image(ψ::AbstractArray{T,N};colormap=:hot,ratio=1) where {T,N}
    convert2image(abs2.(ψ);colormap=colormap,ratio=ratio)
end

"""
    visualize(ψ::AbstractArray{T,2}; colormap=:hot, ratio=1) where T
    visualize(ψ::AbstractArray{T,3}; colormap=:hot, ratio=1) where T
    visualize(ψ::AbstractArray{T,4}; colormap=:hot, ratio=1) where T

Vizualize the beam described by array `ψ`. 

The colors are given by `colormap`. The full catalogue can be found at the [ColorSchemes.jl documentation](https://juliagraphics.github.io/ColorSchemes.jl/stable/catalogue/).

The image is rescaled by `ratio`.

When using the 3D Array signature, the third dimension is interpreted as defining different images, which are displayed in a row.

When using the 4D Array signature, the third and fourth dimensions are interpreted as defining different images, which are displayed in a matrix.
"""
function visualize(ψ::AbstractArray{T,2}; colormap=:hot,ratio=1) where T
    convert2image(normalize(Array(ψ)),colormap=colormap,ratio=ratio)
end

function visualize(ψ::AbstractArray{T,3}; colormap=:hot,ratio=1,normalize_by_first=false) where T
    convert2image(hcat(eachslice(normalize(Array(ψ),normalize_by_first=normalize_by_first),dims=3)...),colormap=colormap,ratio=ratio)
end

function visualize(ψ::AbstractArray{T,4}; colormap=:hot,ratio=1,normalize_by_first=false) where T
    vcat(visualize.( eachslice(Array(ψ),dims=4),colormap=colormap,ratio=ratio,normalize_by_first=normalize_by_first )...)
end


"""
    show_animation(ψs::AbstractArray{T,3}; colormap=:hot,ratio=1,fps=16) where T 

Open a gif showing a visualization of `ψs`with a framerate of `fps`. It displays properly on VSCode, but not on Jupyter.

The colors are given by `colormap`. The full catalogue can be found at the [ColorSchemes.jl documentation](https://juliagraphics.github.io/ColorSchemes.jl/stable/catalogue/).

The image is rescaled by `ratio`.

`normalize_by_first` defines if the intensities should be normalized by the first image.
"""
function show_animation(ψs::AbstractArray{T,3}; colormap=:hot,ratio=1,fps=16) where T 
    Images.gif([convert2image(normalize(Array(ψ)),colormap=colormap,ratio=ratio) for ψ in eachslice(ψs,dims=3)],fps=fps)
end

"""
    save_animation(ψs::AbstractArray{T,3}, path; 
    colormap=:hot,
    ratio=1,
    fps=16, 
    encoder_options = nothing) where T

Save an animation of `ψs` at `path` with a framerate of `fps`.

The colors are given by `colormap`. The full catalogue can be found at the [ColorSchemes.jl documentation](https://juliagraphics.github.io/ColorSchemes.jl/stable/catalogue/).

The image is rescaled by `ratio`.

Follow the [VideoIO.jl documentation](https://juliaio.github.io/VideoIO.jl/stable/writing/) for more information on `encoder_options`.
"""
function save_animation(ψs::AbstractArray{T,3}, path; colormap=:hot,ratio=1,fps=16, encoder_options = nothing) where T
    imgstack = [RGB{N0f8}.(convert2image(normalize(Array(ψ)),colormap=colormap,ratio=ratio)) for ψ in eachslice(ψs,dims=3)]
    if isnothing(encoder_options)
        VideoIO.save(path, imgstack, framerate=fps)
    else
        VideoIO.save(path, imgstack, framerate=fps, encoder_options = encoder_options)
    end    
end

#=
ImageView is based on GTK, which isn't working properly on Windows.
Therefore, I'm disabiling this functionality for now.

"""
    interactive_visualization(ψs::AbstractArray{T,3}; 
        colormap=:hot,
        ratio=1,
        normalize_by_first=false) where T 

Open a window showing an interactive visualization of `ψs`.

The colors are given by `colormap`. The full catalogue can be found at the [ColorSchemes.jl documentation](https://juliagraphics.github.io/ColorSchemes.jl/stable/catalogue/).

The image is rescaled by `ratio`.

`normalize_by_first` defines if the intensities should be normalized by the first image.
"""
function interactive_visualization(ψs::AbstractArray{T,3}; colormap=:hot,ratio=1,normalize_by_first=false) where T 
    imshow(convert2image(normalize(Array(ψs),normalize_by_first=normalize_by_first),colormap=colormap,ratio=ratio))
end
=#
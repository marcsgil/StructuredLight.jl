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

function visualize(ψ::AbstractArray{T,2}; colormap=:hot,ratio=1) where T
    convert2image(normalize(Array(ψ)),colormap=colormap,ratio=ratio)
end

function visualize(ψ::AbstractArray{T,3}; colormap=:hot,ratio=1,normalize_by_first=false) where T
    convert2image(hcat(eachslice(normalize(Array(ψ),normalize_by_first=normalize_by_first),dims=3)...),colormap=colormap,ratio=ratio)
end

function visualize(ψ::AbstractArray{T,4}; colormap=:hot,ratio=1,normalize_by_first=false) where T
    vcat(visualize.( eachslice(Array(ψ),dims=4),colormap=colormap,ratio=ratio,normalize_by_first=normalize_by_first )...)
end

function show_animation(ψs::AbstractArray{T,3}; colormap=:hot,ratio=1,fps=16) where T 
    Images.gif([convert2image(normalize(Array(ψ)),colormap=colormap,ratio=ratio) for ψ in eachslice(ψs,dims=3)],fps=fps)
end

function save_animation(ψs::AbstractArray{T,3}, path; colormap=:hot,ratio=1,fps=16) where T
    imgstack = [RGB{N0f8}.(convert2image(normalize(Array(ψ)),colormap=colormap,ratio=ratio)) for ψ in eachslice(ψs,dims=3)]
    #encoder_options = (crf=23, preset="medium")
    VideoIO.save(path, imgstack, framerate=fps)
end
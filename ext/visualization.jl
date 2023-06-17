function StructuredLight.visualize(ψ::CuArray{T,2}; colormap=:hot,ratio=1,range=extrema(ψ)) where T
    StructuredLight.visualize(Array(ψ); colormap,ratio,range)
end

function StructuredLight.visualize(ψ::CuArray; colormap=:hot,ratio=1,normalize_by_first=false)
    StructuredLight.visualize(Array(ψ); colormap,ratio,normalize_by_first)
end

function StructuredLight.show_animation(ψs::CuArray{T,3}; colormap=:hot,ratio=1,fps=16,normalize_by_first=false) where T 
    StructuredLight.show_animation(Array(ψs); colormap,ratio,fps,normalize_by_first)
end

function StructuredLight.save_animation(ψs::CuArray{T,3}, path; colormap=:hot,ratio=1,fps=16, encoder_options = nothing) where T
    StructuredLight.save_animation(Array(ψs), path; colormap,ratio,fps, encoder_options)  
end
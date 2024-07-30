function get_array_size(x, y, ::Number)
    (length(x), length(y))
end

function get_array_size(x, y, z)
    (length(x), length(y), length(z))
end

function drop_singleton_dims(A, x, y, ::Number)
    reshape(A, length(x), length(y))
end

function drop_singleton_dims(A, x, ::Number, z)
    reshape(A, length(x), length(z))
end

function drop_singleton_dims(A, x::Number, y, z)
    reshape(A, length(x), length(y))
end

function drop_singleton_dims(A, ::Number, ::Number, z)
    reshape(A, length(z))
end

function drop_singleton_dims(A, x, ::Number, z)
    reshape(A, length(x), length(z))
end

function drop_singleton_dims(A, x::Number, y, z)
    reshape(A, length(x), length(y))
end
function get_array_size(x, y, ::Number)
    (length(x), length(y))
end

function get_array_size(x, y, z)
    (length(x), length(y), length(z))
end
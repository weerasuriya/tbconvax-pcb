# Not otherwise specified.
export yf

"""
Year finder
"""
function yf(y)
    ((y - 1900) * 2) .+ [1, 2]
end

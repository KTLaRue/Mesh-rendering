module GfxBase

export Vec3, Vec2

using StaticArrays

# some helpful type aliases:
#const Color = RGB{Float32} # RGB color represented with floats
const Vec3 = SVector{3, Float64} # 3-vector of floats
const Vec2 = SVector{2, Float64} # 2-vector of floats

end # module GfxBase

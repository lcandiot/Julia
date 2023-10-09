using Metal, Test

function vadd(a, b, c)
    idx, idy = thread_position_in_grid_2d()
    c[idx,idy] = a[idx,idy] + b[idx,idy]
    return
end

dims = (30,30)
a = round.(rand(Float32, dims) * 100)
b = round.(rand(Float32, dims) * 100)
c = similar(a)

d_a = MtlArray(a)
d_b = MtlArray(b)
d_c = MtlArray(c)

len = prod(dims)
@metal threads=len vadd(d_a, d_b, d_c)
c = Array(d_c)
@test a+b ≈ c
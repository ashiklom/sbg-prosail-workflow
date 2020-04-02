import ArchGDAL
const AG = ArchGDAL

using RTM
using Images, Colors
using ProgressMeter

using Plots

include("utils.jl")

# fname = "../data/aviris/f130503t01p00r14_refl/f130503t01p00r14rdn_refl_img_corr"
fname = "../data/aviris-ng/ang20170706t203356_rfl_v2p9/ang20170706t203356_corr_v2p9_img"
# fname = "../data/aviris-ng/ang20170706t204753_rfl_v2p9/ang20170706t204753_corr_v2p9_img"

rgb_preview(fname)

# Read wavelengths from header file
fname_hdr = fname * ".hdr"
waves = read_waves(fname_hdr)

# Read a subset of the raster
xs = 400:550
ys = 5270:5400

aviris = AG.read(fname)
data_raw = getallbands(aviris, xs, ys)

# Visualize it
gr()
plot(waves, data_raw[1,1,:]; ylim = (0, 0.3))
# Find the atmospheric window
Plots.vline!([1340, 1450])
Plots.vline!([1800, 1960])

# Remove incorrect values
data = convert(Array{Union{Missing, Float32}}, data_raw)
data[(data .< 0) .| (data .> 1)] .= missing

# Also remove atmospheric absorption features, values < 400, and other noise
atm1 = [findnearest(x, waves) for x in [1348, 1450]]
atm2 = [findnearest(x, waves) for x in [1800, 1960]]
atm3 = findnearest(2460, waves)
data[:,:,atm1[1]:atm1[2]] .= missing
data[:,:,atm2[1]:atm2[2]] .= missing
data[:,:,atm3:end] .= missing
# Everything below 400 nm
data[:,:,1:5] .= missing

avinds = findall(.!ismissing.(data[1,1,:]))
pwaves = round.(waves[avinds])
pinds = Int.(pwaves .- 399)

pixel = data[1,1,:]
plot(waves, pixel)


plot(prospect4(N, Cab, Cw, Cm))

plot(RTM.hapke_soil(0.3))
plot!(RTM.hapke_soil(0.5))
plot!(RTM.hapke_soil(0.8))

function aviris_p4s(N, Cab, Cw, Cm, LAI, psoil, inds = pinds)
    result = prospect4_4sailh(N, Cab, Cw, Cm, LAI, psoil)
    refl = result[:,1]
    plot(refl[:,1])
end

b1 = AG.read(aviris, 1)
# b1[b1 .< 0] .= 0
# b1 = b1 / 10000
b1 = convert(Array{Union{Missing, Float32}}, b1)
b1[b1 .< 0] .= missing
heatmap(b1)
# plot(Gray.(b1))

function getbands(img, xs, ys)
    nx = length(xs)
    ny = length(ys)
    nb = AG.nraster(img)
    # NOTE: X and Y are flipped?
    result = Array{Float32, 3}(undef, nx, ny, nb)
    @showprogress for i in 1:nb
        r = AG.read(aviris, i, ys, xs)
        result[:,:,i] = convert(Array{Float32, 2}, r) / 10000.0
    end
    result
end

xs = 300:350
ys = 4000:4050

# xs = 500:750
# ys = 5000:6000

aviris_sub = getbands(aviris, xs, ys)
plot(aviris_sub[2,2,:])

img = aviris_sub[:,:,[5, 20, 64]]

# Look at an RGB composite
imgc = colorview(RGB, permuteddimsview(img, (3, 1, 2)))
plot(imgc)

# Look at some spectra
sp1 = aviris_sub[200, 800:850, :]
plot(permuteddimsview(sp1, (2, 1)))

# Calculate some stats
img_mean = mean(aviris_sub, dims = (1,2))[1,1,:]
img_sd = std(aviris_sub, dims = (1,2))[1,1,:]
plot(img_mean)
plot!(img_mean .+ img_sd)
plot!(img_mean .- img_sd)

# Flatten the image
img_flat = mapslices(x -> vcat(x...), aviris_sub; dims = (1, 2))[:,1,:];

# b1 = AG.read(aviris, 2, 1:1000, 1:1000);
b1 = AG.read(aviris, 1);
b1f = convert(Array{Float32, 2}, b1) / 10000.;
plot(Gray.(b1f))
is_zero = b1f .!= 0;
plot(Gray.(is_zero))

histogram(vcat(b1f...))

plot(Gray.(sqrt(b1f)))

max(b1f...)
b1f_img = colorview()
plot(b1f)

plot(RGB.(rand(4, 4)))

win = AG.windows()

refl_s1 = AG.read(aviris, 1, 400:500, 400:500);
refl_s1f = convert(Array{Float32, 2}, refl_s1)

# Take samples

refl_small = aviris[1:100,1:100];

# refl = AG.read(aviris)

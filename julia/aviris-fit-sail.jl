import ArchGDAL
const AG = ArchGDAL

using RTM, ProgressMeter, Plots
include("julia/utils.jl")

# fname = "../data/aviris/f130503t01p00r14_refl/f130503t01p00r14rdn_refl_img_corr"
fname = "data/aviris-ng/ang20170706t203356_rfl_v2p9/ang20170706t203356_corr_v2p9_img"
# fname = "../data/aviris-ng/ang20170706t204753_rfl_v2p9/ang20170706t204753_corr_v2p9_img"

# Read wavelengths from header file
fname_hdr = fname * ".hdr"
waves = read_waves(fname_hdr)

# Read a subset of the raster
xs = 400:550
ys = 5270:5400

aviris = AG.read(fname)
data_raw = getallbands(aviris, xs, ys)

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

# Run PROSAIL at these same wavelengths
N = 1.4; Cab = 40; Cw = 0.01; Cm = 0.01
psoil = 0.3; LAI = 3

plot(prospect4_4sailh(N, Cab, Cw, Cm, LAI, psoil)[pinds,:])
plot!(pixel[avinds])

const obs_c = pixel[avinds]
const pinds_c = pinds

function lsq_likelihood(N, Cab, Cw, Cm, LAI, psoil)
    mod = prospect4_4sailh(N, Cab, Cw, Cm, LAI, psoil)[pinds_c,1]
    sum((mod .- obs_c).^2)
end

lsq_likelihood(1.4, 40, 0.01, 0.01, 3, 0.5)

using JuMP, Ipopt

model = Model(Ipopt.Optimizer)
@variable(model, 10 >= N >= 1)
@variable(model, 200 >= Cab >= 0)
@variable(model, 100 >= Car >= 0)
@variable(model, 0.5 >= Cw >= 0)
@variable(model, 0.5 >= Cm >= 0)
@variable(model, 10 >= LAI >= 0)
@variable(model, 0 <= psoil <= 1)

register(model, :lsq_likelihood, 6, lsq_likelihood, autodiff = true)
@NLobjective(model, Min, lsq_likelihood(N, Cab, Cw, Cm, LAI, psoil))
optimize!(model)

# Look at the fit
fit_spec = prospect4_4sailh(value(N), value(Cab), value(Cw), value(Cm),
                            value(LAI), value(psoil))[pinds_c,1]
plot(pwaves, obs_c)
plot!(pwaves, fit_spec)
plot(pwaves, fit_spec .- obs_c)
hline!([0])

# Now, try fitting multiple pixels
const aviris_window = data[1:10,1:10,avinds]

prosail_av(N, Cab, Cw, Cm, LAI, psoil) = prospect4_4sailh(N, Cab, Cw, Cm, LAI, psoil)[pinds_c,1]
Cabi = rand(Float32, 10, 10) .* 30 .+ 40
testimg = prosail_av.(1.4, Cabi, 0.01, 0.01, 3, 0.5)
mod_array = testimg

# Nonlinear optimizaiton doesn't allow vector variables. So we cheat using
# splatting.
function img_like(x...)
    N = x[1:100]
    Cab = x[101:200]
    Cw = x[201:300]
    Cm = x[301:400]
    LAI = x[401:500]
    psoil = x[501:600]
    idx = CartesianIndices(aviris_window[:,:,1])
    ss = 0.0
    for i in 1:length(idx)
        mod = prosail_av(N[i], Cab[i], Cw[i], Cm[i], LAI[i], psoil[i])
        ss += sum((mod .- aviris_window[idx[i][1], idx[i][2], :]) .^ 2)
    end
    ss
end

# Test the likelihood function
img_like(vcat(
    fill(1.4, 100),
    vec(Cabi),
    fill(0.01, 100),
    fill(0.01, 100),
    fill(3, 100),
    rand(100)
)...)

# Define the optimization model
model = Model(Ipopt.Optimizer)
@variable(model, x[1:600])
set_lower_bound.(x[1:100], 1)
set_lower_bound.(x[101:600], 0)
set_upper_bound.(x[1:100], 10)
set_upper_bound.(x[101:200], 200)
set_upper_bound.(x[201:400], 0.5)
set_upper_bound.(x[401:500], 10)
set_upper_bound.(x[501:600], 1)

# @variable(model, 1 <= x[1:100] <= 10)
# @variable(model, 0 <= x[101:200] <= 200)
# @variable(model, 0 <= x[201:300] <= 0.5)
# @variable(model, 0 <= x[301:400] <= 0.5)
# @variable(model, 0 <= x[401:500] <= 10)
# @variable(model, 0 <= x[501:600] <= 1)
register(model, :like, 600, img_like, autodiff = true)
@NLobjective(model, Min, like(x...))
optimize!(model)

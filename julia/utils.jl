"Read wavelengths from an ENVI header file."
function read_waves(filename)
    header = readlines(filename)
    iwave = map(header) do s occursin(r"^wavelength *=", s) end
    wave_string = header[iwave][1]
    wave_s2 = match(r"{(.*)}", wave_string).captures[1]
    wave_split = split(wave_s2, ',')
    waves = map(x -> parse(Float32, x), wave_split)
    waves
end

"Find index of nearest value. Useful for identifying specific wavelengths."
findnearest(target, vector = waves) = findmin(abs.(target .- vector))[2]

"Preview an image in a PyPlot buffer"
function rgb_preview(img, waves)
    red = findnearest(690, waves)
    green = findnearest(550, waves)
    blue = findnearest(400, waves)
    bands = [ArchGDAL.read(img, band) for band in [red, green, blue]]
    imgprep = map(bands) do b
        b[b .< 0] .= 0
        # Stretch to quantile
        b .= b ./ quantile(vec(b), 0.99)
        b
    end
    img2 = cat(imgprep...; dims = 3)
    imgc = colorview(RGB, permuteddimsview(img2, (3, 2, 1)))
    pyplot()
    plot(imgc)
end

function rgb_preview(imgfile::String)
    img = ArchGDAL.read(imgfile)
    fname_hdr = imgfile * ".hdr"
    waves = read_waves(fname_hdr)
    rgb_preview(img, waves)
end

function getallbands(img, xs, ys)
    nx = length(xs)
    ny = length(ys)
    nb = ArchGDAL.nraster(img)
    # NOTE: X and Y are flipped?
    result = Array{Float32, 3}(undef, nx, ny, nb)
    @showprogress for i in 1:nb
        r = ArchGDAL.read(aviris, i, ys, xs)
        result[:,:,i] = convert(Array{Float32, 2}, r)
    end
    result
end

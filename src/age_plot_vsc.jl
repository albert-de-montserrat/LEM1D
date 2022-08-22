
using CairoMakie, ColorSchemes, HDF5, DelimitedFiles
# include GLMakie for interactive plot, but it won't save to pdf

include("C:\\Users\\crosetto\\Documents\\TerracesJulia\\LEM1D-main\\src\\read_txt_files.jl")

sea_level_file = "Spratt2016-450" # take the las X years of the sea lvl 

if sea_level_file == "Spratt2016-800"
    starting_time = 800e3
elseif sea_level_file == "Spratt2016-450"
    starting_time = 450e3
elseif sea_level_file == "Bintanja3"
    starting_time = 3e6
elseif sea_level_file == "Bintanja6" || "Bintanja3x2"
    starting_time = 6e6
end
sea_lvl_curve, sea_age = get_sea_level(sea_level_file)
# sea_age                 = sea_age[end:-1:1]# reverse age array -> 0 = oldest age
id_time = sea_age .>= maximum(sea_age) .- starting_time
sea_age = sea_age[id_time]
sea_age = sea_age .- minimum(sea_age)
sea_lvl_curve = sea_lvl_curve[id_time]

pathio = "C:\\Users\\crosetto\\Documents\\TerracesJulia\\LEM1D-main\\output\\"
#pathio = "E:\\PROJECTS\\fwd modelling\\OUTPUTS_TR_Arg-Cor\\09-13__TR_arg_spr450mod2_Uriode\\"

model = "TRspratt450cor_U1.7-1.6_slope-0.1405_hwb5.0_beta1.0e-5_dx5.h5"

fname = string(pathio, model)

fid = h5open(fname, "r");

T = fid["Terrace"]
x = read(T, "x")
z = read(T, "z")
A = fid["Age"]
age = read(A, "age")

fig = Figure(; resolution=(800, 450), backgroundcolor=RGBf0(1, 1, 1))

ax1 = fig[1, 1] = Axis(fig)
l1 = lines!(
    ax1,
    x ./ 1e3,
    z;
    color=(age[end] .- age) ./ 1e3,
    colormap=:linear_bgyw_15_100_c68_n256,
    rev=true,
    linewidth=3,
) #:tableau_hue_circle
xlims!(ax1, 10, 22)
ylims!(ax1, -400, 1600)
# Colorbar(fig[1,2], l1, label="ka")

ax2 = fig[1, 2] = Axis(fig)
l2 = lines!(
    ax2,
    -sea_lvl_curve,
    sea_age ./ 1e3;
    color=sea_age ./ 1e3,
    colormap=:linear_bgyw_15_100_c68_n256,
    rev=true,
    linewidth=3,
)
xlims!(ax2, extrema(-sea_lvl_curve))
ylims!(ax2, extrema(age / 1e3))
fig
save("TRspratt450cor_U1.7-1.6_slope-0.1405_hwb5.0_beta1.0e-5_dx5.pdf", fig)

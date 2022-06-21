### A Pluto.jl notebook ###
# v0.19.8

using Markdown
using InteractiveUtils

# ╔═╡ 6b400190-a550-4e32-baa6-d54b66d7676e
# ╠═╡ show_logs = false
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
		Pkg.add("Statistics")
		Pkg.add("StatsBase")
		Pkg.add("ImageMorphology")
		Pkg.add("ImageFiltering")
		Pkg.add("CSV")
		Pkg.add("DataFrames")
		Pkg.add("GLM")
		Pkg.add(url="https://github.com/JuliaHealth/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
		Pkg.add(url="https://github.com/Dale-Black/PhantomSegmentation.jl")
		Pkg.add(url="https://github.com/Dale-Black/CalciumScoring.jl")
		Pkg.add("CairoMakie")
		Pkg.add("HypothesisTests")
		Pkg.add("Colors")
		Pkg.add("MLJBase")
	end
	
	using PlutoUI
	using Statistics
	using StatsBase: quantile!, rmsd
	using ImageMorphology
	using ImageFiltering
	using CSV
	using DataFrames
	using GLM
	using DICOM
	using DICOMUtils
	using PhantomSegmentation
	using CalciumScoring
	using CairoMakie
	using HypothesisTests
	using Colors
	using MLJBase
end

# ╔═╡ 38cd4959-e50a-4724-849d-9a246dd4414e
TableOfContents()

# ╔═╡ 787dda66-c0c9-4e2f-8646-c1c5ac95feee
root_path = "/Users/daleblack/Google Drive/dev/MolloiLab/CAC-stanford-data/data/output/";

# ╔═╡ 5017031f-2e8d-44f9-925d-a14376050342
md"""
## Integrated
"""

# ╔═╡ 61a077c9-6bc2-450f-a58f-6622ca46a592
md"""
#### 1 Point Calibration
"""

# ╔═╡ 1a7206ea-4be5-4236-954e-a927f9cba993
df_i1 = CSV.read(string(root_path, "integrated.csv"), DataFrame);

# ╔═╡ 9e71a277-0a70-483d-ab5f-2269ba244ca4
let
	df = df_i1
	gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	global model_i1
	model_i1 = lm(@formula(Y ~ X), data)
	global r2_1
	r2_1 = GLM.r2(model_i1)
	global rms_values1
	rms_values1 = [
		rms(data[!, :X], data[!, :Y]),
		rmsd(data[!, :Y], GLM.predict(model_i1))
	]
end

# ╔═╡ 83f1c3d3-f566-4c69-8eaf-8c363286e031
begin
	newX1 = DataFrame(X=collect(1:1000));
	pred_i1 = GLM.predict(model_i1, newX1)
end

# ╔═╡ 1fe58e55-eafb-491a-8402-b3b07ac3ca2f
co1 = coef(model_i1)

# ╔═╡ 37284bae-a062-4620-ad53-66de09ad9322
md"""
#### 3 Point Calibration
"""

# ╔═╡ 19713b20-4c52-488b-9f0d-13d9e6599645
df_i2 = CSV.read(string(root_path, "integrated3cal.csv"), DataFrame);

# ╔═╡ 9f9c0080-c311-442f-b4d8-1efd0552845b
let
	df = df_i2
	gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	global model_i2
	model_i2 = lm(@formula(Y ~ X), data)
	global r2_2
	r2_2 = GLM.r2(model_i2)
	global rms_values2
	rms_values2 = [
		rms(data[!, :X], data[!, :Y]),
		rmsd(data[!, :Y], GLM.predict(model_i2))
	]
end

# ╔═╡ 97c05609-c9a0-41e6-9f21-2d8dfac412aa
begin
	newX2 = DataFrame(X=collect(1:1000));
	pred_i2 = GLM.predict(model_i2, newX2)
end

# ╔═╡ 736e7c7d-d53d-4706-ba8f-c83c6fde3d14
co2 = coef(model_i2)

# ╔═╡ d017319a-34aa-4d81-b230-203a43c66602
let
	f = Figure()

	##-- A --##
	axtop = Axis(f[1, 1])
	
	df = df_i1
	scatter!(axtop, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
	errorbars!(axtop, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], rms(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large]))
	scatter!(axtop, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
	errorbars!(axtop, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], rms(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium]))
	scatter!(axtop, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
	errorbars!(axtop, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], rms(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small]))
	lines!(axtop, [-1000, 1000], [-1000, 1000],)
	lines!(axtop, collect(1:1000), pred_i1, linestyle=:dashdot)
	Textbox(
		f[1, 1], 
		placeholder = "y = $(trunc(co1[2]; digits=3))x + $(trunc(co1[1]; digits=3)) \nr = $(trunc(r2_1; digits=3)) \nRMSE: $(trunc(rms_values1[1]; digits=3)) \nRMSD: $(trunc(rms_values1[2]; digits=3))", 
		tellheight = false,
        tellwidth = false,
		boxcolor=:white,
		halign=:left,
		valign=:top,
		textsize=12
	)

	xlims!(axtop, low=0, high=200)
	ylims!(axtop, low=0, high=200)
	axtop.xticks = [0, 50, 100, 150, 200]
	axtop.yticks = [0, 50, 100, 150, 200]
	axtop.xlabel = "Known Mass (mg)"
	axtop.ylabel = "Calculated Mass (mg)"
	axtop.title = "Integrated (1 Point Calibration)"
	hidedecorations!(axtop, ticklabels=false, ticks=false, label=false)
	
	##-- B --##
	axtopright = Axis(f[2, 1])
	
	df3 = df_i2
	sc1=scatter!(axtopright, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large])
	errorbars!(axtopright, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], rms(df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large]))
	sc2=scatter!(axtopright, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium])
	errorbars!(axtopright, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], rms(df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium]))
	sc3=scatter!(axtopright, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], color=:red)
	errorbars!(axtopright, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], rms(df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small]))
	ln1=lines!(axtopright, [-1000, 1000], [-1000, 1000])
	ln2=lines!(axtopright, collect(1:1000), pred_i2, linestyle=:dashdot)
	Textbox(
		f[2, 1], 
		placeholder = "y = $(trunc(co2[2]; digits=3))x + $(trunc(co2[1]; digits=3)) \nr = $(trunc(r2_2; digits=3)) \nRMSE: $(trunc(rms_values2[1]; digits=3)) \nRMSD: $(trunc(rms_values2[2]; digits=3))", 
		tellheight = false,
        tellwidth = false,
		boxcolor=:white,
		halign=:left,
		valign=:top,
		textsize=12
	)
	
	xlims!(axtopright, low=0, high=200)
	ylims!(axtopright, low=0, high=200)
	axtopright.xticks = [0, 50, 100, 150, 200]
	axtopright.yticks = [0, 50, 100, 150, 200]
	axtopright.xlabel = "Known Mass (mg)"
	axtopright.ylabel = "Calculated Mass (mg)"
	axtopright.title = "Integrated (3 Point Calibration)"
	hidedecorations!(axtopright, ticklabels=false, ticks=false, label=false)
	
	#-- LABELS --##

	f[1:2, 2] = Legend(f, [sc1, sc2, sc3, ln1, ln2], ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"], framevisible = false)

	
	for (label, layout) in zip(["A", "B"], [f[1,1], f[2,1]])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end

	save("/Users/daleblack/Google Drive/Research/2022-AAPM/linear_reg_norm.png", f)
	f
end

# ╔═╡ 8ff85b4b-96c2-4804-8284-5d64bc690b13
md"""
## Agatson vs Integrated
"""

# ╔═╡ 0df839b7-6d1a-4388-a612-f54167debbf3
df_a = CSV.read(string(root_path, "agatston.csv"), DataFrame);

# ╔═╡ 10fbfb54-27c7-45bf-adf0-1e2501fcaa08
let
	df = df_a
	gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	global model_a
	model_a = lm(@formula(Y ~ X), data)
	global r2a
	r2a = GLM.r2(model_a)
	global rms_valuesa
	rms_valuesa = [
		rms(data[!, :X], data[!, :Y]),
		rmsd(data[!, :Y], GLM.predict(model_a))
	]
end

# ╔═╡ 119b44d7-3d82-4954-b16d-6ffa9c7e2d7b
begin
	newX3 = DataFrame(X=collect(1:1000));
	pred_a = GLM.predict(model_a, newX3)
end

# ╔═╡ 911c757b-272e-4b00-97b8-d62d4c90a998
co3 = coef(model_a)

# ╔═╡ f8815a69-75cb-42be-9ad8-bfebffc8d68c
let
	f = Figure()

	##-- A --##
	axtop = Axis(f[1, 1])
	
	df = df_a
	scatter!(axtop, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
	errorbars!(axtop, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], rms(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large]))
	scatter!(axtop, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
	errorbars!(axtop, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], rms(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium]))
	scatter!(axtop, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
	errorbars!(axtop, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], rms(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small]))
	lines!(axtop, [-1000, 1000], [-1000, 1000],)
	lines!(axtop, collect(1:1000), pred_a, linestyle=:dashdot)
	Textbox(
		f[1, 1], 
		placeholder = "y = $(trunc(co3[2]; digits=3))x + $(trunc(co3[1]; digits=3)) \nr = $(trunc(r2a; digits=3)) \nRMSE: $(trunc(rms_valuesa[1]; digits=3)) \nRMSD: $(trunc(rms_valuesa[2]; digits=3))", 
		tellheight = false,
        tellwidth = false,
		boxcolor=:white,
		halign=:left,
		valign=:top,
		textsize=12
	)

	xlims!(axtop, low=0, high=200)
	ylims!(axtop, low=0, high=200)
	axtop.xticks = [0, 50, 100, 150, 200]
	axtop.yticks = [0, 50, 100, 150, 200]
	axtop.xlabel = "Known Mass (mg)"
	axtop.ylabel = "Calculated Mass (mg)"
	axtop.title = "Agatston"
	hidedecorations!(axtop, ticklabels=false, ticks=false, label=false)
	
	##-- B --##
	axtopright = Axis(f[2, 1])
	
	df3 = df_i2
	sc1=scatter!(axtopright, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large])
	errorbars!(axtopright, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], rms(df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large]))
	sc2=scatter!(axtopright, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium])
	errorbars!(axtopright, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], rms(df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium]))
	sc3=scatter!(axtopright, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], color=:red)
	errorbars!(axtopright, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], rms(df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small]))
	ln1=lines!(axtopright, [-1000, 1000], [-1000, 1000])
	ln2=lines!(axtopright, collect(1:1000), pred_i2, linestyle=:dashdot)
	Textbox(
		f[2, 1], 
		placeholder = "y = $(trunc(co2[2]; digits=3))x + $(trunc(co2[1]; digits=3)) \nr = $(trunc(r2_2; digits=3)) \nRMSE: $(trunc(rms_values2[1]; digits=3)) \nRMSD: $(trunc(rms_values2[2]; digits=3))", 
		tellheight = false,
        tellwidth = false,
		boxcolor=:white,
		halign=:left,
		valign=:top,
		textsize=12
	)
	
	xlims!(axtopright, low=0, high=200)
	ylims!(axtopright, low=0, high=200)
	axtopright.xticks = [0, 50, 100, 150, 200]
	axtopright.yticks = [0, 50, 100, 150, 200]
	axtopright.xlabel = "Known Mass (mg)"
	axtopright.ylabel = "Calculated Mass (mg)"
	axtopright.title = "Integrated (3 Point Calibration)"
	hidedecorations!(axtopright, ticklabels=false, ticks=false, label=false)
	
	#-- LABELS --##

	f[1:2, 2] = Legend(f, [sc1, sc2, sc3, ln1, ln2], ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"], framevisible = false)

	
	for (label, layout) in zip(["A", "B"], [f[1,1], f[2,1]])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end

	save("/Users/daleblack/Google Drive/Research/2022-AAPM/linear_reg_norm.png", f)
	f
end


# ╔═╡ 196fc7a5-b5d3-444f-8d6e-196541e035a9
md"""
## Small Inserts
"""

# ╔═╡ fbdd0e52-0c92-45d8-a5f4-e191eb387523
md"""
#### Agatston
"""

# ╔═╡ 63837a92-18bd-443e-8b31-e3f8ee7b1234
let
	df = df_a
	gt_array = vec(df[!, :ground_truth_mass_small])
	calc_array = vec(df[!, :calculated_mass_small])
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	global model_a_small
	model_a_small = lm(@formula(Y ~ X), data)
	global r2a_s
	r2a_s = GLM.r2(model_a)
	global rms_valuesa_small
	rms_valuesa_small = [
		rms(data[!, :X], data[!, :Y]),
		rmsd(data[!, :Y], GLM.predict(model_a_small))
	]
end

# ╔═╡ 88440b79-bf4e-4c3b-8669-5ab5a193b9e6
begin
	newX4 = DataFrame(X=collect(-1000:1000));
	pred_a_small = GLM.predict(model_a_small, newX4)
end

# ╔═╡ 6b7df142-c10c-4573-a62a-5366eceb7c68
co_a_small = coef(model_a_small)

# ╔═╡ 878ffe43-2c66-45e6-bdc6-2265dcbcdfbe
md"""
#### Integrated
"""

# ╔═╡ b94b3990-9653-4fb3-a44d-52cc11aebf3a
let
	df = df_i2
	gt_array = vec(df[!, :ground_truth_mass_small])
	calc_array = vec(df[!, :calculated_mass_small])
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	global model_i_small
	model_i_small = lm(@formula(Y ~ X), data)
	global r2i_s
	r2i_s = GLM.r2(model_i_small)
	global rms_valuesi_small
	rms_valuesi_small = [
		rms(data[!, :X], data[!, :Y]),
		rmsd(data[!, :Y], GLM.predict(model_i_small))
	]
end

# ╔═╡ 9d197c96-ae44-43e8-9f92-4cef7b77310c
begin
	newX5 = DataFrame(X=collect(-1000:1000));
	pred_i_small = GLM.predict(model_i_small, newX5)
end

# ╔═╡ 208e1303-b344-41b7-a476-a72909132ea0
co_i_small = coef(model_i_small)

# ╔═╡ 22fab6ad-6b31-45e1-9582-c534339f5a26
# f8815a69-75cb-42be-9ad8-bfebffc8d68c
let
	f = Figure()

	##-- A --##
	axtop = Axis(f[1, 1])
	
	df = df_a
	# scatter!(axtop, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
	# errorbars!(axtop, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], rms(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large]))
	# scatter!(axtop, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
	# errorbars!(axtop, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], rms(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium]))
	scatter!(axtop, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
	errorbars!(axtop, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], rms(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small]))
	lines!(axtop, collect(-1000:1000), collect(-1000:1000))
	lines!(axtop, collect(-1000:1000), pred_a_small, linestyle=:dashdot)
	Textbox(
		f[1, 1], 
		placeholder = "y = $(trunc(co_a_small[2]; digits=3))x + $(trunc(co_a_small[1]; digits=3)) \nr = $(trunc(r2a_s; digits=3)) \nRMSE: $(trunc(rms_valuesa_small[1]; digits=3)) \nRMSD: $(trunc(rms_valuesa_small[2]; digits=3))", 
		tellheight = false,
        tellwidth = false,
		boxcolor=:white,
		halign=:left,
		valign=:top,
		textsize=12
	)

	xlims!(axtop, low=-2, high=3)
	ylims!(axtop, low=-2, high=3)
	# axtop.xticks = [0, 50, 100, 150, 200]
	# axtop.yticks = [0, 50, 100, 150, 200]
	axtop.xlabel = "Known Mass (mg)"
	axtop.ylabel = "Calculated Mass (mg)"
	axtop.title = "Agatston"
	hidedecorations!(axtop, ticklabels=false, ticks=false, label=false)
	
	##-- B --##
	axtopright = Axis(f[2, 1])
	
	df3 = df_i2
	# sc1=scatter!(axtopright, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large])
	# errorbars!(axtopright, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], rms(df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large]))
	# sc2=scatter!(axtopright, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium])
	# errorbars!(axtopright, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], rms(df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium]))
	sc3=scatter!(axtopright, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], color=:red)
	errorbars!(axtopright, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], rms(df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small]))
	ln1=lines!(axtopright, [-1000, 1000], [-1000, 1000])
	ln2=lines!(axtopright, collect(-1000:1000), pred_i_small, linestyle=:dashdot)
	Textbox(
		f[2, 1], 
		placeholder = "y = $(trunc(co_i_small[2]; digits=3))x + $(trunc(co_i_small[1]; digits=3)) \nr = $(trunc(r2i_s; digits=3)) \nRMSE: $(trunc(rms_valuesi_small[1]; digits=3)) \nRMSD: $(trunc(rms_valuesi_small[2]; digits=3))", 
		tellheight = false,
        tellwidth = false,
		boxcolor=:white,
		halign=:left,
		valign=:top,
		textsize=12
	)
	
	xlims!(axtopright, low=-2, high=3)
	ylims!(axtopright, low=-2, high=3)
	# axtopright.xticks = [0, 50, 100, 150, 200]
	# axtopright.yticks = [0, 50, 100, 150, 200]
	axtopright.xlabel = "Known Mass (mg)"
	axtopright.ylabel = "Calculated Mass (mg)"
	axtopright.title = "Integrated (3 Point Calibration)"
	hidedecorations!(axtopright, ticklabels=false, ticks=false, label=false)
	
	#-- LABELS --##

	f[1:2, 2] = Legend(f, [sc3, ln1, ln2], ["Small Inserts", "Unity", "Fitted Line"], framevisible = false)

	
	for (label, layout) in zip(["A", "B"], [f[1,1], f[2,1]])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end

	save("/Users/daleblack/Google Drive/Research/2022-AAPM/linear_reg_norm.png", f)
	f
end



# ╔═╡ Cell order:
# ╠═6b400190-a550-4e32-baa6-d54b66d7676e
# ╠═38cd4959-e50a-4724-849d-9a246dd4414e
# ╠═787dda66-c0c9-4e2f-8646-c1c5ac95feee
# ╟─5017031f-2e8d-44f9-925d-a14376050342
# ╟─d017319a-34aa-4d81-b230-203a43c66602
# ╟─61a077c9-6bc2-450f-a58f-6622ca46a592
# ╠═1a7206ea-4be5-4236-954e-a927f9cba993
# ╠═9e71a277-0a70-483d-ab5f-2269ba244ca4
# ╠═83f1c3d3-f566-4c69-8eaf-8c363286e031
# ╠═1fe58e55-eafb-491a-8402-b3b07ac3ca2f
# ╟─37284bae-a062-4620-ad53-66de09ad9322
# ╠═19713b20-4c52-488b-9f0d-13d9e6599645
# ╠═9f9c0080-c311-442f-b4d8-1efd0552845b
# ╠═97c05609-c9a0-41e6-9f21-2d8dfac412aa
# ╠═736e7c7d-d53d-4706-ba8f-c83c6fde3d14
# ╟─8ff85b4b-96c2-4804-8284-5d64bc690b13
# ╟─f8815a69-75cb-42be-9ad8-bfebffc8d68c
# ╠═0df839b7-6d1a-4388-a612-f54167debbf3
# ╠═10fbfb54-27c7-45bf-adf0-1e2501fcaa08
# ╠═119b44d7-3d82-4954-b16d-6ffa9c7e2d7b
# ╠═911c757b-272e-4b00-97b8-d62d4c90a998
# ╟─196fc7a5-b5d3-444f-8d6e-196541e035a9
# ╟─22fab6ad-6b31-45e1-9582-c534339f5a26
# ╟─fbdd0e52-0c92-45d8-a5f4-e191eb387523
# ╠═63837a92-18bd-443e-8b31-e3f8ee7b1234
# ╠═88440b79-bf4e-4c3b-8669-5ab5a193b9e6
# ╠═6b7df142-c10c-4573-a62a-5366eceb7c68
# ╟─878ffe43-2c66-45e6-bdc6-2265dcbcdfbe
# ╠═b94b3990-9653-4fb3-a44d-52cc11aebf3a
# ╠═9d197c96-ae44-43e8-9f92-4cef7b77310c
# ╠═208e1303-b344-41b7-a476-a72909132ea0

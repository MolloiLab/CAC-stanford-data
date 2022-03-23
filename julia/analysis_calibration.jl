### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 6b400190-a550-4e32-baa6-d54b66d7676e
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
end

# ╔═╡ 38cd4959-e50a-4724-849d-9a246dd4414e
TableOfContents()

# ╔═╡ c625ccbb-66ae-4a17-b72e-cbc1e20e00a9
md"""
# CANON
"""

# ╔═╡ 5017031f-2e8d-44f9-925d-a14376050342
md"""
## Calibration 1
For this dataset, the calibration approach used just the mean intensity of the 200 mg/cc calibration rod. This intensity was then used to calculate a "false" volume for each insert. Then using the known density of the calibration rod, a "true" mass was calculated
"""

# ╔═╡ 787dda66-c0c9-4e2f-8646-c1c5ac95feee
canon_path1 = "/Users/daleblack/Google Drive/dev/MolloiLab/CAC-stanford-data/data/output/Canon_Aquilion_One_Vision4";

# ╔═╡ 1a7206ea-4be5-4236-954e-a927f9cba993
df_i_canon1 = CSV.read(string(canon_path1, "/full.csv"), DataFrame);

# ╔═╡ 955d74f0-2ee9-40e9-9275-01555e13daa2
begin
	fInt1 = Figure()
	axInt1 = Axis(fInt1[1, 1])

	scatter!(df_i_canon1[!, :ground_truth_mass_large], df_i_canon1[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df_i_canon1[!, :ground_truth_mass_medium], df_i_canon1[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df_i_canon1[!, :ground_truth_mass_small], df_i_canon1[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity")
	
	axInt1.title = "IHU Mass vs Known Mass"
	axInt1.ylabel = "Calculated Mass (mg)"
	axInt1.xlabel = "Known Mass (mg)"
	axInt1.xticks = [0, 25, 50, 75, 100]
	axInt1.yticks = [0, 25, 50, 75, 100]

	xlims!(axInt1, 0, 100)
	ylims!(axInt1, 0, 100)
	
	fInt1[1, 2] = Legend(fInt1, axInt1, framevisible = false)
	fInt1
end

# ╔═╡ 8d1452f9-032c-41df-830e-ab6459e7c5ae
begin
	fInt1_small = Figure()
	axInt1_small = Axis(fInt1_small[1, 1])

	scatter!(df_i_canon1[!, :ground_truth_mass_small], df_i_canon1[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 1], [0, 1], label="Unity")
	
	axInt1_small.title = "IHU Mass vs Known Mass (Small)"
	axInt1_small.ylabel = "Calculated Mass (mg)"
	axInt1_small.xlabel = "Known Mass (mg)"
	
	fInt1_small[1, 2] = Legend(fInt1_small, axInt1_small, framevisible = false)
	fInt1_small
end

# ╔═╡ 067b9f6b-6b57-40c0-ae5f-0465400f77bb
md"""
### RMSD measurements (large, medium, and small inserts)
"""

# ╔═╡ 7f4d43df-f51e-49c2-b47e-629e3ccedc91
int1_rms_large_c = rmsd(df_i_canon1[!, :calculated_mass_large], df_i_canon1[!, :ground_truth_mass_large])

# ╔═╡ 39886e68-c4ea-4b96-96fb-ba46eeb54ecd
int1_rms_medium_c = rmsd(df_i_canon1[!, :calculated_mass_medium], df_i_canon1[!, :ground_truth_mass_medium])

# ╔═╡ fdad40dc-5e3c-4768-837f-7dbf4fca07ac
int1_rms_small_c = rmsd(df_i_canon1[!, :calculated_mass_small], df_i_canon1[!, :ground_truth_mass_small])

# ╔═╡ 02bb14eb-a269-4d4b-8151-db0d37c0c269
md"""
## Calibration 2
For this dataset, the calibration approach used the mean intensity of the 200 mg/cc calibration rod then calculated a line of best fit assuming the intercept is (0, 0) for a 0 mg/cc rod. This line of best fit was then used to calculate the "true" object intensity of each insert assuming we already know the density of the insert of interest. This allows for the calculation of a "true" volume and mass directly, but requires that we artificially know the density of each insert beforehand
"""

# ╔═╡ 41ef4353-b874-4848-b266-0a5c0e79ce9f
canon_path2 = "/Users/daleblack/Google Drive/dev/MolloiLab/CAC-stanford-data/data/output/Canon_Aquilion_One_Vision3";

# ╔═╡ 4eeb61ce-9e72-433d-829a-1651cc20b196
df_i_canon2 = CSV.read(string(canon_path2, "/full.csv"), DataFrame);

# ╔═╡ 9124a9c9-4553-49a2-8fbe-f9882a8b2c69
begin
	fInt2 = Figure()
	axInt2 = Axis(fInt2[1, 1])

	scatter!(df_i_canon2[!, :ground_truth_mass_large], df_i_canon2[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df_i_canon2[!, :ground_truth_mass_medium], df_i_canon2[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df_i_canon2[!, :ground_truth_mass_small], df_i_canon2[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity")
	
	axInt2.title = "IHU Mass vs Known Mass"
	axInt2.ylabel = "Calculated Mass (mg)"
	axInt2.xlabel = "Known Mass (mg)"
	axInt2.xticks = [0, 25, 50, 75, 100]
	axInt2.yticks = [0, 25, 50, 75, 100]

	xlims!(axInt2, 0, 100)
	ylims!(axInt2, 0, 100)
	
	fInt2[1, 2] = Legend(fInt2, axInt2, framevisible = false)
	fInt2
end

# ╔═╡ 246b0e98-8325-4031-9884-751a39454e40
begin
	fInt2_small = Figure()
	axInt2_small = Axis(fInt2_small[1, 1])
	
	scatter!(df_i_canon2[!, :ground_truth_mass_small], df_i_canon2[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 1], [0, 1], label="Unity")
	
	axInt2_small.title = "IHU Mass vs Known Mass (Small)"
	axInt2_small.ylabel = "Calculated Mass (mg)"
	axInt2_small.xlabel = "Known Mass (mg)"
	
	fInt2_small[1, 2] = Legend(fInt2_small, axInt2_small, framevisible = false)
	fInt2_small
end

# ╔═╡ 8c4745d4-3dae-46ca-9673-9ab1cb84deb8
md"""
### RMSD measurements (large, medium, and small inserts)
"""

# ╔═╡ 7015f1cf-ee36-4cbe-8daa-804de3e198ec
int2_rms_large_c = rmsd(df_i_canon2[!, :calculated_mass_large], df_i_canon2[!, :ground_truth_mass_large])

# ╔═╡ 65e2dca0-1835-4c1e-b492-1c3fe210663d
int2_rms_medium_c = rmsd(df_i_canon2[!, :calculated_mass_medium], df_i_canon2[!, :ground_truth_mass_medium])

# ╔═╡ 2a30236f-1c78-40f8-80d2-8437171c93a1
int2_rms_small_c = rmsd(df_i_canon2[!, :calculated_mass_small], df_i_canon2[!, :ground_truth_mass_small])

# ╔═╡ 565fab58-ef10-43e8-8ba8-e2e792ef872f
md"""
## Agatston Mass Scoring
"""

# ╔═╡ 09bfb23a-db24-42c8-8c34-8234da9c0d56
canon_path_agat = "/Users/daleblack/Google Drive/dev/MolloiLab/CAC-stanford-data/data/output/Canon_Aquilion_One_Vision_agatston_mass_score_script";

# ╔═╡ 5c1bea23-972b-4867-a424-9119cbc5edbc
df_a_canon = CSV.read(string(canon_path_agat, "/output.csv"), DataFrame);

# ╔═╡ 65131928-ad9a-4649-9600-e663ab5a7a28
begin
	fAgat = Figure()
	axAgat = Axis(fAgat[1, 1])

	scatter!(df_a_canon[!, :ground_truth_mass_large], df_a_canon[!, :calculated_mass_large], label="Mass: Large Inserts")
	scatter!(df_a_canon[!, :ground_truth_mass_medium], df_a_canon[!, :calculated_mass_medium], label="Mass: Medium Inserts")
	scatter!(df_a_canon[!, :ground_truth_mass_small], df_a_canon[!, :calculated_mass_small], label="Mass: Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity Line")
	
	axAgat.title = "AS Mass vs Known Mass"
	axAgat.ylabel = "Calculated Mass (mg)"
	axAgat.xlabel = "Known Mass (mg)"
	axAgat.xticks = [0, 25, 50, 75, 100]
	axAgat.yticks = [0, 25, 50, 75, 100]

	xlims!(axAgat, 0, 100)
	ylims!(axAgat, 0, 100)
	
	fAgat[1, 2] = Legend(fAgat, axAgat, framevisible = false)

	fAgat
end

# ╔═╡ c4bc4b99-7adc-48a2-a3b9-114d6c5a4e3a
begin
	fAgat_small = Figure()
	axAgat_small = Axis(fAgat_small[1, 1])

	scatter!(df_a_canon[!, :ground_truth_mass_small], df_a_canon[!, :calculated_mass_small], label="Mass: Small Inserts", color=:red)
	lines!([0, 1], [0, 1], label="Unity Line")
	
	axAgat_small.title = "AS Mass vs Known Mass (Small)"
	axAgat_small.ylabel = "Calculated Mass (mg)"
	axAgat_small.xlabel = "Known Mass (mg)"
	
	fAgat_small[1, 2] = Legend(fAgat_small, axAgat_small, framevisible = false)

	fAgat_small
end

# ╔═╡ c6cd7f7d-d77b-4136-9485-913f962be85a
agat_rms_large_c = rmsd(df_a_canon[!, :calculated_mass_large], df_a_canon[!, :ground_truth_mass_large])

# ╔═╡ 7f8756e1-36e0-444a-a9a9-e20f763298af
agat_rms_medium_c = rmsd(df_a_canon[!, :calculated_mass_medium], df_a_canon[!, :ground_truth_mass_medium])

# ╔═╡ ec955132-42a3-40cd-a031-480ee6e95545
agat_rms_small_c = rmsd(df_a_canon[!, :calculated_mass_small], df_a_canon[!, :ground_truth_mass_small])

# ╔═╡ Cell order:
# ╟─6b400190-a550-4e32-baa6-d54b66d7676e
# ╟─38cd4959-e50a-4724-849d-9a246dd4414e
# ╟─c625ccbb-66ae-4a17-b72e-cbc1e20e00a9
# ╟─5017031f-2e8d-44f9-925d-a14376050342
# ╠═787dda66-c0c9-4e2f-8646-c1c5ac95feee
# ╠═1a7206ea-4be5-4236-954e-a927f9cba993
# ╟─955d74f0-2ee9-40e9-9275-01555e13daa2
# ╟─8d1452f9-032c-41df-830e-ab6459e7c5ae
# ╟─067b9f6b-6b57-40c0-ae5f-0465400f77bb
# ╠═7f4d43df-f51e-49c2-b47e-629e3ccedc91
# ╠═39886e68-c4ea-4b96-96fb-ba46eeb54ecd
# ╠═fdad40dc-5e3c-4768-837f-7dbf4fca07ac
# ╟─02bb14eb-a269-4d4b-8151-db0d37c0c269
# ╠═41ef4353-b874-4848-b266-0a5c0e79ce9f
# ╠═4eeb61ce-9e72-433d-829a-1651cc20b196
# ╟─9124a9c9-4553-49a2-8fbe-f9882a8b2c69
# ╟─246b0e98-8325-4031-9884-751a39454e40
# ╟─8c4745d4-3dae-46ca-9673-9ab1cb84deb8
# ╠═7015f1cf-ee36-4cbe-8daa-804de3e198ec
# ╠═65e2dca0-1835-4c1e-b492-1c3fe210663d
# ╠═2a30236f-1c78-40f8-80d2-8437171c93a1
# ╟─565fab58-ef10-43e8-8ba8-e2e792ef872f
# ╠═09bfb23a-db24-42c8-8c34-8234da9c0d56
# ╠═5c1bea23-972b-4867-a424-9119cbc5edbc
# ╟─65131928-ad9a-4649-9600-e663ab5a7a28
# ╟─c4bc4b99-7adc-48a2-a3b9-114d6c5a4e3a
# ╠═c6cd7f7d-d77b-4136-9485-913f962be85a
# ╠═7f8756e1-36e0-444a-a9a9-e20f763298af
# ╠═ec955132-42a3-40cd-a031-480ee6e95545

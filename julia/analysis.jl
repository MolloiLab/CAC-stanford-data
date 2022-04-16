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

# ╔═╡ 5017031f-2e8d-44f9-925d-a14376050342
md"""
## Integrated Mass Scoring
"""

# ╔═╡ 787dda66-c0c9-4e2f-8646-c1c5ac95feee
canon_path1 = "/Users/daleblack/Google Drive/dev/MolloiLab/CAC-stanford-data/data/output/Canon_Aquilion_One_Vision";

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
	lines!([0, 1.2], [0, 1.2], label="Unity")
	
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

# ╔═╡ 243c0c2b-e84d-47eb-aaf8-5fb79d76d04d
md"""
### RMSD measurements (large, medium, and small inserts)
"""

# ╔═╡ c6cd7f7d-d77b-4136-9485-913f962be85a
agat_rms_large_c = rmsd(df_a_canon[!, :calculated_mass_large], df_a_canon[!, :ground_truth_mass_large])

# ╔═╡ 7f8756e1-36e0-444a-a9a9-e20f763298af
agat_rms_medium_c = rmsd(df_a_canon[!, :calculated_mass_medium], df_a_canon[!, :ground_truth_mass_medium])

# ╔═╡ ec955132-42a3-40cd-a031-480ee6e95545
agat_rms_small_c = rmsd(df_a_canon[!, :calculated_mass_small], df_a_canon[!, :ground_truth_mass_small])

# ╔═╡ Cell order:
# ╠═6b400190-a550-4e32-baa6-d54b66d7676e
# ╠═38cd4959-e50a-4724-849d-9a246dd4414e
# ╟─5017031f-2e8d-44f9-925d-a14376050342
# ╠═787dda66-c0c9-4e2f-8646-c1c5ac95feee
# ╠═1a7206ea-4be5-4236-954e-a927f9cba993
# ╟─955d74f0-2ee9-40e9-9275-01555e13daa2
# ╟─8d1452f9-032c-41df-830e-ab6459e7c5ae
# ╟─067b9f6b-6b57-40c0-ae5f-0465400f77bb
# ╠═7f4d43df-f51e-49c2-b47e-629e3ccedc91
# ╠═39886e68-c4ea-4b96-96fb-ba46eeb54ecd
# ╠═fdad40dc-5e3c-4768-837f-7dbf4fca07ac
# ╟─565fab58-ef10-43e8-8ba8-e2e792ef872f
# ╠═09bfb23a-db24-42c8-8c34-8234da9c0d56
# ╠═5c1bea23-972b-4867-a424-9119cbc5edbc
# ╟─65131928-ad9a-4649-9600-e663ab5a7a28
# ╟─c4bc4b99-7adc-48a2-a3b9-114d6c5a4e3a
# ╟─243c0c2b-e84d-47eb-aaf8-5fb79d76d04d
# ╠═c6cd7f7d-d77b-4136-9485-913f962be85a
# ╠═7f8756e1-36e0-444a-a9a9-e20f763298af
# ╠═ec955132-42a3-40cd-a031-480ee6e95545

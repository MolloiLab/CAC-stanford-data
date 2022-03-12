### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 9d7e1f1b-6a0c-44f2-994e-bed095cef076
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
		Pkg.add("CairoMakie")
		Pkg.add("Statistics")
		Pkg.add("StatsBase")
		Pkg.add("ImageMorphology")
		Pkg.add("ImageFiltering")
		Pkg.add("CSV")
		Pkg.add("DataFrames")
		Pkg.add("GLM")
		Pkg.add(url="https://github.com/JuliaHealth/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
		Pkg.add(url="https://github.com/Dale-Black/Phantoms.jl")
		Pkg.add(url="https://github.com/Dale-Black/CalciumScoring.jl")
	end
	
	using PlutoUI
	using CairoMakie
	using Statistics
	using StatsBase: quantile!
	using ImageMorphology
	using ImageFiltering
	using CSV
	using DataFrames
	using GLM
	using DICOM
	using DICOMUtils
	using Phantoms
	using CalciumScoring
end

# ╔═╡ 5fe79c24-a25d-11ec-2b6a-0f08622c9bf0
function mass_mean_def(array, calcium_threshold)
    # Calculation of mean value within masked image
    if maximum(array) == 0
        mass_mean = 0
    else
        array = Int.(array .> calcium_threshold)
		mass_mean = sum(array) / length(array .!= 0)	
	end
	return mass_mean
end

# ╔═╡ Cell order:
# ╠═9d7e1f1b-6a0c-44f2-994e-bed095cef076
# ╠═5fe79c24-a25d-11ec-2b6a-0f08622c9bf0

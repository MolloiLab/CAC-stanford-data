### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 52c184d2-9a5e-11ec-0e10-976cdfa5f253
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
		Pkg.add("CairoMakie")
		Pkg.add("MAT")
		Pkg.add(url="https://github.com/JuliaHealth/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
	end

	using PlutoUI
	using CairoMakie
	using MAT
	using DICOM
	using DICOMUtils
end

# ╔═╡ 5abcc458-7733-4ce4-a56f-bb7cb8882540
TableOfContents()

# ╔═╡ 1417f262-16f4-41f7-bed6-c895fb49fa80
file = 3

# ╔═╡ 3c5bfe99-3098-43b6-9d0e-0b3ed5884dea
file_inserts = 6

# ╔═╡ 2ae4e5ff-3cfd-4c80-a59a-e91c52480ecf
md"""
## Convert .mat to Array
"""

# ╔═╡ 0fab4688-9f2c-4df6-9cff-d3f9086e7a4b
begin
	path = "/Users/daleblack/Desktop/s120.mat"
	vars1 = matread(path)
	array1 = vars1["r120"]
	array1 = Int16.(round.(array1))
end;

# ╔═╡ a9fcf3de-496f-4ce6-8c52-5682f69cdceb
heatmap(transpose(array1), colormap=:grays)

# ╔═╡ 506a5e6b-94ed-4502-af57-0092eba7d817
begin
	path2 = "/Users/daleblack/Desktop/r120.mat";
	vars2 = matread(path2);
	array2 = vars2["s120"]
	array2 = Int16.(round.(array2))
end;

# ╔═╡ 65c27fd3-3959-4d87-a120-2be6a221ca87
heatmap(transpose(array2), colormap=:grays)

# ╔═╡ 42bdac01-ab21-43ea-928e-04231f0428ba
md"""
## Create DICOM rod image
"""

# ╔═╡ 7ef29728-2fd5-4c70-ae42-4d35ce6b45da
dcm_path = "/Users/daleblack/Google Drive/Datasets/Canon_Aquilion_One_Vision/Large_rep1/96E1EB4F";

# ╔═╡ a70027f1-5fbf-4e9d-8d09-4dd8718fc081
begin
	dcm = dcm_parse(dcm_path)
	dcm[tag"Pixel Data"] = array1
	dcm[tag"Instance Number"] = file
	dcm[tag"Rows"] = size(array1, 1)
	dcm[tag"Columns"] = size(array1, 2)
end;

# ╔═╡ 4062b5d4-9fe7-41f5-92e3-fb19cb6a1274
begin
	output_root_rod = "/Users/daleblack/Google Drive/Datasets/Simulated/rod"
	output_path_rod = string(output_root_rod, "/", file, ".dcm")
end

# ╔═╡ 042bc6e8-edd2-4d56-b3e4-a0ec6ee94ffe
dcm_write(output_path_rod, dcm)

# ╔═╡ 68ab4f72-1d54-444d-b8a8-7d9b96108f36
md"""
## Create DICOM inserts image
"""

# ╔═╡ 6b558de5-751d-4fda-b6ab-c90ca4f5a6e6
begin
	dcm2 = dcm_parse(dcm_path)
	dcm2[tag"Pixel Data"] = array2
	dcm2[tag"Instance Number"] = file_inserts
	dcm2[tag"Rows"] = size(array2, 1)
	dcm2[tag"Columns"] = size(array2, 2)
end;

# ╔═╡ 3371dd16-6a49-430f-afad-7d94c9bc63ad
begin
	output_root_inserts = "/Users/daleblack/Google Drive/Datasets/Simulated/inserts"
	output_path_inserts = string(output_root_inserts, "/", file_inserts, ".dcm")
end;

# ╔═╡ b822d070-167e-4ce4-a876-57ba142e2751
dcm_write(output_path_inserts, dcm2)

# ╔═╡ f78c8cb1-6f5a-439a-873f-68c71c95516f
md"""
## Check DICOM image(s)
"""

# ╔═╡ 0ca2e61d-83be-42b4-beb1-359e8eb12793
md"""
### Inserts
"""

# ╔═╡ e8000d5d-1da4-4fb1-a253-de0838a50fcb
dcmdir_inserts = dcmdir_parse(output_root_inserts);

# ╔═╡ 5e1bb447-9100-4460-a0db-b90eea55b7ce
vol_inserts = load_dcm_array(dcmdir_inserts);

# ╔═╡ 0288226a-4c70-492b-b8ab-15701069f594
@bind a PlutoUI.Slider(1:size(vol_inserts, 3); default=1, show_value=true)

# ╔═╡ dcc8ae65-713c-495a-bc2b-61bd514915f7
heatmap(transpose(vol_inserts[:, :, a]), colormap=:grays)

# ╔═╡ 1a1ef5d1-ba40-4f96-bc48-707759ce96f4
md"""
### Rod
"""

# ╔═╡ 1c1ad710-da52-4f5d-835d-853ab539eeab
dcmdir_rod = dcmdir_parse(output_root_rod);

# ╔═╡ 1773556b-66e8-4786-a003-14130ad06086
vol_rod = load_dcm_array(dcmdir_rod);

# ╔═╡ a6454ba1-2082-4f93-9383-9b07c8b30798
@bind b PlutoUI.Slider(1:size(vol_rod, 3); default=1, show_value=true)

# ╔═╡ 21aaf77f-354d-43d4-83ee-743760e698f0
heatmap(transpose(vol_rod[:, :, b]), colormap=:grays)

# ╔═╡ 3467f536-fbac-4542-8da2-ab6396de541f
md"""
## Combined
"""

# ╔═╡ 719a8d5c-ef3a-4859-9fc9-7c02c911c71d
# begin
# 	pth_root = "/Users/daleblack/Google Drive/Datasets/Simulated/"
# 	pth = string(pth_root, "combined")
# 	if ~isdir(pth)
# 		mkdir(pth)
# 	end
	
# 	pth_inserts = string(pth_root, inserts)
# 	files_inserts = readdir(pth_inserts)
# 	for i in files_inserts
# 		string(pth_inserts, )
# end

# ╔═╡ 6d32567c-6d6b-4a96-8ce1-a296eab482d4
output_root_combined = "/Users/daleblack/Google Drive/Datasets/Simulated/combined"

# ╔═╡ 93d48ee0-83e6-46cd-8567-16273bc269f9
dcmdir_combined = dcmdir_parse(output_root_combined);

# ╔═╡ 81f0b564-52ff-4078-85a7-765375ad99f0
vol_combined = load_dcm_array(dcmdir_combined);

# ╔═╡ a79b50b9-7efb-49b8-8970-ecce09a284d4
@bind c PlutoUI.Slider(1:size(vol_combined, 3); default=1, show_value=true)

# ╔═╡ ec615641-95ee-43f8-9433-587ca30f3551
heatmap(transpose(vol_combined[:, :, c]), colormap=:grays)

# ╔═╡ Cell order:
# ╠═52c184d2-9a5e-11ec-0e10-976cdfa5f253
# ╠═5abcc458-7733-4ce4-a56f-bb7cb8882540
# ╠═1417f262-16f4-41f7-bed6-c895fb49fa80
# ╠═3c5bfe99-3098-43b6-9d0e-0b3ed5884dea
# ╟─2ae4e5ff-3cfd-4c80-a59a-e91c52480ecf
# ╠═0fab4688-9f2c-4df6-9cff-d3f9086e7a4b
# ╠═a9fcf3de-496f-4ce6-8c52-5682f69cdceb
# ╠═506a5e6b-94ed-4502-af57-0092eba7d817
# ╠═65c27fd3-3959-4d87-a120-2be6a221ca87
# ╟─42bdac01-ab21-43ea-928e-04231f0428ba
# ╠═7ef29728-2fd5-4c70-ae42-4d35ce6b45da
# ╠═a70027f1-5fbf-4e9d-8d09-4dd8718fc081
# ╠═4062b5d4-9fe7-41f5-92e3-fb19cb6a1274
# ╠═042bc6e8-edd2-4d56-b3e4-a0ec6ee94ffe
# ╟─68ab4f72-1d54-444d-b8a8-7d9b96108f36
# ╠═6b558de5-751d-4fda-b6ab-c90ca4f5a6e6
# ╠═3371dd16-6a49-430f-afad-7d94c9bc63ad
# ╠═b822d070-167e-4ce4-a876-57ba142e2751
# ╟─f78c8cb1-6f5a-439a-873f-68c71c95516f
# ╟─0ca2e61d-83be-42b4-beb1-359e8eb12793
# ╠═e8000d5d-1da4-4fb1-a253-de0838a50fcb
# ╠═5e1bb447-9100-4460-a0db-b90eea55b7ce
# ╟─0288226a-4c70-492b-b8ab-15701069f594
# ╠═dcc8ae65-713c-495a-bc2b-61bd514915f7
# ╟─1a1ef5d1-ba40-4f96-bc48-707759ce96f4
# ╠═1c1ad710-da52-4f5d-835d-853ab539eeab
# ╠═1773556b-66e8-4786-a003-14130ad06086
# ╟─a6454ba1-2082-4f93-9383-9b07c8b30798
# ╠═21aaf77f-354d-43d4-83ee-743760e698f0
# ╟─3467f536-fbac-4542-8da2-ab6396de541f
# ╠═719a8d5c-ef3a-4859-9fc9-7c02c911c71d
# ╠═6d32567c-6d6b-4a96-8ce1-a296eab482d4
# ╠═93d48ee0-83e6-46cd-8567-16273bc269f9
# ╠═81f0b564-52ff-4078-85a7-765375ad99f0
# ╟─a79b50b9-7efb-49b8-8970-ecce09a284d4
# ╠═ec615641-95ee-43f8-9433-587ca30f3551

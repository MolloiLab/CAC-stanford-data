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
# file = 3

# ╔═╡ 3c5bfe99-3098-43b6-9d0e-0b3ed5884dea
file_inserts = 3

# ╔═╡ 2ae4e5ff-3cfd-4c80-a59a-e91c52480ecf
md"""
## Convert .mat to Array
"""

# ╔═╡ 90f1c6e8-dc60-40ac-bc09-f4ad3d98482e
path1 = "/Users/daleblack/Desktop/r120.mat";

# ╔═╡ aece3048-7d4e-4abf-8b9c-79da131fef58
vars = matread(path1);

# ╔═╡ 0fab4688-9f2c-4df6-9cff-d3f9086e7a4b
# begin
# 	array1 = vars["r120"]
# 	array1 = Int16.(round.(array1))
# 	# array2[2] = -899
# end;

# ╔═╡ a9fcf3de-496f-4ce6-8c52-5682f69cdceb
# heatmap(transpose(array1), colormap=:grays)

# ╔═╡ 506a5e6b-94ed-4502-af57-0092eba7d817
begin
	array2 = vars["s120"]
	array2 = Int16.(round.(array2))
	# array2[2] = -899
end;

# ╔═╡ 65c27fd3-3959-4d87-a120-2be6a221ca87
heatmap(transpose(array2), colormap=:grays)

# ╔═╡ 42bdac01-ab21-43ea-928e-04231f0428ba
md"""
## Create DICOM rod image
"""

# ╔═╡ 7ef29728-2fd5-4c70-ae42-4d35ce6b45da
dcm_path = "/Users/daleblack/Google Drive/Datasets/Canon_Aquilion_One_Vision/Large_rep1/96E1EB4F";

# ╔═╡ 2fec4eed-1320-4ed9-bb00-8e808152fc38
# dcm = dcm_parse(dcm_path);

# ╔═╡ fb9331a9-9fe6-47cc-9e44-f1fc96a974e3
# dcm[tag"Pixel Data"] = array1;

# ╔═╡ 4b33c612-1474-4ff1-b70c-65a3a8c9a4c3
# dcm[tag"Instance Number"] = file;

# ╔═╡ 15e25dff-7a95-4c4e-a06d-1bea2f745721
# dcm[tag"Rows"] = size(array1, 1);

# ╔═╡ 86572058-be5c-481d-9e43-58cb47cb0ec2
# dcm[tag"Columns"] = size(array1, 2);

# ╔═╡ 4062b5d4-9fe7-41f5-92e3-fb19cb6a1274
output_root = "/Users/daleblack/Google Drive/Datasets/Simulated";

# ╔═╡ 4e91f3b7-70dd-45ca-9ffe-fbe13b902208
# output_path = string(output_root, "/", file, ".dcm");

# ╔═╡ 042bc6e8-edd2-4d56-b3e4-a0ec6ee94ffe
# dcm_write(output_path, dcm)

# ╔═╡ 68ab4f72-1d54-444d-b8a8-7d9b96108f36
md"""
## Create DICOM inserts image
"""

# ╔═╡ 6b558de5-751d-4fda-b6ab-c90ca4f5a6e6
dcm2 = dcm_parse(dcm_path);

# ╔═╡ 935a78d3-8621-4f27-a653-721a06742800
dcm2[tag"Pixel Data"] = array2;

# ╔═╡ cae51130-b763-48f7-acdb-9c3b613cb758
dcm2[tag"Instance Number"] = file_inserts;

# ╔═╡ 38454e0b-fa6d-42af-aad7-99505a82bd8c
dcm2[tag"Rows"] = size(array2, 1);

# ╔═╡ 848d05ca-119e-4791-a535-08d6d1b92033
dcm2[tag"Columns"] = size(array2, 2);

# ╔═╡ d2d42f63-6ab3-4c28-98f4-8dc50caf8acb
output_path2 = string(output_root, "/", file_inserts, ".dcm")

# ╔═╡ b822d070-167e-4ce4-a876-57ba142e2751
dcm_write(output_path2, dcm2)

# ╔═╡ f78c8cb1-6f5a-439a-873f-68c71c95516f
md"""
## Check DICOM image(s)
"""

# ╔═╡ b5041f5a-96e9-48ef-a0c0-71c365a12934
# new_dcm = dcm_parse(output_path);

# ╔═╡ 207b9300-90f5-4656-a89e-60475c8c88a2
# new_array = new_dcm[tag"Pixel Data"];

# ╔═╡ 7013d713-0f5d-40d7-be8d-85e8fcdc35c3
# heatmap(transpose(new_array), colormap=:grays)

# ╔═╡ e8000d5d-1da4-4fb1-a253-de0838a50fcb
dcmdir = dcmdir_parse(output_root);

# ╔═╡ 5e1bb447-9100-4460-a0db-b90eea55b7ce
vol = load_dcm_array(dcmdir);

# ╔═╡ 0288226a-4c70-492b-b8ab-15701069f594
@bind a PlutoUI.Slider(1:size(vol, 3); default=1, show_value=true)

# ╔═╡ dcc8ae65-713c-495a-bc2b-61bd514915f7
heatmap(transpose(vol[:, :, a]), colormap=:grays)

# ╔═╡ Cell order:
# ╠═52c184d2-9a5e-11ec-0e10-976cdfa5f253
# ╠═5abcc458-7733-4ce4-a56f-bb7cb8882540
# ╠═1417f262-16f4-41f7-bed6-c895fb49fa80
# ╠═3c5bfe99-3098-43b6-9d0e-0b3ed5884dea
# ╟─2ae4e5ff-3cfd-4c80-a59a-e91c52480ecf
# ╠═90f1c6e8-dc60-40ac-bc09-f4ad3d98482e
# ╠═aece3048-7d4e-4abf-8b9c-79da131fef58
# ╟─0fab4688-9f2c-4df6-9cff-d3f9086e7a4b
# ╟─a9fcf3de-496f-4ce6-8c52-5682f69cdceb
# ╠═506a5e6b-94ed-4502-af57-0092eba7d817
# ╠═65c27fd3-3959-4d87-a120-2be6a221ca87
# ╟─42bdac01-ab21-43ea-928e-04231f0428ba
# ╠═7ef29728-2fd5-4c70-ae42-4d35ce6b45da
# ╟─2fec4eed-1320-4ed9-bb00-8e808152fc38
# ╟─fb9331a9-9fe6-47cc-9e44-f1fc96a974e3
# ╟─4b33c612-1474-4ff1-b70c-65a3a8c9a4c3
# ╟─15e25dff-7a95-4c4e-a06d-1bea2f745721
# ╟─86572058-be5c-481d-9e43-58cb47cb0ec2
# ╠═4062b5d4-9fe7-41f5-92e3-fb19cb6a1274
# ╠═4e91f3b7-70dd-45ca-9ffe-fbe13b902208
# ╟─042bc6e8-edd2-4d56-b3e4-a0ec6ee94ffe
# ╟─68ab4f72-1d54-444d-b8a8-7d9b96108f36
# ╠═6b558de5-751d-4fda-b6ab-c90ca4f5a6e6
# ╠═935a78d3-8621-4f27-a653-721a06742800
# ╠═cae51130-b763-48f7-acdb-9c3b613cb758
# ╠═38454e0b-fa6d-42af-aad7-99505a82bd8c
# ╠═848d05ca-119e-4791-a535-08d6d1b92033
# ╠═d2d42f63-6ab3-4c28-98f4-8dc50caf8acb
# ╠═b822d070-167e-4ce4-a876-57ba142e2751
# ╟─f78c8cb1-6f5a-439a-873f-68c71c95516f
# ╠═b5041f5a-96e9-48ef-a0c0-71c365a12934
# ╠═207b9300-90f5-4656-a89e-60475c8c88a2
# ╠═7013d713-0f5d-40d7-be8d-85e8fcdc35c3
# ╠═e8000d5d-1da4-4fb1-a253-de0838a50fcb
# ╠═5e1bb447-9100-4460-a0db-b90eea55b7ce
# ╟─0288226a-4c70-492b-b8ab-15701069f594
# ╠═dcc8ae65-713c-495a-bc2b-61bd514915f7

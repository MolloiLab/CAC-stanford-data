### A Pluto.jl notebook ###
# v0.17.5

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

# ╔═╡ d21e17f4-701d-11ec-3a26-b76aa97ddfe7
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
		Pkg.add("CairoMakie")
		Pkg.add("Statistics")
		Pkg.add("ImageMorphology")
		Pkg.add(url="https://github.com/JuliaHealth/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
		Pkg.add(url="https://github.com/Dale-Black/Phantoms.jl")
	end
	
	using PlutoUI
	using CairoMakie
	using Statistics
	using ImageMorphology
	using DICOM
	using DICOMUtils
	using Phantoms
end

# ╔═╡ 0652892a-36e0-47ed-ad4f-01ee2542c244
TableOfContents()

# ╔═╡ 01a90e28-55e0-470b-a2f7-ba521c43b074
path = string(cd(pwd, "..") , "/", "data/Large_rep1")

# ╔═╡ 5cec884d-a590-4f39-8f52-074be0c133b3
dcms = dcmdir_parse(path);

# ╔═╡ b3da9520-40f2-43eb-9711-81189fc295bd
dcm_array = load_dcm_array(dcms);

# ╔═╡ 92147219-1fc8-4759-8002-9088a4c50461
header = dcms[1].meta;

# ╔═╡ af38f083-2f60-4d76-b9c0-044e6ad2881e
md"""
## Segment Heart
"""

# ╔═╡ 7a718604-3ec9-4029-8aeb-6efd7fcdb8e9
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ b8dfb742-d3cb-4df9-a68b-80f16ab461b7
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 0de35ee3-e9ab-4b1c-8eb2-549e4fb3191b
heatmap(masked_array[:, :, a], colormap=:grays)

# ╔═╡ c38e0ae7-c055-4443-883b-89e5474a5595
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 15]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig
end

# ╔═╡ fd82d1ef-9099-4028-bb91-d360de732804
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig2
end

# ╔═╡ dc884e2b-dada-4ee9-9451-8ee0f227d3f7
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 23]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig3
end

# ╔═╡ bf741808-b25b-474c-8c33-a00595a9fc54
md"""
## Segment Calcium Rod
"""

# ╔═╡ ecdaa368-6fb0-427a-9ba3-d165ec60058b
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header);

# ╔═╡ 5a285155-5eb5-4c9f-aba2-5c987fc9f141
@bind b PlutoUI.Slider(1:size(calcium_image, 3), default=10, show_value=true)

# ╔═╡ 8b25a31f-6173-4fe2-9f69-9aa0fa861d7c
heatmap(calcium_image[:, :, b], colormap=:grays)

# ╔═╡ f38ee659-eff3-4eb3-b1c2-339aef74c234
md"""
## Segment Calcium Inserts
"""

# ╔═╡ c4221824-889e-429b-b02f-23b7656bd3ab
mask_L_LD, mask_M_LD, mask_S_LD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_HD, mask_M_HD, mask_S_HD = mask_inserts(
            dcm_array, masked_array, header, slice_CCI, center_insert
);

# ╔═╡ 66407b07-7c99-47d1-8021-06810bae6d39
masks = mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD;

# ╔═╡ 2db8c8cc-d286-4edd-ad6b-1fc5541ab0c4
@bind c PlutoUI.Slider(1:size(masks, 3), default=10, show_value=true)

# ╔═╡ 9fd42cf0-b264-4461-90f9-62abf82d29fb
heatmap(masks, colormap=:grays)

# ╔═╡ fcf84719-d2a4-4cf9-a12f-85f7eeb20a4e
md"""
## Overlay Mask Calcium Inserts
"""

# ╔═╡ ac3d2be8-6f4a-44d3-aa01-8a10e34c8add
arr_L_HD = masked_array[:, :, 23:28] .* mask_L_HD;

# ╔═╡ 933cf514-deaf-4c3c-98b8-0ee6e3894fc7
@bind d PlutoUI.Slider(1:size(arr_L_HD, 3), default=2, show_value=true)

# ╔═╡ 9fd03d96-4ba5-4686-b664-44c8bc4754b2
heatmap(arr_L_HD[:, :, d], colormap=:grays)

# ╔═╡ 33223415-a88a-4304-8af0-5ef063d14d40
@bind e PlutoUI.Slider(1:size(dcm_array, 3), default=10, show_value=true)

# ╔═╡ e941a84c-f386-4770-affa-df60bf5dea0b
heatmap(dcm_array[:, :, e], colormap=:grays)

# ╔═╡ 3b1bdc18-cf14-4811-8cba-e20f65bf3188
maximum(masked_array)

# ╔═╡ 8f54662a-efb1-4a2c-9305-53fcfdcd0af0
md"""
## Integrated Calcium Scoring
"""

# ╔═╡ 30cca982-ef9a-4c67-ae53-82d39ff24a55
md"""
### Find Calibration Line
"""

# ╔═╡ 3818e468-db0a-4bca-bbcd-ba4732e0e9df
begin
	arr_L_HD_cal = masked_array[:, :, 25] .* mask_L_HD
	core_HD = erode(erode(erode(arr_L_HD_cal)))
end;

# ╔═╡ 692b4e3d-7aec-4b73-8128-085e5b23d1f5
heatmap(core_HD, colormap=:grays)

# ╔═╡ 9f148b9a-fc4e-4d0d-baed-c37f958afc3e
begin
	arr_L_MD_cal = masked_array[:, :, 25] .* mask_L_MD
	core_MD = erode(erode(erode(arr_L_MD_cal)))
end;

# ╔═╡ 5bdc9b96-8d13-4c55-be09-139682913c9c
heatmap(core_MD, colormap=:grays)

# ╔═╡ 9041e3b6-ccde-4155-bee2-9aab694f5ced
begin
	arr_L_LD_cal = masked_array[:, :, 25] .* mask_L_LD
	core_LD = erode(erode(erode(arr_L_LD_cal)))
end;

# ╔═╡ 3244f883-acec-4dbf-83ad-0c8dfb6326e0
heatmap(core_LD, colormap=:grays)

# ╔═╡ 519e3ed1-8f9a-4ccd-bd6d-115a91aa30d2
cal_array = [mean(core_LD), mean(core_MD), mean(core_HD)]

# ╔═╡ 85f65257-8c39-4836-8517-2109cf55663b
hist(core_HD, bins=100)

# ╔═╡ 61b908cb-41ab-4bd6-9db8-dce88834e43c
scatter(cal_array)

# ╔═╡ 3c5e85b2-fd11-4f70-9e6c-cc26fbcce5b3
md"""
### Large, High Density Insert
"""

# ╔═╡ 547bd984-10e5-4063-85a1-20262f0f1edf


# ╔═╡ Cell order:
# ╠═d21e17f4-701d-11ec-3a26-b76aa97ddfe7
# ╠═0652892a-36e0-47ed-ad4f-01ee2542c244
# ╠═01a90e28-55e0-470b-a2f7-ba521c43b074
# ╠═5cec884d-a590-4f39-8f52-074be0c133b3
# ╠═b3da9520-40f2-43eb-9711-81189fc295bd
# ╠═92147219-1fc8-4759-8002-9088a4c50461
# ╟─af38f083-2f60-4d76-b9c0-044e6ad2881e
# ╠═7a718604-3ec9-4029-8aeb-6efd7fcdb8e9
# ╟─b8dfb742-d3cb-4df9-a68b-80f16ab461b7
# ╠═0de35ee3-e9ab-4b1c-8eb2-549e4fb3191b
# ╠═c38e0ae7-c055-4443-883b-89e5474a5595
# ╠═fd82d1ef-9099-4028-bb91-d360de732804
# ╠═dc884e2b-dada-4ee9-9451-8ee0f227d3f7
# ╟─bf741808-b25b-474c-8c33-a00595a9fc54
# ╠═ecdaa368-6fb0-427a-9ba3-d165ec60058b
# ╟─5a285155-5eb5-4c9f-aba2-5c987fc9f141
# ╠═8b25a31f-6173-4fe2-9f69-9aa0fa861d7c
# ╟─f38ee659-eff3-4eb3-b1c2-339aef74c234
# ╠═c4221824-889e-429b-b02f-23b7656bd3ab
# ╠═66407b07-7c99-47d1-8021-06810bae6d39
# ╟─2db8c8cc-d286-4edd-ad6b-1fc5541ab0c4
# ╠═9fd42cf0-b264-4461-90f9-62abf82d29fb
# ╟─fcf84719-d2a4-4cf9-a12f-85f7eeb20a4e
# ╠═ac3d2be8-6f4a-44d3-aa01-8a10e34c8add
# ╟─933cf514-deaf-4c3c-98b8-0ee6e3894fc7
# ╠═9fd03d96-4ba5-4686-b664-44c8bc4754b2
# ╟─33223415-a88a-4304-8af0-5ef063d14d40
# ╠═e941a84c-f386-4770-affa-df60bf5dea0b
# ╠═3b1bdc18-cf14-4811-8cba-e20f65bf3188
# ╟─8f54662a-efb1-4a2c-9305-53fcfdcd0af0
# ╟─30cca982-ef9a-4c67-ae53-82d39ff24a55
# ╠═3818e468-db0a-4bca-bbcd-ba4732e0e9df
# ╠═692b4e3d-7aec-4b73-8128-085e5b23d1f5
# ╠═9f148b9a-fc4e-4d0d-baed-c37f958afc3e
# ╠═5bdc9b96-8d13-4c55-be09-139682913c9c
# ╠═9041e3b6-ccde-4155-bee2-9aab694f5ced
# ╠═3244f883-acec-4dbf-83ad-0c8dfb6326e0
# ╠═519e3ed1-8f9a-4ccd-bd6d-115a91aa30d2
# ╠═85f65257-8c39-4836-8517-2109cf55663b
# ╠═61b908cb-41ab-4bd6-9db8-dce88834e43c
# ╟─3c5e85b2-fd11-4f70-9e6c-cc26fbcce5b3
# ╠═547bd984-10e5-4063-85a1-20262f0f1edf

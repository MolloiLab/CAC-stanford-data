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
	using ImageMorphology
	using DataFrames
	using GLM
	using DICOM
	using DICOMUtils
	using Phantoms
	using CalciumScoring
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

# ╔═╡ 78dce318-e832-40cc-b250-6dac843da9c7
pixel_size = DICOMUtils.get_pixel_size(header)

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

# ╔═╡ 7b320fab-daa4-4c6b-92ca-083e54fbfa34
c_img = calcium_image[:, :, 12];

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

# ╔═╡ 8f54662a-efb1-4a2c-9305-53fcfdcd0af0
md"""
## Calibration Prep
"""

# ╔═╡ de19275a-5b35-4e56-afd5-ea0c7e13344b
m_arr = masked_array[:, :, 25];

# ╔═╡ 6582ef5a-f3a2-444c-b4a4-804f3eef0ca9
md"""
### Intensity High Density
"""

# ╔═╡ 311d8c79-6fc3-4f1f-8635-8ae45609e490
core_L_HD = Bool.(erode(erode(mask_L_HD)));

# ╔═╡ 2c3a05dc-8a2f-488c-9318-4f40637af04e
mean_L_HD = mean(m_arr[core_L_HD])

# ╔═╡ 1922f9c6-aec3-45e4-a2fd-48532bdc2d97
# md"""
# #### Visualize High Density
# """
# ╔═╡ 3818e468-db0a-4bca-bbcd-ba4732e0e9df
begin
	arr_L_HD_cal = masked_array[:, :, 25] .* mask_L_HD
	core_HD = erode(arr_L_HD_cal)
end;

# ╔═╡ 692b4e3d-7aec-4b73-8128-085e5b23d1f5
# begin
# 	arr_L_HD_cal = m_arr .* core_L_HD
# 	heatmap(arr_L_HD_cal, colormap=:grays)
# end

# ╔═╡ da54fc49-8f69-4ff8-bbb7-800183fc3966
md"""
### Intensity Medium Density
"""

# ╔═╡ 2eb05244-eec6-46b0-99ad-cbbc8dcb9e2f
core_L_MD = Bool.(erode(erode(mask_L_MD)));

# ╔═╡ 36892f7c-e435-4e4d-a13f-584b58422a9c
mean_L_MD = mean(m_arr[core_L_MD])

# ╔═╡ b3a0a168-fc34-4153-af56-313b97eec14f
# md"""
# #### Visualize Medium Density
# """

# ╔═╡ 9e14592b-cc17-4097-8f7b-15b3fa0acba7
# begin
# 	arr_L_MD_cal = m_arr .* core_L_MD
# 	heatmap(arr_L_MD_cal, colormap=:grays)
# end

# ╔═╡ 3435f781-c842-4398-8980-b8a9572a8003
md"""
### Intensity Low Density
"""

# ╔═╡ d3248239-f499-44d2-bdc2-9ce448b87a20
core_L_LD = Bool.(erode(erode(mask_L_LD)));

# ╔═╡ 8e126adc-7eb5-4fdf-9ea0-d0b8166d7bf6
mean_L_LD = mean(m_arr[core_L_LD])

# ╔═╡ 4e0a91b8-27a3-4eb3-aa2f-edae2306d84d
# md"""
# #### Visualize Low Density
# """

# ╔═╡ 2c0c2372-675a-4dc8-9ee6-d82fa656c639
# begin
# 	arr_L_LD_cal = m_arr .* core_L_LD
# 	heatmap(arr_L_LD_cal, colormap=:grays)
# end

# ╔═╡ 4964fa87-188b-4eda-9876-7047e92b2040
md"""
### Intensity No Calcium
"""

# ╔═╡ 73fe3538-4bff-4acf-8219-8c130ef0ab9b
ring = Bool.((dilate(dilate(dilate((dilate(mask_L_LD)))))) - dilate(dilate(dilate((mask_L_LD)))));

# ╔═╡ f193ff12-f020-4348-9f80-5708c1289c1a
mean_no_CA = mean(m_arr[ring])

# ╔═╡ aafe96bc-6ecc-44dc-9842-926bf93a6536
# md"""
# #### Visualize No Calcium
# """

# ╔═╡ f84b8041-f82a-49af-8362-6771f70e2027
# begin
# 	arr_no_CA_cal = m_arr .* ring
# 	heatmap(arr_no_CA_cal, colormap=:grays)
# end

# ╔═╡ 701d01b4-26d6-4f65-8f0a-f6b8dfb24980
md"""
### Calibration Line
"""

# ╔═╡ 2ebc7519-efa2-4faa-bac1-bc8d0d283929
density_array = [0, 200, 400, 800] # g/cc

# ╔═╡ 8d0beb3d-f72b-4d16-9047-166262f8c4e5
intensity_array = [mean_no_CA, mean_L_LD, mean_L_MD, mean_L_HD] # HU

# ╔═╡ 59ccd0a1-6946-4a5c-aa5b-bb51b9cbb38c
df = DataFrame(:density => density_array, :intensity => intensity_array)

# ╔═╡ 0ce341dc-3884-4d14-b39d-d02d0e7ba906
linearRegressor = lm(@formula(intensity ~ density), df)

# ╔═╡ e3972297-eb83-41ec-a9d8-c71330b10ed6
linearFit = predict(linearRegressor)

# ╔═╡ a878d863-798e-4cb2-bd44-31e93f1c6776
md"""
We can see from above that the linear regression returns a best fit line with the formula:

```math
y = 1.4046x + 49.0528
```

Which can be solved for ``x`` and then used to calculate the density (``x``) given some measured intensity (``y``)

```math
x = \frac{y - 49.0528}{1.4046}
```

"""

# ╔═╡ 02b73f48-ed97-4b10-bcaf-cd9555d5964c
density(intensity) = (intensity - 49.0528) / 1.4046

# ╔═╡ 33a0382a-dd3a-46ec-8ff2-672da09e7f26
begin
<<<<<<< HEAD
	f = Figure()
	ax1 = Axis(f[1, 1])
	
	scatter!(density_array, intensity_array)
	lines!(density_array, linearFit, color = :red)
	ax1.title = "Calibration Line (Intensity vs Density)"
	ax1.ylabel = "Intensity (HU)"
	ax1.xlabel = "Density (mg/cm^3)"
	
	f
end

# ╔═╡ d4d9ad90-551c-4197-b4d0-b7497c945404
md"""
## Score Large Inserts
"""

# ╔═╡ f874712b-76ff-4a0a-a200-48e6f2d6237e
arr = masked_array[:, :, 23:27];
=======
	arr_L_MD_cal = masked_array[:, :, 25] .* mask_L_MD
	core_MD = erode(arr_L_MD_cal)
end;
>>>>>>> 2a130bea0f54a18e9306e9847a906f524c862d5e

# ╔═╡ 728682bd-096a-4bd5-8bda-9519d0eaedc3
s_bkg = mean_no_CA

# ╔═╡ 7cf0c601-9917-4c12-b2b6-ccff15fab380
@bind w PlutoUI.Slider(1:size(arr, 3), default=2, show_value=true)

# ╔═╡ 4c720e2d-012a-488e-a20c-4fa0a86115d4
heatmap(arr[:, :, w], colormap=:grays)

# ╔═╡ 04bfc7b6-357e-48cb-9994-5f17b5003ef9
md"""
### High Density
"""

# ╔═╡ 71343f21-9d2f-4d34-87c1-be54a9b8cd94
begin
<<<<<<< HEAD
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_HD_3D[:, :, z] = mask_L_HD
	end
=======
	arr_L_LD_cal = masked_array[:, :, 25] .* mask_L_LD
	core_LD = erode(arr_L_LD_cal)
>>>>>>> 2a130bea0f54a18e9306e9847a906f524c862d5e
end;

# ╔═╡ b10c83c5-040c-4fc8-83c4-1c8d0504ab21
s_obj_L_HD = mean(arr[mask_L_HD_3D])

# ╔═╡ c2212ea1-0644-4522-9f42-5182f10a11b3
ρ_L_HD = density(s_obj_L_HD) / 1e6 # g/cm^3 => g/mm^3

<<<<<<< HEAD
# ╔═╡ e72aefc4-3b2b-41b6-9968-88b3f1eed267
alg_L_HD = Integrated(arr[mask_L_HD_3D]);

# ╔═╡ d9b291f5-6e4c-429b-9782-4ce92de25868
V_l_hd = CalciumScoring.score(s_bkg, s_obj_L_HD, pixel_size, alg_L_HD)  # mm^3

# ╔═╡ edeb0ef5-6533-42e1-b3b5-61492af055f4
l_hd = score(s_bkg, s_obj_L_HD, pixel_size, ρ_L_HD, alg_L_HD)
=======
# ╔═╡ 519e3ed1-8f9a-4ccd-bd6d-115a91aa30d2
cal_array = [mean(core_LD[core_LD.>0]), mean(core_MD[core_MD.>0]), mean(core_HD[core_HD.>0])]

# ╔═╡ f8fbfd68-7cc9-4422-baff-e9079cbce314
plot(cal_array,[200,400,800])

# ╔═╡ d38ac6fd-75de-4943-9fec-5fe2592a4513
md"""
### 	Measured HU vs Ground truth
"""

# ╔═╡ 10b5549f-a3db-443d-9f8c-85be8f72a591
begin
	x = [1 227.609
		 1 393.29
		 1 773.667]
	y = [200
		 400
		 800]
	β̂ = (x'*x)\x'*y
	f(a)= β̂[1] + β̂[2]a
	plot(0:800,f)
end
>>>>>>> 2a130bea0f54a18e9306e9847a906f524c862d5e

# ╔═╡ db94547a-6447-43b3-afad-886c88472220
md"""
### Medium Density
"""

# ╔═╡ b33aa5e8-17e4-4fee-96d3-903c6fa2f4a2
begin
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_MD_3D[:, :, z] = mask_L_MD
	end
end

# ╔═╡ 9e7c0c1a-bc4c-431a-9ac8-27a5315b8d4b
s_obj_L_MD = mean(arr[mask_L_MD_3D])

# ╔═╡ a3082016-5b3c-42da-b7a6-ff6bfad17249
ρ_L_MD = density(s_obj_L_MD) / 1e6 # g/cm^3 => g/mm^3

# ╔═╡ 8e2d5b7d-0615-45bf-9376-7e5f1ee446e9
alg_L_MD = Integrated(arr[mask_L_MD_3D]);

# ╔═╡ 314a36ad-5211-4074-aaed-668be0175fd3
l_md = score(s_bkg, s_obj_L_MD, pixel_size, ρ_L_MD, alg_L_MD)

# ╔═╡ 0ea2438e-be61-4445-8f5e-f66b672139f5
md"""
### Low Density
"""

# ╔═╡ 596dbd94-593d-4736-83b8-62914f5fe34a
begin
	mask_L_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_LD_3D[:, :, z] = mask_L_LD
	end
end

# ╔═╡ 0063b41c-024e-44d7-b3e5-086bf611a4d5
s_obj_L_LD = mean(arr[mask_L_LD_3D])

# ╔═╡ f4e82388-be38-40b0-9203-88508a142f19
ρ_L_LD = density(s_obj_L_LD) / 1e6 # mg/cm^3 => g/mm^3

# ╔═╡ fe30dc82-f071-4bcf-b4c2-cb0b3d3636cf
alg_L_LD = Integrated(arr[mask_L_LD_3D]);

# ╔═╡ f9593c7b-0bfb-4b25-a1e3-235f603e7dc3
l_ld = score(s_bkg, s_obj_L_LD, pixel_size, ρ_L_LD, alg_L_LD)

# ╔═╡ 832250cb-0c57-4a89-be14-129f2ef026e0
md"""
## Score Medium Inserts
"""

# ╔═╡ beb32b0f-c25d-455e-92bd-018e5faa361d
md"""
### High Density
"""

# ╔═╡ ff6668c5-31fc-4764-ba0e-ffb237195074
begin
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_HD_3D[:, :, z] = mask_M_HD
	end
end

# ╔═╡ b3ed970f-0709-465f-9dd6-a8b0cd4c1c03
# arr_m_hd_cal = masked_array .* mask_M_HD;

# ╔═╡ f764cda0-5864-4291-a9c1-bba70231a1d3
# @bind z PlutoUI.Slider(1:size(arr_m_hd_cal, 3), default=10, show_value=true)

# ╔═╡ e4efb753-3d52-483e-bdfe-c8646a5102c6
# heatmap(masked_array[:, :, z], colormap=:grays)

# ╔═╡ 2ccd9e69-4e8d-4480-bd59-96d7fb126370
# heatmap(arr_m_hd_cal[:, :, z], colormap=:grays)

# ╔═╡ 76f2b9c3-5774-4022-b742-7472c3332a5b
s_obj_M_HD = mean(arr[mask_M_HD_3D])

# ╔═╡ 0b79299c-4a6d-40ce-936b-6fd4d7f26992
ρ_M_HD = density(s_obj_M_HD) / 100 # g/cm^3 => g/mm^3

# ╔═╡ d7d5ba32-f240-456c-b511-8e494ee011f1
alg_M_HD = Integrated(arr[mask_M_HD_3D]);

# ╔═╡ c0a2ff94-4bcd-4e4c-8c34-3b9ca2e420f7
m_hd = score(s_bkg, s_obj_M_HD, pixel_size, ρ_M_HD, alg_M_HD)

# ╔═╡ 356f5a4e-0c1d-42c5-9e05-184ccceb7f2d


# ╔═╡ f9f9dda8-a5fb-4535-8092-3dbde2be54e7


# ╔═╡ 7c38f7c5-a554-4697-ad5b-548a86decc9d


# ╔═╡ bca93ffa-1aa6-4c2d-b549-9d87b9b1a587


# ╔═╡ 039b2635-620d-415b-8460-283e9f7e5ec3


# ╔═╡ 8ac60d4a-a730-4766-be04-dd818f237bd3


# ╔═╡ d81c30ed-9ef7-4a6c-bb7f-012f21de942f


# ╔═╡ 310a679b-a374-4c3a-be61-40d8b555ee99


# ╔═╡ b6255f23-cc94-4092-b50d-16691706634b


# ╔═╡ 285c0c5d-ed0c-45c8-b74a-55b1aff4fe98


# ╔═╡ 9bef33a3-db70-4097-a859-3e9a96743147


# ╔═╡ 964fb40c-c6d0-4925-baad-426f7a693ea8


# ╔═╡ Cell order:
# ╠═d21e17f4-701d-11ec-3a26-b76aa97ddfe7
# ╠═0652892a-36e0-47ed-ad4f-01ee2542c244
# ╠═01a90e28-55e0-470b-a2f7-ba521c43b074
# ╠═5cec884d-a590-4f39-8f52-074be0c133b3
# ╠═b3da9520-40f2-43eb-9711-81189fc295bd
# ╠═92147219-1fc8-4759-8002-9088a4c50461
# ╠═78dce318-e832-40cc-b250-6dac843da9c7
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
# ╠═7b320fab-daa4-4c6b-92ca-083e54fbfa34
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
# ╟─8f54662a-efb1-4a2c-9305-53fcfdcd0af0
# ╠═de19275a-5b35-4e56-afd5-ea0c7e13344b
# ╟─6582ef5a-f3a2-444c-b4a4-804f3eef0ca9
# ╠═311d8c79-6fc3-4f1f-8635-8ae45609e490
# ╠═2c3a05dc-8a2f-488c-9318-4f40637af04e
# ╠═1922f9c6-aec3-45e4-a2fd-48532bdc2d97
# ╠═692b4e3d-7aec-4b73-8128-085e5b23d1f5
# ╟─da54fc49-8f69-4ff8-bbb7-800183fc3966
# ╠═2eb05244-eec6-46b0-99ad-cbbc8dcb9e2f
# ╠═36892f7c-e435-4e4d-a13f-584b58422a9c
# ╠═b3a0a168-fc34-4153-af56-313b97eec14f
# ╠═9e14592b-cc17-4097-8f7b-15b3fa0acba7
# ╟─3435f781-c842-4398-8980-b8a9572a8003
# ╠═d3248239-f499-44d2-bdc2-9ce448b87a20
# ╠═8e126adc-7eb5-4fdf-9ea0-d0b8166d7bf6
# ╠═4e0a91b8-27a3-4eb3-aa2f-edae2306d84d
# ╠═2c0c2372-675a-4dc8-9ee6-d82fa656c639
# ╟─4964fa87-188b-4eda-9876-7047e92b2040
# ╠═73fe3538-4bff-4acf-8219-8c130ef0ab9b
# ╠═f193ff12-f020-4348-9f80-5708c1289c1a
# ╠═aafe96bc-6ecc-44dc-9842-926bf93a6536
# ╠═f84b8041-f82a-49af-8362-6771f70e2027
# ╟─701d01b4-26d6-4f65-8f0a-f6b8dfb24980
# ╠═2ebc7519-efa2-4faa-bac1-bc8d0d283929
# ╠═8d0beb3d-f72b-4d16-9047-166262f8c4e5
# ╠═59ccd0a1-6946-4a5c-aa5b-bb51b9cbb38c
# ╠═0ce341dc-3884-4d14-b39d-d02d0e7ba906
# ╠═e3972297-eb83-41ec-a9d8-c71330b10ed6
# ╟─a878d863-798e-4cb2-bd44-31e93f1c6776
# ╠═02b73f48-ed97-4b10-bcaf-cd9555d5964c
# ╟─33a0382a-dd3a-46ec-8ff2-672da09e7f26
# ╟─d4d9ad90-551c-4197-b4d0-b7497c945404
# ╠═f874712b-76ff-4a0a-a200-48e6f2d6237e
# ╠═728682bd-096a-4bd5-8bda-9519d0eaedc3
# ╟─7cf0c601-9917-4c12-b2b6-ccff15fab380
# ╠═4c720e2d-012a-488e-a20c-4fa0a86115d4
# ╟─04bfc7b6-357e-48cb-9994-5f17b5003ef9
# ╠═71343f21-9d2f-4d34-87c1-be54a9b8cd94
# ╠═b10c83c5-040c-4fc8-83c4-1c8d0504ab21
# ╠═c2212ea1-0644-4522-9f42-5182f10a11b3
# ╠═e72aefc4-3b2b-41b6-9968-88b3f1eed267
# ╠═d9b291f5-6e4c-429b-9782-4ce92de25868
# ╠═edeb0ef5-6533-42e1-b3b5-61492af055f4
# ╟─db94547a-6447-43b3-afad-886c88472220
# ╠═b33aa5e8-17e4-4fee-96d3-903c6fa2f4a2
# ╠═9e7c0c1a-bc4c-431a-9ac8-27a5315b8d4b
# ╠═a3082016-5b3c-42da-b7a6-ff6bfad17249
# ╠═8e2d5b7d-0615-45bf-9376-7e5f1ee446e9
# ╠═314a36ad-5211-4074-aaed-668be0175fd3
# ╟─0ea2438e-be61-4445-8f5e-f66b672139f5
# ╠═596dbd94-593d-4736-83b8-62914f5fe34a
# ╠═0063b41c-024e-44d7-b3e5-086bf611a4d5
# ╠═f4e82388-be38-40b0-9203-88508a142f19
# ╠═fe30dc82-f071-4bcf-b4c2-cb0b3d3636cf
# ╠═f9593c7b-0bfb-4b25-a1e3-235f603e7dc3
# ╟─832250cb-0c57-4a89-be14-129f2ef026e0
# ╟─beb32b0f-c25d-455e-92bd-018e5faa361d
# ╠═ff6668c5-31fc-4764-ba0e-ffb237195074
# ╠═b3ed970f-0709-465f-9dd6-a8b0cd4c1c03
# ╠═f764cda0-5864-4291-a9c1-bba70231a1d3
# ╠═e4efb753-3d52-483e-bdfe-c8646a5102c6
# ╠═2ccd9e69-4e8d-4480-bd59-96d7fb126370
# ╠═76f2b9c3-5774-4022-b742-7472c3332a5b
# ╠═0b79299c-4a6d-40ce-936b-6fd4d7f26992
# ╠═d7d5ba32-f240-456c-b511-8e494ee011f1
# ╠═c0a2ff94-4bcd-4e4c-8c34-3b9ca2e420f7
# ╠═356f5a4e-0c1d-42c5-9e05-184ccceb7f2d
# ╠═f9f9dda8-a5fb-4535-8092-3dbde2be54e7
# ╠═7c38f7c5-a554-4697-ad5b-548a86decc9d
# ╠═bca93ffa-1aa6-4c2d-b549-9d87b9b1a587
# ╠═039b2635-620d-415b-8460-283e9f7e5ec3
# ╠═8ac60d4a-a730-4766-be04-dd818f237bd3
# ╠═d81c30ed-9ef7-4a6c-bb7f-012f21de942f
# ╠═310a679b-a374-4c3a-be61-40d8b555ee99
# ╠═b6255f23-cc94-4092-b50d-16691706634b
# ╠═285c0c5d-ed0c-45c8-b74a-55b1aff4fe98
# ╠═9bef33a3-db70-4097-a859-3e9a96743147
# ╠═964fb40c-c6d0-4925-baad-426f7a693ea8

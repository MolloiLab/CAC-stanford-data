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
@bind b1 PlutoUI.Slider(1:size(calcium_image, 3), default=10, show_value=true)

# ╔═╡ 8b25a31f-6173-4fe2-9f69-9aa0fa861d7c
heatmap(calcium_image[:, :, b1], colormap=:grays)

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
md"""
#### Visualize High Density
"""

# ╔═╡ 692b4e3d-7aec-4b73-8128-085e5b23d1f5
begin
	arr_L_HD_cal = m_arr .* core_L_HD
	heatmap(arr_L_HD_cal, colormap=:grays)
end

# ╔═╡ da54fc49-8f69-4ff8-bbb7-800183fc3966
md"""
### Intensity Medium Density
"""

# ╔═╡ 2eb05244-eec6-46b0-99ad-cbbc8dcb9e2f
core_L_MD = Bool.(erode(erode(mask_L_MD)));

# ╔═╡ 36892f7c-e435-4e4d-a13f-584b58422a9c
mean_L_MD = mean(m_arr[core_L_MD])

# ╔═╡ b3a0a168-fc34-4153-af56-313b97eec14f
md"""
#### Visualize Medium Density
"""

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
md"""
#### Visualize Low Density
"""

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

# ╔═╡ fabe0fe4-25a6-4f97-b30a-76e9df7820eb
no_ca_arr = masked_array[:, :, 23];

# ╔═╡ cd410554-5b3a-4bef-ab3b-d6875c72bd76
mean_no_CA = mean(no_ca_arr[mask])

# ╔═╡ aafe96bc-6ecc-44dc-9842-926bf93a6536
md"""
#### Visualize No Calcium
"""

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

# ╔═╡ d1b634cc-9f82-4162-bf1f-e971ec263b3f
m = linearRegressor.model.pp.beta0[2]

# ╔═╡ a8094e09-8d2a-4b45-a96e-a48703928a02
b = linearRegressor.model.rr.mu[1]

# ╔═╡ a878d863-798e-4cb2-bd44-31e93f1c6776
md"""
We can see from above that the linear regression returns a best fit line with the formula:

```math
y = mx + b
```

Which can be solved for ``x`` and then used to calculate the density (``x``) given some measured intensity (``y``)

```math
x = \frac{y - b}{m}
```

"""

# ╔═╡ 02b73f48-ed97-4b10-bcaf-cd9555d5964c
density(intensity) = (intensity - b) / m

# ╔═╡ 33a0382a-dd3a-46ec-8ff2-672da09e7f26
begin
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

# ╔═╡ 0ed5b044-0ef7-4a18-a73f-606615d76a7c
arr = masked_array[:, :, 23:27];

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
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_HD_3D[:, :, z] = mask_L_HD
	end
end;

# ╔═╡ b10c83c5-040c-4fc8-83c4-1c8d0504ab21
s_obj_L_HD = mean(arr[mask_L_HD_3D])

# ╔═╡ c2212ea1-0644-4522-9f42-5182f10a11b3
ρ_L_HD = density(s_obj_L_HD) / 1e6 # g/cm^3 => g/mm^3

# ╔═╡ e72aefc4-3b2b-41b6-9968-88b3f1eed267
alg_L_HD = Integrated(arr[mask_L_HD_3D]);

# ╔═╡ edeb0ef5-6533-42e1-b3b5-61492af055f4
l_hd = score(s_bkg, s_obj_L_HD, pixel_size, ρ_L_HD, alg_L_HD)

# ╔═╡ d0b819a6-59ff-4de0-b6c4-2c15d17cb0f1
l_hd_mg = l_hd * 10^3

# ╔═╡ db94547a-6447-43b3-afad-886c88472220
md"""
### Medium Density
"""

# ╔═╡ b33aa5e8-17e4-4fee-96d3-903c6fa2f4a2
begin
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_MD_3D[:, :, z] = (mask_L_MD)
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

# ╔═╡ 1eba82d2-88b5-45e5-94f9-6f7ae4774ac3
l_md_mg = l_md * 10^3

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

# ╔═╡ 1e1a5fb0-8aaf-4a45-be6c-f697cb5980a4
l_ld_mg = l_ld * 10^3

# ╔═╡ 8261e708-c1b0-4a7a-b772-cf2ef7e4e125
md"""
## Score Medium Inserts
"""

# ╔═╡ 707f978a-f8be-4517-b6f7-7287f8e2d75d
md"""
### High Density
"""

# ╔═╡ d92b01c0-311a-4ed8-945a-57a72b238cc9
begin
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_HD_3D[:, :, z] = mask_M_HD
	end
end;

# ╔═╡ 5e28533f-6b31-4c12-9565-1d9241f58170
s_obj_M_HD = mean(arr[mask_M_HD_3D])

# ╔═╡ 9fca1993-b3c8-4be7-aba1-4fb06e0e1381
ρ_M_HD = density(s_obj_M_HD) / 1e6 # g/cm^3 => g/mm^3

# ╔═╡ a31661cd-fbc1-4c85-bf4f-1dd7b4ed14ec
alg_M_HD = Integrated(arr[mask_M_HD_3D]);

# ╔═╡ 03359f0d-dd88-47db-a8df-7e0a3ce23798
m_hd = score(s_bkg, s_obj_M_HD, pixel_size, ρ_M_HD, alg_M_HD)

# ╔═╡ 93462f29-d951-4256-a640-756c2f8ceca8
m_hd_mg = m_hd * 10^3

# ╔═╡ cda73fb2-07ae-4937-b754-198b2f911896
md"""
### Medium Density
"""

# ╔═╡ a5e3376a-a8b7-4a2e-88ee-304e0af272da
begin
	mask_M_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_MD_3D[:, :, z] = mask_M_MD
	end
end;

# ╔═╡ 44eaf5eb-95e4-40fd-8e27-83802d6e08b9
s_obj_M_MD = mean(arr[mask_M_MD_3D])

# ╔═╡ ccfe54fe-868f-41a4-98b1-8ca4bef07daa
ρ_M_MD = density(s_obj_M_MD) / 1e6 # g/cm^3 => g/mm^3

# ╔═╡ eb908148-2c3d-415b-bfe5-baf0d2471e4c
alg_M_MD = Integrated(arr[mask_M_MD_3D]);

# ╔═╡ 49502a40-2e85-4687-9d87-9f9367375de1
m_md = score(s_bkg, s_obj_M_MD, pixel_size, ρ_M_MD, alg_M_MD)

# ╔═╡ 7e4036f1-214f-4e59-a5db-abbbf3d1aa23
m_md_mg = m_md * 10^3

# ╔═╡ 6a8fd67a-62ca-41aa-87d8-6a96eab16a31
md"""
### Low Density
"""

# ╔═╡ 833f2e9b-458a-42d6-897a-4ce93a516ba5
begin
	mask_M_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_LD_3D[:, :, z] = mask_M_LD
	end
end;

# ╔═╡ d980e52c-0258-4d95-a774-b337e46c7944
s_obj_M_LD = mean(arr[mask_M_LD_3D])

# ╔═╡ 33c9fa09-4bf4-4b92-ac0a-2d772ee23baf
ρ_M_LD = density(s_obj_M_LD) / 1e6 # g/cm^3 => g/mm^3

# ╔═╡ e3ace679-5cf9-4ca6-8835-5c2b018121aa
alg_M_LD = Integrated(arr[mask_M_LD_3D]);

# ╔═╡ 804422b8-8a2e-4ec6-a77a-12a7ac627a3e
m_ld = score(s_bkg, s_obj_M_LD, pixel_size, ρ_M_LD, alg_M_LD)

# ╔═╡ 874424a8-a1f2-4710-91d3-6da3d22578f8
m_ld_mg = m_ld * 10^3

# ╔═╡ 2ca0a9e2-cab4-428d-820d-ddbf05cede43
md"""
## Score Small Inserts
"""

# ╔═╡ 6d967ff1-c870-4cba-8c2d-98d2c2347643
md"""
### High Density
"""

# ╔═╡ b941068e-4348-49ce-a562-f75ce93b0f6c
begin
	mask_S_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_HD_3D[:, :, z] = mask_S_HD
	end
end;

# ╔═╡ 8ee64580-75da-4e7f-9d6b-5a4f6c86dba3
s_obj_S_HD = mean(arr[mask_S_HD_3D])

# ╔═╡ fdcf0efc-b347-47e7-aba4-82e6faba1273
ρ_S_HD = density(s_obj_S_HD) / 1e6 # g/cm^3 => g/mm^3

# ╔═╡ cfd3d799-d720-4586-b6eb-5c26cb61346f
alg_S_HD = Integrated(arr[mask_S_HD_3D]);

# ╔═╡ 64d39cb7-1934-4b79-ae64-afb88935f6f1
s_hd = score(s_bkg, s_obj_S_HD, pixel_size, ρ_S_HD, alg_S_HD)

# ╔═╡ 33f13e1c-dcd4-4144-a4ec-99105d1dd4b8
s_hd_mg = s_hd * 10^3

# ╔═╡ 9d77bfbb-4681-408c-b9cb-7836a2cbf1f9
md"""
### Medium Density
"""

# ╔═╡ 75f87fbf-554e-48b4-91a5-cc2c6299b3c9
begin
	mask_S_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_MD_3D[:, :, z] = mask_S_MD
	end
end;

# ╔═╡ d52e5878-84da-4421-bf18-d6ddfc2aabae
s_obj_S_MD = mean(arr[mask_S_MD_3D])

# ╔═╡ f089268a-ef08-42b5-afe6-ba2fa8d6ff63
ρ_S_MD = density(s_obj_S_MD) / 1e6 # g/cm^3 => g/mm^3

# ╔═╡ 47a8effe-fcf4-4917-8d41-07f88783672b
alg_S_MD = Integrated(arr[mask_S_MD_3D]);

# ╔═╡ 27bf9e2a-db01-421c-8eaf-8c4cf9b4065a
s_md = score(s_bkg, s_obj_S_MD, pixel_size, ρ_S_MD, alg_S_MD)

# ╔═╡ 38e81cac-c0b6-47be-9790-4377540f4b23
s_md_mg = s_md * 10^3

# ╔═╡ 83805a01-e3ee-4bef-a5ad-b90bfe6f3a33
md"""
### Low Density
"""

# ╔═╡ 032e9aa4-4312-47d2-be64-4fd61fd87301
begin
	mask_S_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_LD_3D[:, :, z] = mask_S_LD
	end
end;

# ╔═╡ 29400dab-a97a-4a9f-83a5-e166a1f47b44
s_obj_S_LD = mean(arr[mask_S_LD_3D])

# ╔═╡ 66bde869-c238-42d5-9073-39a5371fdb67
ρ_S_LD = density(s_obj_S_LD) / 1e6 # g/cm^3 => g/mm^3

# ╔═╡ 6f797cf4-485c-4199-afe3-0a5b27281dfd
alg_S_LD = Integrated(arr[mask_S_LD_3D]);

# ╔═╡ 47a150a6-0256-4a06-8811-1fb9bfb70e75
s_ld = score(s_bkg, s_obj_S_LD, pixel_size, ρ_S_LD, alg_S_LD)

# ╔═╡ 891980c6-4db6-4bf9-bece-29ab8dbcb548
s_ld_mg = s_ld * 10^3

# ╔═╡ b877b8f5-da63-47c7-aa3f-e989abfd4c76
md"""
## Visualize Results
"""

# ╔═╡ 68b9757d-3345-43e8-8a03-3694c74914c0
md"""
### High Density
"""

# ╔═╡ b6b38c11-17a0-4b78-b1df-0db54f36d5fe
ground_truth_lg = [19.6, 39.3, 78.5] # mg

# ╔═╡ a00e8fc5-574d-467f-8c5e-e971a2705fb6
measured_lg = [l_ld_mg, l_md_mg, l_hd_mg] # mg

# ╔═╡ 7f72f4cb-e4b9-478a-888d-8b542b588926
begin
	f22 = Figure()
	ax22 = Axis(f22[1, 1])
	
	sc22_1 = scatter!(density_array[2:end], ground_truth_lg, color=:blue)
	sc22_2 = scatter!(density_array[2:end], measured_lg, color=:red)
	ax22.title = "Ground Truth vs Measured Mass"
	ax22.ylabel = "Mass (mg)"
	ax22.xlabel = "Density (mg/cm^3)"
	
	Legend(f22[1, 2],
    [sc22_1, sc22_2],
    ["ground truth", "measured"])
	
	f22
end

# ╔═╡ 65c9a7a5-6be5-43f9-b0c7-6ad57925ea6f
md"""
### Medium Density
"""

# ╔═╡ 9823a3c1-887e-47e4-b8d4-ee755aa2ddf5
ground_truth_md = [4.2, 8.5, 17.0] # mg

# ╔═╡ 19311e20-42ed-4fc5-9e8e-bd579ef12531
measured_md = [m_ld_mg, m_md_mg, m_hd_mg] # mg

# ╔═╡ 2c0d56bb-2c5c-46c2-83f4-61733512e87c
begin
	f23 = Figure()
	ax23 = Axis(f23[1, 1])
	
	sc23_1 = scatter!(density_array[2:end], ground_truth_md, color=:blue)
	sc23_2 = scatter!(density_array[2:end], measured_md, color=:orange)
	ax23.title = "Ground Truth vs Measured Mass"
	ax23.ylabel = "Mass (mg)"
	ax23.xlabel = "Density (mg/cm^3)"

	Legend(f23[1, 2],
    [sc23_1, sc23_2],
    ["ground truth", "measured"])
	
	f23
end

# ╔═╡ 5cc86912-a3e4-44b7-880b-f93c26cf8374
md"""
### Low Density
"""

# ╔═╡ 09ef6539-c621-48de-9b1a-5342b5772ce6
ground_truth_ld = [0.2, 0.3, 0.6] # mg

# ╔═╡ fb0d7d07-fedc-443b-846c-cb004900d2be
measured_ld = [s_ld_mg, s_md_mg, s_hd_mg] # mg

# ╔═╡ ccb72aac-3f98-4c00-b176-f7e157bb5bd7
begin
	f24 = Figure()
	ax24 = Axis(f24[1, 1])
	
	sc24_1 = scatter!(density_array[2:end], ground_truth_ld, color=:blue)
	sc24_2 = scatter!(density_array[2:end], measured_ld, color=:orange)
	ax24.title = "Ground Truth vs Measured Mass"
	ax24.ylabel = "Mass (mg)"
	ax24.xlabel = "Density (mg/cm^3)"

	Legend(f24[1, 2],
    [sc24_1, sc24_2],
    ["ground truth", "measured"])
	
	f24
end

# ╔═╡ 9dc08bae-57e6-4477-bff7-008ae6e996a0
md"""
## Volume Measurements
"""

# ╔═╡ c061bf75-94f2-45f6-b5cc-3b3abe66547a
md"""
### Large
"""

# ╔═╡ 0c2bea33-62c1-4b64-b55f-44a887e32153
vol_l_hd = score(s_bkg, s_obj_L_HD, pixel_size, alg_L_HD)

# ╔═╡ 641bc2b3-b599-42a0-8d4d-376e40cacddc
vol_l_md = score(s_bkg, s_obj_L_MD, pixel_size, alg_L_MD)

# ╔═╡ a2c146a3-f470-4e77-8691-b694ac298615
vol_l_ld = score(s_bkg, s_obj_L_LD, pixel_size, alg_L_LD)

# ╔═╡ a829ec5b-7647-434e-b3d9-0037891cecbc
vol_groundtruth_lrg = [98.2, 98.2, 98.2]

# ╔═╡ e5fc8875-10ce-40a2-afbb-d2e5a80c74a0
vol_measured_lrg = [vol_l_ld, vol_l_md, vol_l_hd]

# ╔═╡ 87b56ab0-e196-4c08-bb56-c8fff5556ded
md"""
### Medium
"""

# ╔═╡ 9698c78c-d1a8-4c5c-ae28-3df1207c84de
vol_m_hd = score(s_bkg, s_obj_M_HD, pixel_size, alg_M_HD)

# ╔═╡ 827caa30-ec25-4f6a-adeb-5c739a2af3c3
vol_m_md = score(s_bkg, s_obj_M_MD, pixel_size, alg_M_MD)

# ╔═╡ f2119d87-8a8d-4e10-80aa-e27052ae5f8a
vol_m_ld = score(s_bkg, s_obj_M_LD, pixel_size, alg_M_LD)

# ╔═╡ 6361f784-0c36-4b8a-9607-078011051758
vol_groundtruth_med = [21.2, 21.2, 21.2]

# ╔═╡ 55accaa0-56c6-4334-b481-604347bd6c63
vol_measured_med = [vol_m_ld, vol_m_md, vol_m_hd]

# ╔═╡ 4497f442-13f1-49fa-96e1-36028e013bfc
md"""
### Small
"""

# ╔═╡ e39753c3-962b-436c-976f-13cf3b438f9f
vol_s_hd = score(s_bkg, s_obj_S_HD, pixel_size, alg_S_HD)

# ╔═╡ 866748cf-dbaa-42ec-8af2-e9c69a17acf4
vol_s_md = score(s_bkg, s_obj_S_MD, pixel_size, alg_S_MD)

# ╔═╡ 6b038fa8-a115-4244-aceb-d88f5d8eba4d
vol_s_ld = score(s_bkg, s_obj_S_LD, pixel_size, alg_S_LD)

# ╔═╡ fbddfbea-d94b-4b4c-ae24-775b337fc92a
vol_groundtruth_small = [21.2, 21.2, 21.2]

# ╔═╡ 8423c6af-b54e-429a-a5c6-07d047c26028
vol_measured_small = [vol_s_ld, vol_s_md, vol_s_hd]

# ╔═╡ db8331c0-ccc4-45a2-9dec-f1a79b8e56be
md"""
### Visualize
"""

# ╔═╡ 76f28ce7-15a8-4feb-a502-6077a6b4b17c
begin
	f25 = Figure()
	ax25 = Axis(f25[1, 1])
	
	sc25_1 = scatter!(density_array[2:end], vol_groundtruth_lrg, color=:blue)
	sc25_2 = scatter!(density_array[2:end], vol_measured_lrg, color=:purple)
	ax25.title = "Ground Truth vs Measured Volume (Large Inserts)"
	ax25.ylabel = "Mass (mg)"
	ax25.xlabel = "Density (mg/cm^3)"

	Legend(f25[1, 2],
    [sc25_1, sc25_2],
    ["ground truth", "measured"])
	
	f25
end

# ╔═╡ 332511a5-b41f-4ed3-8fa7-5e086a1556e5
begin
	f26 = Figure()
	ax26 = Axis(f26[1, 1])
	
	sc26_1 = scatter!(density_array[2:end], vol_groundtruth_med, color=:pink)
	sc26_2 = scatter!(density_array[2:end], vol_measured_med, color=:green)
	ax26.title = "Ground Truth vs Measured Volume (Medium Inserts)"
	ax26.ylabel = "Mass (mg)"
	ax26.xlabel = "Density (mg/cm^3)"

	Legend(f26[1, 2],
    [sc26_1, sc26_2],
    ["ground truth", "measured"])
	
	f26
end

# ╔═╡ 3c93ce01-2695-48fd-abad-fd13d77c097a
begin
	f27 = Figure()
	ax27 = Axis(f27[1, 1])
	
	sc27_1 = scatter!(density_array[2:end], vol_groundtruth_small, color=:orange)
	sc27_2 = scatter!(density_array[2:end], vol_measured_small, color=:red)
	ax27.title = "Ground Truth vs Measured Volume (Small Inserts)"
	ax27.ylabel = "Mass (mg)"
	ax27.xlabel = "Density (mg/cm^3)"

	Legend(f27[1, 2],
    [sc27_1, sc27_2],
    ["ground truth", "measured"])
	
	f27
end

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
# ╟─1922f9c6-aec3-45e4-a2fd-48532bdc2d97
# ╠═692b4e3d-7aec-4b73-8128-085e5b23d1f5
# ╟─da54fc49-8f69-4ff8-bbb7-800183fc3966
# ╠═2eb05244-eec6-46b0-99ad-cbbc8dcb9e2f
# ╠═36892f7c-e435-4e4d-a13f-584b58422a9c
# ╟─b3a0a168-fc34-4153-af56-313b97eec14f
# ╠═9e14592b-cc17-4097-8f7b-15b3fa0acba7
# ╟─3435f781-c842-4398-8980-b8a9572a8003
# ╠═d3248239-f499-44d2-bdc2-9ce448b87a20
# ╠═8e126adc-7eb5-4fdf-9ea0-d0b8166d7bf6
# ╟─4e0a91b8-27a3-4eb3-aa2f-edae2306d84d
# ╠═2c0c2372-675a-4dc8-9ee6-d82fa656c639
# ╟─4964fa87-188b-4eda-9876-7047e92b2040
# ╠═73fe3538-4bff-4acf-8219-8c130ef0ab9b
# ╠═fabe0fe4-25a6-4f97-b30a-76e9df7820eb
# ╠═cd410554-5b3a-4bef-ab3b-d6875c72bd76
# ╟─aafe96bc-6ecc-44dc-9842-926bf93a6536
# ╟─701d01b4-26d6-4f65-8f0a-f6b8dfb24980
# ╠═2ebc7519-efa2-4faa-bac1-bc8d0d283929
# ╠═8d0beb3d-f72b-4d16-9047-166262f8c4e5
# ╠═59ccd0a1-6946-4a5c-aa5b-bb51b9cbb38c
# ╠═0ce341dc-3884-4d14-b39d-d02d0e7ba906
# ╠═e3972297-eb83-41ec-a9d8-c71330b10ed6
# ╠═d1b634cc-9f82-4162-bf1f-e971ec263b3f
# ╠═a8094e09-8d2a-4b45-a96e-a48703928a02
# ╟─a878d863-798e-4cb2-bd44-31e93f1c6776
# ╠═02b73f48-ed97-4b10-bcaf-cd9555d5964c
# ╟─33a0382a-dd3a-46ec-8ff2-672da09e7f26
# ╟─d4d9ad90-551c-4197-b4d0-b7497c945404
# ╠═0ed5b044-0ef7-4a18-a73f-606615d76a7c
# ╠═728682bd-096a-4bd5-8bda-9519d0eaedc3
# ╟─7cf0c601-9917-4c12-b2b6-ccff15fab380
# ╠═4c720e2d-012a-488e-a20c-4fa0a86115d4
# ╟─04bfc7b6-357e-48cb-9994-5f17b5003ef9
# ╠═71343f21-9d2f-4d34-87c1-be54a9b8cd94
# ╠═b10c83c5-040c-4fc8-83c4-1c8d0504ab21
# ╠═c2212ea1-0644-4522-9f42-5182f10a11b3
# ╠═e72aefc4-3b2b-41b6-9968-88b3f1eed267
# ╠═edeb0ef5-6533-42e1-b3b5-61492af055f4
# ╠═d0b819a6-59ff-4de0-b6c4-2c15d17cb0f1
# ╟─db94547a-6447-43b3-afad-886c88472220
# ╠═b33aa5e8-17e4-4fee-96d3-903c6fa2f4a2
# ╠═9e7c0c1a-bc4c-431a-9ac8-27a5315b8d4b
# ╠═a3082016-5b3c-42da-b7a6-ff6bfad17249
# ╠═8e2d5b7d-0615-45bf-9376-7e5f1ee446e9
# ╠═314a36ad-5211-4074-aaed-668be0175fd3
# ╠═1eba82d2-88b5-45e5-94f9-6f7ae4774ac3
# ╟─0ea2438e-be61-4445-8f5e-f66b672139f5
# ╠═596dbd94-593d-4736-83b8-62914f5fe34a
# ╠═0063b41c-024e-44d7-b3e5-086bf611a4d5
# ╠═f4e82388-be38-40b0-9203-88508a142f19
# ╠═fe30dc82-f071-4bcf-b4c2-cb0b3d3636cf
# ╠═f9593c7b-0bfb-4b25-a1e3-235f603e7dc3
# ╠═1e1a5fb0-8aaf-4a45-be6c-f697cb5980a4
# ╟─8261e708-c1b0-4a7a-b772-cf2ef7e4e125
# ╟─707f978a-f8be-4517-b6f7-7287f8e2d75d
# ╠═d92b01c0-311a-4ed8-945a-57a72b238cc9
# ╠═5e28533f-6b31-4c12-9565-1d9241f58170
# ╠═9fca1993-b3c8-4be7-aba1-4fb06e0e1381
# ╠═a31661cd-fbc1-4c85-bf4f-1dd7b4ed14ec
# ╠═03359f0d-dd88-47db-a8df-7e0a3ce23798
# ╠═93462f29-d951-4256-a640-756c2f8ceca8
# ╟─cda73fb2-07ae-4937-b754-198b2f911896
# ╠═a5e3376a-a8b7-4a2e-88ee-304e0af272da
# ╠═44eaf5eb-95e4-40fd-8e27-83802d6e08b9
# ╠═ccfe54fe-868f-41a4-98b1-8ca4bef07daa
# ╠═eb908148-2c3d-415b-bfe5-baf0d2471e4c
# ╠═49502a40-2e85-4687-9d87-9f9367375de1
# ╠═7e4036f1-214f-4e59-a5db-abbbf3d1aa23
# ╟─6a8fd67a-62ca-41aa-87d8-6a96eab16a31
# ╠═833f2e9b-458a-42d6-897a-4ce93a516ba5
# ╠═d980e52c-0258-4d95-a774-b337e46c7944
# ╠═33c9fa09-4bf4-4b92-ac0a-2d772ee23baf
# ╠═e3ace679-5cf9-4ca6-8835-5c2b018121aa
# ╠═804422b8-8a2e-4ec6-a77a-12a7ac627a3e
# ╠═874424a8-a1f2-4710-91d3-6da3d22578f8
# ╟─2ca0a9e2-cab4-428d-820d-ddbf05cede43
# ╟─6d967ff1-c870-4cba-8c2d-98d2c2347643
# ╠═b941068e-4348-49ce-a562-f75ce93b0f6c
# ╠═8ee64580-75da-4e7f-9d6b-5a4f6c86dba3
# ╠═fdcf0efc-b347-47e7-aba4-82e6faba1273
# ╠═cfd3d799-d720-4586-b6eb-5c26cb61346f
# ╠═64d39cb7-1934-4b79-ae64-afb88935f6f1
# ╠═33f13e1c-dcd4-4144-a4ec-99105d1dd4b8
# ╟─9d77bfbb-4681-408c-b9cb-7836a2cbf1f9
# ╠═75f87fbf-554e-48b4-91a5-cc2c6299b3c9
# ╠═d52e5878-84da-4421-bf18-d6ddfc2aabae
# ╠═f089268a-ef08-42b5-afe6-ba2fa8d6ff63
# ╠═47a8effe-fcf4-4917-8d41-07f88783672b
# ╠═27bf9e2a-db01-421c-8eaf-8c4cf9b4065a
# ╠═38e81cac-c0b6-47be-9790-4377540f4b23
# ╟─83805a01-e3ee-4bef-a5ad-b90bfe6f3a33
# ╠═032e9aa4-4312-47d2-be64-4fd61fd87301
# ╠═29400dab-a97a-4a9f-83a5-e166a1f47b44
# ╠═66bde869-c238-42d5-9073-39a5371fdb67
# ╠═6f797cf4-485c-4199-afe3-0a5b27281dfd
# ╠═47a150a6-0256-4a06-8811-1fb9bfb70e75
# ╠═891980c6-4db6-4bf9-bece-29ab8dbcb548
# ╟─b877b8f5-da63-47c7-aa3f-e989abfd4c76
# ╟─68b9757d-3345-43e8-8a03-3694c74914c0
# ╠═b6b38c11-17a0-4b78-b1df-0db54f36d5fe
# ╠═a00e8fc5-574d-467f-8c5e-e971a2705fb6
# ╟─7f72f4cb-e4b9-478a-888d-8b542b588926
# ╟─65c9a7a5-6be5-43f9-b0c7-6ad57925ea6f
# ╠═9823a3c1-887e-47e4-b8d4-ee755aa2ddf5
# ╠═19311e20-42ed-4fc5-9e8e-bd579ef12531
# ╟─2c0d56bb-2c5c-46c2-83f4-61733512e87c
# ╟─5cc86912-a3e4-44b7-880b-f93c26cf8374
# ╠═09ef6539-c621-48de-9b1a-5342b5772ce6
# ╠═fb0d7d07-fedc-443b-846c-cb004900d2be
# ╟─ccb72aac-3f98-4c00-b176-f7e157bb5bd7
# ╟─9dc08bae-57e6-4477-bff7-008ae6e996a0
# ╟─c061bf75-94f2-45f6-b5cc-3b3abe66547a
# ╠═0c2bea33-62c1-4b64-b55f-44a887e32153
# ╠═641bc2b3-b599-42a0-8d4d-376e40cacddc
# ╠═a2c146a3-f470-4e77-8691-b694ac298615
# ╠═a829ec5b-7647-434e-b3d9-0037891cecbc
# ╠═e5fc8875-10ce-40a2-afbb-d2e5a80c74a0
# ╟─87b56ab0-e196-4c08-bb56-c8fff5556ded
# ╠═9698c78c-d1a8-4c5c-ae28-3df1207c84de
# ╠═827caa30-ec25-4f6a-adeb-5c739a2af3c3
# ╠═f2119d87-8a8d-4e10-80aa-e27052ae5f8a
# ╠═6361f784-0c36-4b8a-9607-078011051758
# ╠═55accaa0-56c6-4334-b481-604347bd6c63
# ╟─4497f442-13f1-49fa-96e1-36028e013bfc
# ╠═e39753c3-962b-436c-976f-13cf3b438f9f
# ╠═866748cf-dbaa-42ec-8af2-e9c69a17acf4
# ╠═6b038fa8-a115-4244-aceb-d88f5d8eba4d
# ╠═fbddfbea-d94b-4b4c-ae24-775b337fc92a
# ╠═8423c6af-b54e-429a-a5c6-07d047c26028
# ╟─db8331c0-ccc4-45a2-9dec-f1a79b8e56be
# ╟─76f28ce7-15a8-4feb-a502-6077a6b4b17c
# ╟─332511a5-b41f-4ed3-8fa7-5e086a1556e5
# ╟─3c93ce01-2695-48fd-abad-fd13d77c097a

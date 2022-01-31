### A Pluto.jl notebook ###
# v0.17.6

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

# ╔═╡ c3e3a7a8-9d0b-4363-943d-d723f7f4cb37
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
	using CSV
	using DataFrames
	using GLM
	using DICOM
	using DICOMUtils
	using Phantoms
	using CalciumScoring
end

# ╔═╡ ff527e22-8eb6-413e-b45e-7c07d7623f45
TableOfContents()

# ╔═╡ 903cf8bf-1113-4a0f-9e07-122ba61d34cb
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDER`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ 88af5c56-9841-48ff-95fc-cceddec0184f
begin
	SCAN_NUMBER = 10
	VENDER = "Canon_Aquilion_One_Vision"
	BASE_PATH = "/Users/daleblack/Google Drive/Datasets/"
end

# ╔═╡ 8024a262-7958-4d21-bbe0-7b2d95a562f1
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ 7addd762-8754-4ba2-a05e-493f7171847e
root_path = string(BASE_PATH, VENDER)

# ╔═╡ 992dd715-3740-4f8e-b623-6b99f61ebc90
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ cf504678-ae5b-4327-b956-b3dc1ffac6df
pth = dcm_path_list[SCAN_NUMBER]

# ╔═╡ 7d8776eb-6623-4cde-b91d-e2a91e50b377
scan = basename(pth)

# ╔═╡ 60842b6a-72d2-4501-911e-19388c9df4d9
header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ c7dc5fd3-9f08-43f3-86a7-b77a174f76ca
md"""
## Helper Functions
"""

# ╔═╡ 6fb61984-75c8-4d2f-b5ad-d99b6cb382ad
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ a241f5db-3e70-4893-8b19-3f05f1b94a5a
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ 6b403990-eef6-42ca-8cf0-e479c5eca2f8
function overlay_mask_plot(array, mask, var, title::AbstractString)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	indices_lbl = findall(x -> x == zs[var], label_array[:,3])
	
	fig = Figure()
	ax = Makie.Axis(fig[1, 1])
	ax.title = title
	heatmap!(array[:, :, zs[var]], colormap=:grays)
	scatter!(label_array[:, 1][indices_lbl], label_array[:, 2][indices_lbl], markersize=1, color=:red)
	fig
end

# ╔═╡ 4fefaf32-4cb5-4cf8-b94d-4695d13b5797
md"""
## Segment Heart
"""

# ╔═╡ b1d48c95-a6c4-4b6a-8b14-fa34976cad83
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ 8f6f0951-8f5a-483d-81bb-0928f2f765ce
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ ec6bde85-0dec-46c4-9c59-14b3126927a5
heatmap(masked_array[:, :, a], colormap=:grays)

# ╔═╡ 985c7ce3-3414-4fa4-8c11-647fd3cbf2ce
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 15]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig
end

# ╔═╡ 27ed5fe4-775b-40a9-8821-ef8b4283e40c
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig2
end

# ╔═╡ 749acc75-b49e-4a24-b1a2-c8e57bfdcb58
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 23]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig3
end

# ╔═╡ f63c9679-3c11-4a28-8f6d-7c903419930a
md"""
## Segment Calcium Rod
"""

# ╔═╡ 39c633c8-af27-4d92-804d-5b0263c40e91
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header);

# ╔═╡ 2dc7cf6e-2626-496f-84ee-629320f25bc8
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=cal_rod_slice, show_value=true)

# ╔═╡ faa85dd6-303d-478f-ad12-148635f9d02d
heatmap(calcium_image[:, :, c], colormap=:grays)

# ╔═╡ ae76e2c4-a74e-496a-b6a4-17fe91f2ad57
md"""
## Segment Calcium Inserts
"""

# ╔═╡ cbf71130-ddcf-4418-9e3e-7df8fc57c12f
mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
            dcm_array, masked_array, header, slice_CCI, center_insert
);

# ╔═╡ 212ed9e1-b5de-4f0f-bf71-c811908c5d45
masks = mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD;

# ╔═╡ b87fd07d-8cea-4a5a-94dc-6682d66e8bd8
heatmap(masks, colormap=:grays)

# ╔═╡ 5b9d4af0-3fa7-4e33-8888-4e3b7d1985ba
md"""
## Calibration Prep
"""

# ╔═╡ 5484979d-aa3a-4eab-a1eb-2b81a83068f4
m_arr = masked_array[:, :, slice_CCI:slice_CCI+1];

# ╔═╡ 8557a14a-17b2-4039-a110-f35c87ea77b6
md"""
### Intensity High Density
"""

# ╔═╡ eab14ca7-4299-42d1-b0de-f85e6cf0386c
core_L_HD = Bool.(erode(erode(mask_L_HD)));

# ╔═╡ 43a5e31f-2c82-4980-a4eb-ce0ae1a5f448
begin
	mask_HD_cal = Array{Bool}(undef, size(m_arr))
	for z in 1:size(m_arr, 3)
		mask_HD_cal[:, :, z] = core_L_HD
	end
end;

# ╔═╡ e00ef290-7ffa-4a43-81d3-735b4191dc98
hist(m_arr[mask_HD_cal])

# ╔═╡ 431ff05c-804a-4e86-ac77-c6dd5e933d35
mean_L_HD = quantile!(m_arr[mask_HD_cal], 0.70)

# ╔═╡ 781db9e6-f71b-42ff-a060-0d3d11cd8db8
md"""
### Intensity Medium Density
"""

# ╔═╡ a2262988-1f50-4622-a603-e0925ad0d065
core_L_MD = Bool.(erode(erode((mask_L_MD))));

# ╔═╡ d2033181-98a8-4a46-b435-2acf6eb72485
begin
	mask_MD_cal = Array{Bool}(undef, size(m_arr))
	for z in 1:size(m_arr, 3)
		mask_MD_cal[:, :, z] = core_L_MD
	end
end;

# ╔═╡ 1c62ece3-da87-48a6-a2c0-fa536ead7e73
hist(m_arr[mask_MD_cal])

# ╔═╡ ddeda291-edb8-4bc2-b1cf-130ea996c086
mean_L_MD = quantile!(m_arr[mask_MD_cal], 0.70)

# ╔═╡ 7695a3b7-31c6-4a31-9fcb-24901851e269
md"""
### Intensity Low Density
"""

# ╔═╡ c85b86ba-4f5b-4bfb-a26e-36ac185dc794
core_L_LD = Bool.(erode(erode(mask_L_LD)));

# ╔═╡ 3fff7d44-1586-492a-8e1e-dfcda0979106
begin
	mask_LD_cal = Array{Bool}(undef, size(m_arr))
	for z in 1:size(m_arr, 3)
		mask_LD_cal[:, :, z] = core_L_LD
	end
end;

# ╔═╡ d4cb8eba-17cf-4aa3-8ad4-e8ac35811caf
hist(m_arr[mask_LD_cal])

# ╔═╡ 74b034ae-a9e1-4f72-bc75-5860047ef8c5
mean_L_LD = quantile!(m_arr[mask_LD_cal], 0.70)

# ╔═╡ dcf800fc-7dea-4197-b226-504d76569554
md"""
### Calibration Line
"""

# ╔═╡ ffd33235-091b-4154-a3e6-e76ba1d1b1fa
density_array = [0, 200, 400, 800] # mg/cc

# ╔═╡ bcbc253e-9f66-4488-af84-10022e4d765d
intensity_array = [0, mean_L_LD, mean_L_MD, mean_L_HD] # HU

# ╔═╡ c440667a-deb5-4e21-b3ec-6e715cda538b
df = DataFrame(:density => density_array, :intensity => intensity_array)

# ╔═╡ 1378a515-b757-4776-b5aa-864ddc89e94b
linearRegressor = lm(@formula(intensity ~ density), df)

# ╔═╡ 264fe678-9b84-4a2a-9a66-ca2e6dd3b65d
linearFit = predict(linearRegressor)

# ╔═╡ 243b5511-ec59-4a8a-a29e-334ae2fa8985
m = linearRegressor.model.pp.beta0[2]

# ╔═╡ 9a8c3d70-5eb5-46b4-a8a5-239335512196
b = linearRegressor.model.rr.mu[1]

# ╔═╡ e8cc070c-1850-474c-a2fa-d0c12a09f683
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

# ╔═╡ 885ffebe-4805-43e6-a305-c7076c75ec5b
density(intensity) = (intensity - b) / m

# ╔═╡ 26a4d12b-9fcb-4e90-a96e-5ed9c8cfe2ce
intensity(ρ) = m*ρ + b

# ╔═╡ 847262c3-dd0a-4647-81fa-ea7a2a441605
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

# ╔═╡ d9600f56-c238-4392-a262-0a25ae19f9ba
md"""
# Score Large Inserts
"""

# ╔═╡ 158347a5-966f-42f7-ae1b-b41b39d79b8b
arr = masked_array[:, :, slice_CCI-2:slice_CCI+2];

# ╔═╡ 98dcf986-4e0e-45dd-a52a-7886ce868842
single_arr = masked_array[:, :, slice_CCI];

# ╔═╡ 953cd854-108f-486d-ac58-0405b6b785c4
md"""
## High Density
"""

# ╔═╡ 2024648c-45f4-403a-b9e7-11b2ec3673ae
begin
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_HD_3D[:, :, z] = mask_L_HD
	end
end;

# ╔═╡ 606c40bb-b9c9-4ec4-8505-77eba0799ce5
md"""
#### Dilated mask
"""

# ╔═╡ dbf78411-7e46-439e-aaf9-640a88354d39
dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D));

# ╔═╡ 8ea123b6-4eba-41db-abe8-c994f76c9cda
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ ba2fed15-a79e-4b21-a575-5884cce8c366
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ e17bd8ad-b2b1-4020-a7db-1ec47af05b09
md"""
#### Ring (background) mask
"""

# ╔═╡ de607584-2f2e-4e9e-9c22-77cb05596b44
ring_mask_L_HD = dilate(dilate(dilate(dilate(mask_L_HD_3D)))) - dilate(dilate(dilate(mask_L_HD_3D)));

# ╔═╡ 7d655971-44da-42ad-9420-86528540a9b2
@bind g4 overlay_mask_bind(ring_mask_L_HD)

# ╔═╡ 476309ca-a32b-49b8-8845-3fdaf4a2d847
overlay_mask_plot(arr, ring_mask_L_HD, g4, "ring mask")

# ╔═╡ 7e6f8c7f-3f30-4d41-aa78-77c8d623a499
md"""
### Calculations
"""

# ╔═╡ cf725101-4564-4546-bc1f-d073d2f1cc04
begin
	single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
	s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
end

# ╔═╡ 0b502b66-da66-44cf-8531-945e5a466ec1
S_Obj_HD = intensity(800)

# ╔═╡ 3bdc0b37-67d7-49fa-886d-f4b66fe28c82
pixel_size = DICOMUtils.get_pixel_size(header)

# ╔═╡ 8c73a297-5d95-43c9-ba2c-09847fc5fc71
begin
	alg_L_HD = Integrated(arr[mask_L_HD_3D])
	vol_l_hd = score(s_bkg_L_HD, S_Obj_HD, pixel_size, alg_L_HD)
end

# ╔═╡ 89edc6a0-2933-4dba-be63-48d3be00444e
begin
	ρ_HD = 0.8 # mg/mm^3
	mass_l_hd = score(s_bkg_L_HD, S_Obj_HD, pixel_size, ρ_HD, alg_L_HD)
end

# ╔═╡ c5a58e31-5f64-4d2e-b222-d918d6d41b5d
md"""
## Medium Density
"""

# ╔═╡ f7b97123-ac55-49d8-9bd9-45089245c22a
begin
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_MD_3D[:, :, z] = mask_L_MD
	end
end;

# ╔═╡ 45379669-4897-460c-8ba9-f1635f3bcb13
md"""
#### Dilated mask
"""

# ╔═╡ 5c98b68c-01e9-4cb1-954d-6f0192906e6f
dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D));

# ╔═╡ cd0aef7b-fdc9-4ebc-a9d8-23d8655144ef
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ 1a5150b9-bab3-4ec5-92b0-4734917ec58f
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ 73a5aee9-38a4-4581-b654-8a95cdf74108
md"""
#### Ring (background) mask
"""

# ╔═╡ d86fa8f6-8e54-44b3-8188-4c0373d5e633
ring_mask_L_MD = dilate(dilate(dilate(dilate(mask_L_MD_3D)))) - dilate(dilate(dilate(mask_L_MD_3D)));

# ╔═╡ e2a2de30-b253-4e0d-a0a8-edfd6d30ac95
@bind h4 overlay_mask_bind(ring_mask_L_MD)

# ╔═╡ 99c1afe5-b859-410b-9439-4399c8e5f87e
overlay_mask_plot(arr, ring_mask_L_MD, h4, "ring mask")

# ╔═╡ d35c5aae-2abe-4b94-a7e6-427f98d660cd
md"""
### Calculations
"""

# ╔═╡ 5ff89e4c-2787-43ef-a71c-f9572d9b462c
begin
	single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
	s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
end

# ╔═╡ 5a4d9e40-2c68-40b3-a045-618b784062f1
S_Obj_MD = intensity(400)

# ╔═╡ 7774595a-479b-4af5-bf17-666755f49085
begin
	alg_L_MD = Integrated(arr[mask_L_MD_3D])
	vol_l_md = score(s_bkg_L_MD, S_Obj_MD, pixel_size, alg_L_MD)
end

# ╔═╡ 3852dfd5-ae79-40a6-a408-e94e17d252d5
begin
	ρ_MD = 0.4 # mg/mm^3
	mass_l_md = score(s_bkg_L_MD, S_Obj_MD, pixel_size, ρ_MD, alg_L_MD)
end

# ╔═╡ c26d0521-1b5d-4d3d-89ac-37fa9bb0ecf0
md"""
## Low Density
"""

# ╔═╡ bfb9fdb4-505f-4e84-aa93-73184e535738
begin
	mask_L_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_LD_3D[:, :, z] = mask_L_LD
	end
end;

# ╔═╡ 4071f8b0-d8dd-48e4-bb7e-e8000403964e
md"""
#### Dilated mask
"""

# ╔═╡ 6aa87bd4-3f4f-4ff6-8ce2-5dcd032cb3b8
dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D));

# ╔═╡ 0887d091-e543-47c3-8b3c-6d2863d5f6a0
@bind i2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ 1698c092-dbeb-4fee-b589-de0bb01bc3f9
overlay_mask_plot(arr, dilated_mask_L_LD, i2, "dilated mask")

# ╔═╡ 26fe669c-db0e-4d9a-958b-18a4527d51fa
md"""
#### Ring (background) mask
"""

# ╔═╡ daaaddb3-15d9-40ba-b902-fbfaeec0e6df
ring_mask_L_LD = dilate(dilate(dilate(dilate(mask_L_LD_3D)))) - dilate(dilate(dilate(mask_L_LD_3D)));

# ╔═╡ 04333747-b72b-4680-8bfc-c09b355c64ba
@bind i4 overlay_mask_bind(ring_mask_L_LD)

# ╔═╡ 78afde73-cf4a-45cc-9afc-4f59c24b8f9d
overlay_mask_plot(arr, ring_mask_L_LD, i4, "ring mask")

# ╔═╡ ad8e96da-0f82-45c7-a410-849c75d1d1c1
md"""
### Calculations
"""

# ╔═╡ ef7453f4-7c80-41aa-85c7-b7df1235eb91
begin	
	single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
	s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
end

# ╔═╡ 0e431278-f0e2-43c9-a3e9-c196c1ec521d
S_Obj_LD = intensity(200)

# ╔═╡ f668013d-666c-469e-84a9-a25d9cce1c3c
begin
	alg_L_LD = Integrated(arr[mask_L_LD_3D])
	vol_l_ld = score(s_bkg_L_LD, S_Obj_LD, pixel_size, alg_L_LD)
end

# ╔═╡ e1cc4110-a976-438e-9d79-89ff84a8236a
begin
	ρ_LD = 0.2
	mass_l_ld = score(s_bkg_L_LD, S_Obj_LD, pixel_size, ρ_LD, alg_L_LD)
end

# ╔═╡ d9cd6868-4606-454c-bbaf-d480c52f270c
md"""
# Score Medium Inserts
"""

# ╔═╡ 410db4da-4f91-4d59-aad2-3613fc059763
md"""
## High Density
"""

# ╔═╡ 61af8807-029a-4034-947d-4b965fbb72b1
begin
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_HD_3D[:, :, z] = mask_M_HD
	end
end;

# ╔═╡ 6f5d85b1-1bdb-4603-9ab8-001e917d3e3b
md"""
#### Dilated mask
"""

# ╔═╡ cfbc4226-ab14-476f-93bd-ad29d7bd4608
dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ d56fe084-2fd8-43f5-871e-d84104305794
@bind j2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ d3d43075-51b3-493d-af0a-93f9279c8dac
overlay_mask_plot(arr, dilated_mask_M_HD, j2, "dilated mask")

# ╔═╡ b6da3cf0-a9d1-43c2-a7ca-f7bafa4d32e2
md"""
#### Ring (background) mask
"""

# ╔═╡ 3bfd53ca-87a1-4591-9c78-75147d811301
ring_mask_M_HD = dilate(dilate(dilate(dilate(dilate(mask_M_HD_3D))))) - dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 4b3d30e0-7901-41a7-9e55-0dc12ac25a9d
@bind j4 overlay_mask_bind(ring_mask_M_HD)

# ╔═╡ 1a33ca4f-24d7-4da0-926c-fbf5e62c800f
overlay_mask_plot(arr, ring_mask_M_HD, j4, "ring mask")

# ╔═╡ 4de69b74-1cf3-4ab9-9cca-f7e0608dfb7c
md"""
### Calculations
"""

# ╔═╡ 22618c42-04d8-401f-99c9-6d88ac59a75e
begin
	single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
	s_bkg_M_HD = mean(single_arr[single_ring_mask_M_HD])
end

# ╔═╡ e15f9ae6-493e-48c7-93bb-5e2b501fe068
begin
	alg_M_HD = Integrated(arr[mask_M_HD_3D])
	vol_m_hd = score(s_bkg_M_HD, S_Obj_HD, pixel_size, alg_M_HD)
end

# ╔═╡ 8e3718fb-6fac-4610-9872-48cf8e3bb540
mass_m_hd = score(s_bkg_M_HD, S_Obj_HD, pixel_size, ρ_HD, alg_M_HD)

# ╔═╡ 52cf7b88-2a09-4818-a89c-9338d704f0d8
md"""
## Medium Density
"""

# ╔═╡ c16a1dcd-55c3-4dd2-9aa3-03a25a8dcd11
begin
	mask_M_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_MD_3D[:, :, z] = mask_M_MD
	end
end;

# ╔═╡ 7fddea1e-b017-414e-a500-2b10407db00d
md"""
#### Dilated mask
"""

# ╔═╡ b6198218-e922-4cfb-b863-12d65a05b72c
dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))));

# ╔═╡ 36749afc-0a19-4915-aa23-faf0f3281ceb
@bind k2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ e3de7579-9509-4414-97db-a8a0d1d86770
overlay_mask_plot(arr, dilated_mask_M_MD, k2, "dilated mask")

# ╔═╡ 03563144-7dd4-4f6b-b4a4-cad2d85caa2e
md"""
#### Ring (background) mask
"""

# ╔═╡ a0fe3d8c-531a-434e-947e-c2e6743ae97b
ring_mask_M_MD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))));

# ╔═╡ 8de2647a-3b48-430f-bb5b-7ce04b8d1225
@bind k4 overlay_mask_bind(ring_mask_M_MD)

# ╔═╡ 76a6455f-8bc7-45d3-98a5-5cc8cff6e722
overlay_mask_plot(arr, ring_mask_M_MD, k4, "ring mask")

# ╔═╡ 093f944b-99d3-400f-885a-e27a3e89c8dc
md"""
### Calculations
"""

# ╔═╡ 6fe3c145-ea8e-4bdb-ad32-b455018e4b6e
begin
	single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
	s_bkg_M_MD = mean(single_arr[single_ring_mask_M_MD])
end

# ╔═╡ b9d7d5fe-d89d-4ff0-a08f-cdaed2bea4b7
begin
	alg_M_MD = Integrated(arr[mask_M_MD_3D])
	vol_m_md = score(s_bkg_M_MD, S_Obj_MD, pixel_size, alg_M_MD)
end

# ╔═╡ 511795c6-a26e-44fe-86e2-88030f822d72
mass_m_md = score(s_bkg_M_MD, S_Obj_MD, pixel_size, ρ_MD, alg_M_MD)

# ╔═╡ 2aa0f179-9093-46ac-962c-eeff5d56c503
md"""
## Low Density
"""

# ╔═╡ 1eb1e572-e70a-4ba3-bce6-91aac8766c92
begin
	mask_M_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_LD_3D[:, :, z] = mask_M_LD
	end
end;

# ╔═╡ 6785b3a8-9f85-44f9-b87e-43a3adecb143
md"""
#### Dilated mask
"""

# ╔═╡ 6460961c-b1ce-4652-8911-60ee9c28d3bb
dilated_mask_M_LD = dilate(dilate(dilate(dilate(mask_M_LD_3D))));

# ╔═╡ a3ec0af5-b28f-40c3-ade5-1deb31c645c9
@bind l2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ c19f0cfa-e2e3-4c66-a4e0-85828eb5cc66
overlay_mask_plot(arr, dilated_mask_M_LD, l2, "dilated mask")

# ╔═╡ cb7eb2c0-c571-4f04-9abb-267333892665
md"""
#### Ring (background) mask
"""

# ╔═╡ 86fdc408-22cf-446d-ba81-849b8baf93eb
ring_mask_M_LD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ 5a9a49bd-7ffd-45d6-8824-1120b54522bf
@bind l4 overlay_mask_bind(ring_mask_M_LD)

# ╔═╡ 2b99fc27-e464-4b31-a362-0a5ff5a8c709
overlay_mask_plot(arr, ring_mask_M_LD, l4, "ring mask")

# ╔═╡ d3a7a251-56ac-4d81-a0fa-5c3dcb925f6b
md"""
### Calculations
"""

# ╔═╡ da4ad8a3-36de-4345-9c88-21eb15b8b5d0
begin
	single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
	s_bkg_M_LD = mean(single_arr[single_ring_mask_M_LD])
end

# ╔═╡ d6cc9a06-231e-43d3-a6ee-8241e1a0902b
begin
	alg_M_LD = Integrated(arr[mask_M_LD_3D])
	vol_m_ld = score(s_bkg_M_LD, S_Obj_LD, pixel_size, alg_M_LD)
end

# ╔═╡ 978fe827-d2ab-4c7e-a0f9-3a2090d5a3de
mass_m_ld = score(s_bkg_M_LD, S_Obj_LD, pixel_size, ρ_LD, alg_M_LD)

# ╔═╡ 76a5f4cd-3669-4a63-a6ad-8574ceab79f8
md"""
# Score Small Inserts
"""

# ╔═╡ 111cfe3f-b39d-4498-a160-80046d07e764
md"""
## High Density
"""

# ╔═╡ 26bc388e-c11f-4752-ad1c-a89c08f28976
begin
	mask_S_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_HD_3D[:, :, z] = mask_S_HD
	end
end;

# ╔═╡ 538f2c30-6c51-40b8-9ae5-77d6da9fe8d8
md"""
#### Dilated mask
"""

# ╔═╡ 23fc7f16-8d75-4c8a-aefd-435e094f2f42
dilated_mask_S_HD = dilate(dilate(dilate(dilate(mask_S_HD_3D))));

# ╔═╡ 00fd89db-4094-4f72-a71a-a2fd98a1c0b7
@bind m2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ a78431e2-2910-4422-af23-63b6886694e6
overlay_mask_plot(arr, dilated_mask_S_HD, m2, "dilated mask")

# ╔═╡ b48f6d46-8483-4bf8-bc0c-da21489747e9
md"""
#### Ring (background) mask
"""

# ╔═╡ 6527ecb4-5514-41c9-abc5-967884c8be58
ring_mask_S_HD = dilate(dilate(dilate(dilate(dilate(mask_S_HD_3D))))) - dilate(dilate(dilate(dilate(mask_S_HD_3D))));

# ╔═╡ 193a11b5-85f6-47a5-adca-97a9e700e37a
@bind m4 overlay_mask_bind(ring_mask_S_HD)

# ╔═╡ f8d22847-78bc-4608-82f4-e801366743b9
overlay_mask_plot(arr, ring_mask_S_HD, m4, "ring mask")

# ╔═╡ 5fb21bc8-7bb4-4570-aaa9-87fb48e60721
md"""
### Calculations
"""

# ╔═╡ 8eb8a862-b505-4c54-b7b2-d01a779a2d0d
begin
	single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
	s_bkg_S_HD = mean(single_arr[single_ring_mask_S_HD])
end

# ╔═╡ a35ac97c-49ef-4b39-b3b9-3be9a7c46769
begin
	alg_S_HD = Integrated(arr[mask_S_HD_3D])
	vol_s_hd = score(s_bkg_S_HD, S_Obj_HD, pixel_size, alg_S_HD)
	if vol_s_hd < 0
		vol_s_hd = 0
	end
end

# ╔═╡ 3edad2b5-5221-44bc-9606-93f04cd4fb1f
begin
	mass_s_hd = score(s_bkg_S_HD, S_Obj_HD, pixel_size, ρ_HD, alg_S_HD)
	if mass_s_hd < 0
		mass_s_hd = 0
	end
end

# ╔═╡ 55f63e94-d284-4c5d-ba2f-240115ecd5ea
md"""
## Medium Density
"""

# ╔═╡ db16dde6-7753-4380-9f67-0d027da69ab5
begin
	mask_S_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_MD_3D[:, :, z] = mask_S_MD
	end
end;

# ╔═╡ cf9c0f01-fcd8-46a8-83f5-5ec3fae61eeb
md"""
#### Dilated mask
"""

# ╔═╡ 8fd6f1ff-d821-491b-95d9-20558e2d10e9
dilated_mask_S_MD = dilate(dilate(dilate(dilate(mask_S_MD_3D))));

# ╔═╡ fba330f6-5964-48cb-804f-f53e74740213
@bind n2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ d2140b59-4f0a-4a39-bc83-f70c64a12154
overlay_mask_plot(arr, dilated_mask_S_MD, n2, "dilated mask")

# ╔═╡ a189e091-7607-4f1b-8bc7-1503241ef64b
md"""
#### Ring (background) mask
"""

# ╔═╡ 5da622df-370a-4139-ae3f-eceb56c2f705
ring_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D))))) - dilate(dilate(dilate(dilate(mask_S_MD_3D))));

# ╔═╡ c835eb9c-45a9-4759-8e7d-c86638e94479
@bind n4 overlay_mask_bind(ring_mask_S_MD)

# ╔═╡ 7baa6dd9-4718-4e36-8ea8-8f906abceeb0
overlay_mask_plot(arr, ring_mask_S_MD, n4, "ring mask")

# ╔═╡ a52a944b-978a-47fe-9a5d-bab6bd09f39f
md"""
### Calculations
"""

# ╔═╡ 8420fe90-af8d-4c06-b536-a247622d5027
begin
	single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
	s_bkg_S_MD = mean(single_arr[single_ring_mask_S_MD])
end

# ╔═╡ 35d99a2b-f92d-43e9-806c-c73a66cc5a52
begin
	alg_S_MD = Integrated(arr[mask_S_MD_3D])
	vol_s_md = score(s_bkg_S_MD, S_Obj_HD, pixel_size, alg_S_MD)
	if vol_s_md < 0
		vol_s_md = 0
	end
end

# ╔═╡ 06497bc3-0238-423e-80a3-be97fd1644a1
begin
	mass_s_md = score(s_bkg_S_MD, S_Obj_HD, pixel_size, ρ_MD, alg_S_MD)
	if mass_s_md < 0
		mass_s_md = 0
	end
end

# ╔═╡ 022460e2-59f2-45ca-9472-adeddf79590d
md"""
## Low Density
"""

# ╔═╡ 9627200e-97ce-4499-93f1-e51fc3899e48
begin
	mask_S_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_LD_3D[:, :, z] = mask_S_LD
	end
end;

# ╔═╡ c03de923-33b9-4614-ab94-07345d6928a7
md"""
#### Dilated mask
"""

# ╔═╡ d7b9f44e-8b1b-4c0e-8068-2dad890e6690
dilated_mask_S_LD = dilate(dilate(dilate(dilate(mask_S_LD_3D))));

# ╔═╡ f2884ccf-faab-4113-a5ee-a57c35b25908
@bind o2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ 5a9f9185-660e-4ab0-8260-6bdef67f298f
overlay_mask_plot(arr, dilated_mask_S_LD, o2, "dilated mask")

# ╔═╡ 613b92f0-d97b-49a9-b3d0-0e358903df37
md"""
#### Ring (background) mask
"""

# ╔═╡ c5c5d5e3-9e80-41c4-a493-53c642a7e0dd
ring_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D))))) - dilate(dilate(dilate(dilate(mask_S_LD_3D))));

# ╔═╡ 33f1e03d-e538-4a7a-9a7e-ff56954ca126
@bind o4 overlay_mask_bind(ring_mask_S_LD)

# ╔═╡ 47c0944c-39c9-4e5f-b939-45d25fcf4348
overlay_mask_plot(arr, ring_mask_S_LD, o4, "ring mask")

# ╔═╡ ee50b165-3ee3-44f4-8e99-3e529b31e416
md"""
### Calculations
"""

# ╔═╡ 8302bf47-a164-46e2-8c0c-208c1e9695b8
begin
	single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
	s_bkg_S_LD = mean(single_arr[single_ring_mask_S_LD])
end

# ╔═╡ 7f828ef6-c6fe-442b-94d4-f8c6a4e0a2b0
begin
	alg_S_LD = Integrated(arr[mask_S_LD_3D])
	vol_s_ld = score(s_bkg_S_LD, S_Obj_HD, pixel_size, alg_S_LD)
	if vol_s_ld < 0
		vol_s_ld = 0
	end
end

# ╔═╡ 345e211d-72d2-477e-972e-e4a34664ab82
begin
	mass_s_ld = score(s_bkg_S_LD, S_Obj_HD, pixel_size, ρ_LD, alg_S_LD)
	if mass_s_ld < 0
		mass_s_ld = 0
	end
end

# ╔═╡ bb794b5f-d4f8-4329-81c8-b05280e30edc
md"""
# Results
"""

# ╔═╡ 08c88f23-0223-4e69-b562-4f156bb5cc38
md"""
### Volume
"""

# ╔═╡ 709b68e7-ce66-4726-ab40-e0a3f5cb403f
inserts = [
	"Low Density",
	"Medium Density",
	"High Density"
]

# ╔═╡ 0042c385-edf4-4964-a46f-f6bf2fa9c900
ground_truth_volume_large = [
	98.2,
	98.2,
	98.2,
] # mm^3

# ╔═╡ 6a8a87bd-e9ac-46cd-9f9d-2356c6de48bb
calculated_volume_large = [
	vol_l_ld,
	vol_l_md,
	vol_l_hd
]

# ╔═╡ caef1c5e-e77f-4ed0-99ca-148e113052c0
ground_truth_volume_medium = [
	21.2,
	21.2,
	21.2
]

# ╔═╡ 256d4644-bba4-45d6-908a-8505ae291c93
calculated_volume_medium = [
	vol_m_ld,
	vol_m_md,
	vol_m_hd
]

# ╔═╡ 01941804-7ff6-498e-9a00-53467d3bef6b
ground_truth_volume_small = [
	0.8,
	0.8,
	0.8
]

# ╔═╡ d9ecccb6-6684-449d-8386-f391ec22f255
calculated_volume_small = [
	vol_s_ld,
	vol_s_md,
	vol_s_hd
]

# ╔═╡ 173cb04c-5d3d-4bc6-930c-043e7b923ae6
df1 = DataFrame(
	inserts = inserts,
	ground_truth_volume_large = ground_truth_volume_large,
	calculated_volume_large = calculated_volume_large,
	ground_truth_volume_medium = ground_truth_volume_medium,
	calculated_volume_medium = calculated_volume_medium,
	ground_truth_volume_small = ground_truth_volume_small,
	calculated_volume_small = calculated_volume_small
)

# ╔═╡ 0fa9a2e1-8dee-4f91-90c6-d40dbe9d0f23
begin
	fvol2 = Figure()
	axvol2 = Axis(fvol2[1, 1])
	
	scatter!(density_array[2:end], df1[!, :ground_truth_volume_large], label="ground_truth_volume_large")
	scatter!(density_array[2:end], df1[!, :calculated_volume_large], label="calculated_volume_large")
	
	axvol2.title = "Volume Measurements (Large)"
	axvol2.ylabel = "Volume (mm^3)"
	axvol2.xlabel = "Density (mg/cm^3)"

	xlims!(axvol2, 0, 850)
	ylims!(axvol2, 0, 130)
	
	fvol2[1, 2] = Legend(fvol2, axvol2, framevisible = false)
	
	fvol2
end

# ╔═╡ 5316a9f1-4ca5-47bd-a187-d2787c500f84
begin
	fvol3 = Figure()
	axvol3 = Axis(fvol3[1, 1])
	
	scatter!(density_array[2:end], df1[!, :ground_truth_volume_medium], label="ground_truth_volume_medium")
	scatter!(density_array[2:end], df1[!, :calculated_volume_medium], label="calculated_volume_medium")
	
	axvol3.title = "Volume Measurements (Medium)"
	axvol3.ylabel = "Volume (mm^3)"
	axvol3.xlabel = "Density (mg/cm^3)"

	xlims!(axvol3, 0, 850)
	ylims!(axvol3, 0, 50)
	
	fvol3[1, 2] = Legend(fvol3, axvol3, framevisible = false)
	
	fvol3
end

# ╔═╡ 8a0c68e4-e74a-46d4-9ccf-50815e192ee0
begin
	fvol4 = Figure()
	axvol4 = Axis(fvol4[1, 1])
	
	scatter!(density_array[2:end], df1[!, :ground_truth_volume_small], label="ground_truth_volume_small")
	scatter!(density_array[2:end], df1[!, :calculated_volume_small], label="calculated_volume_small")
	
	axvol4.title = "Volume Measurements (Small)"
	axvol4.ylabel = "Volume (mm^3)"
	axvol4.xlabel = "Density (mg/cm^3)"

	xlims!(axvol4, 0, 850)
	ylims!(axvol4, 0, 2.5)
	
	fvol4[1, 2] = Legend(fvol4, axvol4, framevisible = false)
	
	fvol4
end

# ╔═╡ 8c3e0174-00cb-4310-b112-1c39916a637e
md"""
### Mass
"""

# ╔═╡ da58681d-7b87-4a4c-b50b-f45978fae449
ground_truth_mass_large = [
	19.6,
	39.3,
	78.5
] # mg

# ╔═╡ ef0c1a87-4c24-4f4f-9c03-3872f935e8c5
calculated_mass_large = [
	mass_l_ld,
	mass_l_md,
	mass_l_hd
]

# ╔═╡ 908fbe8b-e911-48cd-8a10-fc6c6ed64f16
ground_truth_mass_medium = [
	4.2,
	8.5,
	17.0
]

# ╔═╡ 667aa294-c0e7-45f0-a2f1-eb798d7e6581
calculated_mass_medium = [
	mass_m_ld,
	mass_m_md,
	mass_m_hd
]

# ╔═╡ f470197f-a528-4d8b-a546-9988bc7ed7b9
ground_truth_mass_small = [
	0.2,
	0.3,
	0.6
]

# ╔═╡ 74e12a13-32fd-4f30-a9ec-f44765144708
calculated_mass_small = [
	mass_s_ld,
	mass_s_md,
	mass_s_hd
]

# ╔═╡ 6c576d2b-cf94-4c8d-bfeb-a772ad295d7b
df2 = DataFrame(
	inserts = inserts,
	ground_truth_mass_large = ground_truth_mass_large,
	calculated_mass_large = calculated_mass_large,
	ground_truth_mass_medium = ground_truth_mass_medium,
	calculated_mass_medium = calculated_mass_medium,
	ground_truth_mass_small = ground_truth_mass_small,
	calculated_mass_small = calculated_mass_small
)

# ╔═╡ bcd80f6a-3154-4b73-8e2c-53c57358d5ee
begin
	fmass2 = Figure()
	axmass2 = Axis(fmass2[1, 1])
	
	scatter!(density_array[2:end], df2[!, :ground_truth_mass_large], label="ground_truth_mass_large")
	scatter!(density_array[2:end], df2[!, :calculated_mass_large], label="calculated_mass_large")
	
	axmass2.title = "Mass Measurements (Large)"
	axmass2.ylabel = "Mass (mg)"
	axmass2.xlabel = "Density (mg/cm^3)"

	xlims!(axmass2, 0, 850)
	ylims!(axmass2, 0, 100)
	
	fmass2[1, 2] = Legend(fmass2, axmass2, framevisible = false)
	
	fmass2
end

# ╔═╡ 4fe08534-fe92-4d2d-aa92-4a4503ccb92c
begin
	fmass3 = Figure()
	axmass3 = Axis(fmass3[1, 1])
	
	scatter!(density_array[2:end], df2[!, :ground_truth_mass_medium], label="ground_truth_mass_medium")
	scatter!(density_array[2:end], df2[!, :calculated_mass_medium], label="calculated_mass_medium")
	
	axmass3.title = "Mass Measurements (Medium)"
	axmass3.ylabel = "Mass (mg)"
	axmass3.xlabel = "Density (mg/cm^3)"

	xlims!(axmass3, 0, 850)
	ylims!(axmass3, 0, 25)
	
	fmass3[1, 2] = Legend(fmass3, axmass3, framevisible = false)
	
	fmass3
end

# ╔═╡ f0b92d0e-3730-47f9-8aa5-0f5d542d3c0c
begin
	fmass4 = Figure()
	axmass4 = Axis(fmass4[1, 1])
	
	scatter!(density_array[2:end], df2[!, :ground_truth_mass_small], label="ground_truth_mass_small")
	scatter!(density_array[2:end], df2[!, :calculated_mass_small], label="calculated_mass_small")
	
	axmass4.title = "Mass Measurements (Small)"
	axmass4.ylabel = "Mass (mg)"
	axmass4.xlabel = "Density (mg/cm^3)"

	xlims!(axmass4, 0, 850)
	ylims!(axmass4, 0, 1.5)
	
	fmass4[1, 2] = Legend(fmass4, axmass4, framevisible = false)
	
	fmass4
end

# ╔═╡ a555540c-c6b9-45af-b57b-4f87f6cb7e16
md"""
### Save Results
"""

# ╔═╡ f6f63c5c-d99f-413b-b19d-3e03aad32d67
df_final = leftjoin(df1, df2, on=:inserts)

# ╔═╡ 157b7758-9481-466e-9b75-33c5d7a02521
if ~isdir(string(cd(pwd, "..") , "/data/output/", VENDER, "3"))
	mkdir(string(cd(pwd, "..") , "/data/output/", VENDER, "3"))
end

# ╔═╡ 9e5e5b78-a752-442e-a8ce-d0170aeaa88e
output_path = string(cd(pwd, "..") , "/data/output/", VENDER, "3", "/", scan, ".csv")

# ╔═╡ 5e3d48a5-fdfb-4dab-900e-6a7538d0b5f2
CSV.write(output_path, df_final)

# ╔═╡ Cell order:
# ╠═c3e3a7a8-9d0b-4363-943d-d723f7f4cb37
# ╠═ff527e22-8eb6-413e-b45e-7c07d7623f45
# ╟─903cf8bf-1113-4a0f-9e07-122ba61d34cb
# ╠═88af5c56-9841-48ff-95fc-cceddec0184f
# ╟─8024a262-7958-4d21-bbe0-7b2d95a562f1
# ╠═7addd762-8754-4ba2-a05e-493f7171847e
# ╠═992dd715-3740-4f8e-b623-6b99f61ebc90
# ╠═cf504678-ae5b-4327-b956-b3dc1ffac6df
# ╠═7d8776eb-6623-4cde-b91d-e2a91e50b377
# ╠═60842b6a-72d2-4501-911e-19388c9df4d9
# ╟─c7dc5fd3-9f08-43f3-86a7-b77a174f76ca
# ╟─6fb61984-75c8-4d2f-b5ad-d99b6cb382ad
# ╟─a241f5db-3e70-4893-8b19-3f05f1b94a5a
# ╟─6b403990-eef6-42ca-8cf0-e479c5eca2f8
# ╟─4fefaf32-4cb5-4cf8-b94d-4695d13b5797
# ╠═b1d48c95-a6c4-4b6a-8b14-fa34976cad83
# ╠═8f6f0951-8f5a-483d-81bb-0928f2f765ce
# ╠═ec6bde85-0dec-46c4-9c59-14b3126927a5
# ╠═985c7ce3-3414-4fa4-8c11-647fd3cbf2ce
# ╠═27ed5fe4-775b-40a9-8821-ef8b4283e40c
# ╠═749acc75-b49e-4a24-b1a2-c8e57bfdcb58
# ╟─f63c9679-3c11-4a28-8f6d-7c903419930a
# ╠═39c633c8-af27-4d92-804d-5b0263c40e91
# ╠═2dc7cf6e-2626-496f-84ee-629320f25bc8
# ╠═faa85dd6-303d-478f-ad12-148635f9d02d
# ╟─ae76e2c4-a74e-496a-b6a4-17fe91f2ad57
# ╠═cbf71130-ddcf-4418-9e3e-7df8fc57c12f
# ╠═212ed9e1-b5de-4f0f-bf71-c811908c5d45
# ╠═b87fd07d-8cea-4a5a-94dc-6682d66e8bd8
# ╟─5b9d4af0-3fa7-4e33-8888-4e3b7d1985ba
# ╠═5484979d-aa3a-4eab-a1eb-2b81a83068f4
# ╟─8557a14a-17b2-4039-a110-f35c87ea77b6
# ╠═eab14ca7-4299-42d1-b0de-f85e6cf0386c
# ╠═43a5e31f-2c82-4980-a4eb-ce0ae1a5f448
# ╠═e00ef290-7ffa-4a43-81d3-735b4191dc98
# ╠═431ff05c-804a-4e86-ac77-c6dd5e933d35
# ╟─781db9e6-f71b-42ff-a060-0d3d11cd8db8
# ╠═a2262988-1f50-4622-a603-e0925ad0d065
# ╠═d2033181-98a8-4a46-b435-2acf6eb72485
# ╠═1c62ece3-da87-48a6-a2c0-fa536ead7e73
# ╠═ddeda291-edb8-4bc2-b1cf-130ea996c086
# ╟─7695a3b7-31c6-4a31-9fcb-24901851e269
# ╠═c85b86ba-4f5b-4bfb-a26e-36ac185dc794
# ╠═3fff7d44-1586-492a-8e1e-dfcda0979106
# ╠═d4cb8eba-17cf-4aa3-8ad4-e8ac35811caf
# ╠═74b034ae-a9e1-4f72-bc75-5860047ef8c5
# ╟─dcf800fc-7dea-4197-b226-504d76569554
# ╠═ffd33235-091b-4154-a3e6-e76ba1d1b1fa
# ╠═bcbc253e-9f66-4488-af84-10022e4d765d
# ╠═c440667a-deb5-4e21-b3ec-6e715cda538b
# ╠═1378a515-b757-4776-b5aa-864ddc89e94b
# ╠═264fe678-9b84-4a2a-9a66-ca2e6dd3b65d
# ╠═243b5511-ec59-4a8a-a29e-334ae2fa8985
# ╠═9a8c3d70-5eb5-46b4-a8a5-239335512196
# ╟─e8cc070c-1850-474c-a2fa-d0c12a09f683
# ╠═885ffebe-4805-43e6-a305-c7076c75ec5b
# ╠═26a4d12b-9fcb-4e90-a96e-5ed9c8cfe2ce
# ╠═847262c3-dd0a-4647-81fa-ea7a2a441605
# ╟─d9600f56-c238-4392-a262-0a25ae19f9ba
# ╠═158347a5-966f-42f7-ae1b-b41b39d79b8b
# ╠═98dcf986-4e0e-45dd-a52a-7886ce868842
# ╟─953cd854-108f-486d-ac58-0405b6b785c4
# ╠═2024648c-45f4-403a-b9e7-11b2ec3673ae
# ╟─606c40bb-b9c9-4ec4-8505-77eba0799ce5
# ╠═dbf78411-7e46-439e-aaf9-640a88354d39
# ╠═8ea123b6-4eba-41db-abe8-c994f76c9cda
# ╠═ba2fed15-a79e-4b21-a575-5884cce8c366
# ╟─e17bd8ad-b2b1-4020-a7db-1ec47af05b09
# ╠═de607584-2f2e-4e9e-9c22-77cb05596b44
# ╠═7d655971-44da-42ad-9420-86528540a9b2
# ╠═476309ca-a32b-49b8-8845-3fdaf4a2d847
# ╟─7e6f8c7f-3f30-4d41-aa78-77c8d623a499
# ╠═cf725101-4564-4546-bc1f-d073d2f1cc04
# ╠═0b502b66-da66-44cf-8531-945e5a466ec1
# ╠═3bdc0b37-67d7-49fa-886d-f4b66fe28c82
# ╠═8c73a297-5d95-43c9-ba2c-09847fc5fc71
# ╠═89edc6a0-2933-4dba-be63-48d3be00444e
# ╟─c5a58e31-5f64-4d2e-b222-d918d6d41b5d
# ╠═f7b97123-ac55-49d8-9bd9-45089245c22a
# ╟─45379669-4897-460c-8ba9-f1635f3bcb13
# ╠═5c98b68c-01e9-4cb1-954d-6f0192906e6f
# ╠═cd0aef7b-fdc9-4ebc-a9d8-23d8655144ef
# ╠═1a5150b9-bab3-4ec5-92b0-4734917ec58f
# ╟─73a5aee9-38a4-4581-b654-8a95cdf74108
# ╠═d86fa8f6-8e54-44b3-8188-4c0373d5e633
# ╠═e2a2de30-b253-4e0d-a0a8-edfd6d30ac95
# ╠═99c1afe5-b859-410b-9439-4399c8e5f87e
# ╟─d35c5aae-2abe-4b94-a7e6-427f98d660cd
# ╠═5ff89e4c-2787-43ef-a71c-f9572d9b462c
# ╠═5a4d9e40-2c68-40b3-a045-618b784062f1
# ╠═7774595a-479b-4af5-bf17-666755f49085
# ╠═3852dfd5-ae79-40a6-a408-e94e17d252d5
# ╟─c26d0521-1b5d-4d3d-89ac-37fa9bb0ecf0
# ╠═bfb9fdb4-505f-4e84-aa93-73184e535738
# ╟─4071f8b0-d8dd-48e4-bb7e-e8000403964e
# ╠═6aa87bd4-3f4f-4ff6-8ce2-5dcd032cb3b8
# ╠═0887d091-e543-47c3-8b3c-6d2863d5f6a0
# ╠═1698c092-dbeb-4fee-b589-de0bb01bc3f9
# ╟─26fe669c-db0e-4d9a-958b-18a4527d51fa
# ╠═daaaddb3-15d9-40ba-b902-fbfaeec0e6df
# ╠═04333747-b72b-4680-8bfc-c09b355c64ba
# ╠═78afde73-cf4a-45cc-9afc-4f59c24b8f9d
# ╟─ad8e96da-0f82-45c7-a410-849c75d1d1c1
# ╠═ef7453f4-7c80-41aa-85c7-b7df1235eb91
# ╠═0e431278-f0e2-43c9-a3e9-c196c1ec521d
# ╠═f668013d-666c-469e-84a9-a25d9cce1c3c
# ╠═e1cc4110-a976-438e-9d79-89ff84a8236a
# ╟─d9cd6868-4606-454c-bbaf-d480c52f270c
# ╟─410db4da-4f91-4d59-aad2-3613fc059763
# ╠═61af8807-029a-4034-947d-4b965fbb72b1
# ╟─6f5d85b1-1bdb-4603-9ab8-001e917d3e3b
# ╠═cfbc4226-ab14-476f-93bd-ad29d7bd4608
# ╠═d56fe084-2fd8-43f5-871e-d84104305794
# ╠═d3d43075-51b3-493d-af0a-93f9279c8dac
# ╟─b6da3cf0-a9d1-43c2-a7ca-f7bafa4d32e2
# ╠═3bfd53ca-87a1-4591-9c78-75147d811301
# ╠═4b3d30e0-7901-41a7-9e55-0dc12ac25a9d
# ╠═1a33ca4f-24d7-4da0-926c-fbf5e62c800f
# ╟─4de69b74-1cf3-4ab9-9cca-f7e0608dfb7c
# ╠═22618c42-04d8-401f-99c9-6d88ac59a75e
# ╠═e15f9ae6-493e-48c7-93bb-5e2b501fe068
# ╠═8e3718fb-6fac-4610-9872-48cf8e3bb540
# ╟─52cf7b88-2a09-4818-a89c-9338d704f0d8
# ╠═c16a1dcd-55c3-4dd2-9aa3-03a25a8dcd11
# ╟─7fddea1e-b017-414e-a500-2b10407db00d
# ╠═b6198218-e922-4cfb-b863-12d65a05b72c
# ╠═36749afc-0a19-4915-aa23-faf0f3281ceb
# ╠═e3de7579-9509-4414-97db-a8a0d1d86770
# ╟─03563144-7dd4-4f6b-b4a4-cad2d85caa2e
# ╠═a0fe3d8c-531a-434e-947e-c2e6743ae97b
# ╠═8de2647a-3b48-430f-bb5b-7ce04b8d1225
# ╠═76a6455f-8bc7-45d3-98a5-5cc8cff6e722
# ╟─093f944b-99d3-400f-885a-e27a3e89c8dc
# ╠═6fe3c145-ea8e-4bdb-ad32-b455018e4b6e
# ╠═b9d7d5fe-d89d-4ff0-a08f-cdaed2bea4b7
# ╠═511795c6-a26e-44fe-86e2-88030f822d72
# ╟─2aa0f179-9093-46ac-962c-eeff5d56c503
# ╠═1eb1e572-e70a-4ba3-bce6-91aac8766c92
# ╟─6785b3a8-9f85-44f9-b87e-43a3adecb143
# ╠═6460961c-b1ce-4652-8911-60ee9c28d3bb
# ╠═a3ec0af5-b28f-40c3-ade5-1deb31c645c9
# ╠═c19f0cfa-e2e3-4c66-a4e0-85828eb5cc66
# ╟─cb7eb2c0-c571-4f04-9abb-267333892665
# ╠═86fdc408-22cf-446d-ba81-849b8baf93eb
# ╠═5a9a49bd-7ffd-45d6-8824-1120b54522bf
# ╠═2b99fc27-e464-4b31-a362-0a5ff5a8c709
# ╟─d3a7a251-56ac-4d81-a0fa-5c3dcb925f6b
# ╠═da4ad8a3-36de-4345-9c88-21eb15b8b5d0
# ╠═d6cc9a06-231e-43d3-a6ee-8241e1a0902b
# ╠═978fe827-d2ab-4c7e-a0f9-3a2090d5a3de
# ╟─76a5f4cd-3669-4a63-a6ad-8574ceab79f8
# ╟─111cfe3f-b39d-4498-a160-80046d07e764
# ╠═26bc388e-c11f-4752-ad1c-a89c08f28976
# ╟─538f2c30-6c51-40b8-9ae5-77d6da9fe8d8
# ╠═23fc7f16-8d75-4c8a-aefd-435e094f2f42
# ╠═00fd89db-4094-4f72-a71a-a2fd98a1c0b7
# ╠═a78431e2-2910-4422-af23-63b6886694e6
# ╟─b48f6d46-8483-4bf8-bc0c-da21489747e9
# ╠═6527ecb4-5514-41c9-abc5-967884c8be58
# ╠═193a11b5-85f6-47a5-adca-97a9e700e37a
# ╠═f8d22847-78bc-4608-82f4-e801366743b9
# ╟─5fb21bc8-7bb4-4570-aaa9-87fb48e60721
# ╠═8eb8a862-b505-4c54-b7b2-d01a779a2d0d
# ╠═a35ac97c-49ef-4b39-b3b9-3be9a7c46769
# ╠═3edad2b5-5221-44bc-9606-93f04cd4fb1f
# ╟─55f63e94-d284-4c5d-ba2f-240115ecd5ea
# ╠═db16dde6-7753-4380-9f67-0d027da69ab5
# ╟─cf9c0f01-fcd8-46a8-83f5-5ec3fae61eeb
# ╠═8fd6f1ff-d821-491b-95d9-20558e2d10e9
# ╠═fba330f6-5964-48cb-804f-f53e74740213
# ╠═d2140b59-4f0a-4a39-bc83-f70c64a12154
# ╟─a189e091-7607-4f1b-8bc7-1503241ef64b
# ╠═5da622df-370a-4139-ae3f-eceb56c2f705
# ╠═c835eb9c-45a9-4759-8e7d-c86638e94479
# ╠═7baa6dd9-4718-4e36-8ea8-8f906abceeb0
# ╟─a52a944b-978a-47fe-9a5d-bab6bd09f39f
# ╠═8420fe90-af8d-4c06-b536-a247622d5027
# ╠═35d99a2b-f92d-43e9-806c-c73a66cc5a52
# ╠═06497bc3-0238-423e-80a3-be97fd1644a1
# ╟─022460e2-59f2-45ca-9472-adeddf79590d
# ╠═9627200e-97ce-4499-93f1-e51fc3899e48
# ╟─c03de923-33b9-4614-ab94-07345d6928a7
# ╠═d7b9f44e-8b1b-4c0e-8068-2dad890e6690
# ╠═f2884ccf-faab-4113-a5ee-a57c35b25908
# ╠═5a9f9185-660e-4ab0-8260-6bdef67f298f
# ╟─613b92f0-d97b-49a9-b3d0-0e358903df37
# ╠═c5c5d5e3-9e80-41c4-a493-53c642a7e0dd
# ╠═33f1e03d-e538-4a7a-9a7e-ff56954ca126
# ╠═47c0944c-39c9-4e5f-b939-45d25fcf4348
# ╟─ee50b165-3ee3-44f4-8e99-3e529b31e416
# ╠═8302bf47-a164-46e2-8c0c-208c1e9695b8
# ╠═7f828ef6-c6fe-442b-94d4-f8c6a4e0a2b0
# ╠═345e211d-72d2-477e-972e-e4a34664ab82
# ╟─bb794b5f-d4f8-4329-81c8-b05280e30edc
# ╠═08c88f23-0223-4e69-b562-4f156bb5cc38
# ╠═709b68e7-ce66-4726-ab40-e0a3f5cb403f
# ╠═0042c385-edf4-4964-a46f-f6bf2fa9c900
# ╠═6a8a87bd-e9ac-46cd-9f9d-2356c6de48bb
# ╠═caef1c5e-e77f-4ed0-99ca-148e113052c0
# ╠═256d4644-bba4-45d6-908a-8505ae291c93
# ╠═01941804-7ff6-498e-9a00-53467d3bef6b
# ╠═d9ecccb6-6684-449d-8386-f391ec22f255
# ╠═173cb04c-5d3d-4bc6-930c-043e7b923ae6
# ╠═0fa9a2e1-8dee-4f91-90c6-d40dbe9d0f23
# ╠═5316a9f1-4ca5-47bd-a187-d2787c500f84
# ╠═8a0c68e4-e74a-46d4-9ccf-50815e192ee0
# ╟─8c3e0174-00cb-4310-b112-1c39916a637e
# ╠═da58681d-7b87-4a4c-b50b-f45978fae449
# ╠═ef0c1a87-4c24-4f4f-9c03-3872f935e8c5
# ╠═908fbe8b-e911-48cd-8a10-fc6c6ed64f16
# ╠═667aa294-c0e7-45f0-a2f1-eb798d7e6581
# ╠═f470197f-a528-4d8b-a546-9988bc7ed7b9
# ╠═74e12a13-32fd-4f30-a9ec-f44765144708
# ╠═6c576d2b-cf94-4c8d-bfeb-a772ad295d7b
# ╠═bcd80f6a-3154-4b73-8e2c-53c57358d5ee
# ╠═4fe08534-fe92-4d2d-aa92-4a4503ccb92c
# ╠═f0b92d0e-3730-47f9-8aa5-0f5d542d3c0c
# ╟─a555540c-c6b9-45af-b57b-4f87f6cb7e16
# ╠═f6f63c5c-d99f-413b-b19d-3e03aad32d67
# ╠═157b7758-9481-466e-9b75-33c5d7a02521
# ╠═9e5e5b78-a752-442e-a8ce-d0170aeaa88e
# ╠═5e3d48a5-fdfb-4dab-900e-6a7538d0b5f2

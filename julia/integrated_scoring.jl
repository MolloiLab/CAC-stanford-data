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

# ╔═╡ 52d4c3cb-c0a6-497a-af42-67eb5372cf21
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

# ╔═╡ 9bffa403-b54c-49f4-b4e4-3ec593c3a400
TableOfContents()

# ╔═╡ afa84ebb-9e0d-4f0a-b980-958b0dbddd1a
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDER`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ 857ea37f-58dd-45d4-96e5-8238718d9a3a
begin
	SCAN_NUMBER = 5
	VENDER = "Canon_Aquilion_One_Vision"
	BASE_PATH = "/Users/daleblack/Google Drive/Datasets/"
end

# ╔═╡ aa2e4340-9eb8-42c7-9765-e27859841fa7
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ 08c9aa2d-77c1-45c3-bbe9-889516fbdf7d
root_path = string(BASE_PATH, VENDER)

# ╔═╡ 2b5dfe45-acb2-48ec-b2ec-9b08ad573973
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ a770154a-491b-4e56-95c6-a1b1378eec7b
pth = dcm_path_list[SCAN_NUMBER]

# ╔═╡ 9dca8618-5a6a-4cae-91a3-bf005de45d13
pth

# ╔═╡ 6daad570-d417-49b1-84b2-dc78cf4ed752
scan = basename(pth)

# ╔═╡ aceb0699-06e1-4c6c-a178-0aaa47683930
header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ d1705d56-40ab-4119-adb7-876126981d8b
md"""
## Helper Functions
"""

# ╔═╡ cfea06e8-f022-48cf-9ab4-b9f622193125
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ 9b3f9d9d-ef24-428f-8c2c-f40316974478
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ fc193d1a-a66c-4dca-ac63-deab0d28128a
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

# ╔═╡ f072e933-fa0e-4fd8-b180-50bfe151e4cd
md"""
## Segment Heart
"""

# ╔═╡ 28d931ac-190a-45ac-a4d4-fb16ac3f6a0e
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ 0b60202c-f0a2-4a94-bd50-2ee55d58e3fa
# @bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ a8cb7c21-f388-46bc-843d-d823b4f4c9b8
# heatmap(masked_array[:, :, a], colormap=:grays)

# ╔═╡ ed320a72-d5b9-4837-ae13-7b5745a90c00
# begin
# 	fig = Figure()
	
# 	ax = Makie.Axis(fig[1, 1])
# 	ax.title = "Raw DICOM Array"
# 	heatmap!(transpose(dcm_array[:, :, 15]), colormap=:grays)
# 	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
# 	fig
# end

# ╔═╡ 7516886f-8b71-42d1-a6ce-bb510c08eaa3
# begin
# 	fig2 = Figure()
	
# 	ax2 = Makie.Axis(fig2[1, 1])
# 	ax2.title = "Mask Array"
# 	heatmap!(transpose(mask), colormap=:grays)
# 	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
# 	fig2
# end

# ╔═╡ a1c813fb-5a00-411d-a2cb-0950382e116e
# begin
# 	fig3 = Figure()
	
# 	ax3 = Makie.Axis(fig3[1, 1])
# 	ax3.title = "Masked DICOM Array"
# 	heatmap!(transpose(masked_array[:, :, 23]), colormap=:grays)
# 	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
# 	fig3
# end

# ╔═╡ 250606d1-95d5-4ca1-a632-b4e0f1369d03
md"""
## Segment Calcium Rod
"""

# ╔═╡ 7518de19-9556-4c9d-b0a3-36629118f8d3
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header);

# ╔═╡ d8a9ae7c-65ba-45ca-8a19-69caa28c5118
# @bind c PlutoUI.Slider(1:size(calcium_image, 3), default=cal_rod_slice, show_value=true)

# ╔═╡ 17aebaec-09f1-4db8-a976-69e4e507f809
# heatmap(transpose(calcium_image[:, :, c]), colormap=:grays)

# ╔═╡ c3dc7106-3c24-4c6e-b668-e0a35dc85db4
md"""
## Segment Calcium Inserts
"""

# ╔═╡ e52b6fdf-9e0c-42e4-9d16-d7a2260dcbb5
mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
            dcm_array, masked_array, header, slice_CCI, center_insert
);

# ╔═╡ 0339c27f-ad59-4db0-a383-92440f2d49ee
# masks = mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD;

# ╔═╡ 2fc1d0ff-58bc-4dbf-af7e-a38585b1359d
# heatmap(masks, colormap=:grays)

# ╔═╡ e0298bd6-1501-4969-a357-f37760df5094
md"""
## Calibration Prep
"""

# ╔═╡ 53be7f97-4ca2-4af9-8f21-c0a63750aa1f
array_filtered = abs.(mapwindow(median, calcium_image[:, :, cal_rod_slice], (3, 3)));

# ╔═╡ 53534629-e632-4b9b-8ff7-909a112080b1
bool_arr = array_filtered .> 0;

# ╔═╡ ed979a13-e65c-4d34-8d8a-7ac5537b9af3
bool_arr_erode = (((erode(erode(bool_arr)))));

# ╔═╡ 747a1c2f-0786-4852-8535-1b62a874d55d
heatmap(bool_arr, colormap=:grays)

# ╔═╡ c363f4b1-2419-47f6-a182-f0c5c100e89f
heatmap(bool_arr_erode, colormap=:grays)

# ╔═╡ 929ee2ef-f5bf-4a5b-97c2-78da4f0e7a53
c_img = calcium_image[:, :, cal_rod_slice-1:cal_rod_slice+1];

# ╔═╡ 0554bef9-d8ba-4e65-aa93-7c874f8da732
begin
	mask_cal_3D = Array{Bool}(undef, size(c_img))
	for z in 1:size(c_img, 3)
		mask_cal_3D[:, :, z] = bool_arr_erode
	end
end;

# ╔═╡ 2d1bf37c-7e0b-4f5c-8fbb-58c87837e8fd
hist(c_img[mask_cal_3D])

# ╔═╡ 38a85949-7efa-491b-a6b8-179fd4d57ae2
# cal_insert_mean = mean(c_img[mask_cal_3D])

# ╔═╡ f93b94fa-d63e-404f-93d3-4a18ab31a91d
cal_insert_mean = quantile!(c_img[mask_cal_3D], 0.7)

# ╔═╡ bb5ec211-871f-4bb1-ac40-c23d33fc079d
md"""
# Score Large Inserts
"""

# ╔═╡ cb4db0fb-b87a-4d81-b140-85ca04986f55
arr = masked_array[:, :, slice_CCI-2:slice_CCI+2];

# ╔═╡ 235ce1d6-0860-47e7-87ce-ccb6884eadab
single_arr = masked_array[:, :, slice_CCI];

# ╔═╡ 45e5b78a-5b35-49c6-a94c-c0e5ad9408e1
md"""
## High Density
"""

# ╔═╡ f15d5916-be0b-4553-9508-065a40874182
begin
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_HD_3D[:, :, z] = mask_L_HD
	end
end;

# ╔═╡ b7ff9029-94d8-4018-ab29-69029e04c7dc
md"""
#### Dilated mask
"""

# ╔═╡ 792c5706-5431-42fa-a30b-7bb09516e85c
dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D));

# ╔═╡ 5d5baf0a-4deb-46d8-99ec-eb5f7200529f
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ 24734b8d-547b-4e0d-9804-e2f57f91086b
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ bed3e046-5d6a-4b2d-b486-017bcd96aa53
md"""
#### Ring (background) mask
"""

# ╔═╡ 3a96e981-07e7-48ea-9489-ab784206d339
ring_mask_L_HD = dilate(dilate(dilate(dilate(mask_L_HD_3D)))) - dilate(dilate(dilate(mask_L_HD_3D)));

# ╔═╡ a944304f-92c8-451b-86a0-24a668e90ee6
@bind g4 overlay_mask_bind(ring_mask_L_HD)

# ╔═╡ 1e52a3fc-44d8-4502-a962-09a61242f453
overlay_mask_plot(arr, ring_mask_L_HD, g4, "ring mask")

# ╔═╡ 4017abfe-fc15-4bfb-af62-2b14863eeb93
md"""
### Calculations
"""

# ╔═╡ 2210245f-8422-4eeb-9ff4-8f0b979e3db9
begin
	single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
	s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
end

# ╔═╡ b9e14fe9-bae8-4137-9857-4a6d085c6484
pixel_size = DICOMUtils.get_pixel_size(header)

# ╔═╡ ebb5f180-a086-4f0d-bd58-47756e4e5cf8
begin
	alg_L_HD = Integrated(arr[mask_L_HD_3D])
	ρ = 0.2 # mg/mm^3
	mass_l_hd = score(s_bkg_L_HD, cal_insert_mean, pixel_size, ρ, alg_L_HD)
end

# ╔═╡ 70588b66-5e0f-4886-9c51-6dd6b91dbfcd
md"""
## Medium Density
"""

# ╔═╡ acd856bb-6990-469b-b0ad-355fd2342b34
begin
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_MD_3D[:, :, z] = mask_L_MD
	end
end;

# ╔═╡ 9f8a2b11-f04a-49c1-a0ba-ff1d2c905e3e
md"""
#### Dilated mask
"""

# ╔═╡ ec783672-cf7b-4b37-bb20-e08b81055a69
dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D));

# ╔═╡ fa5cff5e-4109-4b28-9fe1-4523a2fc574e
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ 28c10aae-84b6-437e-91aa-aafc3b5ef4dc
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ 6b0c94c6-905b-4919-a56b-3287525f515d
md"""
#### Ring (background) mask
"""

# ╔═╡ 39f14235-2db6-4e61-9385-2201da2a8001
ring_mask_L_MD = dilate(dilate(dilate(dilate(mask_L_MD_3D)))) - dilate(dilate(dilate(mask_L_MD_3D)));

# ╔═╡ 70202236-ea9e-4997-973f-f814e1e86ea0
@bind h4 overlay_mask_bind(ring_mask_L_MD)

# ╔═╡ dac72065-1e5e-4490-9f3c-3bca88b994e2
overlay_mask_plot(arr, ring_mask_L_MD, h4, "ring mask")

# ╔═╡ c24928c8-aab9-4732-ab5b-651a1d104625
md"""
### Calculations
"""

# ╔═╡ 7c562ede-2091-4719-8e1a-baae16a3ba24
begin
	single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
	s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
end

# ╔═╡ 044a9b05-913c-4bc4-8c7b-cea9d90070cb
begin
	alg_L_MD = Integrated(arr[mask_L_MD_3D])
	mass_l_md = score(s_bkg_L_MD, cal_insert_mean, pixel_size, ρ, alg_L_MD)
end

# ╔═╡ 9d73f08c-8a8f-45de-946e-61ed7fb596a3
md"""
## Low Density
"""

# ╔═╡ fc02ab6b-4ea1-4b9d-9d14-61a8e3c36953
begin
	mask_L_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_LD_3D[:, :, z] = mask_L_LD
	end
end;

# ╔═╡ 6521aacf-09d6-49b1-a409-5baa89ca489c
md"""
#### Dilated mask
"""

# ╔═╡ cd3375a9-2fa4-4ed9-b356-5e2a4ec04d45
dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D));

# ╔═╡ c06f2713-5d83-4f80-b92c-495c4f8d76bf
@bind i2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ a622973c-6fff-4b1c-8e1b-7b47aff454be
overlay_mask_plot(arr, dilated_mask_L_LD, i2, "dilated mask")

# ╔═╡ ec5b2ddc-db27-4d50-bf3d-902ed17ffcaa
md"""
#### Ring (background) mask
"""

# ╔═╡ 705e21b0-0e81-4907-8dcb-0cbcdb6e1810
ring_mask_L_LD = dilate(dilate(dilate(dilate(mask_L_LD_3D)))) - dilate(dilate(dilate(mask_L_LD_3D)));

# ╔═╡ 0a6f9321-51a9-4497-a978-d28634268f7f
@bind i4 overlay_mask_bind(ring_mask_L_LD)

# ╔═╡ c9011d84-cdcd-4d06-ba06-6a158a56425f
overlay_mask_plot(arr, ring_mask_L_LD, i4, "ring mask")

# ╔═╡ bca1eb9e-81bf-4647-9198-c37518d1d287
md"""
### Calculations
"""

# ╔═╡ bb88ef2f-1d30-461c-9e52-a5db8c2b677d
begin	
	single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
	s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
end

# ╔═╡ 25315456-d95f-431a-9928-9d0c9d3f945a
begin
	alg_L_LD = Integrated(arr[mask_L_LD_3D])
	mass_l_ld = score(s_bkg_L_LD, cal_insert_mean, pixel_size, ρ, alg_L_LD)
end

# ╔═╡ b3bd32a6-6493-4f2f-becb-7940a07250e1
md"""
# Score Medium Inserts
"""

# ╔═╡ f416a685-688e-4498-a2bc-a36db3636de1
md"""
## High Density
"""

# ╔═╡ c8e0778d-7dd0-4092-b9ff-9f0123a0314b
begin
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_HD_3D[:, :, z] = mask_M_HD
	end
end;

# ╔═╡ c7d201f9-1c01-4023-bec1-982b25501112
md"""
#### Dilated mask
"""

# ╔═╡ 96515861-2923-4cd4-bc3b-d1298dc4e0fd
dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 5641a223-34e1-4075-88b7-3d04906e2813
@bind j2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ a5d65fd4-64e0-4b8c-8d78-b56e3e09d366
overlay_mask_plot(arr, dilated_mask_M_HD, j2, "dilated mask")

# ╔═╡ c2b4c8db-4d5e-48a0-a851-8150773100f1
md"""
#### Ring (background) mask
"""

# ╔═╡ 2a514a2a-f125-4921-8604-e72359077237
ring_mask_M_HD = dilate(dilate(dilate(dilate(dilate(mask_M_HD_3D))))) - dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ af10e1a1-d81f-4fb4-9c18-6e68ad125b46
@bind j4 overlay_mask_bind(ring_mask_M_HD)

# ╔═╡ 50da14d8-b584-47b7-811b-661666345ca1
overlay_mask_plot(arr, ring_mask_M_HD, j4, "ring mask")

# ╔═╡ 62021493-9443-4d36-a7b3-0e12b1c0a3d8
md"""
### Calculations
"""

# ╔═╡ 3aca33f9-8cbf-46b8-87e7-3e502838d947
begin
	single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
	s_bkg_M_HD = mean(single_arr[single_ring_mask_M_HD])
end

# ╔═╡ 01228d7c-525f-4824-b618-3f3ea40d2db3
begin
	alg_M_HD = Integrated(arr[mask_M_HD_3D])
	mass_m_hd = score(s_bkg_M_HD, cal_insert_mean, pixel_size, ρ, alg_M_HD)
end

# ╔═╡ b7203954-eeaf-4e2c-a203-b6f77ce087db
md"""
## Medium Density
"""

# ╔═╡ 5f2f4697-a988-4326-bb17-0f7ad110b2dd
begin
	mask_M_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_MD_3D[:, :, z] = mask_M_MD
	end
end;

# ╔═╡ 82d686b7-21d4-436d-adb9-a4efc1e61fb2
md"""
#### Dilated mask
"""

# ╔═╡ 37cdb93a-0d75-4a54-8981-30e710a10448
dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))));

# ╔═╡ 0fe83072-7bf2-49be-9a66-0d7bb7789b16
@bind k2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ 87cceb3a-5fd6-4c22-9319-9a294f61c9bc
overlay_mask_plot(arr, dilated_mask_M_MD, k2, "dilated mask")

# ╔═╡ 198e185e-dd12-49da-83b8-cb235c61e39a
md"""
#### Ring (background) mask
"""

# ╔═╡ 5d4e4afc-edca-486c-8b27-0f0a8a369558
ring_mask_M_MD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))));

# ╔═╡ d7d61d7e-4dc5-4c97-b420-7b6acb712960
@bind k4 overlay_mask_bind(ring_mask_M_MD)

# ╔═╡ edf2edc0-03b8-4303-9291-9bd66daadc76
overlay_mask_plot(arr, ring_mask_M_MD, k4, "ring mask")

# ╔═╡ ef2801d6-1060-4d00-90d5-095a7d49bde7
md"""
### Calculations
"""

# ╔═╡ 71c54468-d8c3-46c4-9df3-e705f5ef407d
begin
	single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
	s_bkg_M_MD = mean(single_arr[single_ring_mask_M_MD])
end

# ╔═╡ 0bf183d7-0f55-430e-8def-c4c5f4a2e3cf
begin
	alg_M_MD = Integrated(arr[mask_M_MD_3D])
	mass_m_md = score(s_bkg_M_MD, cal_insert_mean, pixel_size, ρ, alg_M_MD)
end

# ╔═╡ 34ba311a-1419-445c-a321-24e52e5a1eca
md"""
## Low Density
"""

# ╔═╡ 1a60a5d0-6f3f-4c42-a8bc-d0b009a1f3a6
begin
	mask_M_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_LD_3D[:, :, z] = mask_M_LD
	end
end;

# ╔═╡ f8f712bb-e6cb-4c1d-a8d3-89ceb88f0ecf
md"""
#### Dilated mask
"""

# ╔═╡ a3514e20-2588-44a4-85ff-5804db9afcdb
dilated_mask_M_LD = dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ 3d8e308a-bb4a-4a50-a9b9-65918d3d32bb
@bind l2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ cf9dfcb3-f7f1-4ddc-8904-b12840e85ecc
overlay_mask_plot(arr, dilated_mask_M_LD, l2, "dilated mask")

# ╔═╡ e78527dd-aef9-422f-9d1d-450f454ccca9
md"""
#### Ring (background) mask
"""

# ╔═╡ ebdcb914-5ea2-4c09-b831-e4d5c69cd7a1
ring_mask_M_LD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ 4ea4823a-e80c-41ba-be7d-6a437f72d7c3
@bind l4 overlay_mask_bind(ring_mask_M_LD)

# ╔═╡ ebef33ce-e26d-45cc-b524-3a6a2897196d
overlay_mask_plot(arr, ring_mask_M_LD, l4, "ring mask")

# ╔═╡ e7c82858-d19a-4c8f-b04a-7c13e1a5734a
md"""
### Calculations
"""

# ╔═╡ 401a2072-502f-4d79-b292-26e896331662
begin
	single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
	s_bkg_M_LD = mean(single_arr[single_ring_mask_M_LD])
end

# ╔═╡ 8e51515d-387e-462f-89bf-5f015c90ef32
begin
	alg_M_LD = Integrated(arr[mask_M_LD_3D])
	mass_m_ld = score(s_bkg_M_LD, cal_insert_mean, pixel_size, ρ, alg_M_LD)
end

# ╔═╡ 96f2567b-f67e-45fb-a4c5-74c2192a991a
md"""
# Score Small Inserts
"""

# ╔═╡ f44b12da-4d72-4e70-8b65-b7d3a1962d09
md"""
## High Density
"""

# ╔═╡ b62dd477-3974-45f4-9994-7dce9ef79be4
begin
	mask_S_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_HD_3D[:, :, z] = mask_S_HD
	end
end;

# ╔═╡ b71a6daf-00d4-4e5d-9770-2b220720a106
md"""
#### Dilated mask
"""

# ╔═╡ 118b69e5-3491-4a8a-a86d-8b94f37eb009
dilated_mask_S_HD = dilate(dilate(dilate(dilate(dilate((mask_S_HD_3D))))));

# ╔═╡ 0085775f-2ea1-491f-a052-69dd45e53c49
@bind m2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ cd09e8ab-0cb0-4995-9664-d46dd1f32e20
overlay_mask_plot(arr, dilated_mask_S_HD, m2, "dilated mask")

# ╔═╡ 99b18e46-466c-40a8-8bc8-882540d96f45
md"""
#### Ring (background) mask
"""

# ╔═╡ d243a2f9-e74a-420b-8049-b275c24b389f
ring_mask_S_HD = dilate(dilate(dilate(dilate(dilate(mask_S_HD_3D))))) - dilate(dilate(dilate(dilate(mask_S_HD_3D))));

# ╔═╡ 9cff778a-c5d5-45bc-a111-ae7e8a8aee7b
@bind m4 overlay_mask_bind(ring_mask_S_HD)

# ╔═╡ 4d362c5d-5c58-412a-b3ed-8c9d4f5eff9a
overlay_mask_plot(arr, ring_mask_S_HD, m4, "ring mask")

# ╔═╡ 489f72ac-9a5c-4cbe-b331-df3ddc5e0f9c
md"""
### Calculations
"""

# ╔═╡ 889137b5-2cce-4f9b-98af-0eccd1fe25b8
begin
	single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
	s_bkg_S_HD = mean(single_arr[single_ring_mask_S_HD])
end

# ╔═╡ 7bc707fb-593b-421e-90af-1ec1b406ddd0
begin
	alg_S_HD = Integrated(arr[mask_S_HD_3D])
	mass_s_hd = score(s_bkg_S_HD, cal_insert_mean, pixel_size, ρ, alg_S_HD)
	if mass_s_hd < 0
		mass_s_hd = 0
	end
	mass_s_hd
end

# ╔═╡ db58ff82-d605-4361-aaac-209f79afbcb3
md"""
## Medium Density
"""

# ╔═╡ 8b4e151b-808a-4e2e-93af-66779c563113
begin
	mask_S_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_MD_3D[:, :, z] = mask_S_MD
	end
end;

# ╔═╡ 4db69ff4-78b6-4838-af94-f3dd1cd86d5e
md"""
#### Dilated mask
"""

# ╔═╡ dbb38752-db73-4051-9ca2-d2eafa705a33
dilated_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D)))));

# ╔═╡ dab3cf21-1f58-4bb3-a07f-af508d27f798
@bind n2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ e052a45f-3f04-4180-b5f4-be97eedcffb3
overlay_mask_plot(arr, dilated_mask_S_MD, n2, "dilated mask")

# ╔═╡ 31cb4900-be86-4f53-ae77-504e7fbd64df
md"""
#### Ring (background) mask
"""

# ╔═╡ 9df7465b-c59a-4128-adb5-99923e74c408
ring_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D))))) - dilate(dilate(dilate(dilate(mask_S_MD_3D))));

# ╔═╡ 4e461d7e-16e0-4464-be58-8aca2bcf337b
@bind n4 overlay_mask_bind(ring_mask_S_MD)

# ╔═╡ 48a89dfc-e522-4969-b3d0-9490c8c0b55d
overlay_mask_plot(arr, ring_mask_S_MD, n4, "ring mask")

# ╔═╡ 43de5ce4-b710-4343-a0d9-acae9e2d8e7d
md"""
### Calculations
"""

# ╔═╡ b10dc7d8-15ef-441d-b62d-c8f81587762d
begin
	single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
	s_bkg_S_MD = mean(single_arr[single_ring_mask_S_MD])
end

# ╔═╡ 4285f70f-66c5-4a32-ac2b-e5bfd07dc4fe
begin
	alg_S_MD = Integrated(arr[mask_S_MD_3D])
	mass_s_md = score(s_bkg_S_MD, cal_insert_mean, pixel_size, ρ, alg_S_MD)
	if mass_s_md < 0
		mass_s_md = 0
	end
	mass_s_md
end

# ╔═╡ 6241d1d5-939e-447e-83a7-9588d5623cd2
md"""
## Low Density
"""

# ╔═╡ af2a614a-94ca-4e6a-acad-a536a3492d85
begin
	mask_S_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_LD_3D[:, :, z] = mask_S_LD
	end
end;

# ╔═╡ 6c5d42d8-629a-4666-b6a4-1cc79e5af861
md"""
#### Dilated mask
"""

# ╔═╡ db62af76-da9b-4738-88fe-a532d8a8e9e0
dilated_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D)))));

# ╔═╡ b31c8d79-762b-4e50-8799-15876c0c3628
@bind o2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ 8dce9b3c-de8c-49df-8eb6-ad65df5cd168
overlay_mask_plot(arr, dilated_mask_S_LD, o2, "dilated mask")

# ╔═╡ 4b2eb2b4-a921-4731-9892-b6f7d2e4d955
md"""
#### Ring (background) mask
"""

# ╔═╡ 4d6984ad-3d1c-4270-a3fa-c2ba7738580a
ring_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D))))) - dilate(dilate(dilate(dilate(mask_S_LD_3D))));

# ╔═╡ e4259e8b-6ff2-420c-884d-03c6e6942589
@bind o4 overlay_mask_bind(ring_mask_S_LD)

# ╔═╡ 99e2d268-6b41-45c1-ac2d-e8cfc8ae27c1
overlay_mask_plot(arr, ring_mask_S_LD, o4, "ring mask")

# ╔═╡ 546405df-c777-41bb-aa1d-decb43778d65
md"""
### Calculations
"""

# ╔═╡ 80a75352-19df-410d-b5ed-af99690d9f3b
begin
	single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
	s_bkg_S_LD = mean(single_arr[single_ring_mask_S_LD])
end

# ╔═╡ 7fe9985a-e942-4a08-8e75-396219264588
begin
	alg_S_LD = Integrated(arr[mask_S_LD_3D])
	mass_s_ld = score(s_bkg_S_LD, cal_insert_mean, pixel_size, ρ, alg_S_LD)
	if mass_s_ld < 0
		mass_s_ld = 0
	end
	mass_s_ld
end

# ╔═╡ 0ac56156-0ebf-4c13-8354-34f6f572bb09
md"""
# Results
"""

# ╔═╡ bbb64dcc-568f-40e1-a3c5-0eaf41ec2b21
density_array = [0, 200, 400, 800]

# ╔═╡ d29d31a6-7160-46f0-8e9a-f5b51d1ce068
inserts = [
	"Low Density",
	"Medium Density",
	"High Density"
]

# ╔═╡ 0605b308-a420-43e4-9901-7e8600d966e6
ground_truth_mass_large = [
	19.6,
	39.3,
	78.5
] # mg

# ╔═╡ 21ba3aae-e956-4f88-8cc2-43bf08dc4cb9
calculated_mass_large = [
	mass_l_ld,
	mass_l_md,
	mass_l_hd
]

# ╔═╡ f1c2c491-7c50-4c5f-85a2-6e5756fca7e5
ground_truth_mass_medium = [
	4.2,
	8.5,
	17.0
]

# ╔═╡ a8622859-82cd-4570-b0f2-246ba56bf7c2
calculated_mass_medium = [
	mass_m_ld,
	mass_m_md,
	mass_m_hd
]

# ╔═╡ c43c83f0-885a-4e46-a2df-8432ee5a2e3e
ground_truth_mass_small = [
	0.2,
	0.3,
	0.6
]

# ╔═╡ c8b01534-801f-4eaf-98da-1a9ba9f27986
calculated_mass_small = [
	mass_s_ld,
	mass_s_md,
	mass_s_hd
]

# ╔═╡ 3098a271-eb06-403e-943e-464a1208ee7d
df = DataFrame(
	inserts = inserts,
	ground_truth_mass_large = ground_truth_mass_large,
	calculated_mass_large = calculated_mass_large,
	ground_truth_mass_medium = ground_truth_mass_medium,
	calculated_mass_medium = calculated_mass_medium,
	ground_truth_mass_small = ground_truth_mass_small,
	calculated_mass_small = calculated_mass_small
)

# ╔═╡ 9f092358-c67d-4b1a-845a-fd83fb40f8d9
begin
	fmass2 = Figure()
	axmass2 = Axis(fmass2[1, 1])
	
	scatter!(density_array[2:end], df[!, :ground_truth_mass_large], label="ground_truth_mass_large")
	scatter!(density_array[2:end], df[!, :calculated_mass_large], label="calculated_mass_large")
	
	axmass2.title = "Mass Measurements (Large)"
	axmass2.ylabel = "Mass (mg)"
	axmass2.xlabel = "Density (mg/cm^3)"

	xlims!(axmass2, 0, 850)
	ylims!(axmass2, 0, 100)
	
	fmass2[1, 2] = Legend(fmass2, axmass2, framevisible = false)
	
	fmass2
end

# ╔═╡ 411f7e69-1431-47d4-b8a9-badd507dbb86
begin
	fmass3 = Figure()
	axmass3 = Axis(fmass3[1, 1])
	
	scatter!(density_array[2:end], df[!, :ground_truth_mass_medium], label="ground_truth_mass_medium")
	scatter!(density_array[2:end], df[!, :calculated_mass_medium], label="calculated_mass_medium")
	
	axmass3.title = "Mass Measurements (Medium)"
	axmass3.ylabel = "Mass (mg)"
	axmass3.xlabel = "Density (mg/cm^3)"

	xlims!(axmass3, 0, 850)
	ylims!(axmass3, 0, 25)
	
	fmass3[1, 2] = Legend(fmass3, axmass3, framevisible = false)
	
	fmass3
end

# ╔═╡ 6a2273ce-3066-407e-83d2-cda8db84bcb6
begin
	fmass4 = Figure()
	axmass4 = Axis(fmass4[1, 1])
	
	scatter!(density_array[2:end], df[!, :ground_truth_mass_small], label="ground_truth_mass_small")
	scatter!(density_array[2:end], df[!, :calculated_mass_small], label="calculated_mass_small")
	
	axmass4.title = "Mass Measurements (Small)"
	axmass4.ylabel = "Mass (mg)"
	axmass4.xlabel = "Density (mg/cm^3)"

	xlims!(axmass4, 0, 850)
	ylims!(axmass4, 0, 1.5)
	
	fmass4[1, 2] = Legend(fmass4, axmass4, framevisible = false)
	
	fmass4
end

# ╔═╡ 85de788e-4ec8-41cb-832c-5c29204ec5d9
md"""
### Save Results
"""

# ╔═╡ b1ec58d3-7646-4b83-bdfa-4c6002e7fa26
if ~isdir(string(cd(pwd, "..") , "/data/output/", VENDER))
	mkdir(string(cd(pwd, "..") , "/data/output/", VENDER))
end

# ╔═╡ 98ff7e44-7457-4e76-929a-5c9b2cb26ad3
output_path = string(cd(pwd, "..") , "/data/output/", VENDER, "/", scan, ".csv")

# ╔═╡ ffbc2df4-4d26-4863-af3e-23e8ee438442
CSV.write(output_path, df)

# ╔═╡ Cell order:
# ╠═52d4c3cb-c0a6-497a-af42-67eb5372cf21
# ╠═9bffa403-b54c-49f4-b4e4-3ec593c3a400
# ╟─afa84ebb-9e0d-4f0a-b980-958b0dbddd1a
# ╠═857ea37f-58dd-45d4-96e5-8238718d9a3a
# ╟─aa2e4340-9eb8-42c7-9765-e27859841fa7
# ╠═08c9aa2d-77c1-45c3-bbe9-889516fbdf7d
# ╠═2b5dfe45-acb2-48ec-b2ec-9b08ad573973
# ╠═a770154a-491b-4e56-95c6-a1b1378eec7b
# ╠═9dca8618-5a6a-4cae-91a3-bf005de45d13
# ╠═6daad570-d417-49b1-84b2-dc78cf4ed752
# ╠═aceb0699-06e1-4c6c-a178-0aaa47683930
# ╟─d1705d56-40ab-4119-adb7-876126981d8b
# ╟─cfea06e8-f022-48cf-9ab4-b9f622193125
# ╟─9b3f9d9d-ef24-428f-8c2c-f40316974478
# ╟─fc193d1a-a66c-4dca-ac63-deab0d28128a
# ╟─f072e933-fa0e-4fd8-b180-50bfe151e4cd
# ╠═28d931ac-190a-45ac-a4d4-fb16ac3f6a0e
# ╠═0b60202c-f0a2-4a94-bd50-2ee55d58e3fa
# ╠═a8cb7c21-f388-46bc-843d-d823b4f4c9b8
# ╠═ed320a72-d5b9-4837-ae13-7b5745a90c00
# ╠═7516886f-8b71-42d1-a6ce-bb510c08eaa3
# ╠═a1c813fb-5a00-411d-a2cb-0950382e116e
# ╟─250606d1-95d5-4ca1-a632-b4e0f1369d03
# ╠═7518de19-9556-4c9d-b0a3-36629118f8d3
# ╠═d8a9ae7c-65ba-45ca-8a19-69caa28c5118
# ╠═17aebaec-09f1-4db8-a976-69e4e507f809
# ╟─c3dc7106-3c24-4c6e-b668-e0a35dc85db4
# ╠═e52b6fdf-9e0c-42e4-9d16-d7a2260dcbb5
# ╠═0339c27f-ad59-4db0-a383-92440f2d49ee
# ╠═2fc1d0ff-58bc-4dbf-af7e-a38585b1359d
# ╟─e0298bd6-1501-4969-a357-f37760df5094
# ╠═53be7f97-4ca2-4af9-8f21-c0a63750aa1f
# ╠═53534629-e632-4b9b-8ff7-909a112080b1
# ╠═ed979a13-e65c-4d34-8d8a-7ac5537b9af3
# ╠═747a1c2f-0786-4852-8535-1b62a874d55d
# ╠═c363f4b1-2419-47f6-a182-f0c5c100e89f
# ╠═0554bef9-d8ba-4e65-aa93-7c874f8da732
# ╠═929ee2ef-f5bf-4a5b-97c2-78da4f0e7a53
# ╠═2d1bf37c-7e0b-4f5c-8fbb-58c87837e8fd
# ╠═38a85949-7efa-491b-a6b8-179fd4d57ae2
# ╠═f93b94fa-d63e-404f-93d3-4a18ab31a91d
# ╟─bb5ec211-871f-4bb1-ac40-c23d33fc079d
# ╠═cb4db0fb-b87a-4d81-b140-85ca04986f55
# ╠═235ce1d6-0860-47e7-87ce-ccb6884eadab
# ╟─45e5b78a-5b35-49c6-a94c-c0e5ad9408e1
# ╠═f15d5916-be0b-4553-9508-065a40874182
# ╟─b7ff9029-94d8-4018-ab29-69029e04c7dc
# ╠═792c5706-5431-42fa-a30b-7bb09516e85c
# ╟─5d5baf0a-4deb-46d8-99ec-eb5f7200529f
# ╠═24734b8d-547b-4e0d-9804-e2f57f91086b
# ╟─bed3e046-5d6a-4b2d-b486-017bcd96aa53
# ╠═3a96e981-07e7-48ea-9489-ab784206d339
# ╟─a944304f-92c8-451b-86a0-24a668e90ee6
# ╠═1e52a3fc-44d8-4502-a962-09a61242f453
# ╟─4017abfe-fc15-4bfb-af62-2b14863eeb93
# ╠═2210245f-8422-4eeb-9ff4-8f0b979e3db9
# ╠═b9e14fe9-bae8-4137-9857-4a6d085c6484
# ╠═ebb5f180-a086-4f0d-bd58-47756e4e5cf8
# ╟─70588b66-5e0f-4886-9c51-6dd6b91dbfcd
# ╠═acd856bb-6990-469b-b0ad-355fd2342b34
# ╟─9f8a2b11-f04a-49c1-a0ba-ff1d2c905e3e
# ╠═ec783672-cf7b-4b37-bb20-e08b81055a69
# ╟─fa5cff5e-4109-4b28-9fe1-4523a2fc574e
# ╠═28c10aae-84b6-437e-91aa-aafc3b5ef4dc
# ╟─6b0c94c6-905b-4919-a56b-3287525f515d
# ╠═39f14235-2db6-4e61-9385-2201da2a8001
# ╟─70202236-ea9e-4997-973f-f814e1e86ea0
# ╠═dac72065-1e5e-4490-9f3c-3bca88b994e2
# ╟─c24928c8-aab9-4732-ab5b-651a1d104625
# ╠═7c562ede-2091-4719-8e1a-baae16a3ba24
# ╠═044a9b05-913c-4bc4-8c7b-cea9d90070cb
# ╟─9d73f08c-8a8f-45de-946e-61ed7fb596a3
# ╠═fc02ab6b-4ea1-4b9d-9d14-61a8e3c36953
# ╟─6521aacf-09d6-49b1-a409-5baa89ca489c
# ╠═cd3375a9-2fa4-4ed9-b356-5e2a4ec04d45
# ╟─c06f2713-5d83-4f80-b92c-495c4f8d76bf
# ╠═a622973c-6fff-4b1c-8e1b-7b47aff454be
# ╟─ec5b2ddc-db27-4d50-bf3d-902ed17ffcaa
# ╠═705e21b0-0e81-4907-8dcb-0cbcdb6e1810
# ╟─0a6f9321-51a9-4497-a978-d28634268f7f
# ╠═c9011d84-cdcd-4d06-ba06-6a158a56425f
# ╟─bca1eb9e-81bf-4647-9198-c37518d1d287
# ╠═bb88ef2f-1d30-461c-9e52-a5db8c2b677d
# ╠═25315456-d95f-431a-9928-9d0c9d3f945a
# ╟─b3bd32a6-6493-4f2f-becb-7940a07250e1
# ╟─f416a685-688e-4498-a2bc-a36db3636de1
# ╠═c8e0778d-7dd0-4092-b9ff-9f0123a0314b
# ╟─c7d201f9-1c01-4023-bec1-982b25501112
# ╠═96515861-2923-4cd4-bc3b-d1298dc4e0fd
# ╟─5641a223-34e1-4075-88b7-3d04906e2813
# ╠═a5d65fd4-64e0-4b8c-8d78-b56e3e09d366
# ╟─c2b4c8db-4d5e-48a0-a851-8150773100f1
# ╠═2a514a2a-f125-4921-8604-e72359077237
# ╟─af10e1a1-d81f-4fb4-9c18-6e68ad125b46
# ╠═50da14d8-b584-47b7-811b-661666345ca1
# ╟─62021493-9443-4d36-a7b3-0e12b1c0a3d8
# ╠═3aca33f9-8cbf-46b8-87e7-3e502838d947
# ╠═01228d7c-525f-4824-b618-3f3ea40d2db3
# ╟─b7203954-eeaf-4e2c-a203-b6f77ce087db
# ╠═5f2f4697-a988-4326-bb17-0f7ad110b2dd
# ╟─82d686b7-21d4-436d-adb9-a4efc1e61fb2
# ╠═37cdb93a-0d75-4a54-8981-30e710a10448
# ╟─0fe83072-7bf2-49be-9a66-0d7bb7789b16
# ╠═87cceb3a-5fd6-4c22-9319-9a294f61c9bc
# ╟─198e185e-dd12-49da-83b8-cb235c61e39a
# ╠═5d4e4afc-edca-486c-8b27-0f0a8a369558
# ╟─d7d61d7e-4dc5-4c97-b420-7b6acb712960
# ╠═edf2edc0-03b8-4303-9291-9bd66daadc76
# ╟─ef2801d6-1060-4d00-90d5-095a7d49bde7
# ╠═71c54468-d8c3-46c4-9df3-e705f5ef407d
# ╠═0bf183d7-0f55-430e-8def-c4c5f4a2e3cf
# ╟─34ba311a-1419-445c-a321-24e52e5a1eca
# ╠═1a60a5d0-6f3f-4c42-a8bc-d0b009a1f3a6
# ╟─f8f712bb-e6cb-4c1d-a8d3-89ceb88f0ecf
# ╠═a3514e20-2588-44a4-85ff-5804db9afcdb
# ╟─3d8e308a-bb4a-4a50-a9b9-65918d3d32bb
# ╠═cf9dfcb3-f7f1-4ddc-8904-b12840e85ecc
# ╟─e78527dd-aef9-422f-9d1d-450f454ccca9
# ╠═ebdcb914-5ea2-4c09-b831-e4d5c69cd7a1
# ╟─4ea4823a-e80c-41ba-be7d-6a437f72d7c3
# ╠═ebef33ce-e26d-45cc-b524-3a6a2897196d
# ╟─e7c82858-d19a-4c8f-b04a-7c13e1a5734a
# ╠═401a2072-502f-4d79-b292-26e896331662
# ╠═8e51515d-387e-462f-89bf-5f015c90ef32
# ╟─96f2567b-f67e-45fb-a4c5-74c2192a991a
# ╟─f44b12da-4d72-4e70-8b65-b7d3a1962d09
# ╠═b62dd477-3974-45f4-9994-7dce9ef79be4
# ╟─b71a6daf-00d4-4e5d-9770-2b220720a106
# ╠═118b69e5-3491-4a8a-a86d-8b94f37eb009
# ╟─0085775f-2ea1-491f-a052-69dd45e53c49
# ╠═cd09e8ab-0cb0-4995-9664-d46dd1f32e20
# ╟─99b18e46-466c-40a8-8bc8-882540d96f45
# ╠═d243a2f9-e74a-420b-8049-b275c24b389f
# ╟─9cff778a-c5d5-45bc-a111-ae7e8a8aee7b
# ╠═4d362c5d-5c58-412a-b3ed-8c9d4f5eff9a
# ╟─489f72ac-9a5c-4cbe-b331-df3ddc5e0f9c
# ╠═889137b5-2cce-4f9b-98af-0eccd1fe25b8
# ╠═7bc707fb-593b-421e-90af-1ec1b406ddd0
# ╟─db58ff82-d605-4361-aaac-209f79afbcb3
# ╠═8b4e151b-808a-4e2e-93af-66779c563113
# ╟─4db69ff4-78b6-4838-af94-f3dd1cd86d5e
# ╠═dbb38752-db73-4051-9ca2-d2eafa705a33
# ╟─dab3cf21-1f58-4bb3-a07f-af508d27f798
# ╠═e052a45f-3f04-4180-b5f4-be97eedcffb3
# ╟─31cb4900-be86-4f53-ae77-504e7fbd64df
# ╠═9df7465b-c59a-4128-adb5-99923e74c408
# ╟─4e461d7e-16e0-4464-be58-8aca2bcf337b
# ╠═48a89dfc-e522-4969-b3d0-9490c8c0b55d
# ╟─43de5ce4-b710-4343-a0d9-acae9e2d8e7d
# ╠═b10dc7d8-15ef-441d-b62d-c8f81587762d
# ╠═4285f70f-66c5-4a32-ac2b-e5bfd07dc4fe
# ╟─6241d1d5-939e-447e-83a7-9588d5623cd2
# ╠═af2a614a-94ca-4e6a-acad-a536a3492d85
# ╟─6c5d42d8-629a-4666-b6a4-1cc79e5af861
# ╠═db62af76-da9b-4738-88fe-a532d8a8e9e0
# ╟─b31c8d79-762b-4e50-8799-15876c0c3628
# ╠═8dce9b3c-de8c-49df-8eb6-ad65df5cd168
# ╟─4b2eb2b4-a921-4731-9892-b6f7d2e4d955
# ╠═4d6984ad-3d1c-4270-a3fa-c2ba7738580a
# ╟─e4259e8b-6ff2-420c-884d-03c6e6942589
# ╠═99e2d268-6b41-45c1-ac2d-e8cfc8ae27c1
# ╟─546405df-c777-41bb-aa1d-decb43778d65
# ╠═80a75352-19df-410d-b5ed-af99690d9f3b
# ╠═7fe9985a-e942-4a08-8e75-396219264588
# ╟─0ac56156-0ebf-4c13-8354-34f6f572bb09
# ╟─bbb64dcc-568f-40e1-a3c5-0eaf41ec2b21
# ╟─d29d31a6-7160-46f0-8e9a-f5b51d1ce068
# ╟─0605b308-a420-43e4-9901-7e8600d966e6
# ╟─21ba3aae-e956-4f88-8cc2-43bf08dc4cb9
# ╟─f1c2c491-7c50-4c5f-85a2-6e5756fca7e5
# ╟─a8622859-82cd-4570-b0f2-246ba56bf7c2
# ╟─c43c83f0-885a-4e46-a2df-8432ee5a2e3e
# ╟─c8b01534-801f-4eaf-98da-1a9ba9f27986
# ╠═3098a271-eb06-403e-943e-464a1208ee7d
# ╟─9f092358-c67d-4b1a-845a-fd83fb40f8d9
# ╟─411f7e69-1431-47d4-b8a9-badd507dbb86
# ╟─6a2273ce-3066-407e-83d2-cda8db84bcb6
# ╟─85de788e-4ec8-41cb-832c-5c29204ec5d9
# ╠═b1ec58d3-7646-4b83-bdfa-4c6002e7fa26
# ╠═98ff7e44-7457-4e76-929a-5c9b2cb26ad3
# ╠═ffbc2df4-4d26-4863-af3e-23e8ee438442

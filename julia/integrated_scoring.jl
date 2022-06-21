### A Pluto.jl notebook ###
# v0.19.8

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

# ╔═╡ 3c162d31-407a-4c50-a8a7-30153192f207
# ╠═╡ show_logs = false
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
		Pkg.add(url="https://github.com/Dale-Black/PhantomSegmentation.jl")
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
	using PhantomSegmentation
	using CalciumScoring
end

# ╔═╡ 0aeb13cf-6299-473a-908e-5b0be1b23c86
TableOfContents()

# ╔═╡ c4bba726-9f22-4f76-a4ae-966efd13afcf
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDER`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ 80ea14ec-0128-4899-bcae-5d75d6654a89
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ cb6af1f1-9e54-4172-9e65-bcf534840716
begin
	SCAN_NUMBER = 4
	VENDER = "Canon_Aquilion_One_Vision"
	BASE_PATH = "/Users/daleblack/Google Drive/Datasets/"
end

# ╔═╡ 47504ab4-a52f-418d-b9de-3468304c9578
root_path = string(BASE_PATH, VENDER)

# ╔═╡ d41f610c-928b-4a36-a7f8-ac50ff81f6db
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ 7efcfa5c-8d99-4e26-9538-efa1d2bfeea6
pth = dcm_path_list[SCAN_NUMBER]

# ╔═╡ 017121d9-a3bd-481a-ba89-7cd261f0ffe7
pth

# ╔═╡ 49fe69bc-7c19-4dc6-86df-ceaccbdd65f2
scan = basename(pth)

# ╔═╡ e95b1cbc-82a1-47db-8026-e2773f0d0cff
header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ a4f1ebf0-6f99-42fe-96f3-765caf5f207d
md"""
## Helper Functions
"""

# ╔═╡ 8e874d2a-95bc-414b-91cd-1e6df8220c18
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ 1fec01bc-445a-4961-b518-f25a0de737d3
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ 4494460e-4ffc-4317-9a30-21ac8e005bea
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
	scatter!(label_array[:, 1][indices_lbl], label_array[:, 2][indices_lbl], markersize=0.6, color=:red)
	fig
end

# ╔═╡ 9d4357e3-5ee5-4bf9-8879-473a4da59812
md"""
## Segment Heart
"""

# ╔═╡ 4a947165-c565-4c6d-b63c-8cc714d0d039
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ f1339241-169d-40c1-becc-8f230c12e835
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 641043a6-6d00-4f20-96ad-4b588bda2893
heatmap(masked_array[:, :, a], colormap=:grays)

# ╔═╡ ede1723d-ac8b-43a5-a6d8-6a92fcbaa2fa
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 25]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	fig
end

# ╔═╡ 13671de8-9397-422f-81a4-57781e881178
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	fig2
end

# ╔═╡ c07994fb-87c8-4a61-8dd9-7fe50b087b13
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 23]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	fig3
end

# ╔═╡ 2988d4aa-a2f2-45a5-9e96-fb9c233e9430
md"""
## Segment Calcium Rod
"""

# ╔═╡ 06c46a2f-9f03-46ec-88ed-6219bd2a5f8d
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header);

# ╔═╡ 66c50557-49d0-4d9f-beae-404b6784f079
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=cal_rod_slice, show_value=true)

# ╔═╡ 8806c9d6-c416-47be-981e-a03d35d45602
heatmap(transpose(calcium_image[:, :, c]), colormap=:grays)

# ╔═╡ f30cd88e-c6ca-495d-8ef9-6f029182b5a0
md"""
## Segment Calcium Inserts
"""

# ╔═╡ 7d67fa1f-74de-4a18-9b0d-eb46fa98f420
angle_factor = -4

# ╔═╡ 3bb994fa-e4e8-4e98-bf1b-5e896f21e96e
mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
		dcm_array, masked_array, header, slice_CCI, center_insert; angle_factor=angle_factor);

# ╔═╡ 0aece8c7-c8f7-49ca-b68e-b2fda92060c0
masks = mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD;

# ╔═╡ 5eb2cb0e-5353-4bd2-86ac-62f17be96f4c
heatmap(masks, colormap=:grays)

# ╔═╡ a488da37-5fd7-4d23-a958-d2fe3f1249b1
md"""
## Calibration Prep
"""

# ╔═╡ d7352bc6-2432-4926-bb70-c791a0050e6b
array_filtered = abs.(mapwindow(median, calcium_image[:, :, cal_rod_slice], (3, 3)));

# ╔═╡ e69626d2-a4eb-46f6-b3ca-491d427b7f70
bool_arr = array_filtered .> 0;

# ╔═╡ 88efbea6-138f-4885-9afb-543648e20152
bool_arr_erode = (((erode(erode(bool_arr)))));

# ╔═╡ a50ca0be-2613-498a-bcb8-5bdf10b4eeed
heatmap(bool_arr, colormap=:grays)

# ╔═╡ a2ea95bc-8840-4fbe-a29d-3618c629e9dc
heatmap(bool_arr_erode, colormap=:grays)

# ╔═╡ a574b69a-023e-4f1c-a09d-2f75e2eba8ec
md"""
### 1 Point Calibration
"""

# ╔═╡ 29a2117d-0462-4aac-8ff0-060f1c6ba848
c_img = calcium_image[:, :, cal_rod_slice-1:cal_rod_slice+1];

# ╔═╡ 4c5b5a15-a456-4f67-977d-55a727f6ac7f
begin
	mask_cal_3D = Array{Bool}(undef, size(c_img))
	for z in 1:size(c_img, 3)
		mask_cal_3D[:, :, z] = bool_arr_erode
	end
end;

# ╔═╡ 3f17d44f-59fb-4cc3-93d6-df025d92f1bf
cal_insert_mean = mean(c_img[mask_cal_3D])

# ╔═╡ 0362be1a-514b-4133-8c79-b2f554d7ed9f
hist(c_img[mask_cal_3D])

# ╔═╡ 952a3a84-ce2d-4628-8c84-a2fc8f017a96
begin
	intensity_array = [0, cal_insert_mean]
	density_array_cal = [0, 200]
	df_cal = DataFrame(:density => density_array_cal, :intensity => intensity_array)
	linearRegressor = lm(@formula(intensity ~ density), df_cal)
	linearFit = predict(linearRegressor)
	m = linearRegressor.model.pp.beta0[2]
	b = linearRegressor.model.rr.mu[1]
end

# ╔═╡ 67caa8b2-7fb7-49a0-b515-48a5da35f989
begin
	f = Figure()
	ax1 = Axis(f[1, 1])
	
	scatter!(density_array_cal, intensity_array)
	lines!(density_array_cal, linearFit, color = :red)
	ax1.title = "Calibration Line (Intensity vs Density)"
	ax1.ylabel = "Intensity (HU)"
	ax1.xlabel = "Density (mg/cm^3)"
	
	f
end

# ╔═╡ 18fbaeac-73d2-4880-8309-c262e1425896
md"""
### 3 Point Calibration
"""

# ╔═╡ 52193537-2f1c-4054-8909-612454805202
begin
	arr = masked_array[:, :, slice_CCI-2:slice_CCI+2]
	single_arr = masked_array[:, :, slice_CCI]
	pixel_size = DICOMUtils.get_pixel_size(header)
end

# ╔═╡ 3280a5f0-76b1-4e66-829f-75b0e7f2802f
begin
	eroded_mask_L_HD1d = erode(erode(mask_L_HD))
	high_density_cal1d = mean(single_arr[eroded_mask_L_HD1d])
end

# ╔═╡ 606b52f7-9cd5-4ef8-8af5-5585e12a4688
begin
	eroded_mask_L_MD1d = erode(erode(mask_L_MD))
	med_density_cal1d = mean(single_arr[eroded_mask_L_MD1d])
end

# ╔═╡ 1d752f2a-4478-4391-a2a9-10e98656ad68
begin
	eroded_mask_L_LD1d = erode(erode(mask_L_LD))
	low_density_cal1d = mean(single_arr[eroded_mask_L_LD1d])
end

# ╔═╡ b6425bc8-4cb8-4429-b080-064d392a359f
begin
	intensity_array3 = [0, cal_insert_mean, med_density_cal1d, high_density_cal1d]
	density_array_cal3 = [0, 200, 400, 800]
	df_cal3 = DataFrame(:density => density_array_cal3, :intensity => intensity_array3)
	linearRegressor3 = lm(@formula(intensity ~ density), df_cal3)
	linearFit3 = predict(linearRegressor3)
	m3 = linearRegressor3.model.pp.beta0[2]
	b3 = linearRegressor3.model.rr.mu[1]
end

# ╔═╡ 3cf9a717-e872-41e9-ad3a-22f546a6ebf0
let
	f = Figure()
	ax1 = Axis(f[1, 1])
	
	scatter!(density_array_cal3, intensity_array3)
	lines!(density_array_cal3, linearFit3, color = :red)
	ax1.title = "Calibration Line (Intensity vs Density)"
	ax1.ylabel = "Intensity (HU)"
	ax1.xlabel = "Density (mg/cm^3)"
	
	f
end

# ╔═╡ 422b6329-4b1f-4daa-9a03-4acad9c60d38
begin
	# # 1 Point Calibration
	# density(intensity) = (intensity - b) / m
	# intensity(ρ) = m*ρ + b

	# 3 Point Calibration
	density(intensity) = (intensity - b3) / m3
	intensity(ρ) = m3*ρ + b3
end

# ╔═╡ db84d9bb-e0da-4f81-88c8-c2a20c3ce006
md"""
# Score Large Inserts
"""

# ╔═╡ 2d83f341-24e4-4f65-885b-e5852e1c5d96
md"""
## High Density
"""

# ╔═╡ 2d1f54af-2b00-49da-8c45-fd7cc3277234
begin
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_HD_3D[:, :, z] = mask_L_HD
	end
end;

# ╔═╡ 4d08c2b9-b373-4fc4-ad4f-af7621f8e55f
md"""
#### Dilated mask
"""

# ╔═╡ dca4024a-de02-4f2a-a2eb-b458d8089c7d
dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D));

# ╔═╡ 4da96d00-fc05-4b05-b41e-ff6527c2de36
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ 258d3494-a618-4db6-84a2-189b24bb15db
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ 337109bb-42b2-4f6f-83f6-e0e98e8838c3
md"""
#### Ring (background) mask
"""

# ╔═╡ 6eece123-6547-4145-8c76-414a699823e6
ring_mask_L_HD = dilate(dilate(dilate(dilate(mask_L_HD_3D)))) - dilate(dilate(dilate(mask_L_HD_3D)));

# ╔═╡ 5faa2499-1017-4182-92bc-cb6f18401617
@bind g4 overlay_mask_bind(ring_mask_L_HD)

# ╔═╡ 6af0d533-6eb4-4129-8b12-a5eefe82b2c0
overlay_mask_plot(arr, ring_mask_L_HD, g4, "ring mask")

# ╔═╡ 45557e45-f916-43a9-a2c5-944b154ee264
begin
	single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
	s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
end

# ╔═╡ 1a7c8177-6f86-465e-b0dd-2756ff0db034
S_Obj_HD = intensity(800)

# ╔═╡ 8305356e-9b8e-4ccf-83bd-2258b4bfa2d6
begin
	alg_L_HD = Integrated(arr[mask_L_HD_3D])
	ρ_hd = 0.8 # mg/mm^3
	mass_l_hd = score(s_bkg_L_HD, S_Obj_HD, pixel_size, ρ_hd, alg_L_HD)
end

# ╔═╡ 29d535c8-994d-42e1-ad34-48453a3f92f4
md"""
## Medium Density
"""

# ╔═╡ f8379cf2-ef15-42a3-a82f-63343ef7d9d9
begin
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_MD_3D[:, :, z] = mask_L_MD
	end
end;

# ╔═╡ ec4c1348-6eac-4372-88b9-7323d20a80f9
md"""
#### Dilated mask
"""

# ╔═╡ f67d19c6-9496-480e-8fc2-35931434bb0b
dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D));

# ╔═╡ 1e5bb04d-ecd8-4036-85e9-05ba522e65fd
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ 5b01eec5-b967-4e52-ada2-d7e6f98e8743
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ bd50ca64-f6eb-4827-b8eb-a4186648d3c7
md"""
#### Ring (background) mask
"""

# ╔═╡ 7a920c7e-1859-4237-8499-b00fb530493b
ring_mask_L_MD = dilate(dilate(dilate(dilate(mask_L_MD_3D)))) - dilate(dilate(dilate(mask_L_MD_3D)));

# ╔═╡ 370c4549-69f5-4835-ae55-fc19e61c6c01
@bind h4 overlay_mask_bind(ring_mask_L_MD)

# ╔═╡ 835b7446-c0dc-4c31-a172-d52af2a81b89
overlay_mask_plot(arr, ring_mask_L_MD, h4, "ring mask")

# ╔═╡ d189b76f-d6ca-455a-9432-ba054377557a
begin
	single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
	s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
end

# ╔═╡ 16e21e05-b582-4986-b2e6-f7e40ea77f39
S_Obj_MD = intensity(400)

# ╔═╡ b7a1a6d7-82a6-40aa-afcb-5323a2642c56
begin
	alg_L_MD = Integrated(arr[mask_L_MD_3D])
	ρ_md = 0.4 # mg/mm^3
	mass_l_md = score(s_bkg_L_MD, S_Obj_MD, pixel_size, ρ_md, alg_L_MD)
end

# ╔═╡ a42e982c-16e5-45b3-a9a5-859375b72f05
md"""
## Low Density
"""

# ╔═╡ ebd8554c-8048-4739-a6bb-39cc73f7ee56
begin
	mask_L_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_LD_3D[:, :, z] = mask_L_LD
	end
end;

# ╔═╡ 6d2174b3-16ce-4fe2-9ec1-9c945fa6629e
md"""
#### Dilated mask
"""

# ╔═╡ a7caa2fb-b7a3-4b72-bcb9-07ac9475451b
dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D));

# ╔═╡ ba2c6db8-b246-45ba-b667-e7b1883894db
@bind i2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ a05fb79c-3782-460a-a4ee-54b943198e2e
overlay_mask_plot(arr, dilated_mask_L_LD, i2, "dilated mask")

# ╔═╡ fc0f8bd0-40dc-410b-a980-252179ecd3fc
md"""
#### Ring (background) mask
"""

# ╔═╡ 41313ede-e5d4-4df4-9881-a08c10ac227a
ring_mask_L_LD = dilate(dilate(dilate(dilate(mask_L_LD_3D)))) - dilate(dilate(dilate(mask_L_LD_3D)));

# ╔═╡ e9aa29a4-4480-44da-aeb8-1de14cf382b2
@bind i4 overlay_mask_bind(ring_mask_L_LD)

# ╔═╡ c3638503-8a36-4416-8842-88d4c5bda21c
overlay_mask_plot(arr, ring_mask_L_LD, i4, "ring mask")

# ╔═╡ d1e942c5-8ef5-443d-ab5e-083e0d81973e
begin	
	single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
	s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
end

# ╔═╡ 54ba65d1-e7e8-40cb-b016-43686789e0e6
S_Obj_LD = intensity(200)

# ╔═╡ dc424553-ab65-4013-aa1e-0d8163ac0c6f
begin
	alg_L_LD = Integrated(arr[mask_L_LD_3D])
	ρ_ld = 0.2 # mg/mm^3
	mass_l_ld = score(s_bkg_L_LD, S_Obj_LD, pixel_size, ρ_ld, alg_L_LD)
end

# ╔═╡ 24b943da-9319-4360-9d41-992f41752bb8
md"""
# Score Medium Inserts
"""

# ╔═╡ ea1e7093-51c3-4058-ad9b-9a43e67b019f
md"""
## High Density
"""

# ╔═╡ 24636ba1-6ac6-42d9-9167-823cf31444a3
begin
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_HD_3D[:, :, z] = mask_M_HD
	end
end;

# ╔═╡ fa5e0bef-3991-4b42-921d-60015da8000a
md"""
#### Dilated mask
"""

# ╔═╡ 21f76e13-8748-4eaf-a38e-4fd8a0e9d861
dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 344026d1-d9b5-4dd0-bce5-5b668e9545cc
@bind j2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ f6e68cba-833a-467e-af86-658dcc5bf6c7
overlay_mask_plot(arr, dilated_mask_M_HD, j2, "dilated mask")

# ╔═╡ 33bf73da-96c1-4ea3-a6a2-c022c7b64bb9
md"""
#### Ring (background) mask
"""

# ╔═╡ 31d733e5-ca48-4047-9ba1-91b3e6655647
ring_mask_M_HD = dilate(dilate(dilate(dilate(dilate(mask_M_HD_3D))))) - dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 894f3319-b750-490e-8e45-bc3710efb0ce
@bind j4 overlay_mask_bind(ring_mask_M_HD)

# ╔═╡ 82acb4b7-85ae-42c7-bf42-fe2ae1abdaf1
overlay_mask_plot(arr, ring_mask_M_HD, j4, "ring mask")

# ╔═╡ 3163aaaa-2463-4664-aa6c-6e44988aa337
begin
	single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
	s_bkg_M_HD = mean(single_arr[single_ring_mask_M_HD])
end

# ╔═╡ 8af7d835-c43b-40b7-a934-4f6ce7661535
begin
	alg_M_HD = Integrated(arr[mask_M_HD_3D])
	mass_m_hd = score(s_bkg_M_HD, S_Obj_HD, pixel_size, ρ_hd, alg_M_HD)
end

# ╔═╡ 18eae514-f9f8-4cd3-9626-0d702a8fda80
md"""
## Medium Density
"""

# ╔═╡ bbdb568d-c6d4-4eb7-b2d1-785a7df1ed6c
begin
	mask_M_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_MD_3D[:, :, z] = mask_M_MD
	end
end;

# ╔═╡ 43432b7f-cb36-4fc7-9b2a-97920880004a
md"""
#### Dilated mask
"""

# ╔═╡ 5af97cce-0237-4f41-897a-688f9ed98652
dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))));

# ╔═╡ 5a320923-1966-4a68-b90d-af1da17e79ea
@bind k2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ 439dd6bb-ae38-4cd0-ad89-002c9d91a7c1
overlay_mask_plot(arr, dilated_mask_M_MD, k2, "dilated mask")

# ╔═╡ 63abca9e-2111-4799-be13-ab3f86d5b082
md"""
#### Ring (background) mask
"""

# ╔═╡ 7e7b5ad2-a944-4439-8361-0c88d85b9985
ring_mask_M_MD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))));

# ╔═╡ 495cc420-a24f-41d2-a0cb-116f778c0bff
@bind k4 overlay_mask_bind(ring_mask_M_MD)

# ╔═╡ 2e52ff34-db03-404c-ad98-c8f62df0d181
overlay_mask_plot(arr, ring_mask_M_MD, k4, "ring mask")

# ╔═╡ 0687ccfd-d58a-4fa4-aa00-815101997a75
begin
	single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
	s_bkg_M_MD = mean(single_arr[single_ring_mask_M_MD])
end

# ╔═╡ 980785a4-0c19-4208-89eb-2b2b4e85523d
begin
	alg_M_MD = Integrated(arr[mask_M_MD_3D])
	mass_m_md = score(s_bkg_M_MD, S_Obj_MD, pixel_size, ρ_md, alg_M_MD)
end

# ╔═╡ 48400ead-acd0-4c5b-8ff7-87ee2b07b9d1
md"""
## Low Density
"""

# ╔═╡ 7d68622e-6a35-4213-a74b-5c0736995579
begin
	mask_M_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_LD_3D[:, :, z] = mask_M_LD
	end
end;

# ╔═╡ c8a9335b-c271-4097-be57-455119cd461e
md"""
#### Dilated mask
"""

# ╔═╡ 341c1d34-4b2d-4a10-ab04-da5d0dbd579f
dilated_mask_M_LD = dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ 223efc65-64a0-4d27-b829-3288445967d4
@bind l2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ e2a2251d-9964-4fa3-becc-09d0b12cd7f0
overlay_mask_plot(arr, dilated_mask_M_LD, l2, "dilated mask")

# ╔═╡ f7c88984-fcaf-4173-8f8b-023979aa472f
md"""
#### Ring (background) mask
"""

# ╔═╡ 3f8f91db-d916-4bdc-a678-fdab41a9d05e
ring_mask_M_LD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ b03ba157-8b69-4616-8723-f0396bf184c6
@bind l4 overlay_mask_bind(ring_mask_M_LD)

# ╔═╡ f58e9581-8f7d-477e-8820-9c13e210b385
overlay_mask_plot(arr, ring_mask_M_LD, l4, "ring mask")

# ╔═╡ ceecd152-e1f6-497a-872b-c4de20345499
begin
	single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
	s_bkg_M_LD = mean(single_arr[single_ring_mask_M_LD])
end

# ╔═╡ 1a280788-84f9-4a3d-9084-b1ab9c2ba514
begin
	alg_M_LD = Integrated(arr[mask_M_LD_3D])
	mass_m_ld = score(s_bkg_M_LD, S_Obj_LD, pixel_size, ρ_ld, alg_M_LD)
end

# ╔═╡ 16d56ff5-238e-4916-bac6-56fe70164b54
md"""
# Score Small Inserts
"""

# ╔═╡ 37acb2c8-30e2-441d-90d2-884ea0fd2d49
md"""
## High Density
"""

# ╔═╡ d2d48580-177b-4f67-9bbf-022ce40f4771
begin
	mask_S_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_HD_3D[:, :, z] = mask_S_HD
	end
end;

# ╔═╡ fdbf5185-6a1f-440b-aa81-81529e045b73
md"""
#### Dilated mask
"""

# ╔═╡ 62b45cc8-af83-48a1-bb33-73335a6081bf
dilated_mask_S_HD = dilate(dilate(dilate(dilate(dilate((mask_S_HD_3D))))));

# ╔═╡ b1a6f986-4d86-4c62-a662-552bdf9421f9
@bind m2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ de05e5fd-3972-41e9-bdaa-eba6631d8ecb
overlay_mask_plot(arr, dilated_mask_S_HD, m2, "dilated mask")

# ╔═╡ 339617d5-49c4-4bad-8c03-579438b55799
md"""
#### Ring (background) mask
"""

# ╔═╡ cb248de5-bcfd-4ce8-ba9b-c0c39067bbc3
ring_mask_S_HD = dilate(dilate(dilate(dilate(dilate(mask_S_HD_3D))))) - dilate(dilate(dilate(dilate(mask_S_HD_3D))));

# ╔═╡ b712ebf8-bf99-40c8-b861-19e6dfd968c0
@bind m4 overlay_mask_bind(ring_mask_S_HD)

# ╔═╡ d7bbcda2-37bf-43c2-9a4e-e6027d868d26
overlay_mask_plot(arr, ring_mask_S_HD, m4, "ring mask")

# ╔═╡ 10550a99-76b9-487e-922b-3239e7f811f4
begin
	single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
	s_bkg_S_HD = mean(single_arr[single_ring_mask_S_HD])
end

# ╔═╡ 933791e4-934d-40fd-8928-eccb68273059
begin
	alg_S_HD = Integrated(arr[mask_S_HD_3D])
	mass_s_hd = score(s_bkg_S_HD, S_Obj_HD, pixel_size, ρ_hd, alg_S_HD)
end

# ╔═╡ cdc09f53-4f83-46ff-8124-bdf54cea5097
md"""
## Medium Density
"""

# ╔═╡ a71cac4e-9ade-4ffc-a4f0-b00c9345221b
begin
	mask_S_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_MD_3D[:, :, z] = mask_S_MD
	end
end;

# ╔═╡ d27098ee-43ae-4072-9843-11934f377df6
md"""
#### Dilated mask
"""

# ╔═╡ fbbd86c5-0acf-496a-b6ee-3dcb090435db
dilated_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D)))));

# ╔═╡ 6467393e-8590-45db-924b-0033e9ddb724
@bind n2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ a85ad33c-4d57-417d-b223-c3566f87c9ff
overlay_mask_plot(arr, dilated_mask_S_MD, n2, "dilated mask")

# ╔═╡ f1057c02-9060-4d44-a0d9-6b04245c9ed0
md"""
#### Ring (background) mask
"""

# ╔═╡ c32ec1cb-825b-477b-927d-7d72b1deb99d
ring_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D))))) - dilate(dilate(dilate(dilate(mask_S_MD_3D))));

# ╔═╡ 9b254664-72a6-4030-83b0-27db11c222ba
@bind n4 overlay_mask_bind(ring_mask_S_MD)

# ╔═╡ 60370f6d-51f3-436f-ab32-400de8d98f2b
overlay_mask_plot(arr, ring_mask_S_MD, n4, "ring mask")

# ╔═╡ d62d50a8-dd17-4d25-b353-fdfedb456ff0
begin
	single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
	s_bkg_S_MD = mean(single_arr[single_ring_mask_S_MD])
end

# ╔═╡ a77205bb-eab2-4a2a-90b7-2f4844f039de
begin
	alg_S_MD = Integrated(arr[mask_S_MD_3D])
	mass_s_md = score(s_bkg_S_MD, S_Obj_MD, pixel_size, ρ_md, alg_S_MD)
end

# ╔═╡ b3455451-651e-4f70-87f4-08f9db13d9b5
md"""
## Low Density
"""

# ╔═╡ 7d71e326-ab38-4ad3-9e20-abed47caf9e9
begin
	mask_S_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_LD_3D[:, :, z] = mask_S_LD
	end
end;

# ╔═╡ 35b79339-bb2e-4987-a0f4-583a732c71fd
md"""
#### Dilated mask
"""

# ╔═╡ e0861340-3bef-4dd4-a21c-b0151d0f6d02
dilated_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D)))));

# ╔═╡ e35a7298-79fa-42b3-a2c7-faec6b5c9af4
@bind o2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ 86fac660-6ef4-4041-80d8-61eb6edc2a39
overlay_mask_plot(arr, dilated_mask_S_LD, o2, "dilated mask")

# ╔═╡ 005f8151-6b2b-4476-b92d-e63d0d216325
md"""
#### Ring (background) mask
"""

# ╔═╡ 8c379977-38dd-41e8-8956-05f48102cc67
ring_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D))))) - dilate(dilate(dilate(dilate(mask_S_LD_3D))));

# ╔═╡ b71e334a-8452-453c-a7a4-6d5f6eafd0bc
@bind o4 overlay_mask_bind(ring_mask_S_LD)

# ╔═╡ 13806613-6eca-4b9b-8216-894b25ff1b01
overlay_mask_plot(arr, ring_mask_S_LD, o4, "ring mask")

# ╔═╡ 9a8f53e1-4900-4b64-a09f-5a9b8fd5a0c5
begin
	single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
	s_bkg_S_LD = mean(single_arr[single_ring_mask_S_LD])
end

# ╔═╡ 05233816-9613-4cb6-90f0-bac0f01816c0
begin
	alg_S_LD = Integrated(arr[mask_S_LD_3D])
	mass_s_ld = score(s_bkg_S_LD, S_Obj_LD, pixel_size, ρ_ld, alg_S_LD)
end

# ╔═╡ b7f35164-9e08-49ea-aa40-f968fbbb92c2
md"""
# Results
"""

# ╔═╡ b5e14014-12bf-4004-bccc-2f3c759200b5
density_array = [0, 200, 400, 800]

# ╔═╡ 88edf2c7-3c5a-4b90-9c50-5d448a358ce8
inserts = [
	"Low Density",
	"Medium Density",
	"High Density"
]

# ╔═╡ 786eb037-e39e-41e1-9a63-6106800cb306
ground_truth_mass_large = [
	19.6,
	39.3,
	78.5
] # mg

# ╔═╡ fe51b3b3-31bc-454f-bc09-db5a7e6c6941
calculated_mass_large = [
	mass_l_ld,
	mass_l_md,
	mass_l_hd
]

# ╔═╡ 2867365f-f7fe-460c-b8fd-eb8640f7dc73
ground_truth_mass_medium = [
	4.2,
	8.5,
	17.0
]

# ╔═╡ 61cffd44-a5e8-4223-85ef-eb4f2b22de2b
calculated_mass_medium = [
	mass_m_ld,
	mass_m_md,
	mass_m_hd
]

# ╔═╡ 60324b4b-e880-4ea3-9831-fca775eb22ab
ground_truth_mass_small = [
	0.2,
	0.3,
	0.6
]

# ╔═╡ 7c4fe694-7343-4e87-aed4-8895802e34a9
calculated_mass_small = [
	mass_s_ld,
	mass_s_md,
	mass_s_hd
]

# ╔═╡ 9f12b7bd-ccc2-4108-9171-1902a0c75c30
df = DataFrame(
	scan = scan,
	inserts = inserts,
	ground_truth_mass_large = ground_truth_mass_large,
	calculated_mass_large = calculated_mass_large,
	ground_truth_mass_medium = ground_truth_mass_medium,
	calculated_mass_medium = calculated_mass_medium,
	ground_truth_mass_small = ground_truth_mass_small,
	calculated_mass_small = calculated_mass_small
)

# ╔═╡ fef2fcda-cc64-469e-9f17-946f06a3e4ae
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

# ╔═╡ 4f9c5043-3485-45b8-b493-d62780134251
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

# ╔═╡ 197e0479-f15a-43f9-a7cf-5b5d58b5a901
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

# ╔═╡ b030dcf1-fc55-4e32-be9b-6cba1d3456f0
md"""
### Save Results
"""

# ╔═╡ 5a9c7e99-0e6e-4162-8ff7-bff1935fd318
if ~isdir(string(cd(pwd, "..") , "/data/output/", VENDER))
	mkdir(string(cd(pwd, "..") , "/data/output/", VENDER))
end

# ╔═╡ 7976fbad-ae1f-40e8-961f-2ea68f752caa
output_path = string(cd(pwd, "..") , "/data/output/", VENDER, "/", scan, ".csv")

# ╔═╡ d8aeb68a-697d-4477-8b6c-14712f50a945
# CSV.write(output_path, df)

# ╔═╡ 3d6114a2-1356-4998-be56-800395262bc0
md"""
### Save full df
"""

# ╔═╡ cfca657a-da63-436d-b316-a4e9a63ae8e7
dfs = []

# ╔═╡ a2bf4e79-7f23-45f8-8261-8c10bcd7eb3c
push!(dfs, df)

# ╔═╡ b1a053ad-f18f-487d-8c38-01dbe09fd03a
if length(dfs) == 10
	global new_df = vcat(dfs[1:10]...)
	output_path_new = string(cd(pwd, "..") , "/data/output/", VENDER, "/full.csv")
	CSV.write(output_path_new, new_df)
end

# ╔═╡ Cell order:
# ╠═3c162d31-407a-4c50-a8a7-30153192f207
# ╠═0aeb13cf-6299-473a-908e-5b0be1b23c86
# ╟─c4bba726-9f22-4f76-a4ae-966efd13afcf
# ╟─80ea14ec-0128-4899-bcae-5d75d6654a89
# ╠═cb6af1f1-9e54-4172-9e65-bcf534840716
# ╠═47504ab4-a52f-418d-b9de-3468304c9578
# ╠═d41f610c-928b-4a36-a7f8-ac50ff81f6db
# ╠═7efcfa5c-8d99-4e26-9538-efa1d2bfeea6
# ╠═017121d9-a3bd-481a-ba89-7cd261f0ffe7
# ╠═49fe69bc-7c19-4dc6-86df-ceaccbdd65f2
# ╠═e95b1cbc-82a1-47db-8026-e2773f0d0cff
# ╟─a4f1ebf0-6f99-42fe-96f3-765caf5f207d
# ╟─8e874d2a-95bc-414b-91cd-1e6df8220c18
# ╟─1fec01bc-445a-4961-b518-f25a0de737d3
# ╟─4494460e-4ffc-4317-9a30-21ac8e005bea
# ╟─9d4357e3-5ee5-4bf9-8879-473a4da59812
# ╠═4a947165-c565-4c6d-b63c-8cc714d0d039
# ╟─f1339241-169d-40c1-becc-8f230c12e835
# ╠═641043a6-6d00-4f20-96ad-4b588bda2893
# ╟─ede1723d-ac8b-43a5-a6d8-6a92fcbaa2fa
# ╟─13671de8-9397-422f-81a4-57781e881178
# ╟─c07994fb-87c8-4a61-8dd9-7fe50b087b13
# ╟─2988d4aa-a2f2-45a5-9e96-fb9c233e9430
# ╠═06c46a2f-9f03-46ec-88ed-6219bd2a5f8d
# ╟─66c50557-49d0-4d9f-beae-404b6784f079
# ╠═8806c9d6-c416-47be-981e-a03d35d45602
# ╟─f30cd88e-c6ca-495d-8ef9-6f029182b5a0
# ╠═7d67fa1f-74de-4a18-9b0d-eb46fa98f420
# ╠═3bb994fa-e4e8-4e98-bf1b-5e896f21e96e
# ╠═0aece8c7-c8f7-49ca-b68e-b2fda92060c0
# ╠═5eb2cb0e-5353-4bd2-86ac-62f17be96f4c
# ╟─a488da37-5fd7-4d23-a958-d2fe3f1249b1
# ╠═d7352bc6-2432-4926-bb70-c791a0050e6b
# ╠═e69626d2-a4eb-46f6-b3ca-491d427b7f70
# ╠═88efbea6-138f-4885-9afb-543648e20152
# ╠═a50ca0be-2613-498a-bcb8-5bdf10b4eeed
# ╠═a2ea95bc-8840-4fbe-a29d-3618c629e9dc
# ╟─a574b69a-023e-4f1c-a09d-2f75e2eba8ec
# ╠═4c5b5a15-a456-4f67-977d-55a727f6ac7f
# ╠═29a2117d-0462-4aac-8ff0-060f1c6ba848
# ╠═3f17d44f-59fb-4cc3-93d6-df025d92f1bf
# ╠═0362be1a-514b-4133-8c79-b2f554d7ed9f
# ╠═952a3a84-ce2d-4628-8c84-a2fc8f017a96
# ╠═67caa8b2-7fb7-49a0-b515-48a5da35f989
# ╟─18fbaeac-73d2-4880-8309-c262e1425896
# ╠═52193537-2f1c-4054-8909-612454805202
# ╠═3280a5f0-76b1-4e66-829f-75b0e7f2802f
# ╠═606b52f7-9cd5-4ef8-8af5-5585e12a4688
# ╠═1d752f2a-4478-4391-a2a9-10e98656ad68
# ╠═b6425bc8-4cb8-4429-b080-064d392a359f
# ╟─3cf9a717-e872-41e9-ad3a-22f546a6ebf0
# ╠═422b6329-4b1f-4daa-9a03-4acad9c60d38
# ╟─db84d9bb-e0da-4f81-88c8-c2a20c3ce006
# ╟─2d83f341-24e4-4f65-885b-e5852e1c5d96
# ╠═2d1f54af-2b00-49da-8c45-fd7cc3277234
# ╟─4d08c2b9-b373-4fc4-ad4f-af7621f8e55f
# ╠═dca4024a-de02-4f2a-a2eb-b458d8089c7d
# ╟─4da96d00-fc05-4b05-b41e-ff6527c2de36
# ╠═258d3494-a618-4db6-84a2-189b24bb15db
# ╟─337109bb-42b2-4f6f-83f6-e0e98e8838c3
# ╠═6eece123-6547-4145-8c76-414a699823e6
# ╟─5faa2499-1017-4182-92bc-cb6f18401617
# ╠═6af0d533-6eb4-4129-8b12-a5eefe82b2c0
# ╠═45557e45-f916-43a9-a2c5-944b154ee264
# ╠═1a7c8177-6f86-465e-b0dd-2756ff0db034
# ╠═8305356e-9b8e-4ccf-83bd-2258b4bfa2d6
# ╟─29d535c8-994d-42e1-ad34-48453a3f92f4
# ╠═f8379cf2-ef15-42a3-a82f-63343ef7d9d9
# ╟─ec4c1348-6eac-4372-88b9-7323d20a80f9
# ╠═f67d19c6-9496-480e-8fc2-35931434bb0b
# ╟─1e5bb04d-ecd8-4036-85e9-05ba522e65fd
# ╠═5b01eec5-b967-4e52-ada2-d7e6f98e8743
# ╟─bd50ca64-f6eb-4827-b8eb-a4186648d3c7
# ╠═7a920c7e-1859-4237-8499-b00fb530493b
# ╟─370c4549-69f5-4835-ae55-fc19e61c6c01
# ╠═835b7446-c0dc-4c31-a172-d52af2a81b89
# ╠═d189b76f-d6ca-455a-9432-ba054377557a
# ╠═16e21e05-b582-4986-b2e6-f7e40ea77f39
# ╠═b7a1a6d7-82a6-40aa-afcb-5323a2642c56
# ╟─a42e982c-16e5-45b3-a9a5-859375b72f05
# ╠═ebd8554c-8048-4739-a6bb-39cc73f7ee56
# ╟─6d2174b3-16ce-4fe2-9ec1-9c945fa6629e
# ╠═a7caa2fb-b7a3-4b72-bcb9-07ac9475451b
# ╟─ba2c6db8-b246-45ba-b667-e7b1883894db
# ╠═a05fb79c-3782-460a-a4ee-54b943198e2e
# ╟─fc0f8bd0-40dc-410b-a980-252179ecd3fc
# ╠═41313ede-e5d4-4df4-9881-a08c10ac227a
# ╟─e9aa29a4-4480-44da-aeb8-1de14cf382b2
# ╠═c3638503-8a36-4416-8842-88d4c5bda21c
# ╠═d1e942c5-8ef5-443d-ab5e-083e0d81973e
# ╠═54ba65d1-e7e8-40cb-b016-43686789e0e6
# ╠═dc424553-ab65-4013-aa1e-0d8163ac0c6f
# ╟─24b943da-9319-4360-9d41-992f41752bb8
# ╟─ea1e7093-51c3-4058-ad9b-9a43e67b019f
# ╠═24636ba1-6ac6-42d9-9167-823cf31444a3
# ╟─fa5e0bef-3991-4b42-921d-60015da8000a
# ╠═21f76e13-8748-4eaf-a38e-4fd8a0e9d861
# ╟─344026d1-d9b5-4dd0-bce5-5b668e9545cc
# ╠═f6e68cba-833a-467e-af86-658dcc5bf6c7
# ╟─33bf73da-96c1-4ea3-a6a2-c022c7b64bb9
# ╠═31d733e5-ca48-4047-9ba1-91b3e6655647
# ╟─894f3319-b750-490e-8e45-bc3710efb0ce
# ╠═82acb4b7-85ae-42c7-bf42-fe2ae1abdaf1
# ╠═3163aaaa-2463-4664-aa6c-6e44988aa337
# ╠═8af7d835-c43b-40b7-a934-4f6ce7661535
# ╟─18eae514-f9f8-4cd3-9626-0d702a8fda80
# ╠═bbdb568d-c6d4-4eb7-b2d1-785a7df1ed6c
# ╟─43432b7f-cb36-4fc7-9b2a-97920880004a
# ╠═5af97cce-0237-4f41-897a-688f9ed98652
# ╟─5a320923-1966-4a68-b90d-af1da17e79ea
# ╠═439dd6bb-ae38-4cd0-ad89-002c9d91a7c1
# ╟─63abca9e-2111-4799-be13-ab3f86d5b082
# ╠═7e7b5ad2-a944-4439-8361-0c88d85b9985
# ╟─495cc420-a24f-41d2-a0cb-116f778c0bff
# ╠═2e52ff34-db03-404c-ad98-c8f62df0d181
# ╠═0687ccfd-d58a-4fa4-aa00-815101997a75
# ╠═980785a4-0c19-4208-89eb-2b2b4e85523d
# ╟─48400ead-acd0-4c5b-8ff7-87ee2b07b9d1
# ╠═7d68622e-6a35-4213-a74b-5c0736995579
# ╟─c8a9335b-c271-4097-be57-455119cd461e
# ╠═341c1d34-4b2d-4a10-ab04-da5d0dbd579f
# ╟─223efc65-64a0-4d27-b829-3288445967d4
# ╠═e2a2251d-9964-4fa3-becc-09d0b12cd7f0
# ╟─f7c88984-fcaf-4173-8f8b-023979aa472f
# ╠═3f8f91db-d916-4bdc-a678-fdab41a9d05e
# ╟─b03ba157-8b69-4616-8723-f0396bf184c6
# ╠═f58e9581-8f7d-477e-8820-9c13e210b385
# ╠═ceecd152-e1f6-497a-872b-c4de20345499
# ╠═1a280788-84f9-4a3d-9084-b1ab9c2ba514
# ╟─16d56ff5-238e-4916-bac6-56fe70164b54
# ╟─37acb2c8-30e2-441d-90d2-884ea0fd2d49
# ╠═d2d48580-177b-4f67-9bbf-022ce40f4771
# ╟─fdbf5185-6a1f-440b-aa81-81529e045b73
# ╠═62b45cc8-af83-48a1-bb33-73335a6081bf
# ╟─b1a6f986-4d86-4c62-a662-552bdf9421f9
# ╠═de05e5fd-3972-41e9-bdaa-eba6631d8ecb
# ╟─339617d5-49c4-4bad-8c03-579438b55799
# ╠═cb248de5-bcfd-4ce8-ba9b-c0c39067bbc3
# ╟─b712ebf8-bf99-40c8-b861-19e6dfd968c0
# ╠═d7bbcda2-37bf-43c2-9a4e-e6027d868d26
# ╠═10550a99-76b9-487e-922b-3239e7f811f4
# ╠═933791e4-934d-40fd-8928-eccb68273059
# ╟─cdc09f53-4f83-46ff-8124-bdf54cea5097
# ╠═a71cac4e-9ade-4ffc-a4f0-b00c9345221b
# ╟─d27098ee-43ae-4072-9843-11934f377df6
# ╠═fbbd86c5-0acf-496a-b6ee-3dcb090435db
# ╟─6467393e-8590-45db-924b-0033e9ddb724
# ╠═a85ad33c-4d57-417d-b223-c3566f87c9ff
# ╟─f1057c02-9060-4d44-a0d9-6b04245c9ed0
# ╠═c32ec1cb-825b-477b-927d-7d72b1deb99d
# ╟─9b254664-72a6-4030-83b0-27db11c222ba
# ╠═60370f6d-51f3-436f-ab32-400de8d98f2b
# ╠═d62d50a8-dd17-4d25-b353-fdfedb456ff0
# ╠═a77205bb-eab2-4a2a-90b7-2f4844f039de
# ╟─b3455451-651e-4f70-87f4-08f9db13d9b5
# ╠═7d71e326-ab38-4ad3-9e20-abed47caf9e9
# ╟─35b79339-bb2e-4987-a0f4-583a732c71fd
# ╠═e0861340-3bef-4dd4-a21c-b0151d0f6d02
# ╟─e35a7298-79fa-42b3-a2c7-faec6b5c9af4
# ╠═86fac660-6ef4-4041-80d8-61eb6edc2a39
# ╟─005f8151-6b2b-4476-b92d-e63d0d216325
# ╠═8c379977-38dd-41e8-8956-05f48102cc67
# ╟─b71e334a-8452-453c-a7a4-6d5f6eafd0bc
# ╠═13806613-6eca-4b9b-8216-894b25ff1b01
# ╠═9a8f53e1-4900-4b64-a09f-5a9b8fd5a0c5
# ╠═05233816-9613-4cb6-90f0-bac0f01816c0
# ╟─b7f35164-9e08-49ea-aa40-f968fbbb92c2
# ╠═b5e14014-12bf-4004-bccc-2f3c759200b5
# ╠═88edf2c7-3c5a-4b90-9c50-5d448a358ce8
# ╠═786eb037-e39e-41e1-9a63-6106800cb306
# ╠═fe51b3b3-31bc-454f-bc09-db5a7e6c6941
# ╠═2867365f-f7fe-460c-b8fd-eb8640f7dc73
# ╠═61cffd44-a5e8-4223-85ef-eb4f2b22de2b
# ╠═60324b4b-e880-4ea3-9831-fca775eb22ab
# ╠═7c4fe694-7343-4e87-aed4-8895802e34a9
# ╠═9f12b7bd-ccc2-4108-9171-1902a0c75c30
# ╟─fef2fcda-cc64-469e-9f17-946f06a3e4ae
# ╟─4f9c5043-3485-45b8-b493-d62780134251
# ╟─197e0479-f15a-43f9-a7cf-5b5d58b5a901
# ╟─b030dcf1-fc55-4e32-be9b-6cba1d3456f0
# ╠═5a9c7e99-0e6e-4162-8ff7-bff1935fd318
# ╠═7976fbad-ae1f-40e8-961f-2ea68f752caa
# ╠═d8aeb68a-697d-4477-8b6c-14712f50a945
# ╟─3d6114a2-1356-4998-be56-800395262bc0
# ╠═cfca657a-da63-436d-b316-a4e9a63ae8e7
# ╠═a2bf4e79-7f23-45f8-8261-8c10bcd7eb3c
# ╠═b1a053ad-f18f-487d-8c38-01dbe09fd03a

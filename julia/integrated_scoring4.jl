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

# ╔═╡ 8f0a533f-fdb6-4c1f-bc76-6d1b0f15192c
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

# ╔═╡ 9cff9635-f4f4-4d59-81a6-92390c09526c
TableOfContents()

# ╔═╡ f5090527-3372-48ea-8239-f51bf20fbba7
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDER`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ 383f5377-0435-4ea4-96d9-6413fcce1ee0
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ 34cdbadf-edfa-4708-be18-03e40f6e00fe
begin
	SCAN_NUMBER = 10
	VENDER = "Canon_Aquilion_One_Vision"
	BASE_PATH = "/Users/daleblack/Google Drive/Datasets/"
end

# ╔═╡ 8feec48a-9ea2-4aa5-96a1-544ef298a3ac
root_path = string(BASE_PATH, VENDER)

# ╔═╡ 1c7739dc-6615-4aab-b9fa-a4664c858070
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ 752cfdd4-f2cf-45fb-8f70-959d145264f6
pth = dcm_path_list[SCAN_NUMBER]

# ╔═╡ 96fc1811-ab0f-41e6-8857-31b850a03199
pth

# ╔═╡ ff90198d-14e0-4404-813f-40dd35cc1f47
scan = basename(pth)

# ╔═╡ b58a6353-8574-44bc-8312-603b9b30d672
header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ 60cc7eea-8314-4f41-bf5d-8375a2c996e0
md"""
## Helper Functions
"""

# ╔═╡ 555426d7-4584-41eb-a577-53a5a506e97a
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ aa688492-bafb-4623-af4b-a744ca59335e
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ d16e857c-654a-48fa-81dd-b53c7092cc86
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

# ╔═╡ d6b6104a-df14-4cb0-96e6-e07be85255ad
md"""
## Segment Heart
"""

# ╔═╡ 91243a6c-ce6f-4fc1-b28c-792e8c71dd53
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ 343815f2-6b33-4408-81eb-d4b1fddf60ac
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 127344d7-38d6-4338-9ae8-edd92090b5cf
heatmap(masked_array[:, :, a], colormap=:grays)

# ╔═╡ fc3bbc09-e993-4f1d-838d-72041cb0f579
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 15]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig
end

# ╔═╡ 865f42c6-ca98-498a-89a3-d4c9986237d3
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig2
end

# ╔═╡ 0a2563dd-d104-4f1b-b178-f77ca30e96b3
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 23]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig3
end

# ╔═╡ db231239-5a21-402f-8eb5-b2e44e03b312
md"""
## Segment Calcium Rod
"""

# ╔═╡ ef56f1b2-150f-4608-944d-b5785a46e77c
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header);

# ╔═╡ 252d8c09-fcef-435e-9785-ec7ef607f108
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=cal_rod_slice, show_value=true)

# ╔═╡ cf330459-0792-48fb-9346-f87b3019a0f6
heatmap(transpose(calcium_image[:, :, c]), colormap=:grays)

# ╔═╡ eb4494c0-003d-4967-85fa-b20d8e7f9763
md"""
## Segment Calcium Inserts
"""

# ╔═╡ 23e10eee-6022-4b31-bbc2-9ddbcbda0b52
md"""
## Calibration Prep
"""

# ╔═╡ ed57792e-beb4-45af-a710-e0f6268ddac3
array_filtered = abs.(mapwindow(median, calcium_image[:, :, cal_rod_slice], (3, 3)));

# ╔═╡ f7c31a75-39c7-4be7-acbe-c36db4c252cb
bool_arr = array_filtered .> 0;

# ╔═╡ 9fcdac1c-4307-4f27-a924-b8887406fcd7
bool_arr_erode = (((erode(erode(bool_arr)))));

# ╔═╡ 46112e23-ce69-4130-832d-a26e7d5195f9
heatmap(bool_arr, colormap=:grays)

# ╔═╡ 45f4a72e-a74a-4db2-b144-576b3f22b1b0
heatmap(bool_arr_erode, colormap=:grays)

# ╔═╡ 84efd8af-873e-448d-b70d-1832e824fa3f
c_img = calcium_image[:, :, cal_rod_slice-1:cal_rod_slice+1];

# ╔═╡ e326a760-5795-4892-8dc7-88bfef009a52
begin
	mask_cal_3D = Array{Bool}(undef, size(c_img))
	for z in 1:size(c_img, 3)
		mask_cal_3D[:, :, z] = bool_arr_erode
	end
end;

# ╔═╡ 1666956f-1417-45ae-945c-7f58b98cfe9e
# hist(c_img[mask_cal_3D])

# ╔═╡ ab09d0f1-f072-4f77-9889-a2beb852dc48
# cal_insert_mean = mean(c_img[mask_cal_3D])

# ╔═╡ f12cf6c5-9e82-4572-8321-fc848f1c866e
cal_insert_mean = quantile!(c_img[mask_cal_3D], 0.7)

# ╔═╡ 9f5b80aa-afb7-4c0c-829f-30171360a87e
md"""
# Score Large Inserts
"""

# ╔═╡ 3e7253c2-b6df-4dab-8405-08ed93c03f10
arr = masked_array[:, :, slice_CCI-2:slice_CCI+2];

# ╔═╡ 1f40f903-4c8b-48bc-bb8a-fa16bc6fc297
single_arr = masked_array[:, :, slice_CCI];

# ╔═╡ c7f433f7-c12f-41eb-be35-97997d5142b7
md"""
## High Density
"""

# ╔═╡ b0f14232-0be1-4a0e-b336-3cc3a8aaff6d
md"""
#### Dilated mask
"""

# ╔═╡ 3834b0e8-85c4-4cea-bb74-6350bb3e9e7d
md"""
#### Ring (background) mask
"""

# ╔═╡ c30c8fd9-93b5-404b-a1c4-5cdeda48845f
pixel_size = DICOMUtils.get_pixel_size(header)

# ╔═╡ acfcb465-10e7-417a-9fe3-4231c86c5293
md"""
## Medium Density
"""

# ╔═╡ fa2a052f-452a-4325-a428-8db343d560ec
md"""
#### Dilated mask
"""

# ╔═╡ 7741dcb3-e433-4612-aaf5-eabb2d8aa3aa
md"""
#### Ring (background) mask
"""

# ╔═╡ 7ad75365-cb12-4db1-8354-956540204ee8
md"""
## Low Density
"""

# ╔═╡ 2d04c8c8-b652-410b-8357-b123b488c579
md"""
#### Dilated mask
"""

# ╔═╡ eb213bf2-bbd4-48dc-9834-ca4d3b650bfd
md"""
#### Ring (background) mask
"""

# ╔═╡ 0ca2ffa3-01a2-4484-a324-5d55e5609d10
md"""
# Score Medium Inserts
"""

# ╔═╡ 5779ee2f-5821-4475-aaba-c7aac7bb1d5e
md"""
## High Density
"""

# ╔═╡ dca948cb-f537-49b6-bc69-771cae751e94
md"""
#### Dilated mask
"""

# ╔═╡ 1cbd7214-11e9-4266-8fa1-7363d46b2de2
md"""
#### Ring (background) mask
"""

# ╔═╡ 0cda89b7-c4c3-4b8c-b2f5-332a90f92831
md"""
## Medium Density
"""

# ╔═╡ 39019b75-c629-4e04-970b-2813af45cdcd
md"""
#### Dilated mask
"""

# ╔═╡ 231e9c7b-9d7c-48e5-a39e-b60e7c846889
md"""
#### Ring (background) mask
"""

# ╔═╡ 27e50773-d747-4451-8546-16d5c014f37a
md"""
## Low Density
"""

# ╔═╡ 7a5c44ea-baf4-4b93-924e-f8624fe051ff
md"""
#### Dilated mask
"""

# ╔═╡ e15c09f2-0a07-4af6-9c17-44ba50393383
md"""
#### Ring (background) mask
"""

# ╔═╡ 4ad89eca-ce15-43bd-bcd2-dd36ea759c7c
md"""
# Score Small Inserts
"""

# ╔═╡ 70597737-ace4-4676-836f-44fdab22a798
md"""
## High Density
"""

# ╔═╡ fa1c77f9-5f01-47b5-aef7-d885ddae8029
md"""
#### Dilated mask
"""

# ╔═╡ 6f078eda-1fa3-4c9b-81d5-10909a9c9387
angle_factor = -4

# ╔═╡ 53f711d7-7ec2-41ea-b609-3146c01fa030
mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
		dcm_array, masked_array, header, slice_CCI, center_insert; angle_factor=angle_factor);

# ╔═╡ 95948083-8209-426e-ac12-2d7e925a6c30
masks = mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD;

# ╔═╡ bf3853ea-0b3b-476f-b9b5-b969cbcc9827
heatmap(masks, colormap=:grays)

# ╔═╡ 32085083-902c-444b-9384-b79769c6acd2
begin
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_HD_3D[:, :, z] = mask_L_HD
	end
end;

# ╔═╡ 0e31733e-a36f-4d70-970e-7b44accb2b99
dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D));

# ╔═╡ 0d2e8726-4829-46ce-b85b-fc36b5a651dc
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ 75830069-109b-474c-9c4a-448ee6b19728
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ ccb3a71e-8669-4371-b37d-1d2e38166239
ring_mask_L_HD = dilate(dilate(dilate(dilate(mask_L_HD_3D)))) - dilate(dilate(dilate(mask_L_HD_3D)));

# ╔═╡ afeb1614-5476-45db-aca2-bf3dbf9de6da
@bind g4 overlay_mask_bind(ring_mask_L_HD)

# ╔═╡ 73e0c981-e3b0-4a33-810d-04e6767428de
overlay_mask_plot(arr, ring_mask_L_HD, g4, "ring mask")

# ╔═╡ 6542d34b-89f8-42d9-9448-7d0503d35477
begin
	single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
	s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
end

# ╔═╡ 8fbd3e9c-0195-4483-92c6-9a47e6e26a9b
begin
	alg_L_HD = Integrated(arr[mask_L_HD_3D])
	ρ = 0.2 # mg/mm^3
	mass_l_hd = score(s_bkg_L_HD, cal_insert_mean, pixel_size, ρ, alg_L_HD)
end

# ╔═╡ 6e794c8a-0341-4a0c-9d96-192f5b77a8ab
begin
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_MD_3D[:, :, z] = mask_L_MD
	end
end;

# ╔═╡ d9987a6f-b8d8-4ace-8b40-234d5dde2619
dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D));

# ╔═╡ 6facfb53-1b86-4d66-9f86-0e6e1186c415
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ 338b6552-5845-4257-aa74-5fb3f8933239
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ 23839ac8-1f79-4621-afd9-9ae3ff572c62
ring_mask_L_MD = dilate(dilate(dilate(dilate(mask_L_MD_3D)))) - dilate(dilate(dilate(mask_L_MD_3D)));

# ╔═╡ dc71ee51-82fb-4367-ba3f-2ef6599c6808
@bind h4 overlay_mask_bind(ring_mask_L_MD)

# ╔═╡ 11a32a9f-d66a-46f2-bafe-dd032b9db1b0
overlay_mask_plot(arr, ring_mask_L_MD, h4, "ring mask")

# ╔═╡ e88c694a-f3ae-4114-8a2c-1d20bba7f1bd
begin
	single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
	s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
end

# ╔═╡ c2aba951-9e13-40de-8a80-3e38ed800a01
begin
	alg_L_MD = Integrated(arr[mask_L_MD_3D])
	mass_l_md = score(s_bkg_L_MD, cal_insert_mean, pixel_size, ρ, alg_L_MD)
end

# ╔═╡ 5be373b4-7d52-446d-a39c-3aed2aeb762c
begin
	mask_L_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_LD_3D[:, :, z] = mask_L_LD
	end
end;

# ╔═╡ f02480d0-a373-4a32-bdf2-c40055ed2ab9
dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D));

# ╔═╡ 7ac3b056-a9a8-44c0-a8c4-e5e89062f28c
@bind i2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ 03a74421-7912-49af-92d6-56ba38b8b7ad
overlay_mask_plot(arr, dilated_mask_L_LD, i2, "dilated mask")

# ╔═╡ 12e7a90b-27df-4cdb-99aa-c9aa8e1aec05
ring_mask_L_LD = dilate(dilate(dilate(dilate(mask_L_LD_3D)))) - dilate(dilate(dilate(mask_L_LD_3D)));

# ╔═╡ fee22a56-26bb-4dd8-a500-e973a117dbeb
@bind i4 overlay_mask_bind(ring_mask_L_LD)

# ╔═╡ ce557f39-f2d0-4b48-9e58-f232f3056b9b
overlay_mask_plot(arr, ring_mask_L_LD, i4, "ring mask")

# ╔═╡ c1809d80-7c49-43e4-bff2-d6e14a7d9f97
begin	
	single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
	s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
end

# ╔═╡ fa89795c-a6f5-4d9e-8d33-7cf42f4530d6
begin
	alg_L_LD = Integrated(arr[mask_L_LD_3D])
	mass_l_ld = score(s_bkg_L_LD, cal_insert_mean, pixel_size, ρ, alg_L_LD)
end

# ╔═╡ ab8473c4-85a9-4c66-83b1-6b8822afc855
begin
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_HD_3D[:, :, z] = mask_M_HD
	end
end;

# ╔═╡ 821240b4-19e6-4214-8b1e-cfe69c574ae4
dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 9df3801d-8984-4eb9-bc82-46e9e7c458e2
@bind j2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ d78f713a-1c6c-48c5-b1f8-ec8d89a3931c
overlay_mask_plot(arr, dilated_mask_M_HD, j2, "dilated mask")

# ╔═╡ 625402e5-5191-493a-b428-a4ee442a327d
ring_mask_M_HD = dilate(dilate(dilate(dilate(dilate(mask_M_HD_3D))))) - dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 49ea60a3-61df-4b8d-97f0-5ad4c238b4ff
@bind j4 overlay_mask_bind(ring_mask_M_HD)

# ╔═╡ 828585e7-bc8d-40f4-a76a-624f2513db55
overlay_mask_plot(arr, ring_mask_M_HD, j4, "ring mask")

# ╔═╡ dba099fe-7f13-4f61-b5ff-a4d65b7a19e2
begin
	single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
	s_bkg_M_HD = mean(single_arr[single_ring_mask_M_HD])
end

# ╔═╡ 55740593-5f41-441e-b3b8-dcbfb83eceb5
begin
	alg_M_HD = Integrated(arr[mask_M_HD_3D])
	mass_m_hd = score(s_bkg_M_HD, cal_insert_mean, pixel_size, ρ, alg_M_HD)
end

# ╔═╡ 8f0af863-52ee-4e6d-8d70-e50432720501
begin
	mask_M_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_MD_3D[:, :, z] = mask_M_MD
	end
end;

# ╔═╡ 412c307e-e323-4832-9815-98e475f4b5b3
dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))));

# ╔═╡ 1a9af931-7f46-4c71-a753-8eb0b9bfc64b
@bind k2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ 5d63180c-930c-464e-a8d0-3b2ea572f7cc
overlay_mask_plot(arr, dilated_mask_M_MD, k2, "dilated mask")

# ╔═╡ a79df362-4832-4f98-8bf9-e410df17bf71
ring_mask_M_MD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))));

# ╔═╡ ddcde058-5c8b-40d9-8244-4fe99ec7c9a6
@bind k4 overlay_mask_bind(ring_mask_M_MD)

# ╔═╡ ed112bd8-4af6-4522-a24b-47b39f3a7865
overlay_mask_plot(arr, ring_mask_M_MD, k4, "ring mask")

# ╔═╡ 5c777b1a-6b33-4b69-8797-6931a9d07a44
begin
	single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
	s_bkg_M_MD = mean(single_arr[single_ring_mask_M_MD])
end

# ╔═╡ 0d22d7e4-8519-4652-b794-83756671aa87
begin
	alg_M_MD = Integrated(arr[mask_M_MD_3D])
	mass_m_md = score(s_bkg_M_MD, cal_insert_mean, pixel_size, ρ, alg_M_MD)
end

# ╔═╡ 4168a9e9-6d44-401b-85f1-d80d6a2eeff8
begin
	mask_M_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_LD_3D[:, :, z] = mask_M_LD
	end
end;

# ╔═╡ fb8dfcff-5c15-4c30-b843-fa6db9e34e91
dilated_mask_M_LD = dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ 3c6f8e72-a699-491e-b87d-0e0ee385d470
@bind l2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ 02c01e40-d271-40eb-945f-60212dc1758e
overlay_mask_plot(arr, dilated_mask_M_LD, l2, "dilated mask")

# ╔═╡ 76f15b22-a7f6-4d8d-9911-1c072e493589
ring_mask_M_LD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ 3723dee8-a78b-4c17-944d-c393c7ab92ad
@bind l4 overlay_mask_bind(ring_mask_M_LD)

# ╔═╡ 32956aa9-4338-43f2-bba7-d520716609ec
overlay_mask_plot(arr, ring_mask_M_LD, l4, "ring mask")

# ╔═╡ 3fe8e26e-0152-4869-b8ba-7e67d1e5d101
begin
	single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
	s_bkg_M_LD = mean(single_arr[single_ring_mask_M_LD])
end

# ╔═╡ c4252a6c-878a-493f-abc1-b8c56cda2319
begin
	alg_M_LD = Integrated(arr[mask_M_LD_3D])
	mass_m_ld = score(s_bkg_M_LD, cal_insert_mean, pixel_size, ρ, alg_M_LD)
end

# ╔═╡ d5a8ba59-f167-48f3-954f-1b5a661456d5
begin
	mask_S_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_HD_3D[:, :, z] = mask_S_HD
	end
end;

# ╔═╡ 799cec3a-953d-4fe2-bc06-351f5d6f4890
dilated_mask_S_HD = dilate(dilate(dilate(dilate(dilate((mask_S_HD_3D))))));

# ╔═╡ c2c9bd13-f56b-44e7-97ed-be2e06d63ba3
@bind m2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ fac87764-f188-4058-ae48-7aa36ac3b3bf
overlay_mask_plot(arr, dilated_mask_S_HD, m2, "dilated mask")

# ╔═╡ c8bff4a0-a2a1-4e19-9d35-aceb00b365ff
md"""
#### Ring (background) mask
"""

# ╔═╡ ad67881f-9085-46dc-a645-c3f111738c94
ring_mask_S_HD = dilate(dilate(dilate(dilate(dilate(mask_S_HD_3D))))) - dilate(dilate(dilate(dilate(mask_S_HD_3D))));

# ╔═╡ 0710d8ef-aac6-4d6d-983a-2d650502c180
@bind m4 overlay_mask_bind(ring_mask_S_HD)

# ╔═╡ ba1db1c1-47c8-4c25-b650-67cb8c5c4e1b
overlay_mask_plot(arr, ring_mask_S_HD, m4, "ring mask")

# ╔═╡ ff7731cc-5933-4312-919e-e4703254a622
begin
	single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
	s_bkg_S_HD = mean(single_arr[single_ring_mask_S_HD])
end

# ╔═╡ 54e24ee8-a0f8-452a-ac2d-1379174fb56d
begin
	alg_S_HD = Integrated(arr[mask_S_HD_3D])
	mass_s_hd = score(s_bkg_S_HD, cal_insert_mean, pixel_size, ρ, alg_S_HD)
	# if mass_s_hd < 0
	# 	mass_s_hd = 0
	# end
	# mass_s_hd
end

# ╔═╡ 82bbbd01-7a15-438c-b317-22b9511e8895
md"""
## Medium Density
"""

# ╔═╡ 61025c8f-d189-4ad0-bbe1-bb24e13083d4
begin
	mask_S_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_MD_3D[:, :, z] = mask_S_MD
	end
end;

# ╔═╡ 4ae034d0-b3fa-46ca-acf7-b93062bdb5a2
md"""
#### Dilated mask
"""

# ╔═╡ 6e544d17-f156-43c2-9255-58bd6b0c1888
dilated_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D)))));

# ╔═╡ 50333385-e60d-4dd2-92e1-bcd2d9007355
@bind n2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ c3f5fc15-60f5-457f-8353-9e10926d5a7d
overlay_mask_plot(arr, dilated_mask_S_MD, n2, "dilated mask")

# ╔═╡ f66110c6-f501-4a55-8719-fc752b838de3
md"""
#### Ring (background) mask
"""

# ╔═╡ 2427ac81-bbb4-49ec-9897-d4fa59e4259b
ring_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D))))) - dilate(dilate(dilate(dilate(mask_S_MD_3D))));

# ╔═╡ 65722ef0-fbe9-4259-93f2-f5265e6ac3ba
@bind n4 overlay_mask_bind(ring_mask_S_MD)

# ╔═╡ 5c64cf04-7406-46cd-8650-05b8d7ab0964
overlay_mask_plot(arr, ring_mask_S_MD, n4, "ring mask")

# ╔═╡ 02a578df-5c60-47ad-8e75-967794f32668
begin
	single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
	s_bkg_S_MD = mean(single_arr[single_ring_mask_S_MD])
end

# ╔═╡ e80d2548-83e8-40be-be6c-7a0156cc61a8
begin
	alg_S_MD = Integrated(arr[mask_S_MD_3D])
	mass_s_md = score(s_bkg_S_MD, cal_insert_mean, pixel_size, ρ, alg_S_MD)
	# if mass_s_md < 0
	# 	mass_s_md = 0
	# end
	# mass_s_md
end

# ╔═╡ 15bee139-54fd-43ff-a761-8464f1d9d15b
md"""
## Low Density
"""

# ╔═╡ 706a167f-1cff-46ba-9142-5385c68a5f5f
begin
	mask_S_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_LD_3D[:, :, z] = mask_S_LD
	end
end;

# ╔═╡ abbdd826-41be-481c-87a6-b09864fc9d4f
md"""
#### Dilated mask
"""

# ╔═╡ 332b0566-9e06-4d6b-b606-68d2f65bddec
dilated_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D)))));

# ╔═╡ a3728f4f-fb4e-4056-9276-c55dd7f202c8
@bind o2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ 9c48c322-279f-4ec9-9fa7-f7fc0bfdc8ad
overlay_mask_plot(arr, dilated_mask_S_LD, o2, "dilated mask")

# ╔═╡ 97ae0d1b-0903-4096-bf62-660c0733d094
md"""
#### Ring (background) mask
"""

# ╔═╡ 01596aa1-cc2b-4c5c-86cb-8b1e624518ba
ring_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D))))) - dilate(dilate(dilate(dilate(mask_S_LD_3D))));

# ╔═╡ 627d6568-c714-418b-8640-0a58057dfbce
@bind o4 overlay_mask_bind(ring_mask_S_LD)

# ╔═╡ bf2fd448-0413-4fea-b742-cf3e3594cb4a
overlay_mask_plot(arr, ring_mask_S_LD, o4, "ring mask")

# ╔═╡ 9e061274-32f9-40d3-9501-bf2ef3df8d03
begin
	single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
	s_bkg_S_LD = mean(single_arr[single_ring_mask_S_LD])
end

# ╔═╡ a9c9dff8-8c55-4cb0-9059-69c76df0dbfb
begin
	alg_S_LD = Integrated(arr[mask_S_LD_3D])
	mass_s_ld = score(s_bkg_S_LD, cal_insert_mean, pixel_size, ρ, alg_S_LD)
	# if mass_s_ld < 0
	# 	mass_s_ld = 0
	# end
	# mass_s_ld
end

# ╔═╡ b9d9c618-45e2-4538-92ec-450bc96a2bbb
md"""
# Results
"""

# ╔═╡ c8a07078-ae5f-423f-abde-6ac495cf32b1
density_array = [0, 200, 400, 800]

# ╔═╡ 2c71a78a-0eb3-46d8-a625-9c886ca71450
inserts = [
	"Low Density",
	"Medium Density",
	"High Density"
]

# ╔═╡ ee64e181-80d8-4076-a2a5-42d2b1e98214
ground_truth_mass_large = [
	19.6,
	39.3,
	78.5
] # mg

# ╔═╡ 1b6e52cc-fdc1-48ce-bf0a-a3f729eeaf99
calculated_mass_large = [
	mass_l_ld,
	mass_l_md,
	mass_l_hd
]

# ╔═╡ 91739e5e-1e27-4f53-94ab-b63b2dd65dce
ground_truth_mass_medium = [
	4.2,
	8.5,
	17.0
]

# ╔═╡ 26380cc4-947f-4f22-9fc5-ac748abdac4a
calculated_mass_medium = [
	mass_m_ld,
	mass_m_md,
	mass_m_hd
]

# ╔═╡ 3933fc50-4f7e-4e71-9911-b9013377543e
ground_truth_mass_small = [
	0.2,
	0.3,
	0.6
]

# ╔═╡ 700d8ca0-3f41-498a-bd22-42a0324bfe55
calculated_mass_small = [
	mass_s_ld,
	mass_s_md,
	mass_s_hd
]

# ╔═╡ 64e72733-a646-415d-a4cd-c2497ed6cd33
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

# ╔═╡ 98edd7ca-a6c2-45b9-9f59-0a084316028f
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

# ╔═╡ 406c4f4f-aea9-45dd-8947-b4c503ec66ad
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

# ╔═╡ 9c95ca40-8ace-4228-8216-bea0cd99c084
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

# ╔═╡ 0bae0ce4-c0be-46ca-821f-a34bb7cd7e87
md"""
### Save Results
"""

# ╔═╡ 01b4f8ad-96af-4cfc-9b29-be3aa02f7b84
if ~isdir(string(cd(pwd, "..") , "/data/output/", VENDER, "4"))
	mkdir(string(cd(pwd, "..") , "/data/output/", VENDER, "4"))
end

# ╔═╡ 20724345-3b87-4298-ab0f-22357e8f5169
output_path = string(cd(pwd, "..") , "/data/output/", VENDER, "4", "/", scan, ".csv")

# ╔═╡ 68c98c61-11fc-4658-bed7-3491b46d063d
# CSV.write(output_path, df)

# ╔═╡ 7d5d27e8-15a3-419f-93dc-37134740ee29
md"""
### Save full df
"""

# ╔═╡ 38f890c1-6e5d-48c7-aec1-e3971ce2bc2b
dfs = []

# ╔═╡ e9eb7218-54d1-4131-a3c9-cc19e1104166
push!(dfs, df)

# ╔═╡ 1dc38db1-ad2b-436d-9ac5-4065bedacaf3
if length(dfs) == 10
	global new_df = vcat(dfs[1:10]...)
	output_path_new = string(cd(pwd, "..") , "/data/output/", VENDER, "4", "/full.csv")
	CSV.write(output_path_new, new_df)
end

# ╔═╡ Cell order:
# ╠═8f0a533f-fdb6-4c1f-bc76-6d1b0f15192c
# ╠═9cff9635-f4f4-4d59-81a6-92390c09526c
# ╟─f5090527-3372-48ea-8239-f51bf20fbba7
# ╟─383f5377-0435-4ea4-96d9-6413fcce1ee0
# ╠═34cdbadf-edfa-4708-be18-03e40f6e00fe
# ╠═8feec48a-9ea2-4aa5-96a1-544ef298a3ac
# ╠═1c7739dc-6615-4aab-b9fa-a4664c858070
# ╠═752cfdd4-f2cf-45fb-8f70-959d145264f6
# ╠═96fc1811-ab0f-41e6-8857-31b850a03199
# ╠═ff90198d-14e0-4404-813f-40dd35cc1f47
# ╠═b58a6353-8574-44bc-8312-603b9b30d672
# ╟─60cc7eea-8314-4f41-bf5d-8375a2c996e0
# ╟─555426d7-4584-41eb-a577-53a5a506e97a
# ╟─aa688492-bafb-4623-af4b-a744ca59335e
# ╠═d16e857c-654a-48fa-81dd-b53c7092cc86
# ╟─d6b6104a-df14-4cb0-96e6-e07be85255ad
# ╠═91243a6c-ce6f-4fc1-b28c-792e8c71dd53
# ╟─343815f2-6b33-4408-81eb-d4b1fddf60ac
# ╠═127344d7-38d6-4338-9ae8-edd92090b5cf
# ╟─fc3bbc09-e993-4f1d-838d-72041cb0f579
# ╟─865f42c6-ca98-498a-89a3-d4c9986237d3
# ╟─0a2563dd-d104-4f1b-b178-f77ca30e96b3
# ╟─db231239-5a21-402f-8eb5-b2e44e03b312
# ╠═ef56f1b2-150f-4608-944d-b5785a46e77c
# ╟─252d8c09-fcef-435e-9785-ec7ef607f108
# ╠═cf330459-0792-48fb-9346-f87b3019a0f6
# ╟─eb4494c0-003d-4967-85fa-b20d8e7f9763
# ╠═53f711d7-7ec2-41ea-b609-3146c01fa030
# ╠═95948083-8209-426e-ac12-2d7e925a6c30
# ╠═bf3853ea-0b3b-476f-b9b5-b969cbcc9827
# ╟─23e10eee-6022-4b31-bbc2-9ddbcbda0b52
# ╠═ed57792e-beb4-45af-a710-e0f6268ddac3
# ╠═f7c31a75-39c7-4be7-acbe-c36db4c252cb
# ╠═9fcdac1c-4307-4f27-a924-b8887406fcd7
# ╠═46112e23-ce69-4130-832d-a26e7d5195f9
# ╠═45f4a72e-a74a-4db2-b144-576b3f22b1b0
# ╠═e326a760-5795-4892-8dc7-88bfef009a52
# ╠═84efd8af-873e-448d-b70d-1832e824fa3f
# ╠═1666956f-1417-45ae-945c-7f58b98cfe9e
# ╠═ab09d0f1-f072-4f77-9889-a2beb852dc48
# ╠═f12cf6c5-9e82-4572-8321-fc848f1c866e
# ╟─9f5b80aa-afb7-4c0c-829f-30171360a87e
# ╠═3e7253c2-b6df-4dab-8405-08ed93c03f10
# ╠═1f40f903-4c8b-48bc-bb8a-fa16bc6fc297
# ╟─c7f433f7-c12f-41eb-be35-97997d5142b7
# ╠═32085083-902c-444b-9384-b79769c6acd2
# ╟─b0f14232-0be1-4a0e-b336-3cc3a8aaff6d
# ╠═0e31733e-a36f-4d70-970e-7b44accb2b99
# ╟─0d2e8726-4829-46ce-b85b-fc36b5a651dc
# ╠═75830069-109b-474c-9c4a-448ee6b19728
# ╟─3834b0e8-85c4-4cea-bb74-6350bb3e9e7d
# ╠═ccb3a71e-8669-4371-b37d-1d2e38166239
# ╟─afeb1614-5476-45db-aca2-bf3dbf9de6da
# ╠═73e0c981-e3b0-4a33-810d-04e6767428de
# ╠═6542d34b-89f8-42d9-9448-7d0503d35477
# ╠═c30c8fd9-93b5-404b-a1c4-5cdeda48845f
# ╠═8fbd3e9c-0195-4483-92c6-9a47e6e26a9b
# ╟─acfcb465-10e7-417a-9fe3-4231c86c5293
# ╠═6e794c8a-0341-4a0c-9d96-192f5b77a8ab
# ╟─fa2a052f-452a-4325-a428-8db343d560ec
# ╠═d9987a6f-b8d8-4ace-8b40-234d5dde2619
# ╟─6facfb53-1b86-4d66-9f86-0e6e1186c415
# ╠═338b6552-5845-4257-aa74-5fb3f8933239
# ╟─7741dcb3-e433-4612-aaf5-eabb2d8aa3aa
# ╠═23839ac8-1f79-4621-afd9-9ae3ff572c62
# ╟─dc71ee51-82fb-4367-ba3f-2ef6599c6808
# ╠═11a32a9f-d66a-46f2-bafe-dd032b9db1b0
# ╠═e88c694a-f3ae-4114-8a2c-1d20bba7f1bd
# ╠═c2aba951-9e13-40de-8a80-3e38ed800a01
# ╟─7ad75365-cb12-4db1-8354-956540204ee8
# ╠═5be373b4-7d52-446d-a39c-3aed2aeb762c
# ╟─2d04c8c8-b652-410b-8357-b123b488c579
# ╠═f02480d0-a373-4a32-bdf2-c40055ed2ab9
# ╟─7ac3b056-a9a8-44c0-a8c4-e5e89062f28c
# ╠═03a74421-7912-49af-92d6-56ba38b8b7ad
# ╟─eb213bf2-bbd4-48dc-9834-ca4d3b650bfd
# ╠═12e7a90b-27df-4cdb-99aa-c9aa8e1aec05
# ╠═fee22a56-26bb-4dd8-a500-e973a117dbeb
# ╠═ce557f39-f2d0-4b48-9e58-f232f3056b9b
# ╠═c1809d80-7c49-43e4-bff2-d6e14a7d9f97
# ╠═fa89795c-a6f5-4d9e-8d33-7cf42f4530d6
# ╟─0ca2ffa3-01a2-4484-a324-5d55e5609d10
# ╟─5779ee2f-5821-4475-aaba-c7aac7bb1d5e
# ╠═ab8473c4-85a9-4c66-83b1-6b8822afc855
# ╟─dca948cb-f537-49b6-bc69-771cae751e94
# ╠═821240b4-19e6-4214-8b1e-cfe69c574ae4
# ╟─9df3801d-8984-4eb9-bc82-46e9e7c458e2
# ╠═d78f713a-1c6c-48c5-b1f8-ec8d89a3931c
# ╟─1cbd7214-11e9-4266-8fa1-7363d46b2de2
# ╠═625402e5-5191-493a-b428-a4ee442a327d
# ╟─49ea60a3-61df-4b8d-97f0-5ad4c238b4ff
# ╠═828585e7-bc8d-40f4-a76a-624f2513db55
# ╠═dba099fe-7f13-4f61-b5ff-a4d65b7a19e2
# ╠═55740593-5f41-441e-b3b8-dcbfb83eceb5
# ╟─0cda89b7-c4c3-4b8c-b2f5-332a90f92831
# ╠═8f0af863-52ee-4e6d-8d70-e50432720501
# ╟─39019b75-c629-4e04-970b-2813af45cdcd
# ╠═412c307e-e323-4832-9815-98e475f4b5b3
# ╟─1a9af931-7f46-4c71-a753-8eb0b9bfc64b
# ╠═5d63180c-930c-464e-a8d0-3b2ea572f7cc
# ╟─231e9c7b-9d7c-48e5-a39e-b60e7c846889
# ╠═a79df362-4832-4f98-8bf9-e410df17bf71
# ╟─ddcde058-5c8b-40d9-8244-4fe99ec7c9a6
# ╠═ed112bd8-4af6-4522-a24b-47b39f3a7865
# ╠═5c777b1a-6b33-4b69-8797-6931a9d07a44
# ╠═0d22d7e4-8519-4652-b794-83756671aa87
# ╟─27e50773-d747-4451-8546-16d5c014f37a
# ╠═4168a9e9-6d44-401b-85f1-d80d6a2eeff8
# ╟─7a5c44ea-baf4-4b93-924e-f8624fe051ff
# ╠═fb8dfcff-5c15-4c30-b843-fa6db9e34e91
# ╟─3c6f8e72-a699-491e-b87d-0e0ee385d470
# ╠═02c01e40-d271-40eb-945f-60212dc1758e
# ╟─e15c09f2-0a07-4af6-9c17-44ba50393383
# ╠═76f15b22-a7f6-4d8d-9911-1c072e493589
# ╟─3723dee8-a78b-4c17-944d-c393c7ab92ad
# ╠═32956aa9-4338-43f2-bba7-d520716609ec
# ╠═3fe8e26e-0152-4869-b8ba-7e67d1e5d101
# ╠═c4252a6c-878a-493f-abc1-b8c56cda2319
# ╟─4ad89eca-ce15-43bd-bcd2-dd36ea759c7c
# ╟─70597737-ace4-4676-836f-44fdab22a798
# ╠═d5a8ba59-f167-48f3-954f-1b5a661456d5
# ╟─fa1c77f9-5f01-47b5-aef7-d885ddae8029
# ╠═799cec3a-953d-4fe2-bc06-351f5d6f4890
# ╟─c2c9bd13-f56b-44e7-97ed-be2e06d63ba3
# ╠═6f078eda-1fa3-4c9b-81d5-10909a9c9387
# ╠═fac87764-f188-4058-ae48-7aa36ac3b3bf
# ╟─c8bff4a0-a2a1-4e19-9d35-aceb00b365ff
# ╠═ad67881f-9085-46dc-a645-c3f111738c94
# ╟─0710d8ef-aac6-4d6d-983a-2d650502c180
# ╠═ba1db1c1-47c8-4c25-b650-67cb8c5c4e1b
# ╠═ff7731cc-5933-4312-919e-e4703254a622
# ╠═54e24ee8-a0f8-452a-ac2d-1379174fb56d
# ╟─82bbbd01-7a15-438c-b317-22b9511e8895
# ╠═61025c8f-d189-4ad0-bbe1-bb24e13083d4
# ╟─4ae034d0-b3fa-46ca-acf7-b93062bdb5a2
# ╠═6e544d17-f156-43c2-9255-58bd6b0c1888
# ╟─50333385-e60d-4dd2-92e1-bcd2d9007355
# ╠═c3f5fc15-60f5-457f-8353-9e10926d5a7d
# ╟─f66110c6-f501-4a55-8719-fc752b838de3
# ╠═2427ac81-bbb4-49ec-9897-d4fa59e4259b
# ╟─65722ef0-fbe9-4259-93f2-f5265e6ac3ba
# ╠═5c64cf04-7406-46cd-8650-05b8d7ab0964
# ╠═02a578df-5c60-47ad-8e75-967794f32668
# ╠═e80d2548-83e8-40be-be6c-7a0156cc61a8
# ╟─15bee139-54fd-43ff-a761-8464f1d9d15b
# ╠═706a167f-1cff-46ba-9142-5385c68a5f5f
# ╟─abbdd826-41be-481c-87a6-b09864fc9d4f
# ╠═332b0566-9e06-4d6b-b606-68d2f65bddec
# ╟─a3728f4f-fb4e-4056-9276-c55dd7f202c8
# ╠═9c48c322-279f-4ec9-9fa7-f7fc0bfdc8ad
# ╟─97ae0d1b-0903-4096-bf62-660c0733d094
# ╠═01596aa1-cc2b-4c5c-86cb-8b1e624518ba
# ╟─627d6568-c714-418b-8640-0a58057dfbce
# ╠═bf2fd448-0413-4fea-b742-cf3e3594cb4a
# ╠═9e061274-32f9-40d3-9501-bf2ef3df8d03
# ╠═a9c9dff8-8c55-4cb0-9059-69c76df0dbfb
# ╟─b9d9c618-45e2-4538-92ec-450bc96a2bbb
# ╠═c8a07078-ae5f-423f-abde-6ac495cf32b1
# ╠═2c71a78a-0eb3-46d8-a625-9c886ca71450
# ╠═ee64e181-80d8-4076-a2a5-42d2b1e98214
# ╠═1b6e52cc-fdc1-48ce-bf0a-a3f729eeaf99
# ╠═91739e5e-1e27-4f53-94ab-b63b2dd65dce
# ╠═26380cc4-947f-4f22-9fc5-ac748abdac4a
# ╠═3933fc50-4f7e-4e71-9911-b9013377543e
# ╠═700d8ca0-3f41-498a-bd22-42a0324bfe55
# ╠═64e72733-a646-415d-a4cd-c2497ed6cd33
# ╟─98edd7ca-a6c2-45b9-9f59-0a084316028f
# ╟─406c4f4f-aea9-45dd-8947-b4c503ec66ad
# ╠═9c95ca40-8ace-4228-8216-bea0cd99c084
# ╟─0bae0ce4-c0be-46ca-821f-a34bb7cd7e87
# ╠═01b4f8ad-96af-4cfc-9b29-be3aa02f7b84
# ╠═20724345-3b87-4298-ab0f-22357e8f5169
# ╠═68c98c61-11fc-4658-bed7-3491b46d063d
# ╟─7d5d27e8-15a3-419f-93dc-37134740ee29
# ╠═38f890c1-6e5d-48c7-aec1-e3971ce2bc2b
# ╠═e9eb7218-54d1-4131-a3c9-cc19e1104166
# ╠═1dc38db1-ad2b-436d-9ac5-4065bedacaf3

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

# ╔═╡ 3b42fc55-b563-4f46-ba14-caa099271238
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

# ╔═╡ 50f59532-7ce3-44ce-b5dc-0dd1d5d89f4e
TableOfContents()

# ╔═╡ 54c425e4-fc32-41b9-9a33-8238e9c8bdd5
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDER`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ 4dfeb579-9766-4346-bdfb-938e86755a77
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ f35f89a3-ad1e-43a8-abdb-922e00079199
begin
	SCAN_NUMBER = 5
	VENDER = "Canon_Aquilion_One_Vision"
	BASE_PATH = "/Users/daleblack/Google Drive/Datasets/"
end

# ╔═╡ 93f06e6c-2db5-429d-a5e9-fea11fbfa73e
root_path = string(BASE_PATH, VENDER)

# ╔═╡ 745bc0e3-25e3-49fc-b569-6cc8127e4816
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ fce70883-71fd-4f91-810a-d9b528bdceec
pth = dcm_path_list[SCAN_NUMBER]

# ╔═╡ c24c521d-fd4e-4ccb-9668-c47ccfc12489
pth

# ╔═╡ e64079a2-e22e-402e-87fe-3b5cd181b353
scan = basename(pth)

# ╔═╡ 907eaf46-9098-4c7c-95b4-be04420530f1
header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ a3772a66-58ca-4c28-a199-8bccb2599959
md"""
## Helper Functions
"""

# ╔═╡ 195ca2b0-349d-4b50-afd1-337d0a3783f5
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ 5d6b1ac8-2eb1-44f7-95f3-3eeef358367c
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ 5b51d2a4-9f8d-44d6-b7c1-8ecd54a76b78
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

# ╔═╡ 6db1b3b8-d00f-4bdf-a918-83aa13e81cd1
md"""
## Segment Heart
"""

# ╔═╡ 63f5f6fc-8abc-4e8d-907a-2fdfe198f24c
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ 13368f20-4cb0-467d-95e3-3209f1bb7a0f
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 56e0d020-3fec-4848-8966-e957a4dd7ef3
heatmap(masked_array[:, :, a], colormap=:grays)

# ╔═╡ b1d864db-caf0-47f0-9b09-7d1974e987b3
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 15]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig
end

# ╔═╡ b0ab4f49-eb99-4489-8aa2-d2fc1578072a
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig2
end

# ╔═╡ abfe9027-6ca1-45d6-8b63-5e5904232078
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 23]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig3
end

# ╔═╡ 61cca665-8138-431b-b460-844b0eb6f93e
md"""
## Segment Calcium Rod
"""

# ╔═╡ 79d660ef-546c-4827-a123-3a1cd8bcc053
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header);

# ╔═╡ c73a2ad7-cd60-4526-9be6-dfff04ee9a24
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=cal_rod_slice, show_value=true)

# ╔═╡ 46c7f573-c507-4783-8fd1-c16228052676
heatmap(transpose(calcium_image[:, :, c]), colormap=:grays)

# ╔═╡ 5ba8b642-a56d-449c-b294-7b8ca99387ce
md"""
## Segment Calcium Inserts
"""

# ╔═╡ 4b247654-ffaf-4f94-83f6-3642df8801ec
angle_factor = -4

# ╔═╡ 9c82cbed-eabc-4bba-8aaf-1dc2d77ad0c6
mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
		dcm_array, masked_array, header, slice_CCI, center_insert; angle_factor=angle_factor);

# ╔═╡ fd606156-edb7-4cc2-8449-edafe012499b
masks = mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD;

# ╔═╡ 078b6c88-8eb9-409c-82c7-8e4b0589df05
heatmap(masks, colormap=:grays)

# ╔═╡ 631197e0-d5ca-474d-962f-4cc39b35b408
md"""
## Calibration Prep
"""

# ╔═╡ 9f546951-40fc-46d6-a12e-109ae7c9f831
array_filtered = abs.(mapwindow(median, calcium_image[:, :, cal_rod_slice], (3, 3)));

# ╔═╡ 8d704a53-533e-4412-bab7-ee691e1ab063
bool_arr = array_filtered .> 0;

# ╔═╡ 4f653034-4f42-406d-a2d5-47366cc943db
bool_arr_erode = (((erode(erode(bool_arr)))));

# ╔═╡ 4dcc9b17-d766-4beb-bbf1-4cc8a89e4c54
heatmap(bool_arr, colormap=:grays)

# ╔═╡ aea20094-5820-4c7e-9c61-aadab359dbe5
heatmap(bool_arr_erode, colormap=:grays)

# ╔═╡ 66611c4b-f537-4a1d-b43f-cae8a1166925
c_img = calcium_image[:, :, cal_rod_slice-1:cal_rod_slice+1];

# ╔═╡ 2cbdfcc5-df27-478a-a0ef-2573b68399a7
begin
	mask_cal_3D = Array{Bool}(undef, size(c_img))
	for z in 1:size(c_img, 3)
		mask_cal_3D[:, :, z] = bool_arr_erode
	end
end;

# ╔═╡ c9b96229-a278-4560-8e34-4924c43d2aa3
# hist(c_img[mask_cal_3D])

# ╔═╡ 73b2992c-ab71-4840-beae-2da7edbd5cf8
cal_insert_mean = mean(c_img[mask_cal_3D])

# ╔═╡ 042eb1d0-2fc4-4a5c-8fb8-f464434267d4
# cal_insert_mean = quantile!(c_img[mask_cal_3D], 0.7)

# ╔═╡ eb9219a7-c8a5-49d4-b9c7-6dca7b560466
md"""
### Calibration Line
"""

# ╔═╡ 3c30cdf4-8dae-4ef2-b6fa-1c76b9b5306c
density_array_calc = [0, 200] # mg/cc

# ╔═╡ bdd5e150-fdb6-4623-ac99-cf5d3f14c9d4
intensity_array = [10, cal_insert_mean] # HU

# ╔═╡ 7debed73-00d5-4800-93eb-eb750c9ec0f6
df_cal = DataFrame(:density => density_array_calc, :intensity => intensity_array)

# ╔═╡ ab2adccb-6ad0-4958-8205-f6bcdcd8654c
linearRegressor = lm(@formula(intensity ~ density), df_cal);

# ╔═╡ 8e970c13-db21-49cd-80bb-26027d227bb9
linearFit = predict(linearRegressor)

# ╔═╡ 23f8007c-845e-442c-a42e-7f8eba592fb1
m = linearRegressor.model.pp.beta0[2]

# ╔═╡ 8f92d792-e05d-4b47-bd91-6d5f78f3ed02
b = linearRegressor.model.rr.mu[1]

# ╔═╡ 392fb68d-f6a7-4f3a-a27d-ee1609fcdacd
density(intensity) = (intensity - b) / m

# ╔═╡ c2537de9-7192-48f1-9694-d42184ab7a8d
intensity(ρ) = m*ρ + b

# ╔═╡ dae3264f-58e3-4b90-a9d0-ed1759e7c000
begin
	f = Figure()
	ax1 = Axis(f[1, 1])
	
	scatter!(density_array_calc, intensity_array)
	lines!(density_array_calc, linearFit, color = :red)
	ax1.title = "Calibration Line (Intensity vs Density)"
	ax1.ylabel = "Intensity (HU)"
	ax1.xlabel = "Density (mg/cm^3)"
	
	f
end

# ╔═╡ b0da69d3-7001-4b09-88d5-2285ad0ead97
md"""
# Score Large Inserts
"""

# ╔═╡ 59d6ec72-b0d1-4ffb-a6ff-2de59de0569b
arr = masked_array[:, :, slice_CCI-2:slice_CCI+2];

# ╔═╡ 4fa1675d-d363-4b2e-aad9-5265b39fb71b
single_arr = masked_array[:, :, slice_CCI];

# ╔═╡ 527dc20d-7a75-46be-b879-2c395f9590ad
md"""
## High Density
"""

# ╔═╡ ece6d9b9-ae58-4cad-9bee-d545ec62dab2
begin
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_HD_3D[:, :, z] = mask_L_HD
	end
end;

# ╔═╡ 2e7fe033-7cc8-48bb-884e-72b461a506ac
md"""
#### Dilated mask
"""

# ╔═╡ 18127ab6-04c0-46f5-81b5-c998584307ad
dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D));

# ╔═╡ 360f811c-1aa0-4673-a750-54ff808f8b91
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ a9ecbfda-d339-4919-9737-0a93f4b4dc47
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ 2e97924d-10f3-4ba8-85a4-7613fdfcae19
md"""
#### Ring (background) mask
"""

# ╔═╡ ee9dfda6-b17b-4842-8acb-1adb8f59cd32
ring_mask_L_HD = dilate(dilate(dilate(dilate(mask_L_HD_3D)))) - dilate(dilate(dilate(mask_L_HD_3D)));

# ╔═╡ 3ae3e83a-6344-495d-91a1-5aac553faea6
@bind g4 overlay_mask_bind(ring_mask_L_HD)

# ╔═╡ 12114e7e-1775-4360-b2bf-f645c7498d6b
overlay_mask_plot(arr, ring_mask_L_HD, g4, "ring mask")

# ╔═╡ 8d5c5136-9535-4f78-9d8c-1db380bcfc42
begin
	single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
	s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
end

# ╔═╡ dcde7971-9cea-4c87-b9fa-f5d8301d3483
pixel_size = DICOMUtils.get_pixel_size(header)

# ╔═╡ b1c06528-40f6-464c-9646-17700715d19e
cal_200 = intensity(200)

# ╔═╡ a6adbc8d-548e-4ea8-8d82-4aeeab3e3c88
ρ = 0.2

# ╔═╡ d1865bd6-f10d-4610-a4bf-88d2df99320b
begin
	alg_L_HD = Integrated(arr[mask_L_HD_3D])
	ρ_hd = 0.8 # mg/mm^3
	mass_l_hd = score(s_bkg_L_HD, cal_200, pixel_size, ρ, alg_L_HD)
end

# ╔═╡ 624b81e5-f4a0-4b1f-9291-e8d89ee3a84f
md"""
## Medium Density
"""

# ╔═╡ 2548d256-2aee-41a4-8ed5-39c04cad649f
begin
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_MD_3D[:, :, z] = mask_L_MD
	end
end;

# ╔═╡ 7644b902-ea63-49a8-8858-cdf4ee80d4c1
md"""
#### Dilated mask
"""

# ╔═╡ df5c5604-84bf-4eef-92ca-15346a559598
dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D));

# ╔═╡ 0258cc49-4858-477e-82a7-1bdd442a308d
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ 42bcfa62-0524-44cb-9be0-86aecd7e6bdb
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ 992558b3-e494-47c9-9f5e-96aa86a20ca8
md"""
#### Ring (background) mask
"""

# ╔═╡ 5a733f1d-1075-487f-b5e1-ee684c3622e1
ring_mask_L_MD = dilate(dilate(dilate(dilate(mask_L_MD_3D)))) - dilate(dilate(dilate(mask_L_MD_3D)));

# ╔═╡ 2de7602f-74a4-40ba-a635-70c0fc0b4f5a
@bind h4 overlay_mask_bind(ring_mask_L_MD)

# ╔═╡ d4c194f8-b35c-46d9-a23f-4a67fff0cc58
overlay_mask_plot(arr, ring_mask_L_MD, h4, "ring mask")

# ╔═╡ 50bd5a54-4933-48b3-9a57-c77015414375
begin
	single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
	s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
end

# ╔═╡ c79abc4f-a06a-4f19-892d-7823dcdabe17
S_Obj_MD = intensity(400)

# ╔═╡ faa89b46-7aa8-499d-99fd-d3db67618105
begin
	alg_L_MD = Integrated(arr[mask_L_MD_3D])
	ρ_md = 0.4 # mg/mm^3
	mass_l_md = score(s_bkg_L_MD, cal_200, pixel_size, ρ, alg_L_MD)
end

# ╔═╡ 406d12e0-184e-4310-bf2b-a58c436a2afd
md"""
## Low Density
"""

# ╔═╡ 25d78998-8467-460f-89dc-00c1bbaf3766
begin
	mask_L_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_LD_3D[:, :, z] = mask_L_LD
	end
end;

# ╔═╡ 5f6b05c4-b603-4f68-977a-8be96ccb5a16
md"""
#### Dilated mask
"""

# ╔═╡ b6b506ad-1e36-412b-b3cc-1d7c8ba36e90
dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D));

# ╔═╡ 01ecfe2a-5308-463e-8ac1-182d100b0cec
@bind i2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ 8f9c00af-029d-439d-af58-e4dee2715e57
overlay_mask_plot(arr, dilated_mask_L_LD, i2, "dilated mask")

# ╔═╡ f0d5665a-0ba6-43a9-b1d8-fd74fae88239
md"""
#### Ring (background) mask
"""

# ╔═╡ f96c7d1f-c217-4a5f-8dc3-3855600da4ca
ring_mask_L_LD = dilate(dilate(dilate(dilate(mask_L_LD_3D)))) - dilate(dilate(dilate(mask_L_LD_3D)));

# ╔═╡ 74367ae7-195b-431c-b79b-9bac76d6f17d
@bind i4 overlay_mask_bind(ring_mask_L_LD)

# ╔═╡ 0dcd632d-a522-4e1a-977e-0a457eb9a10a
overlay_mask_plot(arr, ring_mask_L_LD, i4, "ring mask")

# ╔═╡ b600f7c2-9856-45ef-a7d8-2beea8ec20ff
begin	
	single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
	s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
end

# ╔═╡ 13dc6c29-b339-4f72-bde2-0081a77582f5
S_Obj_LD = intensity(200)

# ╔═╡ b3559d64-98aa-4460-b235-eb8a6469c4b1
begin
	alg_L_LD = Integrated(arr[mask_L_LD_3D])
	ρ_ld = 0.2 # mg/mm^3
	mass_l_ld = score(s_bkg_L_LD, cal_200, pixel_size, ρ, alg_L_LD)
end

# ╔═╡ 105d9dc6-3d7b-42c9-94b2-bf74536ff81d
md"""
# Score Medium Inserts
"""

# ╔═╡ 4c66abca-7151-4f65-8587-5910a63b4cad
md"""
## High Density
"""

# ╔═╡ 1a5246a8-d329-4dd5-9333-a7dba84aedc2
begin
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_HD_3D[:, :, z] = mask_M_HD
	end
end;

# ╔═╡ c0ad23a0-d628-4e78-81c6-4fa2f5c1381b
md"""
#### Dilated mask
"""

# ╔═╡ a016b583-42c8-41bc-b5db-70c45da38f4b
dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ deb74dca-f4d1-475b-822d-2fa929486fe7
@bind j2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ 3b35a899-cfb2-4eef-bd0a-615eca4c50a6
overlay_mask_plot(arr, dilated_mask_M_HD, j2, "dilated mask")

# ╔═╡ 25e94145-46e7-476f-ad7e-a0771aef6622
md"""
#### Ring (background) mask
"""

# ╔═╡ 3ebc3e2a-6471-4e83-b85b-73fff382f1f8
ring_mask_M_HD = dilate(dilate(dilate(dilate(dilate(mask_M_HD_3D))))) - dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ f3f7196b-af43-4b2c-995d-9f8382ab2502
@bind j4 overlay_mask_bind(ring_mask_M_HD)

# ╔═╡ 8019c764-087f-4c5d-8c81-4901ce218ac0
overlay_mask_plot(arr, ring_mask_M_HD, j4, "ring mask")

# ╔═╡ 11441e3a-6cee-42ee-a921-ae49ef272585
begin
	single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
	s_bkg_M_HD = mean(single_arr[single_ring_mask_M_HD])
end

# ╔═╡ 1887a797-e6d4-4479-83ba-b2b1098cff1e
begin
	alg_M_HD = Integrated(arr[mask_M_HD_3D])
	mass_m_hd = score(s_bkg_M_HD, cal_200, pixel_size, ρ, alg_M_HD)
end

# ╔═╡ 6ea9406b-7a9f-4527-b913-d51f85bf9846
md"""
## Medium Density
"""

# ╔═╡ 6fd0f83b-0f3f-4776-8ae8-a7fcd371b942
begin
	mask_M_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_MD_3D[:, :, z] = mask_M_MD
	end
end;

# ╔═╡ afab8ba1-19c1-4c2a-930d-2d79719663f2
md"""
#### Dilated mask
"""

# ╔═╡ 20fbed68-0971-461b-8bd9-b6f996187814
dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))));

# ╔═╡ 5da040bd-84ff-4d99-9fdc-b64f2736a9d9
@bind k2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ 8a1a72d0-9845-4067-ae84-7e58d4136292
overlay_mask_plot(arr, dilated_mask_M_MD, k2, "dilated mask")

# ╔═╡ b596cfd7-a2d2-4bef-84d8-51b0f27f568d
md"""
#### Ring (background) mask
"""

# ╔═╡ 66244c90-0a0d-444e-9750-b8cebcb1e1c3
ring_mask_M_MD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))));

# ╔═╡ db1af3dc-cbb5-40b4-87dc-74d73238d3a0
@bind k4 overlay_mask_bind(ring_mask_M_MD)

# ╔═╡ d296d26b-a0f0-4864-852e-ff27ea8ee130
overlay_mask_plot(arr, ring_mask_M_MD, k4, "ring mask")

# ╔═╡ 52cd73f3-8427-4bd7-bdf5-7c2f05e530bd
begin
	single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
	s_bkg_M_MD = mean(single_arr[single_ring_mask_M_MD])
end

# ╔═╡ 1b5cf6d2-c452-4344-aed6-87717c2de3f6
begin
	alg_M_MD = Integrated(arr[mask_M_MD_3D])
	mass_m_md = score(s_bkg_M_MD, cal_200, pixel_size, ρ, alg_M_MD)
end

# ╔═╡ 54ff7f41-5fed-4718-ab3f-9d2abe054a93
md"""
## Low Density
"""

# ╔═╡ 646562df-dc94-4243-b461-b2ccb796fb28
begin
	mask_M_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_LD_3D[:, :, z] = mask_M_LD
	end
end;

# ╔═╡ 2b45bff2-6a80-4086-9e77-e480d2dd1d6f
md"""
#### Dilated mask
"""

# ╔═╡ c6218d03-9cef-4c7c-b63f-6bbcd6af7276
dilated_mask_M_LD = dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ 6eff82b1-c880-4cd8-89c0-8793281fd449
@bind l2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ f5170e26-f3e4-4a80-9fcc-4ae48323e0d9
overlay_mask_plot(arr, dilated_mask_M_LD, l2, "dilated mask")

# ╔═╡ d6c02531-e3cd-4a98-9491-48a76f797cbf
md"""
#### Ring (background) mask
"""

# ╔═╡ 58e247af-983b-4aee-9022-81434ca41756
ring_mask_M_LD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ a5fa60b1-d54c-4d23-8267-dd316b19c82d
@bind l4 overlay_mask_bind(ring_mask_M_LD)

# ╔═╡ 9911614d-cee2-4591-9031-9c1743e185ae
overlay_mask_plot(arr, ring_mask_M_LD, l4, "ring mask")

# ╔═╡ 459c1127-7e5e-4ec5-90d0-1e35222f331d
begin
	single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
	s_bkg_M_LD = mean(single_arr[single_ring_mask_M_LD])
end

# ╔═╡ c132038c-44f4-43e4-b1a6-f722cd67cea4
begin
	alg_M_LD = Integrated(arr[mask_M_LD_3D])
	mass_m_ld = score(s_bkg_M_LD, cal_200, pixel_size, ρ, alg_M_LD)
end

# ╔═╡ ebaa6376-bcc5-4357-b5ac-7cdb2b57bab3
md"""
# Score Small Inserts
"""

# ╔═╡ 1aacfcee-ec62-439e-aac3-4ff4fd44937d
md"""
## High Density
"""

# ╔═╡ 0a66885f-fa1b-444a-b873-ee8342e3e32b
begin
	mask_S_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_HD_3D[:, :, z] = mask_S_HD
	end
end;

# ╔═╡ e63cf220-785e-470a-913b-45d68b5777e7
md"""
#### Dilated mask
"""

# ╔═╡ 03de486e-26e4-4484-9a83-25f6323a849e
dilated_mask_S_HD = dilate(dilate(dilate(dilate(dilate((mask_S_HD_3D))))));

# ╔═╡ fa4f7d96-87c6-4d3b-b068-c3e50da246d4
@bind m2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ f171da23-935b-45a3-a626-b8bbd9c861e8
overlay_mask_plot(arr, dilated_mask_S_HD, m2, "dilated mask")

# ╔═╡ 708c1f02-3196-4f43-b17b-1fcbb3bc9074
md"""
#### Ring (background) mask
"""

# ╔═╡ 840a2655-58f2-480e-a77b-89ef7170467f
ring_mask_S_HD = dilate(dilate(dilate(dilate(dilate(mask_S_HD_3D))))) - dilate(dilate(dilate(dilate(mask_S_HD_3D))));

# ╔═╡ 98b62443-8895-47dd-a0ed-3cff043ed72a
@bind m4 overlay_mask_bind(ring_mask_S_HD)

# ╔═╡ 4370c423-65dd-4713-a6bf-826e01b433e8
overlay_mask_plot(arr, ring_mask_S_HD, m4, "ring mask")

# ╔═╡ d1f67ff3-fd07-451b-bf9d-ab9149bfb33f
begin
	single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
	s_bkg_S_HD = mean(single_arr[single_ring_mask_S_HD])
end

# ╔═╡ 8ab680b0-cb4d-4ce0-b6b0-d0838de4a667
begin
	alg_S_HD = Integrated(arr[mask_S_HD_3D])
	mass_s_hd = score(s_bkg_S_HD, cal_200, pixel_size, ρ, alg_S_HD)
end

# ╔═╡ aa0ac9cd-a5f2-4585-8e22-d607953789f3
md"""
## Medium Density
"""

# ╔═╡ 1a3dc86c-7abe-4395-bcb3-acb24aa005d5
begin
	mask_S_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_MD_3D[:, :, z] = mask_S_MD
	end
end;

# ╔═╡ e6ab3d73-081b-4df3-a3dd-14e29e1161a5
md"""
#### Dilated mask
"""

# ╔═╡ e5262767-a320-4839-8d7c-1dd18cb33a6e
dilated_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D)))));

# ╔═╡ 2c4af320-737d-44fa-9584-1afc2717afe4
@bind n2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ 726d3a60-5cbf-4c56-93f7-f3c56ed9516a
overlay_mask_plot(arr, dilated_mask_S_MD, n2, "dilated mask")

# ╔═╡ 085a86bf-8f19-4263-bf76-9864c5b73241
md"""
#### Ring (background) mask
"""

# ╔═╡ 2c5dd97b-4779-484f-a6aa-302b31140e58
ring_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D))))) - dilate(dilate(dilate(dilate(mask_S_MD_3D))));

# ╔═╡ d6f3132e-9fd3-4acc-9f40-92dcd72a1eef
@bind n4 overlay_mask_bind(ring_mask_S_MD)

# ╔═╡ dc5f2d6d-64ae-4e0b-9237-2160874e4d61
overlay_mask_plot(arr, ring_mask_S_MD, n4, "ring mask")

# ╔═╡ 30171650-6e45-4452-9bc6-90016158008f
begin
	single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
	s_bkg_S_MD = mean(single_arr[single_ring_mask_S_MD])
end

# ╔═╡ 888bd1ec-d980-4d7e-b8b7-b6ba1e16ad75
begin
	alg_S_MD = Integrated(arr[mask_S_MD_3D])
	mass_s_md = score(s_bkg_S_MD, cal_200, pixel_size, ρ, alg_S_MD)
end

# ╔═╡ 90e3bd47-e010-4e0d-b760-77d6abe0c0a8
md"""
## Low Density
"""

# ╔═╡ d631ab8d-af05-4a9c-806d-9f3f4471ef37
begin
	mask_S_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_LD_3D[:, :, z] = mask_S_LD
	end
end;

# ╔═╡ 402ecc18-781c-4277-b1bf-91b8c2c31d81
md"""
#### Dilated mask
"""

# ╔═╡ ae7c0eda-91bb-4af1-80af-055fd59fd1e8
dilated_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D)))));

# ╔═╡ c025f39c-80a8-4fb5-bdcc-2784d17a2a98
@bind o2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ d60825bb-f882-42d6-b157-205afd4af59c
overlay_mask_plot(arr, dilated_mask_S_LD, o2, "dilated mask")

# ╔═╡ dd8e989b-e45b-4f49-81cc-fcd3020fc96c
md"""
#### Ring (background) mask
"""

# ╔═╡ 1ff96b90-b6bf-41be-8b00-d4cd291cbe68
ring_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D))))) - dilate(dilate(dilate(dilate(mask_S_LD_3D))));

# ╔═╡ f3c63244-7048-48c6-852c-92721454314a
@bind o4 overlay_mask_bind(ring_mask_S_LD)

# ╔═╡ cdb53366-6688-4174-9741-925699920a42
overlay_mask_plot(arr, ring_mask_S_LD, o4, "ring mask")

# ╔═╡ d7da119a-9298-41a3-9856-7cb59490b000
begin
	single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
	s_bkg_S_LD = mean(single_arr[single_ring_mask_S_LD])
end

# ╔═╡ 60f1f3ef-029b-475e-a8de-8c40a4dc934a
begin
	alg_S_LD = Integrated(arr[mask_S_LD_3D])
	mass_s_ld = score(s_bkg_S_LD, cal_200, pixel_size, ρ, alg_S_LD)
end

# ╔═╡ 48ab0573-9516-4ce1-a790-eeb830e9f0f1
md"""
# Results
"""

# ╔═╡ 54bc62b2-d58d-4e97-a454-34aa0089d26d
density_array = [0, 200, 400, 800]

# ╔═╡ 63ab0548-dd9c-40d0-89c2-134bf5026d19
inserts = [
	"Low Density",
	"Medium Density",
	"High Density"
]

# ╔═╡ c7cf8c03-a63b-47e3-9cc7-785ba6b7feb0
ground_truth_mass_large = [
	19.6,
	39.3,
	78.5
] # mg

# ╔═╡ c04c5201-6004-41fb-80ad-b1d0f890e1fc
calculated_mass_large = [
	mass_l_ld,
	mass_l_md,
	mass_l_hd
]

# ╔═╡ dd9e78c4-ec0e-4946-b9f9-08743cbf3e56
ground_truth_mass_medium = [
	4.2,
	8.5,
	17.0
]

# ╔═╡ 29338b99-50df-4ca1-8d9a-d0553f0fb08f
calculated_mass_medium = [
	mass_m_ld,
	mass_m_md,
	mass_m_hd
]

# ╔═╡ 5ac5ae9b-a2e3-4d46-9e90-4ee56eb82a18
ground_truth_mass_small = [
	0.2,
	0.3,
	0.6
]

# ╔═╡ a9a1ed89-8732-48f7-a573-72e40a05e2fd
calculated_mass_small = [
	mass_s_ld,
	mass_s_md,
	mass_s_hd
]

# ╔═╡ b4e6ea6e-77ec-44f2-b4e4-56c9312b434f
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

# ╔═╡ 8c02c70b-2f71-48f6-8513-b6d1d3bb469b
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

# ╔═╡ ee8b679f-d04e-4a23-876a-15b9d1c3946f
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

# ╔═╡ 1f898cf8-db16-47d3-8ec7-7c380b52aa27
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

# ╔═╡ 16712d6e-5b8b-400f-8ec1-de3fad245bf2
md"""
### Save Results
"""

# ╔═╡ 67325f1c-9667-4a03-bb72-9f2f6c28f498
if ~isdir(string(cd(pwd, "..") , "/data/output/", VENDER, "2"))
	mkdir(string(cd(pwd, "..") , "/data/output/", VENDER, "2"))
end

# ╔═╡ 2e5d8342-036c-4cc1-9bbe-4792f244a8aa
output_path = string(cd(pwd, "..") , "/data/output/", VENDER, "2", "/", scan, ".csv")

# ╔═╡ ea6df367-85d3-475e-bc46-c3a7b1ed5358
# CSV.write(output_path, df)

# ╔═╡ b26f05d3-782b-493b-ba79-60f07a2d0ae3
md"""
### Save full df
"""

# ╔═╡ 40e72956-2ee3-4558-82a4-bf6fb8422241
dfs = []

# ╔═╡ 733d3b31-62e7-46ef-88b4-192ade9f1c74
push!(dfs, df)

# ╔═╡ 5183af49-1598-41ea-badf-8c3a8ca5b28b
if length(dfs) == 10
	global new_df = vcat(dfs[1:10]...)
	output_path_new = string(cd(pwd, "..") , "/data/output/", VENDER, "2", "/full.csv")
	CSV.write(output_path_new, new_df)
end

# ╔═╡ Cell order:
# ╠═3b42fc55-b563-4f46-ba14-caa099271238
# ╠═50f59532-7ce3-44ce-b5dc-0dd1d5d89f4e
# ╟─54c425e4-fc32-41b9-9a33-8238e9c8bdd5
# ╟─4dfeb579-9766-4346-bdfb-938e86755a77
# ╠═f35f89a3-ad1e-43a8-abdb-922e00079199
# ╠═93f06e6c-2db5-429d-a5e9-fea11fbfa73e
# ╠═745bc0e3-25e3-49fc-b569-6cc8127e4816
# ╠═fce70883-71fd-4f91-810a-d9b528bdceec
# ╠═c24c521d-fd4e-4ccb-9668-c47ccfc12489
# ╠═e64079a2-e22e-402e-87fe-3b5cd181b353
# ╠═907eaf46-9098-4c7c-95b4-be04420530f1
# ╟─a3772a66-58ca-4c28-a199-8bccb2599959
# ╟─195ca2b0-349d-4b50-afd1-337d0a3783f5
# ╟─5d6b1ac8-2eb1-44f7-95f3-3eeef358367c
# ╟─5b51d2a4-9f8d-44d6-b7c1-8ecd54a76b78
# ╟─6db1b3b8-d00f-4bdf-a918-83aa13e81cd1
# ╠═63f5f6fc-8abc-4e8d-907a-2fdfe198f24c
# ╟─13368f20-4cb0-467d-95e3-3209f1bb7a0f
# ╠═56e0d020-3fec-4848-8966-e957a4dd7ef3
# ╟─b1d864db-caf0-47f0-9b09-7d1974e987b3
# ╟─b0ab4f49-eb99-4489-8aa2-d2fc1578072a
# ╟─abfe9027-6ca1-45d6-8b63-5e5904232078
# ╟─61cca665-8138-431b-b460-844b0eb6f93e
# ╠═79d660ef-546c-4827-a123-3a1cd8bcc053
# ╟─c73a2ad7-cd60-4526-9be6-dfff04ee9a24
# ╠═46c7f573-c507-4783-8fd1-c16228052676
# ╟─5ba8b642-a56d-449c-b294-7b8ca99387ce
# ╠═4b247654-ffaf-4f94-83f6-3642df8801ec
# ╠═9c82cbed-eabc-4bba-8aaf-1dc2d77ad0c6
# ╠═fd606156-edb7-4cc2-8449-edafe012499b
# ╠═078b6c88-8eb9-409c-82c7-8e4b0589df05
# ╟─631197e0-d5ca-474d-962f-4cc39b35b408
# ╠═9f546951-40fc-46d6-a12e-109ae7c9f831
# ╠═8d704a53-533e-4412-bab7-ee691e1ab063
# ╠═4f653034-4f42-406d-a2d5-47366cc943db
# ╠═4dcc9b17-d766-4beb-bbf1-4cc8a89e4c54
# ╠═aea20094-5820-4c7e-9c61-aadab359dbe5
# ╠═66611c4b-f537-4a1d-b43f-cae8a1166925
# ╠═2cbdfcc5-df27-478a-a0ef-2573b68399a7
# ╠═c9b96229-a278-4560-8e34-4924c43d2aa3
# ╠═73b2992c-ab71-4840-beae-2da7edbd5cf8
# ╠═042eb1d0-2fc4-4a5c-8fb8-f464434267d4
# ╟─eb9219a7-c8a5-49d4-b9c7-6dca7b560466
# ╠═3c30cdf4-8dae-4ef2-b6fa-1c76b9b5306c
# ╠═bdd5e150-fdb6-4623-ac99-cf5d3f14c9d4
# ╠═7debed73-00d5-4800-93eb-eb750c9ec0f6
# ╠═ab2adccb-6ad0-4958-8205-f6bcdcd8654c
# ╠═8e970c13-db21-49cd-80bb-26027d227bb9
# ╠═23f8007c-845e-442c-a42e-7f8eba592fb1
# ╠═8f92d792-e05d-4b47-bd91-6d5f78f3ed02
# ╠═392fb68d-f6a7-4f3a-a27d-ee1609fcdacd
# ╠═c2537de9-7192-48f1-9694-d42184ab7a8d
# ╟─dae3264f-58e3-4b90-a9d0-ed1759e7c000
# ╟─b0da69d3-7001-4b09-88d5-2285ad0ead97
# ╠═59d6ec72-b0d1-4ffb-a6ff-2de59de0569b
# ╠═4fa1675d-d363-4b2e-aad9-5265b39fb71b
# ╟─527dc20d-7a75-46be-b879-2c395f9590ad
# ╠═ece6d9b9-ae58-4cad-9bee-d545ec62dab2
# ╟─2e7fe033-7cc8-48bb-884e-72b461a506ac
# ╠═18127ab6-04c0-46f5-81b5-c998584307ad
# ╟─360f811c-1aa0-4673-a750-54ff808f8b91
# ╠═a9ecbfda-d339-4919-9737-0a93f4b4dc47
# ╟─2e97924d-10f3-4ba8-85a4-7613fdfcae19
# ╠═ee9dfda6-b17b-4842-8acb-1adb8f59cd32
# ╟─3ae3e83a-6344-495d-91a1-5aac553faea6
# ╠═12114e7e-1775-4360-b2bf-f645c7498d6b
# ╠═8d5c5136-9535-4f78-9d8c-1db380bcfc42
# ╠═dcde7971-9cea-4c87-b9fa-f5d8301d3483
# ╠═b1c06528-40f6-464c-9646-17700715d19e
# ╠═a6adbc8d-548e-4ea8-8d82-4aeeab3e3c88
# ╠═d1865bd6-f10d-4610-a4bf-88d2df99320b
# ╟─624b81e5-f4a0-4b1f-9291-e8d89ee3a84f
# ╠═2548d256-2aee-41a4-8ed5-39c04cad649f
# ╟─7644b902-ea63-49a8-8858-cdf4ee80d4c1
# ╠═df5c5604-84bf-4eef-92ca-15346a559598
# ╟─0258cc49-4858-477e-82a7-1bdd442a308d
# ╠═42bcfa62-0524-44cb-9be0-86aecd7e6bdb
# ╟─992558b3-e494-47c9-9f5e-96aa86a20ca8
# ╠═5a733f1d-1075-487f-b5e1-ee684c3622e1
# ╟─2de7602f-74a4-40ba-a635-70c0fc0b4f5a
# ╠═d4c194f8-b35c-46d9-a23f-4a67fff0cc58
# ╠═50bd5a54-4933-48b3-9a57-c77015414375
# ╠═c79abc4f-a06a-4f19-892d-7823dcdabe17
# ╠═faa89b46-7aa8-499d-99fd-d3db67618105
# ╟─406d12e0-184e-4310-bf2b-a58c436a2afd
# ╠═25d78998-8467-460f-89dc-00c1bbaf3766
# ╟─5f6b05c4-b603-4f68-977a-8be96ccb5a16
# ╠═b6b506ad-1e36-412b-b3cc-1d7c8ba36e90
# ╟─01ecfe2a-5308-463e-8ac1-182d100b0cec
# ╠═8f9c00af-029d-439d-af58-e4dee2715e57
# ╟─f0d5665a-0ba6-43a9-b1d8-fd74fae88239
# ╠═f96c7d1f-c217-4a5f-8dc3-3855600da4ca
# ╟─74367ae7-195b-431c-b79b-9bac76d6f17d
# ╠═0dcd632d-a522-4e1a-977e-0a457eb9a10a
# ╠═b600f7c2-9856-45ef-a7d8-2beea8ec20ff
# ╠═13dc6c29-b339-4f72-bde2-0081a77582f5
# ╠═b3559d64-98aa-4460-b235-eb8a6469c4b1
# ╟─105d9dc6-3d7b-42c9-94b2-bf74536ff81d
# ╟─4c66abca-7151-4f65-8587-5910a63b4cad
# ╠═1a5246a8-d329-4dd5-9333-a7dba84aedc2
# ╟─c0ad23a0-d628-4e78-81c6-4fa2f5c1381b
# ╠═a016b583-42c8-41bc-b5db-70c45da38f4b
# ╟─deb74dca-f4d1-475b-822d-2fa929486fe7
# ╠═3b35a899-cfb2-4eef-bd0a-615eca4c50a6
# ╟─25e94145-46e7-476f-ad7e-a0771aef6622
# ╠═3ebc3e2a-6471-4e83-b85b-73fff382f1f8
# ╟─f3f7196b-af43-4b2c-995d-9f8382ab2502
# ╠═8019c764-087f-4c5d-8c81-4901ce218ac0
# ╠═11441e3a-6cee-42ee-a921-ae49ef272585
# ╠═1887a797-e6d4-4479-83ba-b2b1098cff1e
# ╟─6ea9406b-7a9f-4527-b913-d51f85bf9846
# ╠═6fd0f83b-0f3f-4776-8ae8-a7fcd371b942
# ╟─afab8ba1-19c1-4c2a-930d-2d79719663f2
# ╠═20fbed68-0971-461b-8bd9-b6f996187814
# ╟─5da040bd-84ff-4d99-9fdc-b64f2736a9d9
# ╠═8a1a72d0-9845-4067-ae84-7e58d4136292
# ╟─b596cfd7-a2d2-4bef-84d8-51b0f27f568d
# ╠═66244c90-0a0d-444e-9750-b8cebcb1e1c3
# ╟─db1af3dc-cbb5-40b4-87dc-74d73238d3a0
# ╠═d296d26b-a0f0-4864-852e-ff27ea8ee130
# ╠═52cd73f3-8427-4bd7-bdf5-7c2f05e530bd
# ╠═1b5cf6d2-c452-4344-aed6-87717c2de3f6
# ╟─54ff7f41-5fed-4718-ab3f-9d2abe054a93
# ╠═646562df-dc94-4243-b461-b2ccb796fb28
# ╟─2b45bff2-6a80-4086-9e77-e480d2dd1d6f
# ╠═c6218d03-9cef-4c7c-b63f-6bbcd6af7276
# ╟─6eff82b1-c880-4cd8-89c0-8793281fd449
# ╠═f5170e26-f3e4-4a80-9fcc-4ae48323e0d9
# ╟─d6c02531-e3cd-4a98-9491-48a76f797cbf
# ╠═58e247af-983b-4aee-9022-81434ca41756
# ╟─a5fa60b1-d54c-4d23-8267-dd316b19c82d
# ╠═9911614d-cee2-4591-9031-9c1743e185ae
# ╠═459c1127-7e5e-4ec5-90d0-1e35222f331d
# ╠═c132038c-44f4-43e4-b1a6-f722cd67cea4
# ╟─ebaa6376-bcc5-4357-b5ac-7cdb2b57bab3
# ╟─1aacfcee-ec62-439e-aac3-4ff4fd44937d
# ╠═0a66885f-fa1b-444a-b873-ee8342e3e32b
# ╟─e63cf220-785e-470a-913b-45d68b5777e7
# ╠═03de486e-26e4-4484-9a83-25f6323a849e
# ╟─fa4f7d96-87c6-4d3b-b068-c3e50da246d4
# ╠═f171da23-935b-45a3-a626-b8bbd9c861e8
# ╟─708c1f02-3196-4f43-b17b-1fcbb3bc9074
# ╠═840a2655-58f2-480e-a77b-89ef7170467f
# ╟─98b62443-8895-47dd-a0ed-3cff043ed72a
# ╠═4370c423-65dd-4713-a6bf-826e01b433e8
# ╠═d1f67ff3-fd07-451b-bf9d-ab9149bfb33f
# ╠═8ab680b0-cb4d-4ce0-b6b0-d0838de4a667
# ╟─aa0ac9cd-a5f2-4585-8e22-d607953789f3
# ╠═1a3dc86c-7abe-4395-bcb3-acb24aa005d5
# ╟─e6ab3d73-081b-4df3-a3dd-14e29e1161a5
# ╠═e5262767-a320-4839-8d7c-1dd18cb33a6e
# ╟─2c4af320-737d-44fa-9584-1afc2717afe4
# ╠═726d3a60-5cbf-4c56-93f7-f3c56ed9516a
# ╟─085a86bf-8f19-4263-bf76-9864c5b73241
# ╠═2c5dd97b-4779-484f-a6aa-302b31140e58
# ╟─d6f3132e-9fd3-4acc-9f40-92dcd72a1eef
# ╠═dc5f2d6d-64ae-4e0b-9237-2160874e4d61
# ╠═30171650-6e45-4452-9bc6-90016158008f
# ╠═888bd1ec-d980-4d7e-b8b7-b6ba1e16ad75
# ╟─90e3bd47-e010-4e0d-b760-77d6abe0c0a8
# ╠═d631ab8d-af05-4a9c-806d-9f3f4471ef37
# ╟─402ecc18-781c-4277-b1bf-91b8c2c31d81
# ╠═ae7c0eda-91bb-4af1-80af-055fd59fd1e8
# ╟─c025f39c-80a8-4fb5-bdcc-2784d17a2a98
# ╠═d60825bb-f882-42d6-b157-205afd4af59c
# ╟─dd8e989b-e45b-4f49-81cc-fcd3020fc96c
# ╠═1ff96b90-b6bf-41be-8b00-d4cd291cbe68
# ╟─f3c63244-7048-48c6-852c-92721454314a
# ╠═cdb53366-6688-4174-9741-925699920a42
# ╠═d7da119a-9298-41a3-9856-7cb59490b000
# ╠═60f1f3ef-029b-475e-a8de-8c40a4dc934a
# ╟─48ab0573-9516-4ce1-a790-eeb830e9f0f1
# ╠═54bc62b2-d58d-4e97-a454-34aa0089d26d
# ╠═63ab0548-dd9c-40d0-89c2-134bf5026d19
# ╠═c7cf8c03-a63b-47e3-9cc7-785ba6b7feb0
# ╠═c04c5201-6004-41fb-80ad-b1d0f890e1fc
# ╠═dd9e78c4-ec0e-4946-b9f9-08743cbf3e56
# ╠═29338b99-50df-4ca1-8d9a-d0553f0fb08f
# ╠═5ac5ae9b-a2e3-4d46-9e90-4ee56eb82a18
# ╠═a9a1ed89-8732-48f7-a573-72e40a05e2fd
# ╠═b4e6ea6e-77ec-44f2-b4e4-56c9312b434f
# ╟─8c02c70b-2f71-48f6-8513-b6d1d3bb469b
# ╟─ee8b679f-d04e-4a23-876a-15b9d1c3946f
# ╟─1f898cf8-db16-47d3-8ec7-7c380b52aa27
# ╟─16712d6e-5b8b-400f-8ec1-de3fad245bf2
# ╠═67325f1c-9667-4a03-bb72-9f2f6c28f498
# ╠═2e5d8342-036c-4cc1-9bbe-4792f244a8aa
# ╠═ea6df367-85d3-475e-bc46-c3a7b1ed5358
# ╟─b26f05d3-782b-493b-ba79-60f07a2d0ae3
# ╠═40e72956-2ee3-4558-82a4-bf6fb8422241
# ╠═733d3b31-62e7-46ef-88b4-192ade9f1c74
# ╠═5183af49-1598-41ea-badf-8c3a8ca5b28b

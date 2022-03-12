### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 07787e48-5cf4-4583-b605-59054de37632
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

# ╔═╡ 59ca3f65-5e88-44e9-98f1-381036237c9b
TableOfContents()

# ╔═╡ 5706836e-df8d-4a03-8610-dde49bd357dc
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDER`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ 05b7448d-7377-4a4a-ab8a-cd9e1b3d47c5


# ╔═╡ 89785633-a679-4012-a79d-c76c61c742bb
begin
	# SCAN_NUMBER = 8
	VENDER = "Canon_Aquilion_One_Vision"
	BASE_PATH = "/Users/daleblack/Google Drive/Datasets/"
end

# ╔═╡ c9ab625d-ee86-4cc2-ac23-bd8b573cf8e2
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ 11ea321b-4be9-420f-9c9c-01b6450580b4
scans = collect(1:10)

# ╔═╡ e7cabee9-2bea-4814-a53f-b96eadc9d41e
for s in scans
	global SCAN_NUMBER = 9
	global root_path = string(BASE_PATH, VENDER)
	global dcm_path_list = dcm_list_builder(root_path)
	global pth = dcm_path_list[SCAN_NUMBER]
	global scan = basename(pth)
	global header
	global dcm_array
	global slice_thick_ori1
	header, dcm_array, slice_thick_ori1 = dcm_reader(pth)
end

# ╔═╡ 0d3eb83a-0a11-4eb1-8e28-9b98e9316437
# root_path = string(BASE_PATH, VENDER)

# ╔═╡ d431e54b-55b4-40ea-8aab-40cb7652374d
# dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ 148eab5f-3263-4f19-a7b5-f289d631d119
# pth = dcm_path_list[SCAN_NUMBER]

# ╔═╡ a0d58e6e-9cf7-4802-bc87-27259a310840
# scan = basename(pth)

# ╔═╡ 45c84e2b-5439-4ad3-b249-248abb0d0afb
# header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ daf9758b-2884-46ea-807f-8e17c96f4fda
md"""
## Helper Functions
"""

# ╔═╡ be2e6d55-3ad0-4749-9668-cfb275470e22
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ 59afc75e-fa7b-4b0b-90bd-92bf9550a30f
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ f55a30b7-bbfc-413e-b1ce-a8a1fe7cc1a8
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

# ╔═╡ 4aeb5503-a295-4f2b-87e1-578230e01cde
md"""
## Segment Heart
"""

# ╔═╡ ed209978-0ef0-4f7b-9b17-89e35fd7f4cd
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ 4222f5fe-78b5-4bf8-a42b-0b32379cda4d
md"""
## Segment Calcium Rod
"""

# ╔═╡ f0582426-8c3e-4f47-b8b0-af18920dc7e3
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header);

# ╔═╡ a5577cd7-8d70-4730-aefc-4d5dcb439734
md"""
## Calibration Prep
"""

# ╔═╡ 6809886f-55fd-4899-b7a0-0f2076422643
array_filtered = abs.(mapwindow(median, calcium_image[:, :, cal_rod_slice], (3, 3)));

# ╔═╡ 64210c98-98fb-4413-879c-9bedd2f5060e
bool_arr = array_filtered .> 0;

# ╔═╡ fdb968c5-3a3b-43eb-84ba-42617dea9915
bool_arr_erode = (((erode(erode(bool_arr)))));

# ╔═╡ f403c1fa-35ab-4526-b105-a8056da88ae0
c_img = calcium_image[:, :, cal_rod_slice-1:cal_rod_slice+1];

# ╔═╡ fd90f126-0912-4026-b024-f7cf3968f896
begin
	mask_cal_3D = Array{Bool}(undef, size(c_img))
	for z in 1:size(c_img, 3)
		mask_cal_3D[:, :, z] = bool_arr_erode
	end
end;

# ╔═╡ ef899480-a460-49c8-abce-2ed81922385b
# hist(c_img[mask_cal_3D])

# ╔═╡ 5f305c11-b050-4d03-a7af-d7b76da5d3a8
cal_insert_mean2 = mean(c_img[mask_cal_3D])

# ╔═╡ 27d61bb3-00c2-45de-af02-0ca8108fed08
cal_insert_mean = quantile!(c_img[mask_cal_3D], 0.7)

# ╔═╡ 810a1af6-2d85-46d5-9cb1-e316225e201b
md"""
## Segment Calcium Inserts
"""

# ╔═╡ 10f3d200-f876-4d66-9cfa-8303036310f4
angles = collect(-10:2:10)

# ╔═╡ 041fdc88-5fbc-4845-8514-622fa43d36ab
md"""
**Check Best Segmentation**
"""

# ╔═╡ 122d6e0b-3b00-4056-b68b-8452c9364941
begin
	RMSE_Dict = Dict()
	for angle in angles
		mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
		dcm_array, masked_array, header, slice_CCI, center_insert; angle_factor=angle)
	
		arr = masked_array[:, :, slice_CCI-3:slice_CCI+3]
		single_arr = masked_array[:, :, slice_CCI]
		pixel_size = DICOMUtils.get_pixel_size(header)
		ρ = 0.2 # mg/mm^3
		# Score Large InsertS
		## High Density
		mask_L_HD_3D = Array{Bool}(undef, size(arr))
		for z in 1:size(arr, 3)
			mask_L_HD_3D[:, :, z] = mask_L_HD
		end
		dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D))
		ring_mask_L_HD = dilate(dilate(dilate(dilate(mask_L_HD_3D)))) - dilate(dilate(dilate(mask_L_HD_3D)))
		single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
		s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
		alg_L_HD = Integrated(arr[mask_L_HD_3D])
		mass_l_hd = score(s_bkg_L_HD, cal_insert_mean, pixel_size, ρ, alg_L_HD)
	
		## Medium Density
		mask_L_MD_3D = Array{Bool}(undef, size(arr))
		for z in 1:size(arr, 3)
			mask_L_MD_3D[:, :, z] = mask_L_MD
		end
		dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D))
		ring_mask_L_MD = dilate(dilate(dilate(dilate(mask_L_MD_3D)))) - dilate(dilate(dilate(mask_L_MD_3D)))
		single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
		s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
		alg_L_MD = Integrated(arr[mask_L_MD_3D])
		mass_l_md = score(s_bkg_L_MD, cal_insert_mean, pixel_size, ρ, alg_L_MD)
	
		## Low Density
		mask_L_LD_3D = Array{Bool}(undef, size(arr))
		for z in 1:size(arr, 3)
			mask_L_LD_3D[:, :, z] = mask_L_LD
		end
		dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D))
		ring_mask_L_LD = dilate(dilate(dilate(dilate(mask_L_LD_3D)))) - dilate(dilate(dilate(mask_L_LD_3D)))
		single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
		s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
		alg_L_LD = Integrated(arr[mask_L_LD_3D])
		mass_l_ld = score(s_bkg_L_LD, cal_insert_mean, pixel_size, ρ, alg_L_LD)
		
		# Score Medium Inserts
		## High Density
		mask_M_HD_3D = Array{Bool}(undef, size(arr))
		for z in 1:size(arr, 3)
			mask_M_HD_3D[:, :, z] = mask_M_HD
		end
		dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))))
		ring_mask_M_HD = dilate(dilate(dilate(dilate(dilate(mask_M_HD_3D))))) - dilate(dilate(dilate(dilate(mask_M_HD_3D))))
		single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
		s_bkg_M_HD = mean(single_arr[single_ring_mask_M_HD])
		alg_M_HD = Integrated(arr[mask_M_HD_3D])
		mass_m_hd = score(s_bkg_M_HD, cal_insert_mean, pixel_size, ρ, alg_M_HD)
		
		## Medium Density
		mask_M_MD_3D = Array{Bool}(undef, size(arr))
		for z in 1:size(arr, 3)
			mask_M_MD_3D[:, :, z] = mask_M_MD
		end
		dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))))
		ring_mask_M_MD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))
		single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
		s_bkg_M_MD = mean(single_arr[single_ring_mask_M_MD])
		alg_M_MD = Integrated(arr[mask_M_MD_3D])
		mass_m_md = score(s_bkg_M_MD, cal_insert_mean, pixel_size, ρ, alg_M_MD)
	
		## Low Density
		mask_M_LD_3D = Array{Bool}(undef, size(arr))
		for z in 1:size(arr, 3)
			mask_M_LD_3D[:, :, z] = mask_M_LD
		end
		dilated_mask_M_LD = dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))
		ring_mask_M_LD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))
		single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
		s_bkg_M_LD = mean(single_arr[single_ring_mask_M_LD])
		alg_M_LD = Integrated(arr[mask_M_LD_3D])
		mass_m_ld = score(s_bkg_M_LD, cal_insert_mean, pixel_size, ρ, alg_M_LD)
		
	
		# Score Small Inserts
		## High Density
		mask_S_HD_3D = Array{Bool}(undef, size(arr))
		for z in 1:size(arr, 3)
			mask_S_HD_3D[:, :, z] = mask_S_HD
		end
		dilated_mask_S_HD = dilate(dilate(dilate(dilate(dilate((mask_S_HD_3D))))))
		ring_mask_S_HD = dilate(dilate(dilate(dilate(dilate(mask_S_HD_3D))))) - dilate(dilate(dilate(dilate(mask_S_HD_3D))))
		single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
		s_bkg_S_HD = mean(single_arr[single_ring_mask_S_HD])
		alg_S_HD = Integrated(arr[mask_S_HD_3D])
		mass_s_hd = score(s_bkg_S_HD, cal_insert_mean, pixel_size, ρ, alg_S_HD)
		if mass_s_hd < 0
			mass_s_hd = 0
		end
		mass_s_hd
	
		## Medium Density
		mask_S_MD_3D = Array{Bool}(undef, size(arr))
		for z in 1:size(arr, 3)
			mask_S_MD_3D[:, :, z] = mask_S_MD
		end
		dilated_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D)))))
		ring_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D))))) - dilate(dilate(dilate(dilate(mask_S_MD_3D))))
		single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
		s_bkg_S_MD = mean(single_arr[single_ring_mask_S_MD])
		alg_S_MD = Integrated(arr[mask_S_MD_3D])
		mass_s_md = score(s_bkg_S_MD, cal_insert_mean, pixel_size, ρ, alg_S_MD)
		if mass_s_md < 0
			mass_s_md = 0
		end
		mass_s_md
		
		## Low Density
		mask_S_LD_3D = Array{Bool}(undef, size(arr))
		for z in 1:size(arr, 3)
			mask_S_LD_3D[:, :, z] = mask_S_LD
		end
		dilated_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D)))))
		ring_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D))))) - dilate(dilate(dilate(dilate(mask_S_LD_3D))));
		
		single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
		s_bkg_S_LD = mean(single_arr[single_ring_mask_S_LD])
		alg_S_LD = Integrated(arr[mask_S_LD_3D])
		mass_s_ld = score(s_bkg_S_LD, cal_insert_mean, pixel_size, ρ, alg_S_LD)
		if mass_s_ld < 0
			mass_s_ld = 0
		end
		mass_s_ld
		
	
		# Results
		density_array = [0, 200, 400, 800]
		inserts = [
			"Low Density",
			"Medium Density",
			"High Density"
		]
		ground_truth_mass_large = [
			19.6,
			39.3,
			78.5
		] # mg
		
		calculated_mass_large = [
			mass_l_ld,
			mass_l_md,
			mass_l_hd
		]
		ground_truth_mass_medium = [
			4.2,
			8.5,
			17.0
		]
		calculated_mass_medium = [
			mass_m_ld,
			mass_m_md,
			mass_m_hd
		]
		ground_truth_mass_small = [
			0.2,
			0.3,
			0.6
		]
		calculated_mass_small = [
			mass_s_ld,
			mass_s_md,
			mass_s_hd
		]
		RMSE_check = sqrt(sum((ground_truth_mass_small .- calculated_mass_small)).^2 / 3)
		num_zero = length(findall(x -> x == 0, calculated_mass_small))
		RMSE_check += num_zero # penalize zero values heavily
		if haskey(RMSE_Dict, "value") == false
			RMSE_Dict["value"] = RMSE_check
			RMSE_Dict["factor"] = angle
		end
		if RMSE_check < RMSE_Dict["value"]
			RMSE_Dict["value"] = RMSE_check
			RMSE_Dict["factor"] = angle
		end
	end
end

# ╔═╡ c7c44834-19a2-4ad9-8c9b-dcaa54f67f61
md"""
### Save Results
"""

# ╔═╡ 3bc13f3a-fdb3-4d09-8024-e6b283d7c6ae
angle_factor = RMSE_Dict["factor"]

# ╔═╡ 1be94916-8045-4186-aee3-0afcfe5b4739
begin
	mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
		dcm_array, masked_array, header, slice_CCI, center_insert; angle_factor=angle_factor)

	arr = masked_array[:, :, slice_CCI-3:slice_CCI+3]
	single_arr = masked_array[:, :, slice_CCI]
	pixel_size = DICOMUtils.get_pixel_size(header)
	ρ = 0.2 # mg/mm^3
	# Score Large InsertS
	## High Density
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_HD_3D[:, :, z] = mask_L_HD
	end
	dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D))
	ring_mask_L_HD = dilate(dilate(dilate(dilate(mask_L_HD_3D)))) - dilate(dilate(dilate(mask_L_HD_3D)))
	single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
	s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
	alg_L_HD = Integrated(arr[mask_L_HD_3D])
	mass_l_hd = score(s_bkg_L_HD, cal_insert_mean, pixel_size, ρ, alg_L_HD)

	## Medium Density
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_MD_3D[:, :, z] = mask_L_MD
	end
	dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D))
	ring_mask_L_MD = dilate(dilate(dilate(dilate(mask_L_MD_3D)))) - dilate(dilate(dilate(mask_L_MD_3D)))
	single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
	s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
	alg_L_MD = Integrated(arr[mask_L_MD_3D])
	mass_l_md = score(s_bkg_L_MD, cal_insert_mean, pixel_size, ρ, alg_L_MD)

	## Low Density
	mask_L_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_LD_3D[:, :, z] = mask_L_LD
	end
	dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D))
	ring_mask_L_LD = dilate(dilate(dilate(dilate(mask_L_LD_3D)))) - dilate(dilate(dilate(mask_L_LD_3D)))
	single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
	s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
	alg_L_LD = Integrated(arr[mask_L_LD_3D])
	mass_l_ld = score(s_bkg_L_LD, cal_insert_mean, pixel_size, ρ, alg_L_LD)
	
	# Score Medium Inserts
	## High Density
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_HD_3D[:, :, z] = mask_M_HD
	end
	dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))))
	ring_mask_M_HD = dilate(dilate(dilate(dilate(dilate(mask_M_HD_3D))))) - dilate(dilate(dilate(dilate(mask_M_HD_3D))))
	single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
	s_bkg_M_HD = mean(single_arr[single_ring_mask_M_HD])
	alg_M_HD = Integrated(arr[mask_M_HD_3D])
	mass_m_hd = score(s_bkg_M_HD, cal_insert_mean, pixel_size, ρ, alg_M_HD)
	
	## Medium Density
	mask_M_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_MD_3D[:, :, z] = mask_M_MD
	end
	dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))))
	ring_mask_M_MD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))
	single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
	s_bkg_M_MD = mean(single_arr[single_ring_mask_M_MD])
	alg_M_MD = Integrated(arr[mask_M_MD_3D])
	mass_m_md = score(s_bkg_M_MD, cal_insert_mean, pixel_size, ρ, alg_M_MD)

	## Low Density
	mask_M_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_LD_3D[:, :, z] = mask_M_LD
	end
	dilated_mask_M_LD = dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))
	ring_mask_M_LD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))
	single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
	s_bkg_M_LD = mean(single_arr[single_ring_mask_M_LD])
	alg_M_LD = Integrated(arr[mask_M_LD_3D])
	mass_m_ld = score(s_bkg_M_LD, cal_insert_mean, pixel_size, ρ, alg_M_LD)
	

	# Score Small Inserts
	## High Density
	mask_S_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_HD_3D[:, :, z] = mask_S_HD
	end
	dilated_mask_S_HD = dilate(dilate(dilate(dilate(dilate((mask_S_HD_3D))))))
	ring_mask_S_HD = dilate(dilate(dilate(dilate(dilate(mask_S_HD_3D))))) - dilate(dilate(dilate(dilate(mask_S_HD_3D))))
	single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
	s_bkg_S_HD = mean(single_arr[single_ring_mask_S_HD])
	alg_S_HD = Integrated(arr[mask_S_HD_3D])
	mass_s_hd = score(s_bkg_S_HD, cal_insert_mean, pixel_size, ρ, alg_S_HD)
	if mass_s_hd < 0
		mass_s_hd = 0
	end
	mass_s_hd

	## Medium Density
	mask_S_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_MD_3D[:, :, z] = mask_S_MD
	end
	dilated_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D)))))
	ring_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D))))) - dilate(dilate(dilate(dilate(mask_S_MD_3D))))
	single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
	s_bkg_S_MD = mean(single_arr[single_ring_mask_S_MD])
	alg_S_MD = Integrated(arr[mask_S_MD_3D])
	mass_s_md = score(s_bkg_S_MD, cal_insert_mean, pixel_size, ρ, alg_S_MD)
	if mass_s_md < 0
		mass_s_md = 0
	end
	mass_s_md
	
	## Low Density
	mask_S_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_LD_3D[:, :, z] = mask_S_LD
	end
	dilated_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D)))))
	ring_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D))))) - dilate(dilate(dilate(dilate(mask_S_LD_3D))));
	
	single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
	s_bkg_S_LD = mean(single_arr[single_ring_mask_S_LD])
	alg_S_LD = Integrated(arr[mask_S_LD_3D])
	mass_s_ld = score(s_bkg_S_LD, cal_insert_mean, pixel_size, ρ, alg_S_LD)
	if mass_s_ld < 0
		mass_s_ld = 0
	end
	mass_s_ld
	

	# Results
	density_array = [0, 200, 400, 800]
	inserts = [
		"Low Density",
		"Medium Density",
		"High Density"
	]
	ground_truth_mass_large = [
		19.6,
		39.3,
		78.5
	] # mg
	
	calculated_mass_large = [
		mass_l_ld,
		mass_l_md,
		mass_l_hd
	]
	ground_truth_mass_medium = [
		4.2,
		8.5,
		17.0
	]
	calculated_mass_medium = [
		mass_m_ld,
		mass_m_md,
		mass_m_hd
	]
	ground_truth_mass_small = [
		0.2,
		0.3,
		0.6
	]
	calculated_mass_small = [
		mass_s_ld,
		mass_s_md,
		mass_s_hd
	]
end

# ╔═╡ d2c84a5c-73b0-4cd8-988c-f7bc1d131745
df = DataFrame(
	inserts = inserts,
	ground_truth_mass_large = ground_truth_mass_large,
	calculated_mass_large = calculated_mass_large,
	ground_truth_mass_medium = ground_truth_mass_medium,
	calculated_mass_medium = calculated_mass_medium,
	ground_truth_mass_small = ground_truth_mass_small,
	calculated_mass_small = calculated_mass_small
)

# ╔═╡ 4e096971-4212-430e-b709-233a6ffe0b0e
if ~isdir(string(cd(pwd, "..") , "/data/output/", VENDER, "5"))
	mkdir(string(cd(pwd, "..") , "/data/output/", VENDER, "5"))
end

# ╔═╡ 77935f25-27da-4ab7-979f-dfb3fe7239af
output_path = string(cd(pwd, "..") , "/data/output/", VENDER, "5", "/", scan, ".csv")

# ╔═╡ bdb10013-585a-4b59-8a0f-d5f311256fee
CSV.write(output_path, df)

# ╔═╡ Cell order:
# ╠═07787e48-5cf4-4583-b605-59054de37632
# ╠═59ca3f65-5e88-44e9-98f1-381036237c9b
# ╟─5706836e-df8d-4a03-8610-dde49bd357dc
# ╠═05b7448d-7377-4a4a-ab8a-cd9e1b3d47c5
# ╠═89785633-a679-4012-a79d-c76c61c742bb
# ╟─c9ab625d-ee86-4cc2-ac23-bd8b573cf8e2
# ╠═11ea321b-4be9-420f-9c9c-01b6450580b4
# ╠═e7cabee9-2bea-4814-a53f-b96eadc9d41e
# ╠═0d3eb83a-0a11-4eb1-8e28-9b98e9316437
# ╠═d431e54b-55b4-40ea-8aab-40cb7652374d
# ╠═148eab5f-3263-4f19-a7b5-f289d631d119
# ╠═a0d58e6e-9cf7-4802-bc87-27259a310840
# ╠═45c84e2b-5439-4ad3-b249-248abb0d0afb
# ╟─daf9758b-2884-46ea-807f-8e17c96f4fda
# ╟─be2e6d55-3ad0-4749-9668-cfb275470e22
# ╟─59afc75e-fa7b-4b0b-90bd-92bf9550a30f
# ╟─f55a30b7-bbfc-413e-b1ce-a8a1fe7cc1a8
# ╟─4aeb5503-a295-4f2b-87e1-578230e01cde
# ╠═ed209978-0ef0-4f7b-9b17-89e35fd7f4cd
# ╟─4222f5fe-78b5-4bf8-a42b-0b32379cda4d
# ╠═f0582426-8c3e-4f47-b8b0-af18920dc7e3
# ╟─a5577cd7-8d70-4730-aefc-4d5dcb439734
# ╠═6809886f-55fd-4899-b7a0-0f2076422643
# ╠═64210c98-98fb-4413-879c-9bedd2f5060e
# ╠═fdb968c5-3a3b-43eb-84ba-42617dea9915
# ╠═fd90f126-0912-4026-b024-f7cf3968f896
# ╠═f403c1fa-35ab-4526-b105-a8056da88ae0
# ╠═ef899480-a460-49c8-abce-2ed81922385b
# ╠═5f305c11-b050-4d03-a7af-d7b76da5d3a8
# ╠═27d61bb3-00c2-45de-af02-0ca8108fed08
# ╟─810a1af6-2d85-46d5-9cb1-e316225e201b
# ╠═10f3d200-f876-4d66-9cfa-8303036310f4
# ╟─041fdc88-5fbc-4845-8514-622fa43d36ab
# ╠═122d6e0b-3b00-4056-b68b-8452c9364941
# ╟─c7c44834-19a2-4ad9-8c9b-dcaa54f67f61
# ╠═3bc13f3a-fdb3-4d09-8024-e6b283d7c6ae
# ╠═1be94916-8045-4186-aee3-0afcfe5b4739
# ╠═d2c84a5c-73b0-4cd8-988c-f7bc1d131745
# ╠═4e096971-4212-430e-b709-233a6ffe0b0e
# ╠═77935f25-27da-4ab7-979f-dfb3fe7239af
# ╠═bdb10013-585a-4b59-8a0f-d5f311256fee

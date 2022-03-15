### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 34a6ee83-5c46-4661-910e-925ea33ef625
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
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

# ╔═╡ c34ae2e3-466f-4280-8786-b0e92dbb2792
TableOfContents()

# ╔═╡ ff0a6e93-9c5e-4ac2-af02-ce253b535c6b
begin
	VENDER = "GE_Revolution"
	BASE_PATH = "/Users/daleblack/Google Drive/Datasets/"
end;

# ╔═╡ 7d8366d8-9cc5-49b3-8db3-1b9db6d12c92
scans = collect(1:10)

# ╔═╡ 60fe5916-a1a6-11ec-1fbe-8952690e9d06
begin
	dfs = []
	for s in scans
		SCAN_NUMBER = s
		root_path = string(BASE_PATH, VENDER)
		dcm_path_list = dcm_list_builder(root_path)
		pth = dcm_path_list[SCAN_NUMBER]
		scan = basename(pth)
		header, dcm_array, slice_thick_ori1 = dcm_reader(pth)
	
		# Segment Heart
		masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2)
	
		# Segment Calcium Rod
		calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header)
	
		# Calibration Prep
		array_filtered = abs.(mapwindow(median, calcium_image[:, :, cal_rod_slice], (3, 3)))
		bool_arr = array_filtered .> 0
		bool_arr_erode = (((erode(erode(bool_arr)))))
		c_img = calcium_image[:, :, cal_rod_slice-1:cal_rod_slice+1]
		mask_cal_3D = Array{Bool}(undef, size(c_img))
		for z in 1:size(c_img, 3)
			mask_cal_3D[:, :, z] = bool_arr_erode
		end
		cal_insert_mean2 = mean(c_img[mask_cal_3D])
		cal_insert_mean = quantile!(c_img[mask_cal_3D], 0.7)
	
		# Segment Calcium Inserts
		# Check Best Segmentation
		angles = collect(-10:2:10)
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
			RMSE_check += sqrt(sum((ground_truth_mass_medium .- calculated_mass_medium)).^2 / 3)
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
	
		# Run Top Segmentation
		angle_factor = RMSE_Dict["factor"]
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
		push!(dfs, df)
	end
end

# ╔═╡ c8b89590-c7c4-4c40-93c4-9fd905eb3328
md"""
# Save Results
"""

# ╔═╡ 47dfa770-47ea-47ee-98e7-87da471ab01c
new_df = vcat(dfs[1:10]...)

# ╔═╡ e730f5ae-ae32-4869-9a59-98933bcba510
if ~isdir(string(cd(pwd, "..") , "/data/output/", VENDER, "_integrated_score_script"))
	mkdir(string(cd(pwd, "..") , "/data/output/", VENDER, "_integrated_score_script"))
end

# ╔═╡ 75e9e00c-9603-46e9-a009-7cdc76f1a2cc
output_path = string(cd(pwd, "..") , "/data/output/", VENDER, "_integrated_score_script", "/output.csv")

# ╔═╡ a0872982-dcbb-4278-beaf-70409400b514
CSV.write(output_path, new_df)

# ╔═╡ Cell order:
# ╠═34a6ee83-5c46-4661-910e-925ea33ef625
# ╠═c34ae2e3-466f-4280-8786-b0e92dbb2792
# ╠═ff0a6e93-9c5e-4ac2-af02-ce253b535c6b
# ╠═7d8366d8-9cc5-49b3-8db3-1b9db6d12c92
# ╠═60fe5916-a1a6-11ec-1fbe-8952690e9d06
# ╟─c8b89590-c7c4-4c40-93c4-9fd905eb3328
# ╠═47dfa770-47ea-47ee-98e7-87da471ab01c
# ╠═e730f5ae-ae32-4869-9a59-98933bcba510
# ╠═75e9e00c-9603-46e9-a009-7cdc76f1a2cc
# ╠═a0872982-dcbb-4278-beaf-70409400b514

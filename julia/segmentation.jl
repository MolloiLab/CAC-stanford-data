### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ 14b6c3e4-48dd-11ec-336c-6918bf024f43
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("Revise")
		Pkg.add("ImageFiltering")
		Pkg.add("FixedPointNumbers")
		Pkg.add("Images")
		Pkg.add("Statistics")
		Pkg.add("CairoMakie")
		Pkg.add(url="https://github.com/Dale-Black/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
	end

	using Revise
	using PlutoUI
	using ImageFiltering
	using FixedPointNumbers
	using Images
	using Statistics
	using CairoMakie
	using DICOM
	using DICOMUtils
end

# ╔═╡ 9506ae45-f63f-4fd6-8467-0597cbde7006
TableOfContents()

# ╔═╡ 89d2d135-7669-40be-bba1-24a324cdf14e
md"""
Do the dcm specific functions last

Compare each function to the corresponding python function and test that the return values are the same
"""

# ╔═╡ 0fbe9f2a-fea0-45fc-b429-85afc5db3341
tst_path = "/Users/daleblack/Google Drive/Datasets/Canon_Aquilion_One_Vision/"

# ╔═╡ d4c3ba02-f060-426d-805b-c6cb6e04c507
dir = readdir(tst_path)

# ╔═╡ 3f313b73-2785-4299-a1e9-bc5d4e1ad422
md"""
## `dcm_reader`

TODO: Last
"""

# ╔═╡ 482d8669-c750-450b-aae8-03afc2a90007
md"""
## `dcm_list_builder`

TODO: Last
"""

# ╔═╡ 1b755e07-926a-42bd-871e-a693bb037c81
md"""
## `findCircle`
"""

# ╔═╡ d2d1d02c-e592-4c73-afa2-80f89dd5534e
function findCircle(point_1, point_2, point_3)
    x1, y1 = point_1
    x2, y2 = point_2
    x3, y3 = point_3
    
    x12 = x1 - x2 
    x13 = x1 - x3  
    y12 = y1 - y2  
    y13 = y1 - y3 
    y31 = y3 - y1  
    y21 = y2 - y1
    x31 = x3 - x1  
    x21 = x2 - x1 
 
    sx13 = x1^2 - x3^2  
    sy13 = y1^2 - y3^2
    sx21 = x2^2 - x1^2  
    sy21 = y2^2 - y1^2  
  
    f = (((sx13) * (x12) + (sy13) * (x12) + (sx21) * (x13) + (sy21) * (x13)) ÷ (2 * ((y31) * (x12) - (y21) * (x13)))) 
              
    g = (((sx13) * (y12) + (sy13) * (y12) + (sx21) * (y13) + (sy21) * (y13)) ÷ (2 * ((x31) * (y12) - (x21) * (y13))))  
  
    # eqn of circle be x^2 + y^2 + 2*g*x + 2*f*y + c = 0 where centre is (h = -g, k = -f)  
    center_insert = [-g,-f]

    return center_insert
end

# ╔═╡ 2d16b7cb-584e-4ff8-84a1-4e5670bd49d3
findCircle([309, 309], [312, 200], [155, 155])

# ╔═╡ a69293d3-6093-4b87-abc6-11f1ec9398c3
md"""
## `dcm_masked`
"""

# ╔═╡ 59ea4e4d-463e-4a01-a436-315b2995d888
dcm_path = string(tst_path, dir[4])

# ╔═╡ f8d08a1d-bbed-46cd-9c3e-7306f385321d
dcm_dir = readdir(dcm_path);

# ╔═╡ 0883a26a-54ba-4670-a19b-a6a0a23f5aee
dcms = dcmdir_parse(dcm_path);

# ╔═╡ 59976bc3-271e-4418-874c-b553b59e730a
dcm_array = load_dcm_array(dcms);

# ╔═╡ 9b02d4cc-e9b3-47aa-8d57-e0befd64a5fc
header = dcms[1].meta;

# ╔═╡ dfcda44d-5b2c-4a50-9fb2-44cb70b5488e
function dcm_masked(
	header; 
	array_used=nothing, 
	radius_val=95, 
	slice_used_center=nothing
	)
	
	pixel_size = 
		try
			header[(0x0028, 0x0030)]
		catch
	        FOV = header[(0x0018, 0x1100)]
	        matrix_size = header[(0x0028, 0x0010)]
	    
	        pixel_size = FOV / matrix_size
		end
    
    radius = (radius_val / 2) / pixel_size
    central_image = copy(array_used[:, :, slice_used_center])
    central_image = Int.(central_image .< -200)
    kern = Int.(round(5 / pixel_size[1]))
    if kern % 2 == 0
        kern += 1
	end
	
    central_image = mapwindow(median, central_image, (kern, kern))
    center = [size(central_image, 1) ÷ 2, size(central_image, 2) ÷ 2]
    a = copy(central_image)
    local point_1
	for index in 1:size(central_image, 2) ÷ 2
		if (central_image[center[1] + index, center[2] + index] == 1 && central_image[center[1] + index, center[2] + index + 5] == 1) 
			point_1 = [center[1] + index, center[2] + index]
			break
		else
			a[center[1] + index, center[2] + index] = 2
		end
	end
    
	local point_2
	for index in 1:size(central_image, 2) ÷ 2
		if (central_image[center[1] + index, center[2] - index] == 1 && central_image[center[1] + index, center[2] - index - 5] == 1) 
			point_2 = [center[1] + index, center[2] - index]
			break
		else
			a[center[1] + index, center[2] - index] = 2
		end
	end
	
	local point_3
	for index in 1:size(central_image, 2) ÷ 2
		if (central_image[center[1] - index, center[2] - index] == 1 && central_image[center[1] - index, center[2] - index - 5] == 1)
			point_3 = [center[1] - index, center[2] - index]
			break
		else
			a[center[1] - index, center[2] - index] = 2
		end
	end
	
    center_insert = findCircle(point_1, point_2, point_3)
	rows, cols = Int(header[(0x0028, 0x0010)]), Int(header[(0x0028, 0x0011)])
    Y, X = collect(1:rows), collect(1:cols)'
    dist_from_center = @. sqrt((X - center_insert[2])^2 + (Y-center_insert[1])^2)

    mask = dist_from_center .<= radius[1]  
    masked_array = zeros(size(array_used))
    for index in 1:size(array_used, 3)
        masked_array[:, :, index] = array_used[:, :, index] .* mask
	end

    return masked_array, center_insert, mask
end

# ╔═╡ a8da3e8d-c624-450c-97f1-531715e94931
rows = Int(header[(0x0028, 0x0010)])

# ╔═╡ 2f983fba-55f9-4bcf-be49-a0e30d97259c
norm = collect(1:rows)

# ╔═╡ d07b9f1a-51f1-4367-ad4f-65f177882733
tran = collect(1:rows)'

# ╔═╡ 21702280-ab33-4a97-b436-4e788ae2af26
size(norm)

# ╔═╡ 85b9ed25-4be4-4a90-b0e9-c46b4ee82c87
size(tran)

# ╔═╡ 57730e9f-05b2-4e3a-bff3-a3a6403a66e0
masked_array, center_insert, mask = dcm_masked(header; array_used=dcm_array, slice_used_center=size(dcm_array, 3)÷2);

# ╔═╡ f9eb1a87-9b95-483f-b4e4-5eb5147b41e1
with_terminal() do
	dcm_masked(header; array_used=dcm_array, slice_used_center=size(dcm_array, 3) ÷ 2)
end

# ╔═╡ 163c0220-0070-47c1-a15d-6d6f947a7dfd
heatmap(dcm_array[:, :, 23], colormap=:grays)

# ╔═╡ 0d230175-04c8-4ba6-b89f-933ac80519f9
heatmap(masked_array[:, :, 23], colormap=:grays)

# ╔═╡ Cell order:
# ╠═14b6c3e4-48dd-11ec-336c-6918bf024f43
# ╠═9506ae45-f63f-4fd6-8467-0597cbde7006
# ╟─89d2d135-7669-40be-bba1-24a324cdf14e
# ╠═0fbe9f2a-fea0-45fc-b429-85afc5db3341
# ╠═d4c3ba02-f060-426d-805b-c6cb6e04c507
# ╟─3f313b73-2785-4299-a1e9-bc5d4e1ad422
# ╟─482d8669-c750-450b-aae8-03afc2a90007
# ╟─1b755e07-926a-42bd-871e-a693bb037c81
# ╠═d2d1d02c-e592-4c73-afa2-80f89dd5534e
# ╠═2d16b7cb-584e-4ff8-84a1-4e5670bd49d3
# ╟─a69293d3-6093-4b87-abc6-11f1ec9398c3
# ╠═59ea4e4d-463e-4a01-a436-315b2995d888
# ╠═f8d08a1d-bbed-46cd-9c3e-7306f385321d
# ╠═0883a26a-54ba-4670-a19b-a6a0a23f5aee
# ╠═59976bc3-271e-4418-874c-b553b59e730a
# ╠═9b02d4cc-e9b3-47aa-8d57-e0befd64a5fc
# ╠═dfcda44d-5b2c-4a50-9fb2-44cb70b5488e
# ╠═a8da3e8d-c624-450c-97f1-531715e94931
# ╠═2f983fba-55f9-4bcf-be49-a0e30d97259c
# ╠═d07b9f1a-51f1-4367-ad4f-65f177882733
# ╠═21702280-ab33-4a97-b436-4e788ae2af26
# ╠═85b9ed25-4be4-4a90-b0e9-c46b4ee82c87
# ╠═57730e9f-05b2-4e3a-bff3-a3a6403a66e0
# ╠═f9eb1a87-9b95-483f-b4e4-5eb5147b41e1
# ╠═163c0220-0070-47c1-a15d-6d6f947a7dfd
# ╠═0d230175-04c8-4ba6-b89f-933ac80519f9

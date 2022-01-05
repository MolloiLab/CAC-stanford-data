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

# ╔═╡ 14b6c3e4-48dd-11ec-336c-6918bf024f43
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("Revise")
		Pkg.add("PlutoUI")
		Pkg.add("ImageFiltering")
		Pkg.add("Images")
		Pkg.add("ImageComponentAnalysis")
		Pkg.add("DataFrames")
		Pkg.add("Statistics")
		Pkg.add("CairoMakie")
		Pkg.add(url="https://github.com/JuliaHealth/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
	end

	using Revise
	using PlutoUI
	using ImageFiltering
	using Images
	using ImageComponentAnalysis
	using DataFrames
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
tst_path = string(cd(pwd, "..") , "/", "data/Large_rep1")

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

# ╔═╡ 39688933-f68f-47b4-a9dc-aff25bdc6123
md"""
## Load DICOMs
"""

# ╔═╡ 0883a26a-54ba-4670-a19b-a6a0a23f5aee
dcms = dcmdir_parse(tst_path);

# ╔═╡ 59976bc3-271e-4418-874c-b553b59e730a
dcm_array = load_dcm_array(dcms);

# ╔═╡ 9b02d4cc-e9b3-47aa-8d57-e0befd64a5fc
header = dcms[1].meta;

# ╔═╡ 237cfd63-3479-468d-954f-2ebd05376fca
md"""
# Whole heart mask
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
  
    # eqn of circle be x^2 + y^2 + 2*g*x + 2*f*y + c = 0 where center is (h = -g, k = -f)  
    center_insert = [-g, -f]

    return center_insert
end

# ╔═╡ 2d16b7cb-584e-4ff8-84a1-4e5670bd49d3
findCircle([309, 309], [312, 200], [155, 155])

# ╔═╡ a69293d3-6093-4b87-abc6-11f1ec9398c3
md"""
## `dcm_masked`
"""

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
    
    radius = (radius_val / 2) / pixel_size[1]
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

# ╔═╡ 57730e9f-05b2-4e3a-bff3-a3a6403a66e0
masked_array, center_insert, mask = dcm_masked(header; array_used=dcm_array, slice_used_center=size(dcm_array, 3)÷2);

# ╔═╡ 19e0f1b3-86c3-4ffd-aa39-dd7e1e509312
center_insert

# ╔═╡ 55665585-8076-4e23-a989-71fc9909a57b
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 15]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig
end

# ╔═╡ 78ffb833-27c5-4df5-9989-47d57aca47fe
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig2
end

# ╔═╡ 1ac4e0fd-79fb-41d1-b2c7-8b5e011c80e2
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 23]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig3
end

# ╔═╡ 5147ce8d-2d44-443b-a6a2-edba166356a8
@bind a2 PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 004f4a1e-9d11-44ef-8622-69c311edb2cb
heatmap(transpose(masked_array[:, :, a2]), colormap=:grays)

# ╔═╡ da57bd3c-0d8b-453c-9c4f-ea846850cbec
md"""
# Calcium rod mask
"""

# ╔═╡ 11e0c502-5b1f-4e60-bb43-a8c04a4e335b
md"""
## `get_indices`
"""

# ╔═╡ 316a0f8f-163d-44e1-87c2-2cfca14cfc48
function get_pixel_size(header)
	pixel_size = 
		try
			header[(0x0028, 0x0030)]
		catch
			FOV = header[(0x0018, 0x1100)]
			matrix_size = header[(0x0028, 0x0010)]
		
			pixel_size = FOV / matrix_size
		end
	return pixel_size
end

# ╔═╡ c4c6841a-88f4-4ca9-9353-0e9daa6d4dff
function get_indices(dcm_array, header; calcium_threshold=130, comp_connect=4)
    array = copy(dcm_array)
    array = Int.(array .> (1.1 * calcium_threshold))

	pixel_size = get_pixel_size(header)
    CCI_5mm_num_pixels = Int(round(π * (5/2)^2 / pixel_size[1]^2))
    cal_rod_num_pixels = Int(round(π * (20/2)^2 / pixel_size[1]^2))
    
	kern = Int.(round(5 / pixel_size[1]))
    if kern % 2 == 0
        kern += 1
	end
    
    slice_dict = Dict()
    large_index = []
    cal_rod_dict = Dict()
    for idx in 1:size(array, 3)
		array_filtered = mapwindow(median, array[:,:,idx], (kern, kern))
        components = ImageComponentAnalysis.label_components(array_filtered)
        count_5mm = 0
		a1 = analyze_components(
			components, BasicMeasurement(area=true, perimeter=true)
		)
		a2 = analyze_components(components, BoundingBox(box_area = true))
		df = leftjoin(a1, a2, on = :l)
		count = 0
		for row in eachrow(df)
			count += 1
			df_area = Int(round(row[:area]))
			if df_area in Int(round(CCI_5mm_num_pixels * 0.6)):Int(round(CCI_5mm_num_pixels * 1.5))
				count_5mm += 1
			elseif df_area in Int(round(cal_rod_num_pixels * 0.7)):Int(round(cal_rod_num_pixels * 1.3))
				indices = row[:box_indices]
				x_point = ((indices[1][end] - indices[1][1]) ÷ 2) + indices[1][1]
				y_point = ((indices[2][end] - indices[2][1]) ÷ 2) + indices[2][1]
				cal_rod_dict[count] = [x_point, y_point]
			end
		end
        if count_5mm > 0 && count_5mm < 4
            slice_dict[idx] = count_5mm
		end
    
        poppable_keys = []
        for key in cal_rod_dict
            start_coordinate = [key[2][1], key[2][2]]
			
            x_right = 0
            while array_filtered[start_coordinate[1], start_coordinate[2] + x_right] == 1
                x_right += 1
			end
            
            x_left = 0
            while array_filtered[start_coordinate[1], start_coordinate[2] - x_left] == 1
                x_left += 1
			end
            
            y_top = 0
            while array_filtered[start_coordinate[1] + y_top, start_coordinate[2]] == 1
                y_top += 1
			end
            
            y_bottom = 0
            while array_filtered[start_coordinate[1] - y_bottom, start_coordinate[2]] == 1
                y_bottom += 1
			end
                
            x_dist = x_right + x_left
            y_dist = y_top + y_bottom
            
            if ((x_dist in round(0.7 * y_dist):round(1.2 * y_dist)) == false) || ((round(0.7 * y_dist) == 0) && (round(1.2 * y_dist) == 0))
                push!(poppable_keys, key)
            else
				nothing
			end
		end
		
        for key in poppable_keys
			pop!(cal_rod_dict)
		end
                
        if length(cal_rod_dict) == 0
			nothing
        else
            append!(large_index, idx)
		end
	end
	return slice_dict, large_index
end

# ╔═╡ 462fd944-d4ce-4d7b-82ad-40342daff354
slice_dict, large_index = get_indices(masked_array, header);

# ╔═╡ 8959750c-291c-4c2a-a6a1-218f8ed82ab1
slice_dict, large_index

# ╔═╡ 68c16f9e-8520-41ef-912a-9e61acf2ad1d
md"""
## `find_edges`
"""

# ╔═╡ f49fe2a8-8429-4fd1-98e3-6e75189a3ac2
function find_edges(dcm_array, slice_dict, large_index)
	flipped_index = Int(round(median(large_index)))
    edge_index = []
    if flipped_index < (size(dcm_array, 3) / 2)
        flipped = -1
        for element in large_index
            if element > (size(dcm_array, 3) / 2)
                append!(edge_index, element)
			end
		end
        if length(edge_index) == 0
            nothing
        else
            for index_edge in minimum(edge_index):size(dcm_array, 3)
                try
                    delete!(slice_dict, index_edge)
				catch
                    nothing
				end
			end
            for element2 in edge_index
				deleteat!(large_index, findall(x->x==element2, large_index))
			end
		end
                
        for element in 1:maximum(large_index)
            try
                delete!(slice_dict, element)
			catch
                nothing
			end
		end
    else
        flipped = 1
        for element in large_index
            if element < (size(dcm_array, 3) / 2)
                append!(edge_index, element)
			end
		end
        if length(edge_index) == 0
            nothing
        else
            for index_edge in 1:maximum(edge_index)
                try
                    delete!(slice_dict, index_edge)
				catch
                    nothing
				end
			end
            for element2 in edge_index
				deleteat!(large_index, findall(x->x==element2, large_index))
			end
		end
        for element in minimum(large_index):size(dcm_array, 3)
            try
                delete!(slice_dict, element)
			catch
                nothing
			end
		end
	end
	return slice_dict, flipped, flipped_index
end

# ╔═╡ 494163b1-cb52-45a9-815e-991633018cd5
slice_dict2, flipped, flipped_index = find_edges(masked_array, slice_dict, large_index)

# ╔═╡ aed2176b-dfe2-4a38-bc8d-2e5ec6e4531f
slice_dict2, flipped, flipped_index

# ╔═╡ 880db52b-e39c-4578-a6de-ce4a6afa038e
md"""
## `poppable_keys`
"""

# ╔═╡ 6613bd97-7d53-4b8b-a69c-12e040f4aa36
function poppable_keys(flipped, flipped_index, header, slice_dict)
	SliceThickness = header[(0x0018,0x0050)]
	poppable_keys = []        
    if flipped == -1
        for key in slice_dict
            if key[1] > (flipped_index + (55 / SliceThickness))
                append!(poppable_keys, key)
			elseif flipped == 1
		        for key in slice_dict
		            if key[1] < (flipped_index - (55 / SliceThickness))
		                append!(poppable_keys, key)
					end
				end
			end
		end
	end
    for key in poppable_keys
		pop!(slice_dict)           
	end
	return slice_dict
end

# ╔═╡ bff84e32-1013-44f2-9eda-eb5f6bac4c5a
pop_keys = poppable_keys(flipped, flipped_index, header, slice_dict)

# ╔═╡ 451dd3c6-f740-47df-a5b9-85749caf04fb
md"""
## `compute_CCI`
"""

# ╔═╡ 1ac3d601-0057-475e-a6ac-9522618f7bb9
function compute_CCI(dcm_array, header, slice_dict; calcium_threshold=130)
	SliceThickness = header[(0x0018,0x0050)]
	max_key, _ = maximum(zip(values(slice_dict), keys(slice_dict)))
    max_keys = []
    for key in slice_dict
        if key[2] == max_key
            append!(max_keys, key[1])
		end
	end
    slice_CCI = Int(floor(median(max_keys)))
    
    array = copy(dcm_array)
    array = Int.(array .> calcium_threshold)
    
    calcium_image = array .* dcm_array
    quality_slice = round(slice_CCI - flipped * (20 / SliceThickness))

    cal_rod_slice = slice_CCI + (flipped * Int(30 / SliceThickness))
    
    return calcium_image, slice_CCI, quality_slice, cal_rod_slice
end

# ╔═╡ 5d407e49-fe05-4748-be36-b232b5ce3ccb
calcium_image1, slice_CCI1, quality_slice1, cal_rod_slice1 = compute_CCI(masked_array, header, slice_dict);

# ╔═╡ de1ac95e-2399-45c1-b337-dcfe800394a2
slice_CCI1

# ╔═╡ a8f6ece5-9f59-4af1-bcf7-7dafa839bae2
@bind a PlutoUI.Slider(1:size(calcium_image1, 3), default=10, show_value=true)

# ╔═╡ 267d1d6a-5c35-480e-8425-f000f7f31c7f
begin
	fig4 = Figure()
	
	ax4 = Makie.Axis(fig4[1, 1])
	ax4.title = "Raw DICOM Array"
	heatmap!(transpose(calcium_image1[:, :, a]), colormap=:grays)
	ax5 = Makie.Axis(fig4[1,2])
	heatmap!(transpose(dcm_array[:, :, a]), colormap=:grays)
	fig4
end

# ╔═╡ 9a77474e-a853-4369-a6f3-a9e76a683a27
md"""
## `CCI_calcium_image`
"""

# ╔═╡ a9f499f8-ca32-4310-82d7-f20b790dc04c
function CCI_calcium_image(dcm_array, header; calcium_threshold=130, comp_connect=4)
    slice_dict, large_index = get_indices(
		dcm_array, header; 
		calcium_threshold=calcium_threshold, comp_connect= comp_connect
	)
    slice_dict, flipped, flipped_index = find_edges(
		dcm_array, slice_dict, large_index
	)
    slice_dict = poppable_keys(flipped, flipped_index, header, slice_dict)
    calcium_image, slice_CCI, quality_slice, cal_rod_slice = compute_CCI(
		dcm_array, header, slice_dict; calcium_threshold=calcium_threshold
	)
	return calcium_image, slice_CCI, quality_slice, cal_rod_slice
end

# ╔═╡ 4eece882-8416-4330-b95a-1cfaeac89aac
calcium_image, slice_CCI, quality_slice, cal_rod_slice = CCI_calcium_image(masked_array, header);

# ╔═╡ 214d485e-df7e-4a45-b52f-10f854d839a1
@bind b PlutoUI.Slider(1:size(calcium_image, 3), default=10, show_value=true)

# ╔═╡ 03622243-d2af-4d1b-a6e1-fd50670e48a8
heatmap(calcium_image[:, :, b], colormap=:grays)

# ╔═╡ ee208b64-3cf3-428b-a5f6-74793a82fd93
md"""
# Calcium inserts mask
"""

# ╔═╡ 5930b9c5-2f85-465b-abe0-223580142d49
md"""
## `angle_calc`
"""

# ╔═╡ 480bd663-d822-4417-b7f0-6679d26a58b4
function angle_calc(side1, side2)
    #Calculate angle between two sides of rectangular triangle
    if side1 == 0
        angle = 0
	elseif side2 == 0
        angle = π / 2
    else
        angle = atan(side1 / side2)
	end
    
    return angle
end

# ╔═╡ f498c292-367d-4441-873e-29bec99af2a6
angle_calc(4, 3)

# ╔═╡ 56ce9966-939e-48cd-9c70-2e441bfaea58
md"""
## `create_circular_mask`
"""

# ╔═╡ 228dd95c-fdab-4a5a-a079-e14238ff5e83
function create_circular_mask(h, w, center_circle, radius_circle)
	Y, X = collect(1:h), collect(1:w)'
    dist_from_center = sqrt.((X .- center_circle[1]).^2 .+ (Y .- center_circle[2]).^2)

    mask = dist_from_center .<= radius_circle
    
    return mask
end

# ╔═╡ 0f572e01-a690-4699-8a8c-6684c27085bb
mask1 = create_circular_mask(40, 40, [20, 20], 1);

# ╔═╡ ce5b9c0a-4f25-4b8f-bf16-17d05b7bfff4
heatmap(mask1, colormap=:grays)

# ╔═╡ 5a8963f1-3ec8-4dbe-b567-df2648219b7f
md"""
## `calc_output`
"""

# ╔═╡ 92fca992-b36d-48f3-803e-f92e6463c291
function calc_output(dcm_array, CCI_slice, calcium_threshold=130, comp_connect=trues(3, 3))
	# Actual scoring for CCI insert
    # First step is to remove slices without calcium from arrays
	PixelSpacing = get_pixel_size(header)
	SliceThickness = header[(0x0018,0x0050)]
    CCI_min = Int(round((CCI_slice - (5 ÷ SliceThickness)) - 1))
    CCI_max = Int(round((CCI_slice + (5 ÷ SliceThickness)) + 1))
    central_CCI = Int(round((CCI_max - CCI_min) / 2))
    
    if CCI_min < 0
        CCI_min = 0
	end
    if CCI_max > size(dcm_array, 3)
        CCI_max = size(dcm_array, 3)
	end
    
    CCI_array = copy(dcm_array[:, :, CCI_min:CCI_max])

	image_kernel = Int(round(3 / PixelSpacing[1]))
    if image_kernel % 2 == 0
        image_kernel += 1
	end
    
    CCI_array_binary = copy(CCI_array)
	CCI_array_binary = Int.(CCI_array_binary .> 1.0*calcium_threshold)
	inp = CCI_array_binary[:, :, central_CCI - 1] + CCI_array_binary[:, :, central_CCI] + CCI_array_binary[:, :, central_CCI + 1]
	connectivity = trues(3, 3)
	components = ImageComponentAnalysis.label_components(inp, comp_connect)
	a1 = analyze_components(components, BasicMeasurement(area=true, perimeter=true))
	a2 = analyze_components(components, BoundingBox(box_area = true))
	df = leftjoin(a1, a2, on = :l)
	centroids = []
	for row in eachrow(df)
		indices = row[:box_indices]
		x_point = ((indices[1][end] - indices[1][1]) ÷ 2) + indices[1][1]
		y_point = ((indices[2][end] - indices[2][1]) ÷ 2) + indices[2][1]
		push!(centroids, (x_point, y_point))
	end
    
    centroids = deleteat!(centroids, 1)
        
	i1 = mapwindow(median, CCI_array_binary[:,:,central_CCI - 1], (image_kernel, image_kernel))
	i2 = mapwindow(median, CCI_array_binary[:,:,central_CCI], (image_kernel, image_kernel))
	i3 = mapwindow(median, CCI_array_binary[:,:,central_CCI + 1], (image_kernel, image_kernel))

	image_for_center = i1 + i2 + i3

	components2 = ImageComponentAnalysis.label_components(image_for_center, comp_connect)
	b1 = analyze_components(components2, BasicMeasurement(area=true, perimeter=true))
	b2 = analyze_components(components2, BoundingBox(box_area = true))
	df2 = leftjoin(b1, b2, on = :l)
	centroids2 = []
	for row in eachrow(df2)
		indices = row[:box_indices]
		x_point = ((indices[1][end] - indices[1][1]) ÷ 2) + indices[1][1]
		y_point = ((indices[2][end] - indices[2][1]) ÷ 2) + indices[2][1]
		push!(centroids2, (y_point, x_point))
	end
	output = length(unique(components2)), components2, df2, centroids2
	return output
end

# ╔═╡ fb320178-0d33-4402-b075-7bf7541f94e4
output = calc_output(masked_array, slice_CCI1);

# ╔═╡ feca7b76-fc10-4aca-adff-4cd323dabc74
md"""
## `center_points`
"""

# ╔═╡ ae9ab557-553c-474a-8865-37807b9270fb
function center_points(output, tmp_center, CCI_slice)
	PixelSpacing = get_pixel_size(header)
	rows, cols = Int(header[(0x0028, 0x0010)]), Int(header[(0x0028, 0x0011)])
    sizes = []
    for row in eachrow(output[3])
		area = row[:area]
		append!(sizes, area)
	end

	centroids = output[4]
    largest = Dict()
    for index in 2:length(centroids)
		x = centroids[index][1]
		y = centroids[index][2]
		dist_loc = sqrt((tmp_center[2] - x)^2 + (tmp_center[1] - y)^2)
        dist_loc *= PixelSpacing[1]
        if dist_loc > 31
            largest[index] = [round(y), round(x)]
        else
            nothing
		end
	end

    max_dict = Dict()
	radius = 2.5 ÷ PixelSpacing[1]
    for key in largest
        tmp_arr = create_circular_mask(rows, cols, (key[2][2], key[1][1]), radius)
        tmp_arr = @. tmp_arr * dcm_array[:,:,CCI_slice] + tmp_arr * dcm_array[:,:,CCI_slice - 1] + tmp_arr * dcm_array[:,:,CCI_slice + 1]
        tmp_arr = @. ifelse(tmp_arr == 0, missing, tmp_arr)
        max_dict[key[1]] = median(skipmissing(tmp_arr))
	end

    large1_index, large1_key = maximum(zip(values(max_dict), keys(max_dict)))
	@show large1_index, large1_key
    pop!(max_dict, large1_key)
    large2_index, large2_key = maximum(zip(values(max_dict), keys(max_dict)))
    pop!(max_dict, large2_key)
    large3_index, large3_key = maximum(zip(values(max_dict), keys(max_dict)))

    center1 = largest[large1_key]
    center2 = largest[large2_key]  
    center3 = largest[large3_key]
	
    center = findCircle(center1, center2, center3)
	return center, center1, center2, center3
end

# ╔═╡ 39ce8635-28a7-4eee-875e-325e85d142d5
center_insert

# ╔═╡ 42e86a48-1438-4edd-ac2d-f307d534e165
center, center1, center2, center3 = center_points(output, center_insert, slice_CCI1)

# ╔═╡ 0f73d3a7-c6fa-49ee-8fdf-357b0fe94ec8
md"""
## `calc_centers`
"""

# ╔═╡ 9ac2f9ed-92bb-4540-add5-422df2a3f6da
function calc_centers(output, header, tmp_center, CCI_slice)
	PixelSpacing = get_pixel_size(header)
	_, center1, center2, center3 = center_points(output, center_insert, CCI_slice)
	cents = center1, center2, center3
    centers = Dict()
    for size_index4 in (center1, center2, center3)
        center_index = size_index4
        side_x = abs(center[1]-center_index[1])
        side_y = abs(center[2]-center_index[2])
        
        angle = angle_calc(side_x, side_y)
        if (center_index[1] < center[1] && center_index[2] < center[2])
			medium_calc = [center_index[1] + (12.5 / PixelSpacing[1]) * sin(angle), (center_index[2] + (12.5 / PixelSpacing[2]) * cos(angle))]
			low_calc = [center_index[1] + (25 / PixelSpacing[1]) * sin(angle), (center_index[2] + (25 / PixelSpacing[2]) * cos(angle))]
		elseif (center_index[1] < center[1] && center_index[2] > center[2])
			medium_calc = [center_index[1] + (12.5 / PixelSpacing[1]) * sin(angle), (center_index[2] - (12.5 / PixelSpacing[2]) * cos(angle))]
			low_calc = [center_index[1] + (25 / PixelSpacing[1]) * sin(angle), (center_index[2] - (25 / PixelSpacing[2]) * cos(angle))] 
		elseif (center_index[1] > center[1] && center_index[2] < center[2])
			medium_calc = [center_index[1] - (12.5 / PixelSpacing[1]) * sin(angle), (center_index[2] + (12.5 / PixelSpacing[2]) * cos(angle))]
			low_calc = [center_index[1] - (25 / PixelSpacing[1]) * sin(angle), (center_index[2] + (25 / PixelSpacing[2]) * cos(angle))]
		elseif (center_index[1] > center[1] && center_index[2] > center[2])
			medium_calc = [center_index[1] - (12.5 / PixelSpacing[1]) * sin(angle), (center_index[2] - (12.5 / PixelSpacing[2]) * cos(angle))]
			low_calc = [center_index[1] - (25 / PixelSpacing[1]) * sin(angle), (center_index[2] - (25 / PixelSpacing[2]) * cos(angle))]
		elseif (side_x == 0 && center_index[2] < center[2])
			medium_calc = [center_index[1], center_index[2] + (12.5 / PixelSpacing[2])]
			low_calc = [center_index[1], center_index[2] + (25 / PixelSpacing[2])]
		elseif (side_x == 0 && center_index[2] > center[2])
			medium_calc = [center_index[1], center_index[2] - (12.5 / PixelSpacing[2])]
			low_calc = [center_index[1], center_index[2] - (25 / PixelSpacing[2])]
		elseif (center_index[1] > center[1] && side_y == 0)
            medium_calc = [center_index[1] - (12.5 / PixelSpacing[1]), center_index[2]]
            low_calc = [center_index[1] - (25 / PixelSpacing[1]), center_index[2]]
		elseif (center_index[1] > center[1] && side_y == 0)
			medium_calc = [center_index[1] + (12.5 / PixelSpacing[1]), center_index[2]]
            low_calc = [(center_index[1] + (25 / PixelSpacing[1])), center_index[1]]
        else
			error("unknown angle")
		end
                
        if size_index4 == center1
            centers[:Large_HD] = Int.(round.(center_index))
            centers[:Medium_HD] = Int.(round.(medium_calc))
            centers[:Small_HD] = Int.(round.(low_calc))
        
		elseif size_index4 == center2
            centers[:Large_MD] = Int.(round.(center_index))
            centers[:Medium_MD] = Int.(round.(medium_calc))
            centers[:Small_MD] = Int.(round.(low_calc))
        
		elseif size_index4 == center3
            centers[:Large_LD] = Int.(round.(center_index))
            centers[:Medium_LD] = Int.(round.(medium_calc))
            centers[:Small_LD] = Int.(round.(low_calc))
        
        else
            nothing
		end
	end
    
    return centers
end

# ╔═╡ 533a77b8-d541-4277-8123-cbfe416fd263
calc_centers(output, header, center_insert, slice_CCI1)

# ╔═╡ 891d81f2-f30e-4365-b0ba-b0cebae9f60f
md"""
## `mask_inserts`
"""

# ╔═╡ efad1a40-b466-4e8e-a4ae-ab3c14efc3f6
function mask_inserts(
	dcm_array, header, CCI_slice, center_insert; 
	calcium_threshold=130, comp_connect=4
	)
	
    output = calc_output(masked_array, CCI_slice)
	center, center1, center2, center3 =  center_points(output, center_insert, CCI_slice)
	cents = center1, center2, center3
	tmp_center = copy(center)
    calc_size_density_VS_AS_MS = calc_centers(output, header, tmp_center, CCI_slice)

	PixelSpacing = get_pixel_size(header)
	rows, cols = Int(header[(0x0028, 0x0010)]), Int(header[(0x0028, 0x0011)])

	lg_hd = [calc_size_density_VS_AS_MS[:Large_HD][2], calc_size_density_VS_AS_MS[:Large_HD][1]]
	lg_md = [calc_size_density_VS_AS_MS[:Large_MD][2], calc_size_density_VS_AS_MS[:Large_MD][1]]
	lg_ld = [calc_size_density_VS_AS_MS[:Large_LD][2], calc_size_density_VS_AS_MS[:Large_LD][1]]
	md_hd = [calc_size_density_VS_AS_MS[:Medium_HD][2], calc_size_density_VS_AS_MS[:Medium_HD][1]]
	md_md = [calc_size_density_VS_AS_MS[:Medium_MD][2], calc_size_density_VS_AS_MS[:Medium_MD][1]]
	md_ld = [calc_size_density_VS_AS_MS[:Medium_LD][2], calc_size_density_VS_AS_MS[:Medium_LD][1]]
	sm_hd = [calc_size_density_VS_AS_MS[:Small_HD][2], calc_size_density_VS_AS_MS[:Small_HD][1]]
	sm_md = [calc_size_density_VS_AS_MS[:Small_MD][2], calc_size_density_VS_AS_MS[:Small_MD][1]]
	sm_ld = [calc_size_density_VS_AS_MS[:Small_LD][2], calc_size_density_VS_AS_MS[:Small_LD][1]]
	
    mask_L_HD = create_circular_mask(rows, cols, lg_hd, ((5 ÷ PixelSpacing[1]) / 2) + 1)
    mask_L_MD = create_circular_mask(cols, rows, lg_md, ((5 ÷ PixelSpacing[1]) / 2) + 1)
    mask_L_LD = create_circular_mask(cols, rows, lg_ld, ((5 ÷ PixelSpacing[1]) / 2) + 1)   
    mask_M_HD = create_circular_mask(cols, rows, md_hd,((3 ÷ PixelSpacing[1]) / 2) + 1)
    mask_M_MD = create_circular_mask(cols, rows, md_md, ((3 ÷ PixelSpacing[1]) / 2) + 1)
    mask_M_LD = create_circular_mask(cols, rows, md_ld, ((3 ÷ PixelSpacing[1]) / 2) + 1) 
    mask_S_HD = create_circular_mask(cols, rows, sm_hd, ((1 ÷ PixelSpacing[1]) / 2) + 1)
    mask_S_MD = create_circular_mask(cols, rows, sm_md, ((1 ÷ PixelSpacing[1]) / 2) + 1)
    mask_S_LD = create_circular_mask(cols, rows, sm_ld, ((1 ÷ PixelSpacing[1]) / 2) + 1) 
    
    
    masks1 = mask_L_HD + mask_M_HD + mask_S_HD
    masks2 = mask_L_MD + mask_M_MD + mask_S_MD
    masks3 = mask_L_LD + mask_M_LD + mask_S_LD

    return mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD
end

# ╔═╡ 9cf4a86b-d41f-4814-b14e-0854176e73d4
mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(masked_array, header, 25, center_insert);

# ╔═╡ 492574dc-5674-462d-9df1-a77e0f445994
masks = mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD;

# ╔═╡ 33b38675-7434-46e5-808c-ad4b0cbb0408
heatmap(transpose(masks), colormap=:grays)

# ╔═╡ Cell order:
# ╠═14b6c3e4-48dd-11ec-336c-6918bf024f43
# ╠═9506ae45-f63f-4fd6-8467-0597cbde7006
# ╟─89d2d135-7669-40be-bba1-24a324cdf14e
# ╠═0fbe9f2a-fea0-45fc-b429-85afc5db3341
# ╠═d4c3ba02-f060-426d-805b-c6cb6e04c507
# ╟─3f313b73-2785-4299-a1e9-bc5d4e1ad422
# ╟─482d8669-c750-450b-aae8-03afc2a90007
# ╟─39688933-f68f-47b4-a9dc-aff25bdc6123
# ╠═0883a26a-54ba-4670-a19b-a6a0a23f5aee
# ╠═59976bc3-271e-4418-874c-b553b59e730a
# ╠═9b02d4cc-e9b3-47aa-8d57-e0befd64a5fc
# ╟─237cfd63-3479-468d-954f-2ebd05376fca
# ╟─1b755e07-926a-42bd-871e-a693bb037c81
# ╠═d2d1d02c-e592-4c73-afa2-80f89dd5534e
# ╠═2d16b7cb-584e-4ff8-84a1-4e5670bd49d3
# ╟─a69293d3-6093-4b87-abc6-11f1ec9398c3
# ╠═dfcda44d-5b2c-4a50-9fb2-44cb70b5488e
# ╠═57730e9f-05b2-4e3a-bff3-a3a6403a66e0
# ╠═19e0f1b3-86c3-4ffd-aa39-dd7e1e509312
# ╠═55665585-8076-4e23-a989-71fc9909a57b
# ╠═78ffb833-27c5-4df5-9989-47d57aca47fe
# ╠═1ac4e0fd-79fb-41d1-b2c7-8b5e011c80e2
# ╠═5147ce8d-2d44-443b-a6a2-edba166356a8
# ╠═004f4a1e-9d11-44ef-8622-69c311edb2cb
# ╟─da57bd3c-0d8b-453c-9c4f-ea846850cbec
# ╟─11e0c502-5b1f-4e60-bb43-a8c04a4e335b
# ╠═316a0f8f-163d-44e1-87c2-2cfca14cfc48
# ╠═c4c6841a-88f4-4ca9-9353-0e9daa6d4dff
# ╠═462fd944-d4ce-4d7b-82ad-40342daff354
# ╠═8959750c-291c-4c2a-a6a1-218f8ed82ab1
# ╟─68c16f9e-8520-41ef-912a-9e61acf2ad1d
# ╠═f49fe2a8-8429-4fd1-98e3-6e75189a3ac2
# ╠═494163b1-cb52-45a9-815e-991633018cd5
# ╠═aed2176b-dfe2-4a38-bc8d-2e5ec6e4531f
# ╟─880db52b-e39c-4578-a6de-ce4a6afa038e
# ╠═6613bd97-7d53-4b8b-a69c-12e040f4aa36
# ╠═bff84e32-1013-44f2-9eda-eb5f6bac4c5a
# ╟─451dd3c6-f740-47df-a5b9-85749caf04fb
# ╠═1ac3d601-0057-475e-a6ac-9522618f7bb9
# ╠═5d407e49-fe05-4748-be36-b232b5ce3ccb
# ╠═de1ac95e-2399-45c1-b337-dcfe800394a2
# ╟─a8f6ece5-9f59-4af1-bcf7-7dafa839bae2
# ╠═267d1d6a-5c35-480e-8425-f000f7f31c7f
# ╟─9a77474e-a853-4369-a6f3-a9e76a683a27
# ╠═a9f499f8-ca32-4310-82d7-f20b790dc04c
# ╠═4eece882-8416-4330-b95a-1cfaeac89aac
# ╠═214d485e-df7e-4a45-b52f-10f854d839a1
# ╠═03622243-d2af-4d1b-a6e1-fd50670e48a8
# ╟─ee208b64-3cf3-428b-a5f6-74793a82fd93
# ╟─5930b9c5-2f85-465b-abe0-223580142d49
# ╠═480bd663-d822-4417-b7f0-6679d26a58b4
# ╠═f498c292-367d-4441-873e-29bec99af2a6
# ╟─56ce9966-939e-48cd-9c70-2e441bfaea58
# ╠═228dd95c-fdab-4a5a-a079-e14238ff5e83
# ╠═0f572e01-a690-4699-8a8c-6684c27085bb
# ╠═ce5b9c0a-4f25-4b8f-bf16-17d05b7bfff4
# ╟─5a8963f1-3ec8-4dbe-b567-df2648219b7f
# ╠═92fca992-b36d-48f3-803e-f92e6463c291
# ╠═fb320178-0d33-4402-b075-7bf7541f94e4
# ╟─feca7b76-fc10-4aca-adff-4cd323dabc74
# ╠═ae9ab557-553c-474a-8865-37807b9270fb
# ╠═39ce8635-28a7-4eee-875e-325e85d142d5
# ╠═42e86a48-1438-4edd-ac2d-f307d534e165
# ╟─0f73d3a7-c6fa-49ee-8fdf-357b0fe94ec8
# ╠═9ac2f9ed-92bb-4540-add5-422df2a3f6da
# ╠═533a77b8-d541-4277-8123-cbfe416fd263
# ╟─891d81f2-f30e-4365-b0ba-b0cebae9f60f
# ╠═efad1a40-b466-4e8e-a4ae-ab3c14efc3f6
# ╠═9cf4a86b-d41f-4814-b14e-0854176e73d4
# ╠═492574dc-5674-462d-9df1-a77e0f445994
# ╠═33b38675-7434-46e5-808c-ad4b0cbb0408

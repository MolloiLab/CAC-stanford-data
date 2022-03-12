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

# ╔═╡ 0c2acaa4-9685-46b5-9880-96b19c06dbc6
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("Revise")
		Pkg.add("PlutoUI")
		Pkg.add("ImageFiltering")
		Pkg.add("Images")
		Pkg.add("ImageMorphology")
		Pkg.add("ImageSegmentation")
		Pkg.add("ImageComponentAnalysis")
		Pkg.add("DataFrames")
		Pkg.add("Statistics")
		Pkg.add("CairoMakie")
		Pkg.add("DataStructures")
		Pkg.add(url="https://github.com/JuliaHealth/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
		Pkg.add(url="https://github.com/Dale-Black/Phantoms.jl")
	end

	using Revise
	using PlutoUI
	using ImageFiltering
	using Images
	using ImageMorphology
	using ImageSegmentation
	using ImageComponentAnalysis
	using DataFrames
	using Statistics
	using CairoMakie
	using DataStructures
	using DICOM
	using DICOMUtils
	using Phantoms
end

# ╔═╡ a1b59a03-1c81-4c8c-81fc-8e03217f3a6f
TableOfContents()

# ╔═╡ e2f49188-76a6-44ce-9e2c-42c7294f4ede
md"""
## Helper functions
"""

# ╔═╡ 1d1ebd11-4d24-4b91-9753-3a1181f63ce0
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ 67a0e55a-9e02-44bc-833f-662230ff2d87
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=25, show_value=true)
end

# ╔═╡ 784504ad-f2cc-4cd7-9aba-2c79d106e4ef
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
	scatter!(label_array[:, 1][indices_lbl], label_array[:, 2][indices_lbl], markersize=0.5, color=:red)
	fig
end

# ╔═╡ 50fe8dda-e030-481a-adf1-0c3f117b35ac
md"""
## `dcm_list_builder`
"""

# ╔═╡ a0dbdc99-d102-436a-9b4d-c00f969708d4
function dcm_list_builder(path)
    dcm_path_list = []
    for (dirpath, dirnames, filenames) in walkdir(path, topdown=true)
        if (dirpath in dcm_path_list) == false
            for filename in filenames
                try
                    tmp_str = string(dirpath, "/", filename)
                    ds = dcm_parse(tmp_str)
                    if (dirpath in dcm_path_list) == false
                        push!(dcm_path_list, dirpath)
					end
				catch
                    nothing
				end
			end
        else
                nothing
		end
	end
    return dcm_path_list
end

# ╔═╡ c4e79b12-3aa9-44b7-84b9-3d464ecc1407
md"""
## `dcm_reader`
"""

# ╔═╡ 406f67d3-cb5b-471f-a7c8-698b84d06090
function dcm_reader(dcm_path)
    dcm_files = []
    for (dirpath, dirnames, filenames) in walkdir(dcm_path, topdown=false)
        for filename in filenames
            try
                if (filename == "DIRFILE") == false   
                    dcm_file = string(dirpath, "/", filename)
                    dcm_parse(dcm_file)
                    push!(dcm_files, dcm_file)
				end
			catch
				nothing
			end
		end
	end

    read_RefDs = true
	local RefDs
    while read_RefDs
        for index in range(1, length(dcm_files))
            try
                RefDs = dcm_parse(dcm_files[index])
                read_RefDs = false
                break
			catch
                nothing
			end
		end
	end

	header = RefDs.meta
	slice_thick_ori = header[(0x0018, 0x0050)]
	rows, cols = Int(header[(0x0028, 0x0010)]), Int(header[(0x0028, 0x0011)])
    
    ConstPixelDims = (rows, cols, length(dcm_files))
    dcm_array = zeros(ConstPixelDims...)

    instances = []    
    for filenameDCM in dcm_files
        try
            ds = dcm_parse(filenameDCM)
			head = ds.meta
			InstanceNumber = head[(0x0020, 0x0013)]
            push!(instances, InstanceNumber)
		catch
            nothing
		end
	end
    
    sort!(instances)

    index = 0
    for filenameDCM in dcm_files
        try
            ds = dcm_parse(filenameDCM)
			head = ds.meta
			InstanceNumber = head[(0x0020, 0x0013)]
			index = findall(x -> x==InstanceNumber, instances)
			pixel_array = head[(0x7fe0, 0x0010)]
            dcm_array[:, :, index] = pixel_array
            index += 1
		catch
            nothing
		end
	end
	
    RescaleSlope = header[(0x0028, 0x1053)]
	RescaleIntercept = header[(0x0028, 0x1052)]
    dcm_array = dcm_array .* RescaleSlope .+ RescaleIntercept
    return RefDs.meta, dcm_array, slice_thick_ori
end

# ╔═╡ 5ec8e16a-bd03-4b82-a3e7-e5464000a5af
md"""
## Load DICOMs
"""

# ╔═╡ d0436ba5-c786-43f5-afe1-27b54a6f2a92
root_path = "/Users/daleblack/Google Drive/Datasets/Canon_Aquilion_One_Vision"

# ╔═╡ 79b38689-0b30-4b71-8439-88f65412a2c2
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ 6f7c7587-2424-4fa2-b722-ab3cbe84dd5a


# ╔═╡ acf36216-755f-4636-b02c-d49cea87fe0f
md"""
# Whole heart mask
"""

# ╔═╡ 948ee148-ee3c-430a-b2d9-172be0ea24f8
md"""
## `find_circle`
"""

# ╔═╡ 4c30f238-89c8-446a-b45f-fdceb5b7313d
function find_circle(point_1, point_2, point_3)
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

# ╔═╡ 42cc247a-b14f-4ef0-b966-c51a46bb5d15
find_circle([309, 309], [312, 200], [155, 155])

# ╔═╡ e798fc03-7a95-44dd-838d-efdf363bfa34
md"""
## `mask_heart`
"""

# ╔═╡ 6557a5ac-8c4c-4060-bf78-3089dceaeff2
"""
    mask_heart(
        header; 
        array_used=nothing, 
        radius_val=95, 
        slice_used_center=nothing
        )
Given a QRM Phantom with a heart insert, this function will create a mask of the whole heart for image-processing purposes.
"""
function mask_heart(header, array_used, slice_used_center; radius_val=95)
    pixel_size = Phantoms.get_pixel_size(header)

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
    for index in 1:(size(central_image, 2) ÷ 2)
        if (
            central_image[center[1] + index, center[2] + index] == 1 &&
            central_image[center[1] + index, center[2] + index + 5] == 1
        )
            point_1 = [center[1] + index, center[2] + index]
            break
        else
            a[center[1] + index, center[2] + index] = 2
        end
    end

    local point_2
    for index in 1:(size(central_image, 2) ÷ 2)
        if (
            central_image[center[1] + index, center[2] - index] == 1 &&
            central_image[center[1] + index, center[2] - index - 5] == 1
        )
            point_2 = [center[1] + index, center[2] - index]
            break
        else
            a[center[1] + index, center[2] - index] = 2
        end
    end

    local point_3
    for index in 1:(size(central_image, 2) ÷ 2)
        if (
            central_image[center[1] - index, center[2] - index] == 1 &&
            central_image[center[1] - index, center[2] - index - 5] == 1
        )
            point_3 = [center[1] - index, center[2] - index]
            break
        else
            a[center[1] - index, center[2] - index] = 2
        end
    end

    center_insert = find_circle(point_1, point_2, point_3)
    rows, cols = Int(header[(0x0028, 0x0010)]), Int(header[(0x0028, 0x0011)])
    Y, X = collect(1:rows), collect(1:cols)'
    dist_from_center = @. sqrt((X - center_insert[2])^2 + (Y - center_insert[1])^2)

    mask = dist_from_center .<= radius[1]
    masked_array = zeros(size(array_used))
    for index in 1:size(array_used, 3)
        masked_array[:, :, index] = array_used[:, :, index] .* mask
    end

    return masked_array, center_insert, mask
end

# ╔═╡ 1f304ccf-5908-454d-9f6b-a414c6d8166b
md"""
# Calcium rod mask
"""

# ╔═╡ 15772146-fc70-4886-a78e-70056ee6b1da
md"""
## `get_calcium_slices`
"""

# ╔═╡ f96252ea-f133-4cb3-8a4a-61870c2d65b8
"""
	get_calcium_slices(dcm_array, header; calcium_threshold=130)

Returns the slices that contain the calcium vessel inserts `slice_dict` 
and the slices that contain the calcium calibration rod insert `large_index`.

The calcium rod slice `large_index` usually omits the first and last slice 
of the rod. So to account for all of the slices containing the calcium rod, 
one would want a range like so: `(large_index - 1):(large_index + 1)`
"""
function get_calcium_slices(dcm_array, header; calcium_threshold=130)
    array = copy(dcm_array)
    array = Int.(array .> (1.1 * calcium_threshold))

    pixel_size = Phantoms.get_pixel_size(header)
    CCI_5mm_num_pixels = Int(round(π * (5 / 2)^2 / pixel_size[1]^2))
    cal_rod_num_pixels = Int(round(π * (20 / 2)^2 / pixel_size[1]^2))

    kern = Int.(round(5 / pixel_size[1]))
    if kern % 2 == 0
        kern += 1
    end

    slice_dict = Dict()
    large_index = []
    cal_rod_dict = Dict()
    for idx in 1:size(array, 3)
        array_filtered = mapwindow(median, array[:, :, idx], (kern, kern))
        components = ImageComponentAnalysis.label_components(array_filtered)
        a1 = analyze_components(components, BasicMeasurement(; area=true, perimeter=true))
        a2 = analyze_components(components, BoundingBox(; box_area=true))
        df = leftjoin(a1, a2; on=:l)

        count_5mm = 0
        count = 0
        for row in eachrow(df)
            count += 1
            df_area = Int(round(row[:area]))

            r1_1 = Int(round(CCI_5mm_num_pixels * 0.6))
            r1_2 = Int(round(CCI_5mm_num_pixels * 1.5))
            r2_1 = Int(round(cal_rod_num_pixels * 0.7))
            r2_2 = Int(round(cal_rod_num_pixels * 1.3))

            if df_area in r1_1:r1_2
                count_5mm += 1
            elseif df_area in r2_1:r2_2
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

            range1 = round(0.7 * y_dist):round(1.2 * y_dist)
            if ((x_dist in range1) == false) ||
                ((round(0.7 * y_dist) == 0) && (round(1.2 * y_dist) == 0))
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

# ╔═╡ 7bc9ec3a-e5d9-4f92-9a7a-67cb143e8d35
md"""
## `get_calcium_center_slices`
"""

# ╔═╡ 307e3875-d9d1-44e1-9106-f43368d3fe01
"""
	get_calcium_center_slices(dcm_array, slice_dict, large_index)

Returns the slices that contain the calcium vessel inserts `slice_dict`
and the center of the calcium calibration rod slice `flipped_index`
"""
function get_calcium_center_slices(dcm_array, slice_dict, large_index)
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
                deleteat!(large_index, findall(x -> x == element2, large_index))
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
                deleteat!(large_index, findall(x -> x == element2, large_index))
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

# ╔═╡ 44063684-19f7-4d3d-bb05-5d223b81e0c6
md"""
## `poppable_keys`
"""

# ╔═╡ 23000f9d-f346-4dc7-a9b2-e48d73d9988f
"""
    poppable_keys(flipped, flipped_index, header, slice_dict)

No idea what the function is
"""
function poppable_keys(flipped, flipped_index, header, slice_dict)
    SliceThickness = header[(0x0018, 0x0050)]
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

# ╔═╡ 01e5c044-4fbc-4c3a-a4a6-bb7efecb4319
md"""
## `compute_CCI`
"""

# ╔═╡ f58c4041-63e3-47b1-a77a-7882397d4223
function compute_CCI(dcm_array, header, slice_dict, flipped; calcium_threshold=130)
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
    quality_slice = Int.(round(slice_CCI - flipped * (20 / SliceThickness)))

    cal_rod_slice = slice_CCI + (flipped * Int(30 / SliceThickness))
    
    return calcium_image, slice_CCI, quality_slice, cal_rod_slice
end

# ╔═╡ 6d619783-36f5-4c03-84df-36eeb47e9eb3
md"""
## `mask_rod`
"""

# ╔═╡ fdd32c44-c731-476f-9bce-8f522a5aa204
"""
	mask_rod(dcm_array, header; calcium_threshold=130)
"""
function mask_rod(dcm_array, header; calcium_threshold=130)
    slice_dict, large_index = get_calcium_slices(
        dcm_array, header; calcium_threshold=calcium_threshold
    )
    slice_dict, flipped, flipped_index = get_calcium_center_slices(
        dcm_array, slice_dict, large_index
    )
    slice_dict = poppable_keys(flipped, flipped_index, header, slice_dict)
    calcium_image, slice_CCI, quality_slice, cal_rod_slice = compute_CCI(
        dcm_array, header, slice_dict, flipped; calcium_threshold=calcium_threshold
    )
    return calcium_image, slice_CCI, quality_slice, cal_rod_slice
end

# ╔═╡ a1b38df1-e2a2-4ab3-8681-b90a3655525c
md"""
# Calcium inserts mask
"""

# ╔═╡ 64afcb42-64f0-404b-9fad-66a872119fff
md"""
## `angle_calc`
"""

# ╔═╡ ae739339-01b9-4486-bc9b-43c9424ee96d
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

# ╔═╡ f0b5027e-fb01-4282-8ae1-863179c68f32
angle_calc(4, 3)

# ╔═╡ 6822d54c-ec32-4a2e-a902-aa06e045ef6d
md"""
## `create_circular_mask`
"""

# ╔═╡ 9de3903b-6cfb-4db4-9be6-05f9fd36ddd5
function create_circular_mask(h, w, center_circle, radius_circle)
	Y, X = collect(1:h), collect(1:w)'
    dist_from_center = sqrt.((X .- center_circle[1]).^2 .+ (Y .- center_circle[2]).^2)

    mask = dist_from_center .<= radius_circle
    
    return mask
end

# ╔═╡ bda52033-230f-4141-b753-53c60e907d07
mask1 = create_circular_mask(40, 40, [20, 20], 1);

# ╔═╡ a0f781d2-8af6-4425-8507-04e36bb320df
heatmap(mask1, colormap=:grays)

# ╔═╡ b8091995-6a67-4a18-8abb-7da0167758ea
md"""
## `calc_output`
"""

# ╔═╡ adb53900-090e-4888-abf0-98c07c41fe1f
"""
	calc_output(dcm_array, CCI_slice, calcium_threshold=130, comp_connect=trues(3, 3))

Calculate the output of a dcm_array
"""
function calc_output(
    dcm_array, header, CCI_slice, calcium_threshold=130, comp_connect=trues(3, 3)
)
    # Actual scoring for CCI insert
    # First step is to remove slices without calcium from arrays
    PixelSpacing = Phantoms.get_pixel_size(header)
    SliceThickness = header[(0x0018, 0x0050)]
    CCI_min = Int((CCI_slice - round(5 / SliceThickness, RoundUp)) - 1)
    CCI_max = Int((CCI_slice + round(5 / SliceThickness, RoundUp)) + 1)
    central_CCI = Int(round(CCI_max - CCI_min) / 2)

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
    CCI_array_binary = Int.(CCI_array_binary .> 1.0 * calcium_threshold)
    inp =
        CCI_array_binary[:, :, central_CCI - 1] +
        CCI_array_binary[:, :, central_CCI] +
        CCI_array_binary[:, :, central_CCI + 1]
    components = ImageComponentAnalysis.label_components(inp, comp_connect)
    a1 = analyze_components(components, BasicMeasurement(; area=true, perimeter=true))
    a2 = analyze_components(components, BoundingBox(; box_area=true))
    df = leftjoin(a1, a2; on=:l)
    centroids = []
    for row in eachrow(df)
        indices = row[:box_indices]
        x_point = ((indices[1][end] - indices[1][1]) ÷ 2) + indices[1][1]
        y_point = ((indices[2][end] - indices[2][1]) ÷ 2) + indices[2][1]
        push!(centroids, (x_point, y_point))
    end

    centroids = deleteat!(centroids, 1)

    i1 = mapwindow(
        median, CCI_array_binary[:, :, central_CCI - 1], (image_kernel, image_kernel)
    )
    i2 = mapwindow(
        median, CCI_array_binary[:, :, central_CCI], (image_kernel, image_kernel)
    )
    i3 = mapwindow(
        median, CCI_array_binary[:, :, central_CCI + 1], (image_kernel, image_kernel)
    )

    image_for_center = i1 + i2 + i3

    components2 = ImageComponentAnalysis.label_components(image_for_center, comp_connect)
    components2 = Int.(components2 .> 0)
    components2 = ImageComponentAnalysis.label_components(components2, comp_connect)

    b1 = analyze_components(components2, BasicMeasurement(; area=true, perimeter=true))
    b2 = analyze_components(components2, BoundingBox(; box_area=true))
    df2 = leftjoin(b1, b2; on=:l)
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

# ╔═╡ 9f3d84bd-c9f3-4130-8b83-da6d668ffde9
md"""
## `center_points`
"""

# ╔═╡ 053559b6-230a-4026-bad9-41598bbdeead
"""
	center_points(output, header, tmp_center, CCI_slice)

Function ...
"""
function center_points(dcm_array, output, header, tmp_center, CCI_slice)
    PixelSpacing = Phantoms.get_pixel_size(header)
    rows, cols = Int(header[(0x0028, 0x0010)]), Int(header[(0x0028, 0x0011)])
    sizes = []
    for row in eachrow(output[3])
        area = row[:area]
        append!(sizes, area)
    end

    centroids = output[4]
    largest = Dict()
    for index in 1:length(centroids)
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
    radius = round(2.5 / PixelSpacing[1], RoundUp)
    for key in largest
        tmp_arr = create_circular_mask(rows, cols, (key[2][2], key[2][1]), radius)
        tmp_arr = @. abs(tmp_arr * dcm_array[:, :, CCI_slice]) +
            abs(tmp_arr * dcm_array[:, :, CCI_slice - 1]) +
            abs(tmp_arr * dcm_array[:, :, CCI_slice + 1])
        tmp_arr = @. ifelse(tmp_arr == 0, missing, tmp_arr)
        max_dict[key[1]] = median(skipmissing(tmp_arr))
    end
    large1_index, large1_key = maximum(zip(values(max_dict), keys(max_dict)))
    pop!(max_dict, large1_key)
    large2_index, large2_key = maximum(zip(values(max_dict), keys(max_dict)))
    pop!(max_dict, large2_key)
    large3_index, large3_key = maximum(zip(values(max_dict), keys(max_dict)))

    center1 = largest[large1_key]
    center2 = largest[large2_key]
    center3 = largest[large3_key]

    center = find_circle(center1, center2, center3)
    return center, center1, center2, center3
end

# ╔═╡ 077991f2-4b33-4f34-8d3d-ce6827ed830a
md"""
## `calc_centers`
"""

# ╔═╡ 223086c6-dfe9-44ba-b21f-7bec3028e10e
md"""
## `mask_inserts`
"""

# ╔═╡ 63863a87-d279-48c8-b957-392372dc85ea
md"""
## Visualize
"""

# ╔═╡ c97d23e4-1c41-4e81-8f6e-61ad7ff0396b
"""
	calc_centers(dcm_array, output, header, tmp_center, CCI_slice)

Function ...
"""
function calc_centers(dcm_array, output, header, tmp_center, CCI_slice)
    PixelSpacing = Phantoms.get_pixel_size(header)
    center, center1, center2, center3 = center_points(
        dcm_array, output, header, tmp_center, CCI_slice
    )
    centers = Dict()
    for size_index4 in (center1, center2, center3)
        center_index = size_index4
        side_x = abs(center[1] - center_index[1]) - 5
        side_y = abs(center[2] - center_index[2]) - 5

        angle = angle_calc(side_x, side_y)
        if (center_index[1] < center[1] && center_index[2] < center[2])
            medium_calc = [
                center_index[1] + (12.5 / PixelSpacing[1]) * sin(angle),
                (center_index[2] + (12.5 / PixelSpacing[2]) * cos(angle)),
            ]
            low_calc = [
                center_index[1] + (25 / PixelSpacing[1]) * sin(angle),
                (center_index[2] + (25 / PixelSpacing[2]) * cos(angle)),
            ]
        elseif (center_index[1] < center[1] && center_index[2] > center[2])
            medium_calc = [
                center_index[1] + (12.5 / PixelSpacing[1]) * sin(angle),
                (center_index[2] - (12.5 / PixelSpacing[2]) * cos(angle)),
            ]
            low_calc = [
                center_index[1] + (25 / PixelSpacing[1]) * sin(angle),
                (center_index[2] - (25 / PixelSpacing[2]) * cos(angle)),
            ]
        elseif (center_index[1] > center[1] && center_index[2] < center[2])
            medium_calc = [
                center_index[1] - (12.5 / PixelSpacing[1]) * sin(angle),
                (center_index[2] + (12.5 / PixelSpacing[2]) * cos(angle)),
            ]
            low_calc = [
                center_index[1] - (25 / PixelSpacing[1]) * sin(angle),
                (center_index[2] + (25 / PixelSpacing[2]) * cos(angle)),
            ]
        elseif (center_index[1] > center[1] && center_index[2] > center[2])
            medium_calc = [
                center_index[1] - (12.5 / PixelSpacing[1]) * sin(angle),
                (center_index[2] - (12.5 / PixelSpacing[2]) * cos(angle)),
            ]
            low_calc = [
                center_index[1] - (25 / PixelSpacing[1]) * sin(angle),
                (center_index[2] - (25 / PixelSpacing[2]) * cos(angle)),
            ]
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

# ╔═╡ e4649f2b-a6da-4038-a15a-6ee0272a8572
"""
	mask_inserts(
		dcm_array, masked_array, header, CCI_slice, center_insert; 
		calcium_threshold=130, comp_connect=trues(3, 3)
		)

Function ...
"""
function mask_inserts(
    dcm_array,
    masked_array,
    header,
    CCI_slice,
    center_insert;
    calcium_threshold=130,
    comp_connect=trues(3, 3),
)
    output = calc_output(masked_array, header, CCI_slice, calcium_threshold, comp_connect)
    insert_centers = calc_centers(dcm_array, output, header, center_insert, CCI_slice)

    PixelSpacing = Phantoms.get_pixel_size(header)
    rows, cols = Int(header[tag"Rows"]), Int(header[tag"Columns"])

    mask_L_HD = create_circular_mask(
        cols, 
		rows, 
		insert_centers[:Large_HD], 
		(round(5 / PixelSpacing[1], RoundUp) / 2) + 1
	)
    mask_L_MD = create_circular_mask(
        cols, 
		rows, 
		insert_centers[:Large_MD], 
		(round(5 / PixelSpacing[1], RoundUp) / 2) + 1
    )
    mask_L_LD = create_circular_mask(
        cols, 
		rows, 
		insert_centers[:Large_LD], 
		(round(5 / PixelSpacing[1], RoundUp) / 2) + 1
    )
    mask_M_HD = create_circular_mask(
        cols,
        rows,
        insert_centers[:Medium_HD],
        (round(3 / PixelSpacing[1], RoundUp) / 2) + 1,
    )
    mask_M_MD = create_circular_mask(
        cols,
        rows,
        insert_centers[:Medium_MD],
        (round(3 / PixelSpacing[1], RoundUp) / 2) + 1,
    )
    mask_M_LD = create_circular_mask(
        cols,
        rows,
        insert_centers[:Medium_LD],
        (round(3 / PixelSpacing[1], RoundUp) / 2) + 1,
    )
    mask_S_HD = create_circular_mask(
        cols, 
		rows, 
		insert_centers[:Small_HD], 
		(round(1 / PixelSpacing[1], RoundUp) / 2) + 1
    )
    mask_S_MD = create_circular_mask(
        cols, 
		rows, 
		insert_centers[:Small_MD], 
		(round(1 / PixelSpacing[1], RoundUp) / 2) + 1
    )
    mask_S_LD = create_circular_mask(
        cols, 
		rows, 
		insert_centers[:Small_LD], 
		(round(1 / PixelSpacing[1], RoundUp) / 2) + 1
    )
 
    return transpose(mask_L_HD), transpose(mask_M_HD), transpose(mask_S_HD), transpose(mask_L_MD), transpose(mask_M_MD), transpose(mask_S_MD), transpose(mask_L_LD), transpose(mask_M_LD), transpose(mask_S_LD)
end

# ╔═╡ 4f49aca7-ae33-4195-8b7e-0c715e6a69d6
header, dcm_array, slice_thick_ori = dcm_reader(dcm_path_list[1]);

# ╔═╡ 2f5b300f-7632-43f0-ba41-c79564924450
heatmap(transpose(dcm_array[:, :, 8]), colormap=:grays)

# ╔═╡ 06a491cd-155b-40d5-835f-c40602188d3b
maximum(dcm_array)

# ╔═╡ 0726bc84-32dc-48f2-a4e9-c0b1658dfccc
minimum(dcm_array)

# ╔═╡ e4fa8e2f-5c09-4e98-b996-2286ff500f42
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ c2cae050-074e-42eb-bf1c-19e4796ec96f
center_insert

# ╔═╡ f3301be4-1296-4e0b-b08e-c3bbd3abf409
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig2
end

# ╔═╡ 296a28ec-f23f-424a-a817-a1534ec781dc
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 23]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig3
end

# ╔═╡ aea5a28e-c234-4715-a84e-ec2d31e61904
@bind a2 PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 24a3a30d-d74b-4e6b-9a1e-ce0a32c620f5
heatmap(transpose(masked_array[:, :, a2]), colormap=:grays)

# ╔═╡ f44a7afb-2f42-405e-a07b-65f6ea654759
center_insert

# ╔═╡ 24bbea22-7a34-46ff-90df-c92c8f507616
masked_array_slices = masked_array[:, :, 25:26];

# ╔═╡ b9c8beaa-03b8-4ee1-880f-7adf747e8eeb
@bind a3 PlutoUI.Slider(1:size(masked_array, 3), default=26, show_value=true)

# ╔═╡ 438f11a7-d5a4-4406-adf5-7f80931cf2d3
heatmap(transpose(masked_array[:, :, a3]), colormap=:grays)

# ╔═╡ bfbf586e-9f6d-41de-8ede-ffa2a9989871
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 15]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig
end

# ╔═╡ 704ed93d-e481-4e5d-9287-89dc3f0e4288
slice_dict, large_index = get_calcium_slices(masked_array, header);

# ╔═╡ 02b85a5a-a59c-4c3f-80a4-c8a884a604e1
slice_dict, large_index

# ╔═╡ 8f76519b-a37e-4982-8be8-6fee4d3ea050
slice_dict2, flipped, flipped_index = get_calcium_center_slices(masked_array, slice_dict, large_index)

# ╔═╡ 0030a80f-f590-4f88-94c1-1329328f6f01
pop_keys = poppable_keys(flipped, flipped_index, header, slice_dict)

# ╔═╡ 39ddffe8-2677-4dc2-aafd-8d2ba80d60f8
calcium_image1, slice_CCI1, quality_slice1, cal_rod_slice1 = compute_CCI(masked_array, header, slice_dict, flipped);

# ╔═╡ e40391a9-aab9-4a8f-883b-52fc30fcf033
slice_CCI1, quality_slice1, cal_rod_slice1

# ╔═╡ b1f15e4c-3459-4ada-b218-d9f97befbf5d
@bind a PlutoUI.Slider(1:size(calcium_image1, 3), default=10, show_value=true)

# ╔═╡ 036e3e9d-c565-48d7-85c8-880eef44ec1a
begin
	fig4 = Figure()
	
	ax4 = Makie.Axis(fig4[1, 1])
	heatmap!(transpose(calcium_image1[:, :, a]), colormap=:grays)
	ax5 = Makie.Axis(fig4[1,2])
	heatmap!(transpose(dcm_array[:, :, a]), colormap=:grays)
	fig4
end

# ╔═╡ 70d00362-4c39-4cec-93a9-379ff789be4c
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header);

# ╔═╡ c8e743e1-ecf2-4626-8f7f-c4f0fdbd8e34
slice_CCI

# ╔═╡ b29bcac4-1cb1-4770-affb-db3343fa2e89
@bind b PlutoUI.Slider(1:size(calcium_image, 3), default=10, show_value=true)

# ╔═╡ 41a56afe-2e85-47dd-9f64-fa19073f5a13
heatmap(transpose(calcium_image[:, :, b]), colormap=:grays)

# ╔═╡ c3df2597-b162-45c9-b64f-017354f1948f
heatmap(transpose(calcium_image[:, :, slice_CCI]), colormap=:grays)

# ╔═╡ 10f16ad5-48ec-42e0-bd81-bb3598fdd492
output = calc_output(masked_array, header, slice_CCI1);

# ╔═╡ 534157ba-5705-4c3a-b43a-df640603bf54
output

# ╔═╡ 2c97a166-9026-4ec7-9c67-b8d42395a135
heatmap(transpose(output[2]))

# ╔═╡ a9b10635-4769-4638-a374-dcb89056e422
begin
	label_map = output[2]
	label_map = cat(label_map, label_map, dims=3)
end;

# ╔═╡ 1e868ff4-957a-4fe3-bd17-7e00db3c4d19
function centerpoint_to_mask(mask, center)
	center_mask = zeros(size(label_map))
	center_mask[CartesianIndex(center[1], center[2]), 1] = 1
	center_mask[CartesianIndex(center[1], center[2]), 2] = 1
	center_mask = dilate(dilate(Bool.(center_mask)))
	return center_mask
end

# ╔═╡ 945977a6-38de-4908-a854-6620ef038abd
center, center1, center2, center3 = center_points(dcm_array, output, header, center_insert, slice_CCI1)

# ╔═╡ 81256b8b-501f-44fb-be0f-7152ecf411d4
center1_mask = centerpoint_to_mask(label_map, center1);

# ╔═╡ 7994035a-0397-4061-9453-8e9e2ba31baa
@bind a10 overlay_mask_bind(center1_mask)

# ╔═╡ 104af49f-351d-4e12-ad13-3e7c362d5d25
overlay_mask_plot(label_map, center1_mask, a10, "masks overlayed")

# ╔═╡ 6658f160-0460-41cf-a54e-e384a599dc74
dict = calc_centers(dcm_array, output, header, center_insert, slice_CCI1)

# ╔═╡ e9684180-e451-49be-8e29-7f74e96c9bed
dict[:Large_HD], dict[:Medium_HD], dict[:Small_HD], dict[:Large_MD], dict[:Medium_MD], dict[:Small_MD]

# ╔═╡ c61ea8d5-aa39-4a8e-ab28-a495184defd9
Medium_MD_mask = centerpoint_to_mask(label_map, dict[:Small_MD]);

# ╔═╡ 96bb9e54-78be-453b-9ff3-8d8bb309f2bf
@bind a11 overlay_mask_bind(Medium_MD_mask)

# ╔═╡ 3e509ca0-af53-4b92-978e-eedee1cdac1a
overlay_mask_plot(label_map, Medium_MD_mask, a11, "masks overlayed")

# ╔═╡ 91ffd750-7627-45a2-93ab-fe4495cbc351
overlay_mask_plot(masked_array_slices, Medium_MD_mask, a11, "masks overlayed")

# ╔═╡ 65d20ce2-7af0-462f-8dcc-47caf963ec18
mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(dcm_array, masked_array, header, slice_CCI, center_insert);

# ╔═╡ 1a266b68-e672-41ba-b75c-8bf735675518
masks = mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD;

# ╔═╡ 563078ac-32ec-4ff7-ae93-f654dab7d8b4
heatmap(transpose(mask_L_HD), colormap=:grays)

# ╔═╡ c2197a05-dbec-4976-bf1b-1b9ddceedba9
begin
	masks_3D = Array{Bool}(undef, size(dcm_array))
	for z in 1:size(dcm_array, 3)
		masks_3D[:, :, z] = masks
	end
end;

# ╔═╡ d9bf252a-6249-46d3-81df-ed72a9e0065a
@bind v1 overlay_mask_bind(masks_3D)

# ╔═╡ 393954fc-4598-49a6-b8ce-24e7c1ed6568
function overlay_mask_plot_t(array, mask, var, title::AbstractString)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	indices_lbl = findall(x -> x == zs[var], label_array[:,3])
	
	fig = Figure()
	ax = Makie.Axis(fig[1, 1])
	ax.title = title
	heatmap!(transpose(array[:, :, zs[var]]), colormap=:grays)
	scatter!(label_array[:, 1][indices_lbl], label_array[:, 2][indices_lbl], markersize=0.5, color=:red)
	fig
end

# ╔═╡ 723917fd-79fc-45b9-9129-68f537ab0b78
begin
	masks_3Dt = Array{Bool}(undef, size(dcm_array))
	for z in 1:size(dcm_array, 3)
		masks_3Dt[:, :, z] = transpose(masks)
	end
end;

# ╔═╡ ce6d8cc8-24c9-4792-9baa-b79e43865343
f = overlay_mask_plot_t(dcm_array, masks_3Dt, v1, "masks overlayed")

# ╔═╡ 8e6a02ce-de72-4a34-bfe2-a136a0638af5
save("physical_figure.png", f, px_per_unit = 2)

# ╔═╡ 96d9a243-1d85-4ef1-b4be-90f25c3e7cea
md"""
# Test
"""

# ╔═╡ Cell order:
# ╠═0c2acaa4-9685-46b5-9880-96b19c06dbc6
# ╠═a1b59a03-1c81-4c8c-81fc-8e03217f3a6f
# ╟─e2f49188-76a6-44ce-9e2c-42c7294f4ede
# ╟─1d1ebd11-4d24-4b91-9753-3a1181f63ce0
# ╟─67a0e55a-9e02-44bc-833f-662230ff2d87
# ╟─784504ad-f2cc-4cd7-9aba-2c79d106e4ef
# ╟─1e868ff4-957a-4fe3-bd17-7e00db3c4d19
# ╟─50fe8dda-e030-481a-adf1-0c3f117b35ac
# ╠═a0dbdc99-d102-436a-9b4d-c00f969708d4
# ╟─c4e79b12-3aa9-44b7-84b9-3d464ecc1407
# ╠═406f67d3-cb5b-471f-a7c8-698b84d06090
# ╟─5ec8e16a-bd03-4b82-a3e7-e5464000a5af
# ╠═d0436ba5-c786-43f5-afe1-27b54a6f2a92
# ╠═79b38689-0b30-4b71-8439-88f65412a2c2
# ╠═6f7c7587-2424-4fa2-b722-ab3cbe84dd5a
# ╠═2f5b300f-7632-43f0-ba41-c79564924450
# ╠═06a491cd-155b-40d5-835f-c40602188d3b
# ╠═0726bc84-32dc-48f2-a4e9-c0b1658dfccc
# ╟─acf36216-755f-4636-b02c-d49cea87fe0f
# ╟─948ee148-ee3c-430a-b2d9-172be0ea24f8
# ╠═4c30f238-89c8-446a-b45f-fdceb5b7313d
# ╠═42cc247a-b14f-4ef0-b966-c51a46bb5d15
# ╟─e798fc03-7a95-44dd-838d-efdf363bfa34
# ╠═6557a5ac-8c4c-4060-bf78-3089dceaeff2
# ╠═e4fa8e2f-5c09-4e98-b996-2286ff500f42
# ╠═c2cae050-074e-42eb-bf1c-19e4796ec96f
# ╟─bfbf586e-9f6d-41de-8ede-ffa2a9989871
# ╟─f3301be4-1296-4e0b-b08e-c3bbd3abf409
# ╟─296a28ec-f23f-424a-a817-a1534ec781dc
# ╟─aea5a28e-c234-4715-a84e-ec2d31e61904
# ╠═24a3a30d-d74b-4e6b-9a1e-ce0a32c620f5
# ╟─1f304ccf-5908-454d-9f6b-a414c6d8166b
# ╟─15772146-fc70-4886-a78e-70056ee6b1da
# ╠═f96252ea-f133-4cb3-8a4a-61870c2d65b8
# ╠═704ed93d-e481-4e5d-9287-89dc3f0e4288
# ╠═02b85a5a-a59c-4c3f-80a4-c8a884a604e1
# ╟─7bc9ec3a-e5d9-4f92-9a7a-67cb143e8d35
# ╠═307e3875-d9d1-44e1-9106-f43368d3fe01
# ╠═8f76519b-a37e-4982-8be8-6fee4d3ea050
# ╟─44063684-19f7-4d3d-bb05-5d223b81e0c6
# ╠═23000f9d-f346-4dc7-a9b2-e48d73d9988f
# ╠═0030a80f-f590-4f88-94c1-1329328f6f01
# ╟─01e5c044-4fbc-4c3a-a4a6-bb7efecb4319
# ╠═f58c4041-63e3-47b1-a77a-7882397d4223
# ╠═39ddffe8-2677-4dc2-aafd-8d2ba80d60f8
# ╠═e40391a9-aab9-4a8f-883b-52fc30fcf033
# ╟─b1f15e4c-3459-4ada-b218-d9f97befbf5d
# ╠═036e3e9d-c565-48d7-85c8-880eef44ec1a
# ╟─6d619783-36f5-4c03-84df-36eeb47e9eb3
# ╠═fdd32c44-c731-476f-9bce-8f522a5aa204
# ╠═70d00362-4c39-4cec-93a9-379ff789be4c
# ╠═c8e743e1-ecf2-4626-8f7f-c4f0fdbd8e34
# ╟─b29bcac4-1cb1-4770-affb-db3343fa2e89
# ╠═41a56afe-2e85-47dd-9f64-fa19073f5a13
# ╠═c3df2597-b162-45c9-b64f-017354f1948f
# ╟─a1b38df1-e2a2-4ab3-8681-b90a3655525c
# ╟─64afcb42-64f0-404b-9fad-66a872119fff
# ╠═ae739339-01b9-4486-bc9b-43c9424ee96d
# ╠═f0b5027e-fb01-4282-8ae1-863179c68f32
# ╟─6822d54c-ec32-4a2e-a902-aa06e045ef6d
# ╠═9de3903b-6cfb-4db4-9be6-05f9fd36ddd5
# ╠═bda52033-230f-4141-b753-53c60e907d07
# ╠═a0f781d2-8af6-4425-8507-04e36bb320df
# ╟─b8091995-6a67-4a18-8abb-7da0167758ea
# ╠═adb53900-090e-4888-abf0-98c07c41fe1f
# ╠═10f16ad5-48ec-42e0-bd81-bb3598fdd492
# ╠═534157ba-5705-4c3a-b43a-df640603bf54
# ╠═2c97a166-9026-4ec7-9c67-b8d42395a135
# ╟─9f3d84bd-c9f3-4130-8b83-da6d668ffde9
# ╠═053559b6-230a-4026-bad9-41598bbdeead
# ╠═f44a7afb-2f42-405e-a07b-65f6ea654759
# ╠═945977a6-38de-4908-a854-6620ef038abd
# ╠═a9b10635-4769-4638-a374-dcb89056e422
# ╠═81256b8b-501f-44fb-be0f-7152ecf411d4
# ╟─7994035a-0397-4061-9453-8e9e2ba31baa
# ╟─104af49f-351d-4e12-ad13-3e7c362d5d25
# ╟─077991f2-4b33-4f34-8d3d-ce6827ed830a
# ╠═6658f160-0460-41cf-a54e-e384a599dc74
# ╠═e9684180-e451-49be-8e29-7f74e96c9bed
# ╠═c61ea8d5-aa39-4a8e-ab28-a495184defd9
# ╠═96bb9e54-78be-453b-9ff3-8d8bb309f2bf
# ╠═3e509ca0-af53-4b92-978e-eedee1cdac1a
# ╠═24bbea22-7a34-46ff-90df-c92c8f507616
# ╠═91ffd750-7627-45a2-93ab-fe4495cbc351
# ╟─223086c6-dfe9-44ba-b21f-7bec3028e10e
# ╠═e4649f2b-a6da-4038-a15a-6ee0272a8572
# ╠═65d20ce2-7af0-462f-8dcc-47caf963ec18
# ╠═1a266b68-e672-41ba-b75c-8bf735675518
# ╠═563078ac-32ec-4ff7-ae93-f654dab7d8b4
# ╠═b9c8beaa-03b8-4ee1-880f-7adf747e8eeb
# ╠═438f11a7-d5a4-4406-adf5-7f80931cf2d3
# ╟─63863a87-d279-48c8-b957-392372dc85ea
# ╠═c2197a05-dbec-4976-bf1b-1b9ddceedba9
# ╠═c97d23e4-1c41-4e81-8f6e-61ad7ff0396b
# ╠═4f49aca7-ae33-4195-8b7e-0c715e6a69d6
# ╠═d9bf252a-6249-46d3-81df-ed72a9e0065a
# ╠═393954fc-4598-49a6-b8ce-24e7c1ed6568
# ╠═723917fd-79fc-45b9-9129-68f537ab0b78
# ╠═ce6d8cc8-24c9-4792-9baa-b79e43865343
# ╠═8e6a02ce-de72-4a34-bfe2-a136a0638af5
# ╠═96d9a243-1d85-4ef1-b4be-90f25c3e7cea

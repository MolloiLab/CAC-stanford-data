### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 14b6c3e4-48dd-11ec-336c-6918bf024f43
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("DICOM")
	end
	
	using PlutoUI
	using DICOM
end

# ╔═╡ 9506ae45-f63f-4fd6-8467-0597cbde7006
TableOfContents()

# ╔═╡ 89d2d135-7669-40be-bba1-24a324cdf14e
md"""
Do the dcm specific functions last

Compare each function to the corresponding python function and test that the return values are the same
"""

# ╔═╡ 3f313b73-2785-4299-a1e9-bc5d4e1ad422
md"""
## `dcm_reader`

TODO: Last
"""

# ╔═╡ 0fbe9f2a-fea0-45fc-b429-85afc5db3341
# tst_path = "/Users/daleblack/Google Drive/Datasets/Canon_Aquilion_One_Vision/"

# ╔═╡ 482d8669-c750-450b-aae8-03afc2a90007
md"""
## `dcm_list_builder`

TODO: Last
"""

# ╔═╡ 9e825e49-f9b3-4897-8dc7-755bb4e99b57
# function dcm_list_builder(path, test_text = "")
#     # function to get list of dcm_files from dcm directory
#     dcm_path_list = []
# 	# dir = readdir(path)
#     for (dirpath, dirnames, filenames) in walkdir(path)
# 		@show dirpath 
# 		@show dirnames
# 		@show filenames
#         if dirpath not in dcm_path_list
#    #          for filename in filenames
#    #              try
#    #                  tmp_str = str(os.path.join(dirpath, filename))
#    #                  ds = pydicom.read_file(tmp_str, stop_before_pixels = True)
#    #                  if dirpath not in dcm_path_list
#    #                      dcm_path_list.append(dirpath)
# 			# 		end
# 			# 	catch
#    #                  pass
# 			# 	end
# 			# end
# 		else
#             pass
# 		end
# 	end
#     return dcm_path_list
# end

# ╔═╡ de7805f3-e670-4fe4-beb9-cfcf6c1bc021
# with_terminal() do
# 	dcm_list_builder(tst_path)
# end

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

# ╔═╡ dfcda44d-5b2c-4a50-9fb2-44cb70b5488e
function dcm_masked(array_used = None, radius_val = 95, slice_used_center = None)
    try
        pixel_size = header.PixelSpacing[0]
	catch
        FOV = header.ReconstructionDiameter
        matrix_size = header.Rows
    
        pixel_size = FOV / matrix_size
	end
    
    radius = (radius_val/2) / pixel_size
    
    central_image = array_used[:,:,slice_used_center].copy()
    
    central_image[central_image > -200] = 0
    central_image[central_image != 0] = 1

    image_kernel = math.ceil(5 / header.PixelSpacing[0])
    if image_kernel % 2 == 0
        image_kernel += 1
	end
    central_image = scipy.signal.medfilt2d(central_image, image_kernel)
    
    # plt.imshow(central_image)
    # plt.show()

    center = [int(array_used.shape[0] / 2), int(array_used.shape[1] / 2)]
    
    a = central_image.copy()
    for index in range(int(array_used.shape[1] / 2))
        if (central_image[center[0] + index, center[1] + index] == 1 &&
            central_image[center[0] + index, center[1] + index + 5] == 1)
            point_1 = [center[0] + index, center[1] + index]
            break
        else
            a[center[0] + index, center[1] + index] = 2
            pass
		end
	end
    
    for index in range(int(array_used.shape[1] / 2))
        if (central_image[center[0] + index, center[1] - index] == 1 &&
            central_image[center[0] + index, center[1] - index - 5] == 1)
            point_2 = [center[0] + index, center[1] - index]
            break
        else
            a[center[0] + index, center[1] - index] = 2
            pass
		end
	end
        
    # x_iter = range(center[0], int(array_used.shape[0]), 1)
    # y_iter = range(center[1], int(array_used.shape[1]), 2)
    # for tmp_x, tmp_y in zip(x_iter, y_iter):
    #     if (central_image[tmp_x, tmp_y] == 1 and central_image[tmp_x + 5, tmp_y + 5] == 1):
    #         point_3 = [tmp_x, tmp_y]
    #         break
    #     else:
    #         a[tmp_x, tmp_y] = 2
    #         pass
        
        
    for index in range(int(array_used.shape[1] / 2))
        if (central_image[center[0] - index, center[1] - index] == 1 &&
            central_image[center[0] - index, center[1] - index - 5] == 1)
            point_3 = [center[0] - index, center[1] - index]
            break
        else
            a[center[0] - index, center[1] - index] = 2
            pass
		end
	end
        
    # plt.imshow(a)
    # plt.show()
    print(point_1, point_2, point_3)
    center_insert = findCircle(point_1, point_2, point_3)

    Y, X = np.ogrid[:header.Rows, :header.Columns]
    dist_from_center = np.sqrt((X - center_insert[1])^2 + (Y-center_insert[0])^2)

    mask = dist_from_center <= radius  
    masked_array = np.zeros_like(array_used)
    for index in range(array_used.shape[2])
        masked_array[:,:,index] = array_used[:,:,index] * mask
	end

    return masked_array, center_insert, mask
end

# ╔═╡ 2a3d8ba6-d4c3-4df9-81c4-c3cc7431b646


# ╔═╡ Cell order:
# ╠═14b6c3e4-48dd-11ec-336c-6918bf024f43
# ╠═9506ae45-f63f-4fd6-8467-0597cbde7006
# ╟─89d2d135-7669-40be-bba1-24a324cdf14e
# ╟─3f313b73-2785-4299-a1e9-bc5d4e1ad422
# ╠═0fbe9f2a-fea0-45fc-b429-85afc5db3341
# ╟─482d8669-c750-450b-aae8-03afc2a90007
# ╠═9e825e49-f9b3-4897-8dc7-755bb4e99b57
# ╠═de7805f3-e670-4fe4-beb9-cfcf6c1bc021
# ╟─1b755e07-926a-42bd-871e-a693bb037c81
# ╠═d2d1d02c-e592-4c73-afa2-80f89dd5534e
# ╠═2d16b7cb-584e-4ff8-84a1-4e5670bd49d3
# ╟─a69293d3-6093-4b87-abc6-11f1ec9398c3
# ╠═dfcda44d-5b2c-4a50-9fb2-44cb70b5488e
# ╠═2a3d8ba6-d4c3-4df9-81c4-c3cc7431b646

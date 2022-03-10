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

# ╔═╡ a9be3ff8-1210-4028-a506-453df7199f1b
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

# ╔═╡ c398de2c-a02d-11ec-2dcf-39af2a0a9ba9


# ╔═╡ 529781bc-5884-4b8e-98b3-a5f480d73b25
TableOfContents()

# ╔═╡ 4b91294f-da07-4751-92c9-02ef714a3b16
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDER`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ 3733919f-339d-4733-9a31-43cc0961f2e0
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ 7743b460-caca-4272-96e9-4150d24c6653
md"""
## Helper Functions
"""

# ╔═╡ 4a3d40d0-e091-49b8-a71c-ac724655449a
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ dd619294-f1a4-4ae3-b21e-35df163fe29f
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ d83ee887-c1a7-4b7f-b854-097e3835af69
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

# ╔═╡ c9a1e22f-b166-4a1e-b84f-856992a36884
md"""
## Segment Heart
"""

# ╔═╡ 4761a8e4-352c-479b-b008-cefa2965ef7b
md"""
## Segment Calcium Rod
"""

# ╔═╡ 42da21d0-cb8c-4af8-902a-877a61568fa2
md"""
## Segment Calcium Inserts
"""

# ╔═╡ e5a00f35-335c-4ad4-be20-1b20f56896c0
md"""
## Calibration Prep
"""

# ╔═╡ 6e4374d3-706d-43b4-9125-e52bb026cb21
md"""
### Calibration Insert
"""

# ╔═╡ 40f47597-407e-4c10-b3be-c9ff724dbc27
md"""
### Calibration Line
"""

# ╔═╡ e8868d3f-aa56-44fb-ab68-bc3c994548f6
density_array = [0, 200, 400, 800]

# ╔═╡ 1d0e87e1-7b99-4212-868d-c9a41c90fe83
density_array_calc = [0, 200] # mg/cc

# ╔═╡ 273e6416-50e8-42c6-b3e9-beb9824b9b60
md"""
We can see from above that the linear regression returns a best fit line with the formula:

```math
y = mx + b
```

Which can be solved for ``x`` and then used to calculate the density (``x``) given some measured intensity (``y``)

```math
x = \frac{y - b}{m}
```

"""

# ╔═╡ c60e47c6-b1bb-408f-8bec-f00d12fbe422
md"""
# Score Large Inserts
"""

# ╔═╡ c4c5a245-4da9-4ffc-942e-86237bc1b5ec
md"""
## High Density
"""

# ╔═╡ 35dad19b-9469-4f0f-9c3d-14dffb7b8c1b
md"""
#### Dilated mask
"""

# ╔═╡ 78eea972-f777-4243-b81e-ea3e219450c0
md"""
#### Ring (background) mask
"""

# ╔═╡ 4510d8d7-1dae-4946-805d-bbe99dff3011
md"""
### Calculations
"""

# ╔═╡ 9943af5f-523a-4292-a5b2-d21eb3307029
md"""
## Medium Density
"""

# ╔═╡ 4a390212-d39d-4de1-b785-b28ac15b9fab
md"""
#### Dilated mask
"""

# ╔═╡ 3cdf12e8-82da-4f41-9a0f-2dc5f0fa1a05
md"""
#### Ring (background) mask
"""

# ╔═╡ dd419a0f-4596-45b3-b71a-e9ccd01b67eb
md"""
### Calculations
"""

# ╔═╡ d513b4f4-c4cb-4147-ab47-14bacb68c9c0
md"""
## Low Density
"""

# ╔═╡ 5df80a8d-da3a-49c5-80bd-a1c6adf599da
md"""
#### Dilated mask
"""

# ╔═╡ d9d4bf44-fbf5-412b-b728-7eeb95c59e56
md"""
#### Ring (background) mask
"""

# ╔═╡ 66919448-f722-46ae-86cc-d8dca63da388
md"""
### Calculations
"""

# ╔═╡ 60f6e3dd-199c-4e78-a8aa-c2ce46f0661f
md"""
# Score Medium Inserts
"""

# ╔═╡ 36f99e43-f511-4656-97c5-7420db121198
md"""
## High Density
"""

# ╔═╡ 8cbfedbf-be74-4be7-9f20-dbd1ef09cf64
md"""
#### Dilated mask
"""

# ╔═╡ 8538186f-0d1b-44d6-aa1b-3f0654aa1521
md"""
#### Ring (background) mask
"""

# ╔═╡ d13b4c3a-2abd-4f3e-8e7b-34a826b2ff30
md"""
### Calculations
"""

# ╔═╡ 5d0280c0-3b0b-4179-9db6-101c98fdfb40
md"""
## Medium Density
"""

# ╔═╡ 9251c1cb-25ab-471e-9e7a-aad542dcd704
md"""
#### Dilated mask
"""

# ╔═╡ e0068fa8-1b98-422d-ae94-2c43d2c0e30d
md"""
#### Ring (background) mask
"""

# ╔═╡ 00ae29a5-c995-4ac6-802a-92566992e795
md"""
### Calculations
"""

# ╔═╡ e22d4160-b807-4a9b-b408-ad4f5d7fe8f9
md"""
## Low Density
"""

# ╔═╡ 483b0382-f170-40b2-b99f-02e66dc78a5a
md"""
#### Dilated mask
"""

# ╔═╡ 4cc98873-9f0d-4c87-945c-b55b28f37d8e
md"""
#### Ring (background) mask
"""

# ╔═╡ ead96638-696d-4fce-a2e3-6dff1e5a01dd
md"""
### Calculations
"""

# ╔═╡ 1a96ce44-7f24-4b4e-b9ab-b033ee27f65b
md"""
# Score Small Inserts
"""

# ╔═╡ 7cdc209f-8f25-4e73-bf3a-0aab35d89a6e
md"""
## High Density
"""

# ╔═╡ 49440be4-83ee-4e3a-adfc-4128b7b1b6ea
md"""
#### Dilated mask
"""

# ╔═╡ 4dac27d1-2a45-4beb-97ad-531f7260742e
md"""
#### Ring (background) mask
"""

# ╔═╡ fee24495-a54c-40c6-9bbc-2f167a2797e7
md"""
### Calculations
"""

# ╔═╡ 808e1b47-8893-4190-9383-1299e9c118a2
md"""
## Medium Density
"""

# ╔═╡ 1b773cd7-e973-4cf9-b11b-909803274ae6
md"""
#### Dilated mask
"""

# ╔═╡ 3cd04a12-c37a-4f2e-b14c-37d61a0d6ddb
md"""
#### Ring (background) mask
"""

# ╔═╡ beda070b-4a56-4c8c-a06b-13dbd3b4a77d
md"""
### Calculations
"""

# ╔═╡ e4622e9a-6159-4cfa-b1e1-23c7717af64f
md"""
## Low Density
"""

# ╔═╡ dafbf407-b3c1-40a6-a7b7-32ae4e98194b
md"""
#### Dilated mask
"""

# ╔═╡ 0610fac4-fd87-4904-bc15-40b0e818016f
md"""
#### Ring (background) mask
"""

# ╔═╡ 6aee61cf-fe45-4a94-9f9e-5cc89843d3e4
md"""
### Calculations
"""

# ╔═╡ 978a5803-90c4-4618-8891-d2c8840068e0
md"""
# Results
"""

# ╔═╡ 1bea042a-8a76-4417-8619-d3f8b555afb5
md"""
### Volume
"""

# ╔═╡ 43ea44f6-d28f-48f8-8d16-0607c5e56e5a
inserts = [
	"Low Density",
	"Medium Density",
	"High Density"
]

# ╔═╡ 99fc6305-a48e-4d87-a19b-1cf8c4aacfba
ground_truth_volume_large = [
	98.2,
	98.2,
	98.2,
] # mm^3

# ╔═╡ 666e0fcc-8b63-4be1-bd3b-398bf19e3522
ground_truth_volume_medium = [
	21.2,
	21.2,
	21.2
]

# ╔═╡ 7302f016-2233-4efa-a0fc-32ac78e501b3
ground_truth_volume_small = [
	0.8,
	0.8,
	0.8
]

# ╔═╡ cb9f592e-85ca-4bd7-b504-6bbdc58a504d
md"""
### Mass
"""

# ╔═╡ 483c152c-aaf3-4c4a-a526-87045288a9b1
ground_truth_mass_large = [
	19.6,
	39.3,
	78.5
] # mg

# ╔═╡ 39907b2e-12c3-4e72-adcc-abaaa91093dc
ground_truth_mass_medium = [
	4.2,
	8.5,
	17.0
]

# ╔═╡ 828a2446-7604-4bf5-b419-a14ab48d7408
ground_truth_mass_small = [
	0.2,
	0.3,
	0.6
]

# ╔═╡ b0fda2d7-dc70-41a2-822a-ad001c2d2897
begin
	SCAN_NUMBER = 7
	VENDER = "Canon_Aquilion_One_Vision"
	BASE_PATH = "/Users/daleblack/Google Drive/Datasets/"
end

# ╔═╡ fd1a63ab-c466-4148-a885-de31bd5561fa
root_path = string(BASE_PATH, VENDER)

# ╔═╡ 8989d500-4fbc-42ba-8195-817fa3b7d7c5
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ 75cc76b2-d60d-4218-8c77-44955c897967
pth = dcm_path_list[SCAN_NUMBER]

# ╔═╡ 20ee21eb-014b-4698-a8e9-24c0530bf763
pth

# ╔═╡ 8953cb87-a45c-455e-9d45-da1e93c538b6
scan = basename(pth)

# ╔═╡ 090c609f-bb07-4fd2-b8de-9882caddee02
header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ a027756b-e806-4231-aab0-a5afd12fa8d8
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ f35a7017-f224-45fc-b8bc-ced378f71468
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ d441b9f9-bdf6-4d31-a49d-667d70d17614
heatmap(masked_array[:, :, a], colormap=:grays)

# ╔═╡ 77b68cbf-bbfa-430f-8b00-b80de824f281
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig2
end

# ╔═╡ 80eab964-fa01-4534-8db5-378aa1845798
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 23]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig3
end

# ╔═╡ 2ce0afbc-884c-492f-b43b-ec865c12d603
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 15]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig
end

# ╔═╡ fefe96d0-b38d-45ba-9642-aec08cad916b
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header);

# ╔═╡ 7c484188-d3c7-443e-b556-d08dbbd71f5d
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=cal_rod_slice, show_value=true)

# ╔═╡ c2faeff9-ff2a-45b6-ba49-8273f662fa08
heatmap(transpose(calcium_image[:, :, c]), colormap=:grays)

# ╔═╡ d4ee3231-440e-4a9f-b7a4-09cffe6b377d
array_filtered = abs.(mapwindow(median, calcium_image[:, :, cal_rod_slice], (3, 3)));

# ╔═╡ 8ad24746-8774-462c-b8ef-92ada042af15
bool_arr = array_filtered .> 0;

# ╔═╡ 7597b33e-3b5e-4f81-99a0-e0b68ad409f3
bool_arr_erode = (((erode(erode(bool_arr)))));

# ╔═╡ 2e21f0d4-d615-47fe-9ca8-c60ada3b7e9a
heatmap(bool_arr_erode, colormap=:grays)

# ╔═╡ 3d5c0a31-0dde-4ca8-820d-7bb99281252d
heatmap(bool_arr, colormap=:grays)

# ╔═╡ 747966ab-7bdc-4239-b0a3-d41cffb144b6
c_img = calcium_image[:, :, cal_rod_slice-1:cal_rod_slice+1];

# ╔═╡ a21b0a2e-4b40-4cd7-8889-a7f47a857375
begin
	mask_cal_3D = Array{Bool}(undef, size(c_img))
	for z in 1:size(c_img, 3)
		mask_cal_3D[:, :, z] = bool_arr_erode
	end
end;

# ╔═╡ 7aaa7a18-c5fb-443e-a73a-a5de12cdc6be
cal_insert_mean = mean(c_img[mask_cal_3D])

# ╔═╡ 9f54f3b1-5465-4b01-8c25-7acbdd8a30aa
intensity_array = [0, cal_insert_mean] # HU

# ╔═╡ f8e90739-7ea5-4c6c-9c39-9b6107560514
df = DataFrame(:density => density_array_calc, :intensity => intensity_array)

# ╔═╡ 59acf7b1-b824-4362-b231-1a88db152696
linearRegressor = lm(@formula(intensity ~ density), df);

# ╔═╡ 3866a0f1-334e-4891-a44b-d26a163ef651
linearFit = predict(linearRegressor)

# ╔═╡ c383e58a-5edb-4b2d-a261-07f676ba2995
m = linearRegressor.model.pp.beta0[2]

# ╔═╡ b20b1ab0-fd7d-4357-a002-8c8eca617913
b = linearRegressor.model.rr.mu[1]

# ╔═╡ 5097ad44-75da-48ae-b6b3-906588ee2018
density(intensity) = (intensity - b) / m

# ╔═╡ 7ff473f5-d225-4fc0-9512-f40f0ebccacb
intensity(ρ) = m*ρ + b

# ╔═╡ 3c38900a-6701-448d-941f-d636446d573a
S_Obj_HD = intensity(800)

# ╔═╡ 6ebf32db-f14a-4f1b-8f38-33a3950bb5b7
S_Obj_MD = intensity(400)

# ╔═╡ 96266ea0-9f1d-4ce5-9422-8077911b4cb0
S_Obj_LD = intensity(200)

# ╔═╡ fd983a24-79a4-49cd-b480-35b29b601663
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

# ╔═╡ 086272cd-15ab-4da8-a821-a9d02c8f1f25
arr = masked_array[:, :, slice_CCI-3:slice_CCI+3];

# ╔═╡ 1e48d983-d021-446b-aa45-e5c2d8027bd3
single_arr = masked_array[:, :, slice_CCI];

# ╔═╡ 107e05a3-9355-4c1e-b622-033338d0ff13
pixel_size = DICOMUtils.get_pixel_size(header)

# ╔═╡ ff709494-23a9-4e0c-9fb2-1beb34f6aaf3
"""
	calc_centers(dcm_array, output, header, tmp_center, CCI_slice)

Function ...
"""
function calc_centers_new(dcm_array, output, header, tmp_center, CCI_slice)
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

# ╔═╡ 42771927-9007-4d1f-bc05-dce17583c6d5
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
    insert_centers = calc_centers_new(dcm_array, output, header, center_insert, CCI_slice)

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

# ╔═╡ c8d9211d-ad92-41d1-b20f-24fb66a848e2
mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
            dcm_array, masked_array, header, slice_CCI, center_insert
);

# ╔═╡ 2ea8bf1b-4e68-4d85-b44d-dc8678534018
masks = mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD;

# ╔═╡ 87105fec-f539-4844-a981-ef93833d9bad
heatmap(masks, colormap=:grays)

# ╔═╡ 509a3ddf-0bf6-4a44-b8e6-e9858d1c1069
begin
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_HD_3D[:, :, z] = mask_L_HD
	end
end;

# ╔═╡ f8c2e5e7-847a-40b8-bb68-2dc765ba59f3
dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D));

# ╔═╡ 7efdb30e-8793-4f88-bb10-52afc04e0061
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ 2e990a74-28de-4f6e-9eeb-6aa3354cbffb
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ f4e63330-828e-409d-9ee1-4746c2fcf32c
ring_mask_L_HD = dilate(dilate(dilate(dilate(mask_L_HD_3D)))) - dilate(dilate(dilate(mask_L_HD_3D)));

# ╔═╡ 62a3b20a-da16-4e37-9326-8d2f76280ebc
@bind g4 overlay_mask_bind(ring_mask_L_HD)

# ╔═╡ b0329f5f-4017-4bcc-91cb-d497f44b29d3
overlay_mask_plot(arr, ring_mask_L_HD, g4, "ring mask")

# ╔═╡ fc60001c-19cb-40e6-b19a-ef69345d5ad5
begin
	single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
	s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
end

# ╔═╡ 0d757443-4bcd-4e06-a67d-85f7bd97b9bb
begin
	alg_L_HD = Integrated(arr[mask_L_HD_3D])
	vol_l_hd = score(s_bkg_L_HD, S_Obj_HD, pixel_size, alg_L_HD)
end

# ╔═╡ 0e970645-065e-47f4-b152-6fe31be96739
begin
	ρ_HD = 0.8 # mg/mm^3
	mass_l_hd = score(s_bkg_L_HD, S_Obj_HD, pixel_size, ρ_HD, alg_L_HD)
end

# ╔═╡ 6fe9e268-f64d-494a-86db-f92d0e94fd4c
begin
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_MD_3D[:, :, z] = mask_L_MD
	end
end;

# ╔═╡ 8c2460f5-580b-451d-98ca-1ca572386ec9
dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D));

# ╔═╡ 06b869dc-b3e7-48e3-9c38-12c79c7ad545
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ 71cfb8dd-5f81-48e5-a12f-4e8667f7e48d
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ c0177d1d-969e-454a-885d-9c3202ee1460
ring_mask_L_MD = dilate(dilate(dilate(dilate(mask_L_MD_3D)))) - dilate(dilate(dilate(mask_L_MD_3D)));

# ╔═╡ 4e0dd5db-99e3-4a2f-9d63-a0a5672ff5e9
@bind h4 overlay_mask_bind(ring_mask_L_MD)

# ╔═╡ 4e0cd0b4-2325-4157-9c7b-2ac02079cb71
overlay_mask_plot(arr, ring_mask_L_MD, h4, "ring mask")

# ╔═╡ 87bcfbb3-5f3a-4189-af44-520caaf9ee3f
begin
	single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
	s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
end

# ╔═╡ 88481f50-29a8-4ae4-a975-97d0f6193520
begin
	alg_L_MD = Integrated(arr[mask_L_MD_3D])
	vol_l_md = score(s_bkg_L_MD, S_Obj_MD, pixel_size, alg_L_MD)
end

# ╔═╡ d357179c-47c2-4147-b091-1f29f672beb0
begin
	ρ_MD = 0.4 # mg/mm^3
	mass_l_md = score(s_bkg_L_MD, S_Obj_MD, pixel_size, ρ_MD, alg_L_MD)
end

# ╔═╡ 20cf4db8-266e-46f5-b206-14b421bc7959
begin
	mask_L_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_LD_3D[:, :, z] = mask_L_LD
	end
end;

# ╔═╡ 929716ac-b027-4f91-baea-728b922ab457
dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D));

# ╔═╡ d9793581-fc5e-4689-ba5e-c7cae99db598
@bind i2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ 98bf3104-139f-43c2-b480-c87d1d3695ee
overlay_mask_plot(arr, dilated_mask_L_LD, i2, "dilated mask")

# ╔═╡ 5fb5183f-d7d7-4177-a17a-514e4b9948ff
ring_mask_L_LD = dilate(dilate(dilate(dilate(mask_L_LD_3D)))) - dilate(dilate(dilate(mask_L_LD_3D)));

# ╔═╡ 5f571808-d497-4ed6-bde9-e49aea6a63ea
@bind i4 overlay_mask_bind(ring_mask_L_LD)

# ╔═╡ 3fab4d83-6b32-42d1-b0e0-b882b0e85eae
overlay_mask_plot(arr, ring_mask_L_LD, i4, "ring mask")

# ╔═╡ 78d005d1-3773-4e1a-833c-a4581e14863a
begin	
	single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
	s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
end

# ╔═╡ c1c70959-7a7c-491e-b73d-0ac772739844
begin
	alg_L_LD = Integrated(arr[mask_L_LD_3D])
	vol_l_ld = score(s_bkg_L_LD, S_Obj_LD, pixel_size, alg_L_LD)
end

# ╔═╡ c28e7cec-bba2-4fa5-bf35-b8f20df9a60a
begin
	ρ_LD = 0.2
	mass_l_ld = score(s_bkg_L_LD, S_Obj_LD, pixel_size, ρ_LD, alg_L_LD)
end

# ╔═╡ 135c3654-9439-40e5-9a17-b0d98b0e8ae9
calculated_mass_large = [
	mass_l_ld,
	mass_l_md,
	mass_l_hd
]

# ╔═╡ 6f61d797-fa44-4a69-8e76-1e599fb34c8d
calculated_volume_large = [
	vol_l_ld,
	vol_l_md,
	vol_l_hd
]

# ╔═╡ 0c3ff619-0912-4574-98fe-adf24a29f346
begin
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_HD_3D[:, :, z] = mask_M_HD
	end
end;

# ╔═╡ 7ea83662-abef-468c-a3ab-5bb003353528
dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ aac94707-f1cc-43a5-bc78-84a7dd22fcc9
@bind j2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ 4902a46d-af52-4e40-a8c3-4d734664e776
overlay_mask_plot(arr, dilated_mask_M_HD, j2, "dilated mask")

# ╔═╡ ba631579-8dc6-462b-a393-78b903964b29
ring_mask_M_HD = dilate(dilate(dilate(dilate(dilate(mask_M_HD_3D))))) - dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 1acd8997-a448-4677-87c0-8c9c11d224a3
@bind j4 overlay_mask_bind(ring_mask_M_HD)

# ╔═╡ ef49a28c-522b-43e9-80ee-97c876c3e008
overlay_mask_plot(arr, ring_mask_M_HD, j4, "ring mask")

# ╔═╡ 1df53e68-7f6a-4b52-bf89-5166b453bfca
begin
	single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
	s_bkg_M_HD = mean(single_arr[single_ring_mask_M_HD])
end

# ╔═╡ c01e76c5-4e03-4697-b9ef-dd068f58886e
begin
	alg_M_HD = Integrated(arr[mask_M_HD_3D])
	vol_m_hd = score(s_bkg_M_HD, S_Obj_HD, pixel_size, alg_M_HD)
end

# ╔═╡ 7a3e945c-49d3-41d4-bb8e-7642695454f5
mass_m_hd = score(s_bkg_M_HD, S_Obj_HD, pixel_size, ρ_HD, alg_M_HD)

# ╔═╡ 438c25ab-5790-45d2-aa04-64f5fbb6a695
begin
	mask_M_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_MD_3D[:, :, z] = mask_M_MD
	end
end;

# ╔═╡ 6dece886-319e-4baa-93d8-1102a9c7f5d8
dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))));

# ╔═╡ 2fcbe86d-6f79-4909-b846-4d5b4587f03e
@bind k2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ 28191883-fb69-4bf0-a0cd-41888879efa7
overlay_mask_plot(arr, dilated_mask_M_MD, k2, "dilated mask")

# ╔═╡ 4f1683df-4821-4e7c-bc47-61e3f608b73b
ring_mask_M_MD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))));

# ╔═╡ 135d83f9-1e11-4c80-92d0-dbd65dab1dc6
@bind k4 overlay_mask_bind(ring_mask_M_MD)

# ╔═╡ 38bd1cd4-a46a-4d06-b95e-97b24b816bbe
overlay_mask_plot(arr, ring_mask_M_MD, k4, "ring mask")

# ╔═╡ dc0ec820-e77e-4310-92b8-6cc41c718139
begin
	single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
	s_bkg_M_MD = mean(single_arr[single_ring_mask_M_MD])
end

# ╔═╡ dda45dfd-f211-492e-9fab-a8dca3f6dc01
begin
	alg_M_MD = Integrated(arr[mask_M_MD_3D])
	vol_m_md = score(s_bkg_M_MD, S_Obj_MD, pixel_size, alg_M_MD)
end

# ╔═╡ 957e4458-a7da-43a5-85e8-f48891f626c8
mass_m_md = score(s_bkg_M_MD, S_Obj_MD, pixel_size, ρ_MD, alg_M_MD)

# ╔═╡ e53644eb-8554-44bd-9d32-fb1417919935
begin
	mask_M_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_LD_3D[:, :, z] = mask_M_LD
	end
end;

# ╔═╡ 7dd79601-9654-462e-869d-6a143557d7f2
dilated_mask_M_LD = dilate(dilate(dilate(dilate(mask_M_LD_3D))));

# ╔═╡ a5edeefb-3bdb-43ea-a44e-f0ada26f0e88
@bind l2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ 96daeb15-58c1-4197-9476-adace4395b37
overlay_mask_plot(arr, dilated_mask_M_LD, l2, "dilated mask")

# ╔═╡ 1697c395-99c9-487f-9a65-4c89c2ec60a8
ring_mask_M_LD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ f4affa5d-345a-4f15-9cc9-6e1fd10d5f3a
@bind l4 overlay_mask_bind(ring_mask_M_LD)

# ╔═╡ 74e88b49-a435-4f88-84cd-836738ca1d17
overlay_mask_plot(arr, ring_mask_M_LD, l4, "ring mask")

# ╔═╡ 6881d31a-edd7-46da-af13-4f40937cd8c4
begin
	single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
	s_bkg_M_LD = mean(single_arr[single_ring_mask_M_LD])
end

# ╔═╡ af7775e1-4d45-413f-92a7-9e0071ba70cf
begin
	alg_M_LD = Integrated(arr[mask_M_LD_3D])
	vol_m_ld = score(s_bkg_M_LD, S_Obj_LD, pixel_size, alg_M_LD)
end

# ╔═╡ 5e67c7c1-692c-4990-b27c-9219ccb3c237
mass_m_ld = score(s_bkg_M_LD, S_Obj_LD, pixel_size, ρ_LD, alg_M_LD)

# ╔═╡ 120cdf3f-1055-4ff2-84da-86b11346c1f4
calculated_mass_medium = [
	mass_m_ld,
	mass_m_md,
	mass_m_hd
]

# ╔═╡ 5e4b92bf-acd4-4fb1-ac2a-e846d1662eab
calculated_volume_medium = [
	vol_m_ld,
	vol_m_md,
	vol_m_hd
]

# ╔═╡ cb674784-a009-4aec-8dc4-c4c64db7437f
begin
	mask_S_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_HD_3D[:, :, z] = mask_S_HD
	end
end;

# ╔═╡ 7107b8de-3b16-4fdf-9d0d-3513ed6abd6b
dilated_mask_S_HD = dilate(dilate(dilate(dilate(mask_S_HD_3D))));

# ╔═╡ cc304c91-b1df-4286-a3a7-a70f06d0500a
@bind m2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ ce3e5bd4-3fa7-452a-896c-f8a105bc9ac9
overlay_mask_plot(arr, dilated_mask_S_HD, m2, "dilated mask")

# ╔═╡ 0457e5a0-c966-4799-8df5-8f1bdc5be0ff
ring_mask_S_HD = dilate(dilate(dilate(dilate(dilate(mask_S_HD_3D))))) - dilate(dilate(dilate(dilate(mask_S_HD_3D))));

# ╔═╡ b0ed79d5-fc9c-4c64-8d1b-cb80930ce138
@bind m4 overlay_mask_bind(ring_mask_S_HD)

# ╔═╡ 1ccefa42-9cda-4f34-831a-80f902837d91
overlay_mask_plot(arr, ring_mask_S_HD, m4, "ring mask")

# ╔═╡ f2048b7f-b8bd-419f-977d-29f73c6ec9bb
begin
	single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
	s_bkg_S_HD = mean(single_arr[single_ring_mask_S_HD])
end

# ╔═╡ 7645bd28-ba16-43d5-9d08-95e7e0ebb9d4
begin
	alg_S_HD = Integrated(arr[mask_S_HD_3D])
	vol_s_hd = score(s_bkg_S_HD, S_Obj_HD, pixel_size, alg_S_HD)
	if vol_s_hd < 0
		vol_s_hd = 0
	end
end

# ╔═╡ dfe3274d-d2f4-46b3-8a51-6f0672b381e3
begin
	mass_s_hd = score(s_bkg_S_HD, S_Obj_HD, pixel_size, ρ_HD, alg_S_HD)
	if mass_s_hd < 0
		mass_s_hd = 0
	end
end

# ╔═╡ 8a8a551f-8f99-47d9-b4ec-2fd8c0fc11eb
begin
	mask_S_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_MD_3D[:, :, z] = mask_S_MD
	end
end;

# ╔═╡ 807d7093-432e-46cd-af9c-57632f7f8fc5
dilated_mask_S_MD = dilate(dilate(dilate(dilate(mask_S_MD_3D))));

# ╔═╡ dd2b13bd-af75-4fbc-9b00-9335f1d432e8
@bind n2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ d92c53b7-f822-42ec-a8bb-6f22db636305
overlay_mask_plot(arr, dilated_mask_S_MD, n2, "dilated mask")

# ╔═╡ 81753c81-dd01-4774-89f3-899265c62e99
ring_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D))))) - dilate(dilate(dilate(dilate(mask_S_MD_3D))));

# ╔═╡ d2a29609-8d60-4c72-92e8-80da5fbfbdb6
@bind n4 overlay_mask_bind(ring_mask_S_MD)

# ╔═╡ 4d470abc-66fb-444a-af23-3eedf3984e02
overlay_mask_plot(arr, ring_mask_S_MD, n4, "ring mask")

# ╔═╡ 9bfe9c23-7c61-4edf-84a5-742a9c3e5fa3
begin
	single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
	s_bkg_S_MD = mean(single_arr[single_ring_mask_S_MD])
end

# ╔═╡ 4ddea5a8-0273-410d-9d49-e514585b4bbc
begin
	alg_S_MD = Integrated(arr[mask_S_MD_3D])
	vol_s_md = score(s_bkg_S_MD, S_Obj_HD, pixel_size, alg_S_MD)
	if vol_s_md < 0
		vol_s_md = 0
	end
end

# ╔═╡ b11a557d-b3a1-4885-8286-bf942061f05f
begin
	mass_s_md = score(s_bkg_S_MD, S_Obj_HD, pixel_size, ρ_MD, alg_S_MD)
	if mass_s_md < 0
		mass_s_md = 0
	end
end

# ╔═╡ 70386c43-e5bb-4a4b-b433-0c2b0078536c
begin
	mask_S_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_LD_3D[:, :, z] = mask_S_LD
	end
end;

# ╔═╡ 3ea5fef8-0e41-41eb-a790-a0fadae7fa26
dilated_mask_S_LD = dilate(dilate(dilate(dilate(mask_S_LD_3D))));

# ╔═╡ 969987c0-bc72-4c12-a196-c57e7fa54b14
@bind o2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ ad5f61e9-4046-46ab-8213-3fac6fddd136
overlay_mask_plot(arr, dilated_mask_S_LD, o2, "dilated mask")

# ╔═╡ a50d94d8-2272-470f-a5ab-30a9cfe2a448
ring_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D))))) - dilate(dilate(dilate(dilate(mask_S_LD_3D))));

# ╔═╡ 1424585e-e971-4d39-bacb-3dfe3b5dfc02
@bind o4 overlay_mask_bind(ring_mask_S_LD)

# ╔═╡ a14fbc73-a3dd-4fd9-86a9-0200ffa0ed77
overlay_mask_plot(arr, ring_mask_S_LD, o4, "ring mask")

# ╔═╡ f57e77e8-7f79-4c95-813d-fcfa60ee2207
begin
	single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
	s_bkg_S_LD = mean(single_arr[single_ring_mask_S_LD])
end

# ╔═╡ 5d57099b-4098-47f7-b687-7ab753d60cbc
begin
	alg_S_LD = Integrated(arr[mask_S_LD_3D])
	vol_s_ld = score(s_bkg_S_LD, S_Obj_HD, pixel_size, alg_S_LD)
	if vol_s_ld < 0
		vol_s_ld = 0
	end
end

# ╔═╡ 6ae87f69-61b9-415d-8aab-69d9ff4d9726
begin
	mass_s_ld = score(s_bkg_S_LD, S_Obj_HD, pixel_size, ρ_LD, alg_S_LD)
	if mass_s_ld < 0
		mass_s_ld = 0
	end
end

# ╔═╡ 5a1aa2d2-20f7-4d29-a134-8cd295fc28c6
calculated_mass_small = [
	mass_s_ld,
	mass_s_md,
	mass_s_hd
]

# ╔═╡ e4a2e60e-5487-4f0c-a3ea-ca2aa238795f
calculated_volume_small = [
	vol_s_ld,
	vol_s_md,
	vol_s_hd
]

# ╔═╡ 584f426f-1b91-468f-94d1-25a3e16ec563
df1 = DataFrame(
	inserts = inserts,
	ground_truth_volume_large = ground_truth_volume_large,
	calculated_volume_large = calculated_volume_large,
	ground_truth_volume_medium = ground_truth_volume_medium,
	calculated_volume_medium = calculated_volume_medium,
	ground_truth_volume_small = ground_truth_volume_small,
	calculated_volume_small = calculated_volume_small
)

# ╔═╡ a078a864-d9e1-426e-b6fc-3133af9f4298
begin
	fvol2 = Figure()
	axvol2 = Axis(fvol2[1, 1])
	
	scatter!(density_array[2:end], df1[!, :ground_truth_volume_large], label="ground_truth_volume_large")
	scatter!(density_array[2:end], df1[!, :calculated_volume_large], label="calculated_volume_large")
	
	axvol2.title = "Volume Measurements (Large)"
	axvol2.ylabel = "Volume (mm^3)"
	axvol2.xlabel = "Density (mg/cm^3)"

	xlims!(axvol2, 0, 850)
	ylims!(axvol2, 0, 200)
	
	fvol2[1, 2] = Legend(fvol2, axvol2, framevisible = false)
	
	fvol2
end

# ╔═╡ e866663b-0a5d-49e9-a808-7f7da4b2123d
begin
	fvol3 = Figure()
	axvol3 = Axis(fvol3[1, 1])
	
	scatter!(density_array[2:end], df1[!, :ground_truth_volume_medium], label="ground_truth_volume_medium")
	scatter!(density_array[2:end], df1[!, :calculated_volume_medium], label="calculated_volume_medium")
	
	axvol3.title = "Volume Measurements (Medium)"
	axvol3.ylabel = "Volume (mm^3)"
	axvol3.xlabel = "Density (mg/cm^3)"

	xlims!(axvol3, 0, 850)
	ylims!(axvol3, 0, 50)
	
	fvol3[1, 2] = Legend(fvol3, axvol3, framevisible = false)
	
	fvol3
end

# ╔═╡ 5249cc07-96a1-4558-906e-ccc1a326cced
begin
	fvol4 = Figure()
	axvol4 = Axis(fvol4[1, 1])
	
	scatter!(density_array[2:end], df1[!, :ground_truth_volume_small], label="ground_truth_volume_small")
	scatter!(density_array[2:end], df1[!, :calculated_volume_small], label="calculated_volume_small")
	
	axvol4.title = "Volume Measurements (Small)"
	axvol4.ylabel = "Volume (mm^3)"
	axvol4.xlabel = "Density (mg/cm^3)"

	xlims!(axvol4, 0, 850)
	ylims!(axvol4, 0, 2.5)
	
	fvol4[1, 2] = Legend(fvol4, axvol4, framevisible = false)
	
	fvol4
end

# ╔═╡ c3410597-a5a7-4903-8c69-b2e3a00b06a0
df2 = DataFrame(
	inserts = inserts,
	ground_truth_mass_large = ground_truth_mass_large,
	calculated_mass_large = calculated_mass_large,
	ground_truth_mass_medium = ground_truth_mass_medium,
	calculated_mass_medium = calculated_mass_medium,
	ground_truth_mass_small = ground_truth_mass_small,
	calculated_mass_small = calculated_mass_small
)

# ╔═╡ e6e6fd48-3a9f-4592-ba42-c1fda46c76ee
begin
	fmass2 = Figure()
	axmass2 = Axis(fmass2[1, 1])
	
	scatter!(density_array[2:end], df2[!, :ground_truth_mass_large], label="ground_truth_mass_large")
	scatter!(density_array[2:end], df2[!, :calculated_mass_large], label="calculated_mass_large")
	
	axmass2.title = "Mass Measurements (Large)"
	axmass2.ylabel = "Mass (mg)"
	axmass2.xlabel = "Density (mg/cm^3)"

	xlims!(axmass2, 0, 850)
	ylims!(axmass2, 0, 100)
	
	fmass2[1, 2] = Legend(fmass2, axmass2, framevisible = false)
	
	fmass2
end

# ╔═╡ a19061de-59e1-435b-a9d6-d2ae394f2d49
begin
	fmass3 = Figure()
	axmass3 = Axis(fmass3[1, 1])
	
	scatter!(density_array[2:end], df2[!, :ground_truth_mass_medium], label="ground_truth_mass_medium")
	scatter!(density_array[2:end], df2[!, :calculated_mass_medium], label="calculated_mass_medium")
	
	axmass3.title = "Mass Measurements (Medium)"
	axmass3.ylabel = "Mass (mg)"
	axmass3.xlabel = "Density (mg/cm^3)"

	xlims!(axmass3, 0, 850)
	ylims!(axmass3, 0, 25)
	
	fmass3[1, 2] = Legend(fmass3, axmass3, framevisible = false)
	
	fmass3
end

# ╔═╡ c0699caa-f893-46e8-9e00-deea883ca4da
begin
	fmass4 = Figure()
	axmass4 = Axis(fmass4[1, 1])
	
	scatter!(density_array[2:end], df2[!, :ground_truth_mass_small], label="ground_truth_mass_small")
	scatter!(density_array[2:end], df2[!, :calculated_mass_small], label="calculated_mass_small")
	
	axmass4.title = "Mass Measurements (Small)"
	axmass4.ylabel = "Mass (mg)"
	axmass4.xlabel = "Density (mg/cm^3)"

	xlims!(axmass4, 0, 850)
	ylims!(axmass4, 0, 1.5)
	
	fmass4[1, 2] = Legend(fmass4, axmass4, framevisible = false)
	
	fmass4
end

# ╔═╡ 62650e9c-5d94-4953-a9be-2b2cbd57a76f
md"""
### Save Results
"""

# ╔═╡ 50d0b8e3-0448-45d5-bb40-ad7146c0b22e
df_final = leftjoin(df1, df2, on=:inserts)

# ╔═╡ 15e136d4-8f9d-4c3c-ad96-bd1b771e0038
if ~isdir(string(cd(pwd, "..") , "/data/output/", VENDER, "3"))
	mkdir(string(cd(pwd, "..") , "/data/output/", VENDER, "3"))
end

# ╔═╡ 6f0f8099-0530-4dde-a89a-7995b2c11eaf
output_path = string(cd(pwd, "..") , "/data/output/", VENDER, "3", "/", scan, ".csv")

# ╔═╡ f925eb9d-f4bd-46f8-af75-3d95def7b44c
CSV.write(output_path, df_final)

# ╔═╡ Cell order:
# ╠═c398de2c-a02d-11ec-2dcf-39af2a0a9ba9
# ╠═a9be3ff8-1210-4028-a506-453df7199f1b
# ╠═529781bc-5884-4b8e-98b3-a5f480d73b25
# ╠═4b91294f-da07-4751-92c9-02ef714a3b16
# ╠═3733919f-339d-4733-9a31-43cc0961f2e0
# ╠═fd1a63ab-c466-4148-a885-de31bd5561fa
# ╠═8989d500-4fbc-42ba-8195-817fa3b7d7c5
# ╠═75cc76b2-d60d-4218-8c77-44955c897967
# ╠═20ee21eb-014b-4698-a8e9-24c0530bf763
# ╠═8953cb87-a45c-455e-9d45-da1e93c538b6
# ╠═090c609f-bb07-4fd2-b8de-9882caddee02
# ╠═7743b460-caca-4272-96e9-4150d24c6653
# ╠═4a3d40d0-e091-49b8-a71c-ac724655449a
# ╠═dd619294-f1a4-4ae3-b21e-35df163fe29f
# ╠═d83ee887-c1a7-4b7f-b854-097e3835af69
# ╠═c9a1e22f-b166-4a1e-b84f-856992a36884
# ╠═a027756b-e806-4231-aab0-a5afd12fa8d8
# ╠═f35a7017-f224-45fc-b8bc-ced378f71468
# ╠═d441b9f9-bdf6-4d31-a49d-667d70d17614
# ╠═2ce0afbc-884c-492f-b43b-ec865c12d603
# ╠═77b68cbf-bbfa-430f-8b00-b80de824f281
# ╠═80eab964-fa01-4534-8db5-378aa1845798
# ╠═4761a8e4-352c-479b-b008-cefa2965ef7b
# ╠═fefe96d0-b38d-45ba-9642-aec08cad916b
# ╠═7c484188-d3c7-443e-b556-d08dbbd71f5d
# ╠═c2faeff9-ff2a-45b6-ba49-8273f662fa08
# ╠═42da21d0-cb8c-4af8-902a-877a61568fa2
# ╠═42771927-9007-4d1f-bc05-dce17583c6d5
# ╠═c8d9211d-ad92-41d1-b20f-24fb66a848e2
# ╠═2ea8bf1b-4e68-4d85-b44d-dc8678534018
# ╠═87105fec-f539-4844-a981-ef93833d9bad
# ╠═e5a00f35-335c-4ad4-be20-1b20f56896c0
# ╠═6e4374d3-706d-43b4-9125-e52bb026cb21
# ╠═d4ee3231-440e-4a9f-b7a4-09cffe6b377d
# ╠═8ad24746-8774-462c-b8ef-92ada042af15
# ╠═7597b33e-3b5e-4f81-99a0-e0b68ad409f3
# ╠═3d5c0a31-0dde-4ca8-820d-7bb99281252d
# ╠═2e21f0d4-d615-47fe-9ca8-c60ada3b7e9a
# ╠═a21b0a2e-4b40-4cd7-8889-a7f47a857375
# ╠═747966ab-7bdc-4239-b0a3-d41cffb144b6
# ╠═7aaa7a18-c5fb-443e-a73a-a5de12cdc6be
# ╠═40f47597-407e-4c10-b3be-c9ff724dbc27
# ╠═e8868d3f-aa56-44fb-ab68-bc3c994548f6
# ╠═1d0e87e1-7b99-4212-868d-c9a41c90fe83
# ╠═9f54f3b1-5465-4b01-8c25-7acbdd8a30aa
# ╠═f8e90739-7ea5-4c6c-9c39-9b6107560514
# ╠═59acf7b1-b824-4362-b231-1a88db152696
# ╠═3866a0f1-334e-4891-a44b-d26a163ef651
# ╠═c383e58a-5edb-4b2d-a261-07f676ba2995
# ╠═b20b1ab0-fd7d-4357-a002-8c8eca617913
# ╠═273e6416-50e8-42c6-b3e9-beb9824b9b60
# ╠═5097ad44-75da-48ae-b6b3-906588ee2018
# ╠═7ff473f5-d225-4fc0-9512-f40f0ebccacb
# ╠═fd983a24-79a4-49cd-b480-35b29b601663
# ╠═c60e47c6-b1bb-408f-8bec-f00d12fbe422
# ╠═086272cd-15ab-4da8-a821-a9d02c8f1f25
# ╠═1e48d983-d021-446b-aa45-e5c2d8027bd3
# ╠═c4c5a245-4da9-4ffc-942e-86237bc1b5ec
# ╠═509a3ddf-0bf6-4a44-b8e6-e9858d1c1069
# ╠═35dad19b-9469-4f0f-9c3d-14dffb7b8c1b
# ╠═f8c2e5e7-847a-40b8-bb68-2dc765ba59f3
# ╠═7efdb30e-8793-4f88-bb10-52afc04e0061
# ╠═2e990a74-28de-4f6e-9eeb-6aa3354cbffb
# ╠═78eea972-f777-4243-b81e-ea3e219450c0
# ╠═f4e63330-828e-409d-9ee1-4746c2fcf32c
# ╠═62a3b20a-da16-4e37-9326-8d2f76280ebc
# ╠═b0329f5f-4017-4bcc-91cb-d497f44b29d3
# ╠═4510d8d7-1dae-4946-805d-bbe99dff3011
# ╠═fc60001c-19cb-40e6-b19a-ef69345d5ad5
# ╠═3c38900a-6701-448d-941f-d636446d573a
# ╠═107e05a3-9355-4c1e-b622-033338d0ff13
# ╠═0d757443-4bcd-4e06-a67d-85f7bd97b9bb
# ╠═0e970645-065e-47f4-b152-6fe31be96739
# ╠═9943af5f-523a-4292-a5b2-d21eb3307029
# ╠═6fe9e268-f64d-494a-86db-f92d0e94fd4c
# ╠═4a390212-d39d-4de1-b785-b28ac15b9fab
# ╠═8c2460f5-580b-451d-98ca-1ca572386ec9
# ╠═06b869dc-b3e7-48e3-9c38-12c79c7ad545
# ╠═71cfb8dd-5f81-48e5-a12f-4e8667f7e48d
# ╠═3cdf12e8-82da-4f41-9a0f-2dc5f0fa1a05
# ╠═c0177d1d-969e-454a-885d-9c3202ee1460
# ╠═4e0dd5db-99e3-4a2f-9d63-a0a5672ff5e9
# ╠═4e0cd0b4-2325-4157-9c7b-2ac02079cb71
# ╠═dd419a0f-4596-45b3-b71a-e9ccd01b67eb
# ╠═87bcfbb3-5f3a-4189-af44-520caaf9ee3f
# ╠═6ebf32db-f14a-4f1b-8f38-33a3950bb5b7
# ╠═88481f50-29a8-4ae4-a975-97d0f6193520
# ╠═d357179c-47c2-4147-b091-1f29f672beb0
# ╠═d513b4f4-c4cb-4147-ab47-14bacb68c9c0
# ╠═20cf4db8-266e-46f5-b206-14b421bc7959
# ╠═5df80a8d-da3a-49c5-80bd-a1c6adf599da
# ╠═929716ac-b027-4f91-baea-728b922ab457
# ╠═d9793581-fc5e-4689-ba5e-c7cae99db598
# ╠═98bf3104-139f-43c2-b480-c87d1d3695ee
# ╠═d9d4bf44-fbf5-412b-b728-7eeb95c59e56
# ╠═5fb5183f-d7d7-4177-a17a-514e4b9948ff
# ╠═5f571808-d497-4ed6-bde9-e49aea6a63ea
# ╠═3fab4d83-6b32-42d1-b0e0-b882b0e85eae
# ╠═66919448-f722-46ae-86cc-d8dca63da388
# ╠═78d005d1-3773-4e1a-833c-a4581e14863a
# ╠═96266ea0-9f1d-4ce5-9422-8077911b4cb0
# ╠═c1c70959-7a7c-491e-b73d-0ac772739844
# ╠═c28e7cec-bba2-4fa5-bf35-b8f20df9a60a
# ╠═60f6e3dd-199c-4e78-a8aa-c2ce46f0661f
# ╠═36f99e43-f511-4656-97c5-7420db121198
# ╠═0c3ff619-0912-4574-98fe-adf24a29f346
# ╠═8cbfedbf-be74-4be7-9f20-dbd1ef09cf64
# ╠═7ea83662-abef-468c-a3ab-5bb003353528
# ╠═aac94707-f1cc-43a5-bc78-84a7dd22fcc9
# ╠═4902a46d-af52-4e40-a8c3-4d734664e776
# ╠═8538186f-0d1b-44d6-aa1b-3f0654aa1521
# ╠═ba631579-8dc6-462b-a393-78b903964b29
# ╠═1acd8997-a448-4677-87c0-8c9c11d224a3
# ╠═ef49a28c-522b-43e9-80ee-97c876c3e008
# ╠═d13b4c3a-2abd-4f3e-8e7b-34a826b2ff30
# ╠═1df53e68-7f6a-4b52-bf89-5166b453bfca
# ╠═c01e76c5-4e03-4697-b9ef-dd068f58886e
# ╠═7a3e945c-49d3-41d4-bb8e-7642695454f5
# ╠═5d0280c0-3b0b-4179-9db6-101c98fdfb40
# ╠═438c25ab-5790-45d2-aa04-64f5fbb6a695
# ╠═9251c1cb-25ab-471e-9e7a-aad542dcd704
# ╠═6dece886-319e-4baa-93d8-1102a9c7f5d8
# ╠═2fcbe86d-6f79-4909-b846-4d5b4587f03e
# ╠═28191883-fb69-4bf0-a0cd-41888879efa7
# ╠═e0068fa8-1b98-422d-ae94-2c43d2c0e30d
# ╠═4f1683df-4821-4e7c-bc47-61e3f608b73b
# ╠═135d83f9-1e11-4c80-92d0-dbd65dab1dc6
# ╠═38bd1cd4-a46a-4d06-b95e-97b24b816bbe
# ╠═00ae29a5-c995-4ac6-802a-92566992e795
# ╠═dc0ec820-e77e-4310-92b8-6cc41c718139
# ╠═dda45dfd-f211-492e-9fab-a8dca3f6dc01
# ╠═957e4458-a7da-43a5-85e8-f48891f626c8
# ╠═e22d4160-b807-4a9b-b408-ad4f5d7fe8f9
# ╠═e53644eb-8554-44bd-9d32-fb1417919935
# ╠═483b0382-f170-40b2-b99f-02e66dc78a5a
# ╠═7dd79601-9654-462e-869d-6a143557d7f2
# ╠═a5edeefb-3bdb-43ea-a44e-f0ada26f0e88
# ╠═96daeb15-58c1-4197-9476-adace4395b37
# ╠═4cc98873-9f0d-4c87-945c-b55b28f37d8e
# ╠═1697c395-99c9-487f-9a65-4c89c2ec60a8
# ╠═f4affa5d-345a-4f15-9cc9-6e1fd10d5f3a
# ╠═74e88b49-a435-4f88-84cd-836738ca1d17
# ╠═ead96638-696d-4fce-a2e3-6dff1e5a01dd
# ╠═6881d31a-edd7-46da-af13-4f40937cd8c4
# ╠═af7775e1-4d45-413f-92a7-9e0071ba70cf
# ╠═5e67c7c1-692c-4990-b27c-9219ccb3c237
# ╠═1a96ce44-7f24-4b4e-b9ab-b033ee27f65b
# ╠═7cdc209f-8f25-4e73-bf3a-0aab35d89a6e
# ╠═cb674784-a009-4aec-8dc4-c4c64db7437f
# ╟─49440be4-83ee-4e3a-adfc-4128b7b1b6ea
# ╠═7107b8de-3b16-4fdf-9d0d-3513ed6abd6b
# ╠═cc304c91-b1df-4286-a3a7-a70f06d0500a
# ╠═ce3e5bd4-3fa7-452a-896c-f8a105bc9ac9
# ╠═4dac27d1-2a45-4beb-97ad-531f7260742e
# ╠═0457e5a0-c966-4799-8df5-8f1bdc5be0ff
# ╠═b0ed79d5-fc9c-4c64-8d1b-cb80930ce138
# ╠═1ccefa42-9cda-4f34-831a-80f902837d91
# ╠═fee24495-a54c-40c6-9bbc-2f167a2797e7
# ╠═f2048b7f-b8bd-419f-977d-29f73c6ec9bb
# ╠═7645bd28-ba16-43d5-9d08-95e7e0ebb9d4
# ╠═dfe3274d-d2f4-46b3-8a51-6f0672b381e3
# ╠═808e1b47-8893-4190-9383-1299e9c118a2
# ╠═8a8a551f-8f99-47d9-b4ec-2fd8c0fc11eb
# ╠═1b773cd7-e973-4cf9-b11b-909803274ae6
# ╠═807d7093-432e-46cd-af9c-57632f7f8fc5
# ╠═dd2b13bd-af75-4fbc-9b00-9335f1d432e8
# ╠═d92c53b7-f822-42ec-a8bb-6f22db636305
# ╠═3cd04a12-c37a-4f2e-b14c-37d61a0d6ddb
# ╠═81753c81-dd01-4774-89f3-899265c62e99
# ╠═d2a29609-8d60-4c72-92e8-80da5fbfbdb6
# ╠═4d470abc-66fb-444a-af23-3eedf3984e02
# ╠═beda070b-4a56-4c8c-a06b-13dbd3b4a77d
# ╠═9bfe9c23-7c61-4edf-84a5-742a9c3e5fa3
# ╠═4ddea5a8-0273-410d-9d49-e514585b4bbc
# ╠═b11a557d-b3a1-4885-8286-bf942061f05f
# ╠═e4622e9a-6159-4cfa-b1e1-23c7717af64f
# ╠═70386c43-e5bb-4a4b-b433-0c2b0078536c
# ╠═dafbf407-b3c1-40a6-a7b7-32ae4e98194b
# ╠═3ea5fef8-0e41-41eb-a790-a0fadae7fa26
# ╠═969987c0-bc72-4c12-a196-c57e7fa54b14
# ╠═ad5f61e9-4046-46ab-8213-3fac6fddd136
# ╠═0610fac4-fd87-4904-bc15-40b0e818016f
# ╠═a50d94d8-2272-470f-a5ab-30a9cfe2a448
# ╠═1424585e-e971-4d39-bacb-3dfe3b5dfc02
# ╠═a14fbc73-a3dd-4fd9-86a9-0200ffa0ed77
# ╠═6aee61cf-fe45-4a94-9f9e-5cc89843d3e4
# ╠═f57e77e8-7f79-4c95-813d-fcfa60ee2207
# ╠═5d57099b-4098-47f7-b687-7ab753d60cbc
# ╠═6ae87f69-61b9-415d-8aab-69d9ff4d9726
# ╠═978a5803-90c4-4618-8891-d2c8840068e0
# ╠═1bea042a-8a76-4417-8619-d3f8b555afb5
# ╠═43ea44f6-d28f-48f8-8d16-0607c5e56e5a
# ╠═99fc6305-a48e-4d87-a19b-1cf8c4aacfba
# ╠═6f61d797-fa44-4a69-8e76-1e599fb34c8d
# ╠═666e0fcc-8b63-4be1-bd3b-398bf19e3522
# ╠═5e4b92bf-acd4-4fb1-ac2a-e846d1662eab
# ╠═7302f016-2233-4efa-a0fc-32ac78e501b3
# ╠═e4a2e60e-5487-4f0c-a3ea-ca2aa238795f
# ╠═584f426f-1b91-468f-94d1-25a3e16ec563
# ╠═a078a864-d9e1-426e-b6fc-3133af9f4298
# ╠═e866663b-0a5d-49e9-a808-7f7da4b2123d
# ╠═5249cc07-96a1-4558-906e-ccc1a326cced
# ╠═cb9f592e-85ca-4bd7-b504-6bbdc58a504d
# ╠═483c152c-aaf3-4c4a-a526-87045288a9b1
# ╠═135c3654-9439-40e5-9a17-b0d98b0e8ae9
# ╠═39907b2e-12c3-4e72-adcc-abaaa91093dc
# ╠═120cdf3f-1055-4ff2-84da-86b11346c1f4
# ╠═828a2446-7604-4bf5-b419-a14ab48d7408
# ╠═5a1aa2d2-20f7-4d29-a134-8cd295fc28c6
# ╠═b0fda2d7-dc70-41a2-822a-ad001c2d2897
# ╠═ff709494-23a9-4e0c-9fb2-1beb34f6aaf3
# ╠═c3410597-a5a7-4903-8c69-b2e3a00b06a0
# ╠═e6e6fd48-3a9f-4592-ba42-c1fda46c76ee
# ╠═a19061de-59e1-435b-a9d6-d2ae394f2d49
# ╠═c0699caa-f893-46e8-9e00-deea883ca4da
# ╠═62650e9c-5d94-4953-a9be-2b2cbd57a76f
# ╠═50d0b8e3-0448-45d5-bb40-ad7146c0b22e
# ╠═15e136d4-8f9d-4c3c-ad96-bd1b771e0038
# ╠═6f0f8099-0530-4dde-a89a-7995b2c11eaf
# ╠═f925eb9d-f4bd-46f8-af75-3d95def7b44c

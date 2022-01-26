### A Pluto.jl notebook ###
# v0.17.6

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

# ╔═╡ c50510a7-e981-4933-a4b1-fab2efcc43a2
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
		Pkg.add("CairoMakie")
		Pkg.add("Statistics")
		Pkg.add("ImageMorphology")
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
	using ImageMorphology
	using CSV
	using DataFrames
	using GLM
	using DICOM
	using DICOMUtils
	using Phantoms
	using CalciumScoring
end

# ╔═╡ 97a264e6-1037-4862-b6bd-2d2486bd419e
TableOfContents()

# ╔═╡ 9ecf691f-8523-4ea9-a635-7482a1a7b697
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDER`, once for every set, and the `SCAN_NUMBER`, once for each scan
"""

# ╔═╡ c72aaf65-8ba2-4c81-9bad-a059b983a833
SCAN_NUMBER = 2

# ╔═╡ e33a903d-d9fb-4960-bc81-128cfbff8af6
VENDER = "Canon_Aquilion_One_Vision"

# ╔═╡ 762cde8a-82c0-4cf0-b2b3-8ed00a4abb3d
BASE_PATH = "/Users/daleblack/Google Drive/Datasets/"

# ╔═╡ 416c3d5f-9953-4a0d-be1c-1049e39edd4b
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ 916be3b8-b7e3-447d-82a2-07e5185443e1
root_path = string(BASE_PATH, VENDER)

# ╔═╡ 4669bb68-6224-4a31-ae50-c2e44744b87e
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ df50149d-ad91-4607-8cc1-6a17da8d8562
pth = dcm_path_list[SCAN_NUMBER]

# ╔═╡ 35c5e1e1-c88d-4202-8f7d-0e6c434ea161
scan = basename(pth)

# ╔═╡ 5cc4cb9b-b0d7-4158-816f-0cbd4d7a7463
header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ 1fae526e-7db7-46e8-aa63-19882cb2e172
pixel_size = DICOMUtils.get_pixel_size(header)

# ╔═╡ fea80b94-e9c5-4ff7-9b91-4e35bf03b5e3
md"""
## Segment Heart
"""

# ╔═╡ 49184d1a-1e84-48c7-8f71-d179fd4671c5
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ 21266f0e-b535-4e72-950e-ba7d125684de
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 536506ce-2031-4c32-a4e3-d6fc6d7da27c
heatmap(masked_array[:, :, a], colormap=:grays)

# ╔═╡ e2de04cd-2f66-4192-abae-95693b1b5176
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 15]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig
end

# ╔═╡ be5f9be6-e10c-48c1-b9fe-ed229a1dff26
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig2
end

# ╔═╡ c0704706-87c2-4078-8345-7e2097639492
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 23]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig3
end

# ╔═╡ 6e597bb0-fa70-46f3-8526-dce3cd4e7b78
md"""
## Segment Calcium Rod
"""

# ╔═╡ 07c0f484-6225-4c6f-ab20-1941bd6abc7a
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header);

# ╔═╡ e9af6401-8f3a-4f86-8d23-02d88dc7ee48
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=cal_rod_slice, show_value=true)

# ╔═╡ 4f795d56-a67b-4b9f-be94-c37890a5c4c4
heatmap(calcium_image[:, :, c], colormap=:grays)

# ╔═╡ 69d3e9a8-68c9-4a5a-8a4a-86c5fd64a08a
md"""
## Segment Calcium Inserts
"""

# ╔═╡ 0508d26f-09a1-457b-8cb0-5693f5c9ba3f
mask_L_LD, mask_M_LD, mask_S_LD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_HD, mask_M_HD, mask_S_HD = mask_inserts(
            dcm_array, masked_array, header, slice_CCI, center_insert
);

# ╔═╡ 544ae005-6661-487f-8662-ea3006ec4069
masks = mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD;

# ╔═╡ 2647d239-58d6-4448-9e3f-66918a10dbad
heatmap(masks, colormap=:grays)

# ╔═╡ 2afdc85a-a8bf-4fd4-bdef-74c39c0fb452
md"""
## Overlay Mask Calcium Inserts
"""

# ╔═╡ c2394467-eb15-4610-9bfd-7b3f88a7fb00
arr_L_HD = masked_array[:, :, slice_CCI-2:slice_CCI+2] .* mask_L_HD;

# ╔═╡ 9b6e3a13-9044-4b7e-89e0-b552af52bd2b
@bind d PlutoUI.Slider(1:size(arr_L_HD, 3), default=2, show_value=true)

# ╔═╡ 0f85bfc3-f0d8-4ae7-b5b4-f72a072a3d55
heatmap(arr_L_HD[:, :, d], colormap=:grays)

# ╔═╡ b8892b4e-3be7-4f82-a8d9-3be96ca5a86d
@bind e PlutoUI.Slider(1:size(dcm_array, 3), default=10, show_value=true)

# ╔═╡ 03c94bf0-cbd1-4473-a242-186d81317a0e
heatmap(dcm_array[:, :, e], colormap=:grays)

# ╔═╡ c92f7ea0-08fb-4fa3-9bb6-8d0030a06754
md"""
## Calibration Prep
"""

# ╔═╡ b1d01296-26a5-4740-8a88-15fc76d1fa35
m_arr = masked_array[:, :, slice_CCI];

# ╔═╡ f93589a5-76c8-403c-9800-1796ca1f6229
md"""
### Intensity High Density
"""

# ╔═╡ da6e7dcb-92a1-4dcd-b99a-338dcef1a2ad
core_L_HD = Bool.(erode((((mask_L_HD)))));

# ╔═╡ bcac79e6-3218-4abe-8f49-f32d2a7ba9b8
mean_L_HD = mean(m_arr[core_L_HD])

# ╔═╡ 126688f5-d5a8-456b-9aee-44f8bffc29a4
md"""
#### Visualize High Density
"""

# ╔═╡ 0b63b3b6-0c43-42d9-a6be-f64139d04b96
begin
	arr_L_HD_cal = m_arr .* core_L_HD
	heatmap(arr_L_HD_cal, colormap=:grays)
end

# ╔═╡ 8fbbf202-176d-4582-acc6-dc2346486561
md"""
### Intensity Medium Density
"""

# ╔═╡ ad828113-99eb-4260-b715-560a0a559df7
core_L_MD = Bool.(erode(erode((mask_L_MD))));

# ╔═╡ 79305a4d-618c-43e9-a0f5-e8e3d57d15c2
mean_L_MD = mean(m_arr[core_L_MD])

# ╔═╡ 2e7a9cfc-b1a2-4d7f-8db7-69c71f537bc9
md"""
### Intensity Low Density
"""

# ╔═╡ f1756d4e-2a47-4061-94ad-e7601d13b10e
core_L_LD = Bool.(erode(erode(mask_L_LD)));

# ╔═╡ 65b074d6-cdd9-47ea-a84e-528b35f3e703
mean_L_LD = mean(m_arr[core_L_LD])

# ╔═╡ e9dd75fb-82c5-4d6f-98b6-c5f7990948ff
md"""
### Intensity No Calcium
"""

# ╔═╡ 900c7113-04aa-4914-860a-79b4bdde4ee5
ring = Bool.((dilate(dilate(dilate((dilate(mask_L_LD)))))) - dilate(dilate(dilate((mask_L_LD)))));

# ╔═╡ 08c99409-696e-483c-a19c-247b0d93a9a6
no_ca_arr = masked_array[:, :, 23];

# ╔═╡ 7a51a933-2328-477b-89bd-ea2ef510bedd
mean_no_CA = mean(no_ca_arr[mask])

# ╔═╡ 8530ab2f-2c28-4643-a933-911654d01aa6
md"""
### Calibration Line
"""

# ╔═╡ cfa9c97f-b010-48e9-b601-5a729eb0edfa
density_array = [0, 200, 400, 800] # mg/cc

# ╔═╡ 066ea0e5-42e0-4d61-8edb-5319f8ef61ea
intensity_array = [mean_no_CA, mean_L_LD, mean_L_MD, mean_L_HD] # HU

# ╔═╡ 920f7fa2-5229-428e-a95e-2d2b7355f062
df = DataFrame(:density => density_array, :intensity => intensity_array)

# ╔═╡ e66bcd52-d434-43e9-8469-0cc33280923b
linearRegressor = lm(@formula(intensity ~ density), df)

# ╔═╡ 1a26f30f-b9dd-4268-a402-6687432e64ff
linearFit = predict(linearRegressor)

# ╔═╡ 996729b5-95ef-4a5f-b10a-a5a581203d2b
m = linearRegressor.model.pp.beta0[2]

# ╔═╡ 678a969b-b9fd-4cd4-8c80-700a04d42d4c
b = linearRegressor.model.rr.mu[1]

# ╔═╡ 2d5152ad-2f05-4d2a-a651-09f9039e3063
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

# ╔═╡ fdba55a9-c04c-49d6-95fe-e05b9d75ba65
density(intensity) = (intensity - b) / m

# ╔═╡ 0c8a6dc2-2273-4377-8043-418f7f9dffa7
intensity(ρ) = m*ρ + b

# ╔═╡ 2f5f42e9-d96c-4ac7-9505-646f0393986c
begin
	f = Figure()
	ax1 = Axis(f[1, 1])
	
	scatter!(density_array, intensity_array)
	lines!(density_array, linearFit, color = :red)
	ax1.title = "Calibration Line (Intensity vs Density)"
	ax1.ylabel = "Intensity (HU)"
	ax1.xlabel = "Density (mg/cm^3)"
	
	f
end

# ╔═╡ ed01b680-a24d-44b9-8ef7-58ff4ab5a5ba
md"""
# Score Large Inserts
"""

# ╔═╡ 546334f3-b411-43f9-8d0f-061aed8f22ca
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ 8d9504f9-9fa8-435c-9d1b-1ae42f98aa13
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ 9097cbfa-cf27-411c-9b83-8efeda2cdee6
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

# ╔═╡ 3146fe0f-bfe6-4cc9-9782-0ad2a097e777
arr = masked_array[:, :, slice_CCI-2:slice_CCI+2];

# ╔═╡ 2694e3e0-f428-470b-a70c-f3b110d68f37
single_arr = masked_array[:, :, slice_CCI];

# ╔═╡ 0483eebd-8106-4444-9037-fb52057fb0b4
md"""
## High Density
"""

# ╔═╡ 99319b7b-639d-4207-852d-0bd54b0537a4
begin
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_HD_3D[:, :, z] = mask_L_HD
	end
end;

# ╔═╡ fc18a1ee-8042-44e1-bc4c-fb56e345d90d
md"""
#### Dilated mask
"""

# ╔═╡ 16057812-b925-4571-8d5b-35254c42683c
dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D));

# ╔═╡ 53d5b9ae-85b5-4f37-93ca-c76edd0d7916
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ 0a6451a3-19d6-4f4a-9e5b-cc2eb52b5b9d
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ 5bfd1bf6-e24c-42e3-b911-9d87f7c3491c
md"""
#### Ring (background) mask
"""

# ╔═╡ 588a129a-fe05-485f-8441-cd99dac82a03
ring_mask_L_HD = dilate(dilate(dilate(dilate(mask_L_HD_3D)))) - dilate(dilate(dilate(mask_L_HD_3D)));

# ╔═╡ 11b5faad-c162-4cf7-9b6f-2d61750b7a1a
@bind g4 overlay_mask_bind(ring_mask_L_HD)

# ╔═╡ 975f5f5c-3ada-4ff2-9e43-974f53d1b329
overlay_mask_plot(arr, ring_mask_L_HD, g4, "ring mask")

# ╔═╡ 820b0f29-867b-448f-9fb9-9ff1c3474422
md"""
### Calculations
"""

# ╔═╡ 7532aa8e-7dbf-46ad-86c2-7f9ac7f2feca
begin
	single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
	s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
end

# ╔═╡ b74ffab8-cdc6-4009-8e63-476471cecce3
S_Obj_HD = intensity(800)

# ╔═╡ af011205-0ba2-4e37-bced-a8fe9fe88c2f
begin
	alg_L_HD = Integrated(arr[mask_L_HD_3D])
	vol_l_hd = score(s_bkg_L_HD, S_Obj_HD, pixel_size, alg_L_HD)
end

# ╔═╡ 0db23983-2ff1-45c0-878b-e6837bed77d6
begin
	ρ_HD = 0.8 # mg/mm^3
	mass_l_hd = score(s_bkg_L_HD, S_Obj_HD, pixel_size, ρ_HD, alg_L_HD)
end

# ╔═╡ 02be0555-01e0-45d5-9711-cb203ef41b4c
md"""
## Medium Density
"""

# ╔═╡ 04f317d4-5994-488c-92f2-63cd2582d946
begin
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_MD_3D[:, :, z] = mask_L_MD
	end
end;

# ╔═╡ 93192cc0-b222-4506-b3f4-dc12f0cff32f
md"""
#### Dilated mask
"""

# ╔═╡ 265b30ee-7ad8-4cc2-851f-7374e1d28da8
dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D));

# ╔═╡ 04dfe430-4cf4-4cca-8fca-279fd00b286f
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ 184c7b7f-6a17-4a8f-9751-4a2f8d1878c1
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ fbb2fd43-4e42-4090-a4e6-1975d3957026
md"""
#### Ring (background) mask
"""

# ╔═╡ c08d31ac-2fbd-43d6-b334-5157ed18b51b
ring_mask_L_MD = dilate(dilate(dilate(dilate(mask_L_MD_3D)))) - dilate(dilate(dilate(mask_L_MD_3D)));

# ╔═╡ 3fc579af-153a-41ff-9ab4-a78df87764e0
@bind h4 overlay_mask_bind(ring_mask_L_MD)

# ╔═╡ 66e4f9f7-0f35-479c-851c-93318cd7715a
overlay_mask_plot(arr, ring_mask_L_MD, h4, "ring mask")

# ╔═╡ bb1950cc-0235-40f8-b1f5-ba4a1ae48e4b
md"""
### Calculations
"""

# ╔═╡ 5bc01566-11ac-443a-8528-f949f9697de8
begin
	single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
	s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
end

# ╔═╡ d421bb7e-9b75-43e1-8d5d-2203364dad3a
S_Obj_MD = intensity(400)

# ╔═╡ 553473b7-9a59-490f-bf47-becbc3f06bd5
begin
	alg_L_MD = Integrated(arr[mask_L_MD_3D])
	vol_l_md = score(s_bkg_L_MD, S_Obj_MD, pixel_size, alg_L_MD)
end

# ╔═╡ b960f3ce-9212-4f8e-97e3-a5ebd8b3c674
begin
	ρ_MD = 0.4 # mg/mm^3
	mass_l_md = score(s_bkg_L_MD, S_Obj_MD, pixel_size, ρ_MD, alg_L_MD)
end

# ╔═╡ 686d7820-e823-4d42-a4e1-5e463411b3e1
md"""
## Low Density
"""

# ╔═╡ 6113970c-9776-4cbc-9a3e-92cab30b5039
begin
	mask_L_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_LD_3D[:, :, z] = mask_L_LD
	end
end;

# ╔═╡ fe7df65d-c15d-4764-83c7-b18f38c4b35c
md"""
#### Dilated mask
"""

# ╔═╡ 8b5fa573-c221-458b-9bd3-c29fabb902c4
dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D));

# ╔═╡ 7b4df8a4-43bb-490a-9da4-0ac229bc21b4
@bind i2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ aa70bb01-446f-4ac2-b460-a8727b977bb3
overlay_mask_plot(arr, dilated_mask_L_LD, i2, "dilated mask")

# ╔═╡ 612e2381-330b-475f-8614-07c0559c365a
md"""
#### Ring (background) mask
"""

# ╔═╡ 41788e2f-c0c2-4404-875f-576646f97a27
ring_mask_L_LD = dilate(dilate(dilate(dilate(mask_L_LD_3D)))) - dilate(dilate(dilate(mask_L_LD_3D)));

# ╔═╡ 0f2444f8-1a8c-4ab6-abbc-5329aec8cf97
@bind i4 overlay_mask_bind(ring_mask_L_LD)

# ╔═╡ e6795bfb-cd7c-49c4-8e25-6eed0c26a476
overlay_mask_plot(arr, ring_mask_L_LD, i4, "ring mask")

# ╔═╡ 76c230bc-cdf4-4d9f-968f-60700e819de4
md"""
### Calculations
"""

# ╔═╡ 0d0c639e-71ff-479d-8163-961a10db47ab
begin	
	single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
	s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
end

# ╔═╡ d5aa0939-064b-430d-9273-dfa00dea29da
S_Obj_LD = intensity(200)

# ╔═╡ c6fb0eca-176a-45b3-a1ad-7e8aaa4803da
begin
	alg_L_LD = Integrated(arr[mask_L_LD_3D])
	vol_l_ld = score(s_bkg_L_LD, S_Obj_LD, pixel_size, alg_L_LD)
end

# ╔═╡ c3811bc1-8c11-4a26-a245-95f15c6a5f57
begin
	ρ_LD = 0.2
	mass_l_ld = score(s_bkg_L_LD, S_Obj_LD, pixel_size, ρ_LD, alg_L_LD)
end

# ╔═╡ afe4d899-11bc-4d69-b050-007c332c9120
md"""
# Score Medium Inserts
"""

# ╔═╡ cbe46322-d086-413d-9c7a-7e8f536831d0
md"""
## High Density
"""

# ╔═╡ 8ee23d24-b1e2-446f-a317-dbf058abd466
begin
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_HD_3D[:, :, z] = mask_M_HD
	end
end;

# ╔═╡ 27d1c5b7-8e89-4375-b555-32498265fd94
md"""
#### Dilated mask
"""

# ╔═╡ 0709e61c-dd72-44ba-86b6-6b57e3a28c17
dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 8b3ab168-9ce5-41c4-b41a-686dda708f8a
@bind j2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ e10ceddc-e096-41f2-b9d0-01c105d9abf9
overlay_mask_plot(arr, dilated_mask_M_HD, j2, "dilated mask")

# ╔═╡ e3644178-61fd-44c7-80c3-014f3b2bfe9d
md"""
#### Ring (background) mask
"""

# ╔═╡ 3eea0698-362c-4832-bbe1-4276c18260b4
ring_mask_M_HD = dilate(dilate(dilate(dilate(dilate(mask_M_HD_3D))))) - dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 8d7881b9-9dfe-4249-a4d9-d81646ef2e29
@bind j4 overlay_mask_bind(ring_mask_M_HD)

# ╔═╡ 4a38ec26-ce88-4998-8730-598868561897
overlay_mask_plot(arr, ring_mask_M_HD, j4, "ring mask")

# ╔═╡ 78e304c7-bb20-47d4-b644-8d268b04c495
md"""
### Calculations
"""

# ╔═╡ 0b9013f9-b3ed-47b7-92c9-255dd2a5492c
begin
	single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
	s_bkg_M_HD = mean(single_arr[single_ring_mask_M_HD])
end

# ╔═╡ 2275fec1-2c74-42a6-b410-f51b0b776d39
begin
	alg_M_HD = Integrated(arr[mask_M_HD_3D])
	vol_m_hd = score(s_bkg_M_HD, S_Obj_HD, pixel_size, alg_M_HD)
end

# ╔═╡ a9c0bf2d-c608-4d49-b903-64e8ca204569
mass_m_hd = score(s_bkg_M_HD, S_Obj_HD, pixel_size, ρ_HD, alg_M_HD)

# ╔═╡ 2ad58faa-6909-449c-9900-d0a57d4a0858
md"""
## Medium Density
"""

# ╔═╡ 07cf7bb6-1b04-47dc-af53-d5cccf7eed23
begin
	mask_M_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_MD_3D[:, :, z] = mask_M_MD
	end
end;

# ╔═╡ d324d896-a54f-46a9-aa29-27e0085f58f8
md"""
#### Dilated mask
"""

# ╔═╡ 3976ee5d-c7ca-4f26-8a00-e0e5e1ddd716
dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))));

# ╔═╡ a5509f31-91af-4fea-8d73-cad563f54d03
@bind k2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ 004fd1d8-fcae-4c75-abf4-b1f0379e3f74
overlay_mask_plot(arr, dilated_mask_M_MD, k2, "dilated mask")

# ╔═╡ b93be9ec-1e1d-4241-81de-f8720c9c19ee
md"""
#### Ring (background) mask
"""

# ╔═╡ 90867401-9abf-4177-84fe-fe98333e8eaf
ring_mask_M_MD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))));

# ╔═╡ 3a48a62d-cf82-45f9-a7ca-328fc0f63b24
@bind k4 overlay_mask_bind(ring_mask_M_MD)

# ╔═╡ 45c9f09d-4d38-448c-aa2d-40fcb9b4d68b
overlay_mask_plot(arr, ring_mask_M_MD, k4, "ring mask")

# ╔═╡ 9608c3dc-cc23-470a-ac7d-0ba47c3a1925
md"""
### Calculations
"""

# ╔═╡ 18c8d7e1-fedd-4c51-aedc-46bc757e112f
begin
	single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
	s_bkg_M_MD = mean(single_arr[single_ring_mask_M_MD])
end

# ╔═╡ f65e4903-4056-4bfb-9f34-155397f0c687
begin
	alg_M_MD = Integrated(arr[mask_M_MD_3D])
	vol_m_md = score(s_bkg_M_MD, S_Obj_MD, pixel_size, alg_M_MD)
end

# ╔═╡ 133f969c-66ad-4aee-b613-92d462ec0db3
mass_m_md = score(s_bkg_M_MD, S_Obj_MD, pixel_size, ρ_MD, alg_M_MD)

# ╔═╡ fcec4b1e-30d2-4234-b7c1-41df1e5c92f2
md"""
## Low Density
"""

# ╔═╡ 5e1aa4a6-e86a-442c-85d3-61a56bb64503
begin
	mask_M_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_LD_3D[:, :, z] = mask_M_LD
	end
end;

# ╔═╡ 90ad8780-9909-455a-b78d-c8bc6a2aed1d
md"""
#### Dilated mask
"""

# ╔═╡ 57b0c3b9-fb93-4a10-9631-95ef6c33a15c
dilated_mask_M_LD = dilate(dilate(dilate(dilate(mask_M_LD_3D))));

# ╔═╡ 5c3f479c-bd24-4062-af54-92ca0beeaa55
@bind l2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ 04c0033e-b9ef-4df1-92f6-78b53ade82b4
overlay_mask_plot(arr, dilated_mask_M_LD, l2, "dilated mask")

# ╔═╡ d97994f4-d171-444a-9b34-371de87e0352
md"""
#### Ring (background) mask
"""

# ╔═╡ f355865c-3496-4c68-8014-f537c444403b
ring_mask_M_LD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ ffe1a1b4-878f-4ba9-bda1-8545c0c1b0f7
@bind l4 overlay_mask_bind(ring_mask_M_LD)

# ╔═╡ 85b8c2a9-0dfa-4c4c-826f-f6bde4f2e0c0
overlay_mask_plot(arr, ring_mask_M_LD, l4, "ring mask")

# ╔═╡ f9667c35-28ae-4562-943b-d19a33c63bb4
md"""
### Calculations
"""

# ╔═╡ d785fc61-fda8-44eb-893f-77dbb2f416c1
begin
	single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
	s_bkg_M_LD = mean(single_arr[single_ring_mask_M_LD])
end

# ╔═╡ fa1d645d-0cdb-4695-a6d9-c66201ad73eb
begin
	alg_M_LD = Integrated(arr[mask_M_LD_3D])
	vol_m_ld = score(s_bkg_M_LD, S_Obj_LD, pixel_size, alg_M_LD)
end

# ╔═╡ cd7dd665-1311-4f28-8e90-019685584060
mass_m_ld = score(s_bkg_M_LD, S_Obj_LD, pixel_size, ρ_LD, alg_M_LD)

# ╔═╡ 8a808485-0be7-41be-b8ad-121ca7e9e675
md"""
# Score Small Inserts
"""

# ╔═╡ 8975cc7f-f705-4d46-96e4-a3d60b54a46a
md"""
## High Density
"""

# ╔═╡ 1afc8cc4-a4a7-40ce-8098-a4ec31c342d2
begin
	mask_S_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_HD_3D[:, :, z] = mask_S_HD
	end
end;

# ╔═╡ 1470eb83-e075-4a0f-aec6-096169ccb0a2
md"""
#### Dilated mask
"""

# ╔═╡ a918affc-de9c-41dd-981b-0102846902c4
dilated_mask_S_HD = dilate(dilate(dilate(dilate(mask_S_HD_3D))));

# ╔═╡ 044f678c-9542-413a-92b9-9adf10af79b3
@bind m2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ b7b57c8f-d431-4109-9572-d9ef7307e47c
overlay_mask_plot(arr, dilated_mask_S_HD, m2, "dilated mask")

# ╔═╡ e5addbfc-1262-499f-b0b1-1e4be3a0ba48
md"""
#### Ring (background) mask
"""

# ╔═╡ f77e1a56-97b5-4c24-88f1-4cffc3bbdafa
ring_mask_S_HD = dilate(dilate(dilate(dilate(dilate(mask_S_HD_3D))))) - dilate(dilate(dilate(dilate(mask_S_HD_3D))));

# ╔═╡ f712abed-a68c-4605-b991-8e977f2391a3
@bind m4 overlay_mask_bind(ring_mask_S_HD)

# ╔═╡ d48901a6-e3a1-468e-b5cc-ec110a7b8b9f
overlay_mask_plot(arr, ring_mask_S_HD, m4, "ring mask")

# ╔═╡ 0c7571ca-de27-49ef-b44f-c21b25681349
md"""
### Calculations
"""

# ╔═╡ 116b862f-3936-45d3-b085-ee6a557fb91a
begin
	single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
	s_bkg_S_HD = mean(single_arr[single_ring_mask_S_HD])
end

# ╔═╡ 4156a33d-e231-465b-91e5-87b853027098
begin
	alg_S_HD = Integrated(arr[mask_S_HD_3D])
	vol_s_hd = score(s_bkg_S_HD, S_Obj_HD, pixel_size, alg_S_HD)
end

# ╔═╡ 6bccdc83-d013-45c0-b939-12ed7dc7f5a8
mass_s_hd = score(s_bkg_S_HD, S_Obj_HD, pixel_size, ρ_HD, alg_S_HD)

# ╔═╡ 1169d8aa-e370-48fe-a9a1-9de5f7dbd716
md"""
## Medium Density
"""

# ╔═╡ 5677bcfe-e320-4ead-a000-4f53b9381df3
begin
	mask_S_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_MD_3D[:, :, z] = mask_S_MD
	end
end;

# ╔═╡ 95ced72e-ccce-412a-b4f6-99e45233a4ca
md"""
#### Dilated mask
"""

# ╔═╡ d386c3a1-f71e-4463-bec3-61f9be22aa5b
dilated_mask_S_MD = dilate(dilate(dilate(dilate(mask_S_MD_3D))));

# ╔═╡ 7207c1dd-bc61-4dc6-9e4b-13612ba32c88
@bind n2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ 960cbd14-3cfd-41fb-b0e0-628e7571b0b0
overlay_mask_plot(arr, dilated_mask_S_MD, n2, "dilated mask")

# ╔═╡ 47c54cf8-6a91-46bb-9607-6cb7d3ba57d8
md"""
#### Ring (background) mask
"""

# ╔═╡ e02d0c8a-45a5-4365-8833-22f481f49f83
ring_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D))))) - dilate(dilate(dilate(dilate(mask_S_MD_3D))));

# ╔═╡ aa35ede5-cafb-4ace-9562-a91a6de11e31
@bind n4 overlay_mask_bind(ring_mask_S_MD)

# ╔═╡ f9793d90-916b-435b-834e-6cae9f29c35a
overlay_mask_plot(arr, ring_mask_S_MD, n4, "ring mask")

# ╔═╡ ca7b1fd0-0747-49e2-a5c4-af56aeaa5dfa
md"""
### Calculations
"""

# ╔═╡ 322111ba-3b91-4260-8475-acfdbcfa16ae
begin
	single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
	s_bkg_S_MD = mean(single_arr[single_ring_mask_S_MD])
end

# ╔═╡ 2de9a556-982c-4fb8-ba27-2587c8fb7f73
begin
	alg_S_MD = Integrated(arr[mask_S_MD_3D])
	vol_s_md = score(s_bkg_S_MD, S_Obj_HD, pixel_size, alg_S_MD)
end

# ╔═╡ 6c653671-fe17-426f-8d1c-91db79be4f9c
mass_s_md = score(s_bkg_S_MD, S_Obj_HD, pixel_size, ρ_MD, alg_S_MD)

# ╔═╡ 896271ac-b7c4-48db-a407-b6c349901d8c
md"""
## Low Density
"""

# ╔═╡ ab266ea4-3aa3-41d4-bd9d-c3ae519d20ca
begin
	mask_S_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_LD_3D[:, :, z] = mask_S_LD
	end
end;

# ╔═╡ 0c970697-fb8c-43f2-9023-050bb371b8c8
md"""
#### Dilated mask
"""

# ╔═╡ ed4ac4ce-fe4f-4f63-a069-d81a7b217117
dilated_mask_S_LD = dilate(dilate(dilate(dilate(mask_S_LD_3D))));

# ╔═╡ abe7ee15-dd92-482f-a673-90aac90e2dd6
@bind o2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ b85d36dd-5a81-42c4-b173-404926e75308
overlay_mask_plot(arr, dilated_mask_S_LD, o2, "dilated mask")

# ╔═╡ 0dbb862f-595b-4a3b-bcbb-e5686ef4c594
md"""
#### Ring (background) mask
"""

# ╔═╡ 73a4767a-7b1d-4d36-9c22-f14917c3909e
ring_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D))))) - dilate(dilate(dilate(dilate(mask_S_LD_3D))));

# ╔═╡ 01377eb9-a5c3-496c-aab8-fda0bfa5c3c7
@bind o4 overlay_mask_bind(ring_mask_S_LD)

# ╔═╡ f878a2b1-ad66-473f-82da-6089b68f220e
overlay_mask_plot(arr, ring_mask_S_LD, o4, "ring mask")

# ╔═╡ 7866b3d0-5b77-41d1-af22-7d52de6d9a4f
md"""
### Calculations
"""

# ╔═╡ eb71f843-fa07-42a1-b0c5-e1bed1cfdde8
begin
	single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
	s_bkg_S_LD = mean(single_arr[single_ring_mask_S_LD])
end

# ╔═╡ 639ed28d-1516-4a28-8ed7-4483e686ef7a
begin
	alg_S_LD = Integrated(arr[mask_S_LD_3D])
	vol_s_ld = score(s_bkg_S_LD, S_Obj_HD, pixel_size, alg_S_LD)
end

# ╔═╡ 97a1b5ec-262c-4d63-8730-7213ff42abef
mass_s_ld = score(s_bkg_S_LD, S_Obj_HD, pixel_size, ρ_LD, alg_S_LD)

# ╔═╡ d021f42d-33a0-4681-9669-a44f491626e3
md"""
# Results
"""

# ╔═╡ 475696f1-7594-4f79-a6a2-952f649fc9df
md"""
### Volume
"""

# ╔═╡ 39d67c99-7170-4581-9e00-7fceaee9626a
inserts = [
	"Low Density",
	"Medium Density",
	"High Density"
]

# ╔═╡ 85bc17ef-9cae-40c6-a0ab-fbaf78aedb36
ground_truth_volume_large = [
	98.2,
	98.2,
	98.2,
] # mm^3

# ╔═╡ f8acffe4-293b-4393-a744-c0c4ba011d99
calculated_volume_large = [
	vol_l_ld,
	vol_l_md,
	vol_l_hd
]

# ╔═╡ e3c6d740-e6cd-4c68-bf4a-762a4aea53d7
ground_truth_volume_medium = [
	21.2,
	21.2,
	21.2
]

# ╔═╡ 1f9963bc-755a-4e16-8a84-ceb30484e5d6
calculated_volume_medium = [
	vol_m_ld,
	vol_m_md,
	vol_m_hd
]

# ╔═╡ 7bcdaa45-b383-4225-b04a-f7a10d0e530e
ground_truth_volume_small = [
	0.8,
	0.8,
	0.8
]

# ╔═╡ 7b42be71-e0f7-439b-a028-8344ef4a18ae
calculated_volume_small = [
	vol_s_ld,
	vol_s_md,
	vol_s_hd
]

# ╔═╡ 8634fef2-5c95-4b45-bbce-7e2fb9d3e4ca
df1 = DataFrame(
	inserts = inserts,
	ground_truth_volume_large = ground_truth_volume_large,
	calculated_volume_large = calculated_volume_large,
	ground_truth_volume_medium = ground_truth_volume_medium,
	calculated_volume_medium = calculated_volume_medium,
	ground_truth_volume_small = ground_truth_volume_small,
	calculated_volume_small = calculated_volume_small
)

# ╔═╡ 7b5413c6-1946-4c89-a983-dbec78b52035
begin
	fvol2 = Figure()
	axvol2 = Axis(fvol2[1, 1])
	
	scatter!(density_array[2:end], df1[!, :ground_truth_volume_large], label="ground_truth_volume_large")
	scatter!(density_array[2:end], df1[!, :calculated_volume_large], label="calculated_volume_large")
	
	axvol2.title = "Volume Measurements (Large)"
	axvol2.ylabel = "Volume (mm^3)"
	axvol2.xlabel = "Density (mg/cm^3)"

	xlims!(axvol2, 0, 850)
	ylims!(axvol2, 0, 130)
	
	fvol2[1, 2] = Legend(fvol2, axvol2, framevisible = false)
	
	fvol2
end

# ╔═╡ a6e0ce7b-cfa0-4eba-8eca-2f4f3f82ca2d
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

# ╔═╡ 4aaf0fdd-ea48-4231-bfe6-53d05c5163e8
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

# ╔═╡ 9d9e8022-990b-4023-91aa-b87bbe20df6c
md"""
### Mass
"""

# ╔═╡ 097a38db-286f-40b2-8690-3e0a2a58b284
ground_truth_mass_large = [
	19.6,
	39.3,
	78.5
] # mg

# ╔═╡ cfd4c13a-6461-40f8-ab88-b2a05eab2513
calculated_mass_large = [
	mass_l_ld,
	mass_l_md,
	mass_l_hd
]

# ╔═╡ 4d50a14c-5993-42d3-8b19-6e583e6fc944
ground_truth_mass_medium = [
	4.2,
	8.5,
	17.0
]

# ╔═╡ 5ca8fc85-00aa-4372-a001-3c450a2be85c
calculated_mass_medium = [
	mass_m_ld,
	mass_m_md,
	mass_m_hd
]

# ╔═╡ 09fb42d8-af42-4301-afaa-aaa742ea5ab5
ground_truth_mass_small = [
	0.2,
	0.3,
	0.6
]

# ╔═╡ 7eb9b991-a469-4114-b4a9-236f098238d5
calculated_mass_small = [
	mass_s_ld,
	mass_s_md,
	mass_s_hd
]

# ╔═╡ a2a487c8-cacd-4e64-b0ea-bd488e50cfec
df2 = DataFrame(
	inserts = inserts,
	ground_truth_mass_large = ground_truth_mass_large,
	calculated_mass_large = calculated_mass_large,
	ground_truth_mass_medium = ground_truth_mass_medium,
	calculated_mass_medium = calculated_mass_medium,
	ground_truth_mass_small = ground_truth_mass_small,
	calculated_mass_small = calculated_mass_small
)

# ╔═╡ 03b6b414-13aa-4185-ae8f-b1903069c266
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

# ╔═╡ b202cace-72c1-4f83-8493-0ee7b968a238
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

# ╔═╡ 32610d65-672d-47a4-bb28-b69d43ac8592
begin
	fmass4 = Figure()
	axmass4 = Axis(fmass4[1, 1])
	
	scatter!(density_array[2:end], df2[!, :ground_truth_mass_small], label="ground_truth_mass_small")
	scatter!(density_array[2:end], df2[!, :calculated_mass_small], label="calculated_mass_small")
	
	axmass4.title = "Mass Measurements (Small)"
	axmass4.ylabel = "Mass (mg)"
	axmass4.xlabel = "Density (mg/cm^3)"

	xlims!(axmass4, 0, 850)
	ylims!(axmass4, 0, 1.2)
	
	fmass4[1, 2] = Legend(fmass4, axmass4, framevisible = false)
	
	fmass4
end

# ╔═╡ 07aa5b5e-b98a-4eec-8f13-70ac034dfd32
md"""
### Save Results
"""

# ╔═╡ c30dd8f3-1cd4-418a-af45-931ce0c705ff
df_final = leftjoin(df1, df2, on=:inserts)

# ╔═╡ 6ace3811-70cf-4732-89a4-706c17bb75ed
output_path = string(cd(pwd, "..") , "/data/output/", VENDER, "/", scan)

# ╔═╡ 6426d783-b13c-4517-800f-3e54f6e7eee7
CSV.write(output_path, df_final)

# ╔═╡ Cell order:
# ╠═c50510a7-e981-4933-a4b1-fab2efcc43a2
# ╠═97a264e6-1037-4862-b6bd-2d2486bd419e
# ╟─9ecf691f-8523-4ea9-a635-7482a1a7b697
# ╠═c72aaf65-8ba2-4c81-9bad-a059b983a833
# ╠═e33a903d-d9fb-4960-bc81-128cfbff8af6
# ╠═762cde8a-82c0-4cf0-b2b3-8ed00a4abb3d
# ╟─416c3d5f-9953-4a0d-be1c-1049e39edd4b
# ╠═916be3b8-b7e3-447d-82a2-07e5185443e1
# ╠═4669bb68-6224-4a31-ae50-c2e44744b87e
# ╠═df50149d-ad91-4607-8cc1-6a17da8d8562
# ╠═35c5e1e1-c88d-4202-8f7d-0e6c434ea161
# ╠═5cc4cb9b-b0d7-4158-816f-0cbd4d7a7463
# ╠═1fae526e-7db7-46e8-aa63-19882cb2e172
# ╟─fea80b94-e9c5-4ff7-9b91-4e35bf03b5e3
# ╠═49184d1a-1e84-48c7-8f71-d179fd4671c5
# ╟─21266f0e-b535-4e72-950e-ba7d125684de
# ╟─536506ce-2031-4c32-a4e3-d6fc6d7da27c
# ╟─e2de04cd-2f66-4192-abae-95693b1b5176
# ╟─be5f9be6-e10c-48c1-b9fe-ed229a1dff26
# ╟─c0704706-87c2-4078-8345-7e2097639492
# ╟─6e597bb0-fa70-46f3-8526-dce3cd4e7b78
# ╠═07c0f484-6225-4c6f-ab20-1941bd6abc7a
# ╠═e9af6401-8f3a-4f86-8d23-02d88dc7ee48
# ╠═4f795d56-a67b-4b9f-be94-c37890a5c4c4
# ╟─69d3e9a8-68c9-4a5a-8a4a-86c5fd64a08a
# ╠═0508d26f-09a1-457b-8cb0-5693f5c9ba3f
# ╠═544ae005-6661-487f-8662-ea3006ec4069
# ╠═2647d239-58d6-4448-9e3f-66918a10dbad
# ╟─2afdc85a-a8bf-4fd4-bdef-74c39c0fb452
# ╠═c2394467-eb15-4610-9bfd-7b3f88a7fb00
# ╟─9b6e3a13-9044-4b7e-89e0-b552af52bd2b
# ╠═0f85bfc3-f0d8-4ae7-b5b4-f72a072a3d55
# ╟─b8892b4e-3be7-4f82-a8d9-3be96ca5a86d
# ╠═03c94bf0-cbd1-4473-a242-186d81317a0e
# ╟─c92f7ea0-08fb-4fa3-9bb6-8d0030a06754
# ╠═b1d01296-26a5-4740-8a88-15fc76d1fa35
# ╟─f93589a5-76c8-403c-9800-1796ca1f6229
# ╠═da6e7dcb-92a1-4dcd-b99a-338dcef1a2ad
# ╠═bcac79e6-3218-4abe-8f49-f32d2a7ba9b8
# ╟─126688f5-d5a8-456b-9aee-44f8bffc29a4
# ╠═0b63b3b6-0c43-42d9-a6be-f64139d04b96
# ╟─8fbbf202-176d-4582-acc6-dc2346486561
# ╠═ad828113-99eb-4260-b715-560a0a559df7
# ╠═79305a4d-618c-43e9-a0f5-e8e3d57d15c2
# ╟─2e7a9cfc-b1a2-4d7f-8db7-69c71f537bc9
# ╠═f1756d4e-2a47-4061-94ad-e7601d13b10e
# ╠═65b074d6-cdd9-47ea-a84e-528b35f3e703
# ╟─e9dd75fb-82c5-4d6f-98b6-c5f7990948ff
# ╠═900c7113-04aa-4914-860a-79b4bdde4ee5
# ╠═08c99409-696e-483c-a19c-247b0d93a9a6
# ╠═7a51a933-2328-477b-89bd-ea2ef510bedd
# ╟─8530ab2f-2c28-4643-a933-911654d01aa6
# ╠═cfa9c97f-b010-48e9-b601-5a729eb0edfa
# ╠═066ea0e5-42e0-4d61-8edb-5319f8ef61ea
# ╠═920f7fa2-5229-428e-a95e-2d2b7355f062
# ╠═e66bcd52-d434-43e9-8469-0cc33280923b
# ╠═1a26f30f-b9dd-4268-a402-6687432e64ff
# ╠═996729b5-95ef-4a5f-b10a-a5a581203d2b
# ╠═678a969b-b9fd-4cd4-8c80-700a04d42d4c
# ╟─2d5152ad-2f05-4d2a-a651-09f9039e3063
# ╠═fdba55a9-c04c-49d6-95fe-e05b9d75ba65
# ╠═0c8a6dc2-2273-4377-8043-418f7f9dffa7
# ╟─2f5f42e9-d96c-4ac7-9505-646f0393986c
# ╟─ed01b680-a24d-44b9-8ef7-58ff4ab5a5ba
# ╠═546334f3-b411-43f9-8d0f-061aed8f22ca
# ╠═8d9504f9-9fa8-435c-9d1b-1ae42f98aa13
# ╠═9097cbfa-cf27-411c-9b83-8efeda2cdee6
# ╠═3146fe0f-bfe6-4cc9-9782-0ad2a097e777
# ╠═2694e3e0-f428-470b-a70c-f3b110d68f37
# ╟─0483eebd-8106-4444-9037-fb52057fb0b4
# ╠═99319b7b-639d-4207-852d-0bd54b0537a4
# ╟─fc18a1ee-8042-44e1-bc4c-fb56e345d90d
# ╠═16057812-b925-4571-8d5b-35254c42683c
# ╟─53d5b9ae-85b5-4f37-93ca-c76edd0d7916
# ╠═0a6451a3-19d6-4f4a-9e5b-cc2eb52b5b9d
# ╟─5bfd1bf6-e24c-42e3-b911-9d87f7c3491c
# ╠═588a129a-fe05-485f-8441-cd99dac82a03
# ╟─11b5faad-c162-4cf7-9b6f-2d61750b7a1a
# ╟─975f5f5c-3ada-4ff2-9e43-974f53d1b329
# ╟─820b0f29-867b-448f-9fb9-9ff1c3474422
# ╠═7532aa8e-7dbf-46ad-86c2-7f9ac7f2feca
# ╠═b74ffab8-cdc6-4009-8e63-476471cecce3
# ╠═af011205-0ba2-4e37-bced-a8fe9fe88c2f
# ╠═0db23983-2ff1-45c0-878b-e6837bed77d6
# ╟─02be0555-01e0-45d5-9711-cb203ef41b4c
# ╠═04f317d4-5994-488c-92f2-63cd2582d946
# ╟─93192cc0-b222-4506-b3f4-dc12f0cff32f
# ╠═265b30ee-7ad8-4cc2-851f-7374e1d28da8
# ╟─04dfe430-4cf4-4cca-8fca-279fd00b286f
# ╟─184c7b7f-6a17-4a8f-9751-4a2f8d1878c1
# ╟─fbb2fd43-4e42-4090-a4e6-1975d3957026
# ╠═c08d31ac-2fbd-43d6-b334-5157ed18b51b
# ╟─3fc579af-153a-41ff-9ab4-a78df87764e0
# ╟─66e4f9f7-0f35-479c-851c-93318cd7715a
# ╟─bb1950cc-0235-40f8-b1f5-ba4a1ae48e4b
# ╠═5bc01566-11ac-443a-8528-f949f9697de8
# ╠═d421bb7e-9b75-43e1-8d5d-2203364dad3a
# ╠═553473b7-9a59-490f-bf47-becbc3f06bd5
# ╠═b960f3ce-9212-4f8e-97e3-a5ebd8b3c674
# ╟─686d7820-e823-4d42-a4e1-5e463411b3e1
# ╠═6113970c-9776-4cbc-9a3e-92cab30b5039
# ╟─fe7df65d-c15d-4764-83c7-b18f38c4b35c
# ╠═8b5fa573-c221-458b-9bd3-c29fabb902c4
# ╟─7b4df8a4-43bb-490a-9da4-0ac229bc21b4
# ╟─aa70bb01-446f-4ac2-b460-a8727b977bb3
# ╟─612e2381-330b-475f-8614-07c0559c365a
# ╠═41788e2f-c0c2-4404-875f-576646f97a27
# ╟─0f2444f8-1a8c-4ab6-abbc-5329aec8cf97
# ╟─e6795bfb-cd7c-49c4-8e25-6eed0c26a476
# ╟─76c230bc-cdf4-4d9f-968f-60700e819de4
# ╠═0d0c639e-71ff-479d-8163-961a10db47ab
# ╠═d5aa0939-064b-430d-9273-dfa00dea29da
# ╠═c6fb0eca-176a-45b3-a1ad-7e8aaa4803da
# ╠═c3811bc1-8c11-4a26-a245-95f15c6a5f57
# ╟─afe4d899-11bc-4d69-b050-007c332c9120
# ╟─cbe46322-d086-413d-9c7a-7e8f536831d0
# ╠═8ee23d24-b1e2-446f-a317-dbf058abd466
# ╟─27d1c5b7-8e89-4375-b555-32498265fd94
# ╠═0709e61c-dd72-44ba-86b6-6b57e3a28c17
# ╟─8b3ab168-9ce5-41c4-b41a-686dda708f8a
# ╟─e10ceddc-e096-41f2-b9d0-01c105d9abf9
# ╟─e3644178-61fd-44c7-80c3-014f3b2bfe9d
# ╠═3eea0698-362c-4832-bbe1-4276c18260b4
# ╟─8d7881b9-9dfe-4249-a4d9-d81646ef2e29
# ╟─4a38ec26-ce88-4998-8730-598868561897
# ╟─78e304c7-bb20-47d4-b644-8d268b04c495
# ╠═0b9013f9-b3ed-47b7-92c9-255dd2a5492c
# ╠═2275fec1-2c74-42a6-b410-f51b0b776d39
# ╠═a9c0bf2d-c608-4d49-b903-64e8ca204569
# ╟─2ad58faa-6909-449c-9900-d0a57d4a0858
# ╠═07cf7bb6-1b04-47dc-af53-d5cccf7eed23
# ╟─d324d896-a54f-46a9-aa29-27e0085f58f8
# ╠═3976ee5d-c7ca-4f26-8a00-e0e5e1ddd716
# ╟─a5509f31-91af-4fea-8d73-cad563f54d03
# ╟─004fd1d8-fcae-4c75-abf4-b1f0379e3f74
# ╟─b93be9ec-1e1d-4241-81de-f8720c9c19ee
# ╠═90867401-9abf-4177-84fe-fe98333e8eaf
# ╟─3a48a62d-cf82-45f9-a7ca-328fc0f63b24
# ╟─45c9f09d-4d38-448c-aa2d-40fcb9b4d68b
# ╟─9608c3dc-cc23-470a-ac7d-0ba47c3a1925
# ╠═18c8d7e1-fedd-4c51-aedc-46bc757e112f
# ╠═f65e4903-4056-4bfb-9f34-155397f0c687
# ╠═133f969c-66ad-4aee-b613-92d462ec0db3
# ╟─fcec4b1e-30d2-4234-b7c1-41df1e5c92f2
# ╠═5e1aa4a6-e86a-442c-85d3-61a56bb64503
# ╟─90ad8780-9909-455a-b78d-c8bc6a2aed1d
# ╠═57b0c3b9-fb93-4a10-9631-95ef6c33a15c
# ╟─5c3f479c-bd24-4062-af54-92ca0beeaa55
# ╟─04c0033e-b9ef-4df1-92f6-78b53ade82b4
# ╟─d97994f4-d171-444a-9b34-371de87e0352
# ╠═f355865c-3496-4c68-8014-f537c444403b
# ╟─ffe1a1b4-878f-4ba9-bda1-8545c0c1b0f7
# ╟─85b8c2a9-0dfa-4c4c-826f-f6bde4f2e0c0
# ╟─f9667c35-28ae-4562-943b-d19a33c63bb4
# ╠═d785fc61-fda8-44eb-893f-77dbb2f416c1
# ╠═fa1d645d-0cdb-4695-a6d9-c66201ad73eb
# ╠═cd7dd665-1311-4f28-8e90-019685584060
# ╟─8a808485-0be7-41be-b8ad-121ca7e9e675
# ╟─8975cc7f-f705-4d46-96e4-a3d60b54a46a
# ╠═1afc8cc4-a4a7-40ce-8098-a4ec31c342d2
# ╟─1470eb83-e075-4a0f-aec6-096169ccb0a2
# ╠═a918affc-de9c-41dd-981b-0102846902c4
# ╟─044f678c-9542-413a-92b9-9adf10af79b3
# ╠═b7b57c8f-d431-4109-9572-d9ef7307e47c
# ╟─e5addbfc-1262-499f-b0b1-1e4be3a0ba48
# ╠═f77e1a56-97b5-4c24-88f1-4cffc3bbdafa
# ╟─f712abed-a68c-4605-b991-8e977f2391a3
# ╠═d48901a6-e3a1-468e-b5cc-ec110a7b8b9f
# ╟─0c7571ca-de27-49ef-b44f-c21b25681349
# ╠═116b862f-3936-45d3-b085-ee6a557fb91a
# ╠═4156a33d-e231-465b-91e5-87b853027098
# ╠═6bccdc83-d013-45c0-b939-12ed7dc7f5a8
# ╟─1169d8aa-e370-48fe-a9a1-9de5f7dbd716
# ╠═5677bcfe-e320-4ead-a000-4f53b9381df3
# ╟─95ced72e-ccce-412a-b4f6-99e45233a4ca
# ╠═d386c3a1-f71e-4463-bec3-61f9be22aa5b
# ╟─7207c1dd-bc61-4dc6-9e4b-13612ba32c88
# ╠═960cbd14-3cfd-41fb-b0e0-628e7571b0b0
# ╟─47c54cf8-6a91-46bb-9607-6cb7d3ba57d8
# ╠═e02d0c8a-45a5-4365-8833-22f481f49f83
# ╟─aa35ede5-cafb-4ace-9562-a91a6de11e31
# ╠═f9793d90-916b-435b-834e-6cae9f29c35a
# ╟─ca7b1fd0-0747-49e2-a5c4-af56aeaa5dfa
# ╠═322111ba-3b91-4260-8475-acfdbcfa16ae
# ╠═2de9a556-982c-4fb8-ba27-2587c8fb7f73
# ╠═6c653671-fe17-426f-8d1c-91db79be4f9c
# ╟─896271ac-b7c4-48db-a407-b6c349901d8c
# ╠═ab266ea4-3aa3-41d4-bd9d-c3ae519d20ca
# ╟─0c970697-fb8c-43f2-9023-050bb371b8c8
# ╠═ed4ac4ce-fe4f-4f63-a069-d81a7b217117
# ╟─abe7ee15-dd92-482f-a673-90aac90e2dd6
# ╠═b85d36dd-5a81-42c4-b173-404926e75308
# ╟─0dbb862f-595b-4a3b-bcbb-e5686ef4c594
# ╠═73a4767a-7b1d-4d36-9c22-f14917c3909e
# ╟─01377eb9-a5c3-496c-aab8-fda0bfa5c3c7
# ╠═f878a2b1-ad66-473f-82da-6089b68f220e
# ╟─7866b3d0-5b77-41d1-af22-7d52de6d9a4f
# ╠═eb71f843-fa07-42a1-b0c5-e1bed1cfdde8
# ╠═639ed28d-1516-4a28-8ed7-4483e686ef7a
# ╠═97a1b5ec-262c-4d63-8730-7213ff42abef
# ╟─d021f42d-33a0-4681-9669-a44f491626e3
# ╟─475696f1-7594-4f79-a6a2-952f649fc9df
# ╠═39d67c99-7170-4581-9e00-7fceaee9626a
# ╠═85bc17ef-9cae-40c6-a0ab-fbaf78aedb36
# ╠═f8acffe4-293b-4393-a744-c0c4ba011d99
# ╠═e3c6d740-e6cd-4c68-bf4a-762a4aea53d7
# ╠═1f9963bc-755a-4e16-8a84-ceb30484e5d6
# ╠═7bcdaa45-b383-4225-b04a-f7a10d0e530e
# ╠═7b42be71-e0f7-439b-a028-8344ef4a18ae
# ╠═8634fef2-5c95-4b45-bbce-7e2fb9d3e4ca
# ╟─7b5413c6-1946-4c89-a983-dbec78b52035
# ╟─a6e0ce7b-cfa0-4eba-8eca-2f4f3f82ca2d
# ╟─4aaf0fdd-ea48-4231-bfe6-53d05c5163e8
# ╟─9d9e8022-990b-4023-91aa-b87bbe20df6c
# ╠═097a38db-286f-40b2-8690-3e0a2a58b284
# ╠═cfd4c13a-6461-40f8-ab88-b2a05eab2513
# ╠═4d50a14c-5993-42d3-8b19-6e583e6fc944
# ╠═5ca8fc85-00aa-4372-a001-3c450a2be85c
# ╠═09fb42d8-af42-4301-afaa-aaa742ea5ab5
# ╠═7eb9b991-a469-4114-b4a9-236f098238d5
# ╠═a2a487c8-cacd-4e64-b0ea-bd488e50cfec
# ╟─03b6b414-13aa-4185-ae8f-b1903069c266
# ╟─b202cace-72c1-4f83-8493-0ee7b968a238
# ╟─32610d65-672d-47a4-bb28-b69d43ac8592
# ╟─07aa5b5e-b98a-4eec-8f13-70ac034dfd32
# ╠═c30dd8f3-1cd4-418a-af45-931ce0c705ff
# ╠═6ace3811-70cf-4732-89a4-706c17bb75ed
# ╠═6426d783-b13c-4517-800f-3e54f6e7eee7

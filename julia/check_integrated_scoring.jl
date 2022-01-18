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
	using DataFrames
	using GLM
	using DICOM
	using DICOMUtils
	using Phantoms
	using CalciumScoring
end

# ╔═╡ 97a264e6-1037-4862-b6bd-2d2486bd419e
TableOfContents()

# ╔═╡ 84926791-1a72-4b9a-add4-7707efa16a08
path = string(cd(pwd, "..") , "/", "data/Large_rep1")

# ╔═╡ 448bb1f2-6548-4bc8-bc7b-99c093833e7e
dcms = dcmdir_parse(path);

# ╔═╡ 2df42240-1fe7-47bb-baee-3150a6408393
dcm_array = load_dcm_array(dcms);

# ╔═╡ c515d4f4-38b1-4490-8201-774bfec91440
header = dcms[1].meta;

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
@bind b1 PlutoUI.Slider(1:size(calcium_image, 3), default=10, show_value=true)

# ╔═╡ 4f795d56-a67b-4b9f-be94-c37890a5c4c4
heatmap(calcium_image[:, :, b1], colormap=:grays)

# ╔═╡ 78793039-df7f-4dbf-9188-155ecc72d944
c_img = calcium_image[:, :, 12];

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

# ╔═╡ 9d0a4171-eb47-4cbf-8f72-4adf041a11d1
@bind c PlutoUI.Slider(1:size(masks, 3), default=10, show_value=true)

# ╔═╡ 2647d239-58d6-4448-9e3f-66918a10dbad
heatmap(masks, colormap=:grays)

# ╔═╡ 2afdc85a-a8bf-4fd4-bdef-74c39c0fb452
md"""
## Overlay Mask Calcium Inserts
"""

# ╔═╡ c2394467-eb15-4610-9bfd-7b3f88a7fb00
arr_L_HD = masked_array[:, :, 23:28] .* mask_L_HD;

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
m_arr = masked_array[:, :, 25];

# ╔═╡ f93589a5-76c8-403c-9800-1796ca1f6229
md"""
### Intensity High Density
"""

# ╔═╡ da6e7dcb-92a1-4dcd-b99a-338dcef1a2ad
core_L_HD = Bool.(erode(erode(erode((mask_L_HD)))));

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
core_L_MD = Bool.(erode(erode(mask_L_MD)));

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
density_array = [0, 200, 400, 800] # g/cc

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
arr = masked_array[:, :, 23:27];

# ╔═╡ 2694e3e0-f428-470b-a70c-f3b110d68f37
single_arr = arr[:, :, 3];

# ╔═╡ 78b0b7fc-ad4a-4121-bd30-2fdadcd95ea0
# @bind f1 PlutoUI.Slider(1:size(arr, 3), default=3, show_value=true)

# ╔═╡ 13f3098d-e6ee-4fac-a429-a588cc8a2bbd
# heatmap(arr[:, :, f1], colormap=:grays)

# ╔═╡ 384ff1c4-91a2-44b8-ac63-756c813c8432
# heatmap(single_arr, colormap=:grays)

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

# ╔═╡ bb6dd528-05b1-48c1-b1f6-f9ef2d2b6777
md"""
#### Original mask
"""

# ╔═╡ 44f53f62-6a13-45b3-a1e9-ceea1576cbf4
@bind g1 overlay_mask_bind(mask_L_HD_3D)

# ╔═╡ 777e64c7-f05d-4af9-9417-268c3ea9f5b4
overlay_mask_plot(arr, mask_L_HD_3D, g1, "original mask")

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

# ╔═╡ 6fb82cc1-94c8-4a68-af6b-d69201a0fc6d
md"""
#### Eroded (object) mask
"""

# ╔═╡ e999b3f6-7d6f-4945-81fd-6b74216f0eec
eroded_mask_L_HD = erode(erode(mask_L_HD_3D));

# ╔═╡ 0608e3a1-8a56-409c-8258-e2919c51be5a
@bind g3 overlay_mask_bind(eroded_mask_L_HD)

# ╔═╡ 1d2fd28a-25f1-4c81-94a4-70317d8172aa
overlay_mask_plot(arr, eroded_mask_L_HD, g3, "eroded mask")

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
### Calculate Volume
"""

# ╔═╡ 7532aa8e-7dbf-46ad-86c2-7f9ac7f2feca
begin
	single_erode_mask_L_HD = eroded_mask_L_HD[:, :, 3]
	s_obj_L_HD = mean(single_arr[single_erode_mask_L_HD])
	
	single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
	s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
end

# ╔═╡ af011205-0ba2-4e37-bced-a8fe9fe88c2f
begin
	alg_L_HD = Integrated(arr[mask_L_HD_3D])
	vol_l_hd = score(s_bkg_L_HD, s_obj_L_HD, pixel_size, alg_L_HD)
end

# ╔═╡ 6667aa1c-c662-4243-8813-dd1ba595bfe8
md"""
### Calculate Mass
"""

# ╔═╡ 0db23983-2ff1-45c0-878b-e6837bed77d6
begin
	ρ_L_HD = density(mean(arr[dilated_mask_L_HD])) / 1e6 # g/cm^3 => g/mm^3
	mass_l_hd = score(s_bkg_L_HD, s_obj_L_HD, pixel_size, ρ_L_HD, alg_L_HD)
end

# ╔═╡ 01a6cf33-0209-4573-8469-d829e1d3f2bd
mass_l_hd_mg = mass_l_hd * 10^3

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

# ╔═╡ 66beab54-960e-4b4c-9295-aac336994686
md"""
#### Original mask
"""

# ╔═╡ 365ab95b-965c-4604-a25b-ff145ae9981a
@bind h1 overlay_mask_bind(mask_L_MD_3D)

# ╔═╡ 16407d15-90e3-4675-8920-669b5ea4c0c7
overlay_mask_plot(arr, mask_L_MD_3D, h1, "original mask")

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

# ╔═╡ e9cccfb2-38a7-4d2e-9d72-80c11ab65e1c
md"""
#### Eroded (object) mask
"""

# ╔═╡ 63f3d8f7-def9-4ee0-9533-ebcb6b17c19c
eroded_mask_L_MD = erode(erode(mask_L_MD_3D));

# ╔═╡ 990b7390-3d09-4b23-8f56-d02ecbbd4657
@bind h3 overlay_mask_bind(eroded_mask_L_MD)

# ╔═╡ 22f86ed8-1d36-4335-9317-a3a241f54ea3
overlay_mask_plot(arr, eroded_mask_L_MD, h3, "eroded mask")

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
### Calculate Volume
"""

# ╔═╡ 5bc01566-11ac-443a-8528-f949f9697de8
begin
	single_erode_mask_L_MD = eroded_mask_L_MD[:, :, 3]
	s_obj_L_MD = mean(single_arr[single_erode_mask_L_MD])
	
	single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
	s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
end

# ╔═╡ 553473b7-9a59-490f-bf47-becbc3f06bd5
begin
	alg_L_MD = Integrated(arr[mask_L_MD_3D])
	vol_l_md = score(s_bkg_L_MD, s_obj_L_MD, pixel_size, alg_L_MD)
end

# ╔═╡ 9e19c72d-1fdc-4bb8-9a34-cd23aa6892e4
md"""
### Calculate mass
"""

# ╔═╡ 917bb053-88da-4446-a303-0a280ee64ffe
begin
	ρ_L_MD = density(mean(arr[dilated_mask_L_MD])) / 1e6 # g/cm^3 => g/mm^3
	mass_l_md = score(s_bkg_L_MD, s_obj_L_MD, pixel_size, ρ_L_MD, alg_L_MD)
end

# ╔═╡ 98ebe591-20a4-4974-9235-eb0b7ee282f5
mass_l_md_mg = mass_l_md * 10^3

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

# ╔═╡ e24cc386-8e54-4b98-b98f-bce410ebcb40
md"""
#### Original mask
"""

# ╔═╡ 4d9a365a-864e-4b16-b7ac-c47ddb1acf60
@bind i1 overlay_mask_bind(mask_L_LD_3D)

# ╔═╡ c20f313d-7160-49e1-88eb-82bf24d97cd5
overlay_mask_plot(arr, mask_L_LD_3D, i1, "original mask")

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

# ╔═╡ 78ae5dbf-a713-4815-89b5-691291623a50
md"""
#### Eroded (object) mask
"""

# ╔═╡ ac2618ee-027e-4118-84b6-852eec2bcfc4
eroded_mask_L_LD = erode(erode(mask_L_LD_3D));

# ╔═╡ 1aa3b586-adfa-4396-b06a-41f61586e078
@bind i3 overlay_mask_bind(eroded_mask_L_LD)

# ╔═╡ af5f278d-f344-4ea9-a1b0-f8e7132dd72f
overlay_mask_plot(arr, eroded_mask_L_LD, i3, "eroded mask")

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
### Calculate Volume
"""

# ╔═╡ 0d0c639e-71ff-479d-8163-961a10db47ab
begin
	single_erode_mask_L_LD = eroded_mask_L_LD[:, :, 3]
	s_obj_L_LD = mean(single_arr[single_erode_mask_L_LD])
	
	single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
	s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
end

# ╔═╡ c6fb0eca-176a-45b3-a1ad-7e8aaa4803da
begin
	alg_L_LD = Integrated(arr[mask_L_LD_3D])
	vol_l_ld = score(s_bkg_L_LD, s_obj_L_LD, pixel_size, alg_L_LD)
end

# ╔═╡ e1cae446-b4fc-4d98-99b6-6c16c20ee53f
md"""
### Calculate mass
"""

# ╔═╡ c3811bc1-8c11-4a26-a245-95f15c6a5f57
begin
	ρ_L_LD = density(mean(arr[dilated_mask_L_LD])) / 1e6 # g/cm^3 => g/mm^3
	mass_l_ld = score(s_bkg_L_LD, s_obj_L_LD, pixel_size, ρ_L_LD, alg_L_LD)
end

# ╔═╡ e2f3fd61-3499-40ac-ad3b-95a486e3843e
mass_l_ld_mg = mass_l_ld * 10^3

# ╔═╡ Cell order:
# ╠═c50510a7-e981-4933-a4b1-fab2efcc43a2
# ╠═97a264e6-1037-4862-b6bd-2d2486bd419e
# ╠═84926791-1a72-4b9a-add4-7707efa16a08
# ╠═448bb1f2-6548-4bc8-bc7b-99c093833e7e
# ╠═2df42240-1fe7-47bb-baee-3150a6408393
# ╠═c515d4f4-38b1-4490-8201-774bfec91440
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
# ╠═78793039-df7f-4dbf-9188-155ecc72d944
# ╟─69d3e9a8-68c9-4a5a-8a4a-86c5fd64a08a
# ╠═0508d26f-09a1-457b-8cb0-5693f5c9ba3f
# ╠═544ae005-6661-487f-8662-ea3006ec4069
# ╠═9d0a4171-eb47-4cbf-8f72-4adf041a11d1
# ╠═2647d239-58d6-4448-9e3f-66918a10dbad
# ╟─2afdc85a-a8bf-4fd4-bdef-74c39c0fb452
# ╠═c2394467-eb15-4610-9bfd-7b3f88a7fb00
# ╠═9b6e3a13-9044-4b7e-89e0-b552af52bd2b
# ╠═0f85bfc3-f0d8-4ae7-b5b4-f72a072a3d55
# ╠═b8892b4e-3be7-4f82-a8d9-3be96ca5a86d
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
# ╠═2f5f42e9-d96c-4ac7-9505-646f0393986c
# ╟─ed01b680-a24d-44b9-8ef7-58ff4ab5a5ba
# ╠═546334f3-b411-43f9-8d0f-061aed8f22ca
# ╠═8d9504f9-9fa8-435c-9d1b-1ae42f98aa13
# ╠═9097cbfa-cf27-411c-9b83-8efeda2cdee6
# ╠═3146fe0f-bfe6-4cc9-9782-0ad2a097e777
# ╠═2694e3e0-f428-470b-a70c-f3b110d68f37
# ╠═78b0b7fc-ad4a-4121-bd30-2fdadcd95ea0
# ╠═13f3098d-e6ee-4fac-a429-a588cc8a2bbd
# ╠═384ff1c4-91a2-44b8-ac63-756c813c8432
# ╟─0483eebd-8106-4444-9037-fb52057fb0b4
# ╠═99319b7b-639d-4207-852d-0bd54b0537a4
# ╟─bb6dd528-05b1-48c1-b1f6-f9ef2d2b6777
# ╟─44f53f62-6a13-45b3-a1e9-ceea1576cbf4
# ╟─777e64c7-f05d-4af9-9417-268c3ea9f5b4
# ╟─fc18a1ee-8042-44e1-bc4c-fb56e345d90d
# ╠═16057812-b925-4571-8d5b-35254c42683c
# ╟─53d5b9ae-85b5-4f37-93ca-c76edd0d7916
# ╟─0a6451a3-19d6-4f4a-9e5b-cc2eb52b5b9d
# ╟─6fb82cc1-94c8-4a68-af6b-d69201a0fc6d
# ╠═e999b3f6-7d6f-4945-81fd-6b74216f0eec
# ╟─0608e3a1-8a56-409c-8258-e2919c51be5a
# ╠═1d2fd28a-25f1-4c81-94a4-70317d8172aa
# ╟─5bfd1bf6-e24c-42e3-b911-9d87f7c3491c
# ╠═588a129a-fe05-485f-8441-cd99dac82a03
# ╟─11b5faad-c162-4cf7-9b6f-2d61750b7a1a
# ╟─975f5f5c-3ada-4ff2-9e43-974f53d1b329
# ╟─820b0f29-867b-448f-9fb9-9ff1c3474422
# ╠═7532aa8e-7dbf-46ad-86c2-7f9ac7f2feca
# ╠═af011205-0ba2-4e37-bced-a8fe9fe88c2f
# ╟─6667aa1c-c662-4243-8813-dd1ba595bfe8
# ╠═0db23983-2ff1-45c0-878b-e6837bed77d6
# ╠═01a6cf33-0209-4573-8469-d829e1d3f2bd
# ╟─02be0555-01e0-45d5-9711-cb203ef41b4c
# ╠═04f317d4-5994-488c-92f2-63cd2582d946
# ╟─66beab54-960e-4b4c-9295-aac336994686
# ╟─365ab95b-965c-4604-a25b-ff145ae9981a
# ╠═16407d15-90e3-4675-8920-669b5ea4c0c7
# ╟─93192cc0-b222-4506-b3f4-dc12f0cff32f
# ╠═265b30ee-7ad8-4cc2-851f-7374e1d28da8
# ╟─04dfe430-4cf4-4cca-8fca-279fd00b286f
# ╟─184c7b7f-6a17-4a8f-9751-4a2f8d1878c1
# ╟─e9cccfb2-38a7-4d2e-9d72-80c11ab65e1c
# ╠═63f3d8f7-def9-4ee0-9533-ebcb6b17c19c
# ╟─990b7390-3d09-4b23-8f56-d02ecbbd4657
# ╟─22f86ed8-1d36-4335-9317-a3a241f54ea3
# ╟─fbb2fd43-4e42-4090-a4e6-1975d3957026
# ╠═c08d31ac-2fbd-43d6-b334-5157ed18b51b
# ╟─3fc579af-153a-41ff-9ab4-a78df87764e0
# ╟─66e4f9f7-0f35-479c-851c-93318cd7715a
# ╟─bb1950cc-0235-40f8-b1f5-ba4a1ae48e4b
# ╠═5bc01566-11ac-443a-8528-f949f9697de8
# ╠═553473b7-9a59-490f-bf47-becbc3f06bd5
# ╟─9e19c72d-1fdc-4bb8-9a34-cd23aa6892e4
# ╠═917bb053-88da-4446-a303-0a280ee64ffe
# ╠═98ebe591-20a4-4974-9235-eb0b7ee282f5
# ╟─686d7820-e823-4d42-a4e1-5e463411b3e1
# ╠═6113970c-9776-4cbc-9a3e-92cab30b5039
# ╟─e24cc386-8e54-4b98-b98f-bce410ebcb40
# ╟─4d9a365a-864e-4b16-b7ac-c47ddb1acf60
# ╟─c20f313d-7160-49e1-88eb-82bf24d97cd5
# ╟─fe7df65d-c15d-4764-83c7-b18f38c4b35c
# ╠═8b5fa573-c221-458b-9bd3-c29fabb902c4
# ╟─7b4df8a4-43bb-490a-9da4-0ac229bc21b4
# ╟─aa70bb01-446f-4ac2-b460-a8727b977bb3
# ╟─78ae5dbf-a713-4815-89b5-691291623a50
# ╠═ac2618ee-027e-4118-84b6-852eec2bcfc4
# ╟─1aa3b586-adfa-4396-b06a-41f61586e078
# ╟─af5f278d-f344-4ea9-a1b0-f8e7132dd72f
# ╟─612e2381-330b-475f-8614-07c0559c365a
# ╠═41788e2f-c0c2-4404-875f-576646f97a27
# ╟─0f2444f8-1a8c-4ab6-abbc-5329aec8cf97
# ╟─e6795bfb-cd7c-49c4-8e25-6eed0c26a476
# ╟─76c230bc-cdf4-4d9f-968f-60700e819de4
# ╠═0d0c639e-71ff-479d-8163-961a10db47ab
# ╠═c6fb0eca-176a-45b3-a1ad-7e8aaa4803da
# ╟─e1cae446-b4fc-4d98-99b6-6c16c20ee53f
# ╠═c3811bc1-8c11-4a26-a245-95f15c6a5f57
# ╠═e2f3fd61-3499-40ac-ad3b-95a486e3843e

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

# ╔═╡ d693d25a-40c8-443d-8cab-264b611fdbb8
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

# ╔═╡ d04784c6-71b9-4c9b-b7a9-267557518ce4
TableOfContents()

# ╔═╡ 10f3154b-79be-4a0f-b1aa-77d751b15f77
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDER`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ 9fcf4eca-1593-494e-95f2-e992c9107ebe
begin
	SCAN_NUMBER = 1
	VENDER = "Canon_Aquilion_One_Vision"
	BASE_PATH = "/Users/daleblack/Google Drive/Datasets/"
end

# ╔═╡ f03483ff-7567-4147-9ec3-722eb93ee498
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ e528161c-e925-4891-9c92-66fcdcda7b3e
root_path = string(BASE_PATH, VENDER)

# ╔═╡ 8103fab5-f51b-4aed-b877-84764adb58ad
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ 879b7ea4-91a7-4423-8e14-14b383c4f1a6
pth = dcm_path_list[SCAN_NUMBER]

# ╔═╡ 8ab76bd2-33dc-4834-bc82-b69f9a4bdec3
scan = basename(pth)

# ╔═╡ 042c7d6f-9427-4dfa-8712-f7d07f2ffd8c
header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ 852c5afc-fdae-4b6f-9b5e-706cb79483ed
md"""
## Helper Functions
"""

# ╔═╡ bc91a832-2c72-4ef5-a0c3-90a9069fb1ea
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ abb3fbbe-b7e1-450e-a020-7dcde6fcabd4
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ 132a7e0f-1748-4178-9db4-4b42f900b715
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

# ╔═╡ f97db31e-7308-47be-9f78-2b745b95e5bd
md"""
## Segment Heart
"""

# ╔═╡ f655ef59-87a9-43af-8938-a760060ee2a9
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ 13a28169-e4a8-4e10-85aa-4f910f8b3e08
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 5752445b-1311-4ada-a78f-78110bb3a1d9
heatmap(masked_array[:, :, a], colormap=:grays)

# ╔═╡ b27f97a9-3cbb-4db7-9c1b-e1d2f28135ef
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 15]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig
end

# ╔═╡ 0789b8c8-1472-4c8e-9458-d28c48217074
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig2
end

# ╔═╡ 55e8902c-44b1-477f-9a5f-d6c8065b9a95
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 23]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig3
end

# ╔═╡ 19e46b98-c5f3-4372-a7bc-63e00e739573
md"""
## Segment Calcium Rod
"""

# ╔═╡ 90463a54-7ae1-44e4-a55c-4406fa359206
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header);

# ╔═╡ 65edf7a9-ee57-4cf5-aa08-03a9dc2f504c
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=cal_rod_slice, show_value=true)

# ╔═╡ c5bba353-1293-46ae-8bdd-f3dbb3437934
heatmap(transpose(calcium_image[:, :, c]), colormap=:grays)

# ╔═╡ 2c257af8-8c92-4a9b-b135-fa86c57fc0a1
md"""
## Segment Calcium Inserts
"""

# ╔═╡ 8e31710e-a3b6-4a8f-a00c-8966900fca15
mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
            dcm_array, masked_array, header, slice_CCI, center_insert
);

# ╔═╡ c22d1c6f-47e4-4797-9e88-b89b84faf386
masks = mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD;

# ╔═╡ 3cead08a-a32b-4c19-bc71-ea4b88f59e93
heatmap(masks, colormap=:grays)

# ╔═╡ 5010860c-112c-4ae5-8044-cf95e207b9c5
md"""
## Calibration Prep
"""

# ╔═╡ 4911c60b-928e-4f31-9bed-374a238cbc38
md"""
### Calibration Insert
"""

# ╔═╡ 66429a9e-86a8-42b1-b042-9aa7bb159e40
array_filtered = abs.(mapwindow(median, calcium_image[:, :, cal_rod_slice], (3, 3)));

# ╔═╡ 7b36195d-0579-43bb-9305-9b62164338bf
bool_arr = array_filtered .> 0;

# ╔═╡ 63f9e598-a8b1-42d2-a958-d369d2bc77ac
bool_arr_erode = (((erode(erode(bool_arr)))));

# ╔═╡ a4bdaca0-ad29-4845-9879-e72462d2204d
heatmap(bool_arr, colormap=:grays)

# ╔═╡ 28c810fc-012a-423e-8e5b-613b2f8ee664
heatmap(bool_arr_erode, colormap=:grays)

# ╔═╡ b7c82ec4-9f9e-47e6-a0ae-7d0cf526bc8f
c_img = calcium_image[:, :, cal_rod_slice-1:cal_rod_slice+1];

# ╔═╡ c1afe8a1-24d9-477e-8f98-126c789b2801
begin
	mask_cal_3D = Array{Bool}(undef, size(c_img))
	for z in 1:size(c_img, 3)
		mask_cal_3D[:, :, z] = bool_arr_erode
	end
end;

# ╔═╡ 4e654494-bbdc-45af-82fd-d26bcb8d106b
cal_insert_mean = mean(c_img[mask_cal_3D])

# ╔═╡ ebdb78c0-3b17-42b8-80cd-6eff134805b0
md"""
### Calibration Line
"""

# ╔═╡ e24af383-13dc-4cd4-92a2-2e35b44e0910
density_array = [0, 200, 400, 800]

# ╔═╡ 57d3fb05-2516-4e9d-9acd-6263dd572f02
density_array_calc = [0, 200] # mg/cc

# ╔═╡ 0220131d-cd4c-4f4a-aa02-726685da39f0
intensity_array = [0, cal_insert_mean] # HU

# ╔═╡ 4f032201-37c5-4823-b314-cc0e8761e196
df = DataFrame(:density => density_array_calc, :intensity => intensity_array)

# ╔═╡ a0954527-e227-40aa-a7ae-5bd2e1a44f48
linearRegressor = lm(@formula(intensity ~ density), df);

# ╔═╡ fcb1b8de-cba2-4989-b71f-5e1fd1591297
linearFit = predict(linearRegressor)

# ╔═╡ a8590648-5378-48dc-b9bd-c46d77e27aa6
m = linearRegressor.model.pp.beta0[2]

# ╔═╡ 8f360e64-11d8-420e-bbb2-a55c06650af1
b = linearRegressor.model.rr.mu[1]

# ╔═╡ a8773915-c2be-4127-b1db-299be8781f7c
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

# ╔═╡ eba83ef8-efa8-4cad-b6dc-aebdca73bee2
density(intensity) = (intensity - b) / m

# ╔═╡ dea725bd-de75-44b8-84b3-f393f66adc5c
intensity(ρ) = m*ρ + b

# ╔═╡ f0c50fef-3bfe-430e-b060-9981b8d417d7
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

# ╔═╡ 5710f522-d980-473b-9c1d-90c5c9ead5af
md"""
# Score Large Inserts
"""

# ╔═╡ b33e0312-63d8-44d2-9ca4-bd07192dcff5
arr = masked_array[:, :, slice_CCI-2:slice_CCI+2];

# ╔═╡ 1b0f7d16-f1cf-4a5f-80c0-4fa3ee8e3e36
single_arr = masked_array[:, :, slice_CCI];

# ╔═╡ 8b95766c-4f50-4923-8384-9aea67860d24
md"""
## High Density
"""

# ╔═╡ f44689fc-41c1-49a8-9198-bf4f46a7cee7
begin
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_HD_3D[:, :, z] = mask_L_HD
	end
end;

# ╔═╡ 17211561-e3c6-4b13-8c4d-75f8a05b284b
md"""
#### Dilated mask
"""

# ╔═╡ 871f1466-0be9-4d7f-8797-8ec54cc8cb7a
dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D));

# ╔═╡ 23963e27-0957-4751-9ca2-b3817f7a485a
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ d1704206-df85-4fdb-8ba0-8cf5784bb1b3
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ c91ea7f7-ffd5-41a1-80f1-952daa00e2d0
md"""
#### Ring (background) mask
"""

# ╔═╡ 3fff8462-1e59-4692-a107-4cfeaa3d435a
ring_mask_L_HD = dilate(dilate(dilate(dilate(mask_L_HD_3D)))) - dilate(dilate(dilate(mask_L_HD_3D)));

# ╔═╡ 483b5682-d35f-4f82-8b45-c679fab9bba9
@bind g4 overlay_mask_bind(ring_mask_L_HD)

# ╔═╡ 9dfa799a-98e3-4c57-9317-17c21fac0f8b
overlay_mask_plot(arr, ring_mask_L_HD, g4, "ring mask")

# ╔═╡ 74b3220b-7b75-4a45-b76f-ede69567d942
md"""
### Calculations
"""

# ╔═╡ 20253e73-f846-4094-9196-8509e7aa207f
begin
	single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
	s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
end

# ╔═╡ 5e7e2d95-59fb-467b-946f-d6a0917ef971
S_Obj_HD = intensity(800)

# ╔═╡ 20d94be3-5244-4762-8f6a-9b7a397c7ef6
pixel_size = DICOMUtils.get_pixel_size(header)

# ╔═╡ c0ffcf34-0ac4-49aa-9687-4a867c3683d4
begin
	alg_L_HD = Integrated(arr[mask_L_HD_3D])
	vol_l_hd = score(s_bkg_L_HD, S_Obj_HD, pixel_size, alg_L_HD)
end

# ╔═╡ 64ecc3de-c80c-4aa3-b8c6-8f539c40936d
begin
	ρ_HD = 0.8 # mg/mm^3
	mass_l_hd = score(s_bkg_L_HD, S_Obj_HD, pixel_size, ρ_HD, alg_L_HD)
end

# ╔═╡ 32f5c372-f092-4e2c-ae9a-e1127ed36209
md"""
## Medium Density
"""

# ╔═╡ efe7eadd-0be4-46c7-8e05-1cf0c78b8c39
begin
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_MD_3D[:, :, z] = mask_L_MD
	end
end;

# ╔═╡ ff0c4b5e-201b-495b-951a-fb283fa8df02
md"""
#### Dilated mask
"""

# ╔═╡ 35d5e53a-e543-4288-a5bc-d41bb0b0dd85
dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D));

# ╔═╡ 550c3434-bb91-4721-bfa4-b3f7fddab95f
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ f80dfb77-0a01-42a8-9d8e-b4f5b78a10bb
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ e6d277f1-7cb8-4887-829a-a991cbfa32bc
md"""
#### Ring (background) mask
"""

# ╔═╡ 7ae02410-b40e-414a-ad44-e1d62b48281f
ring_mask_L_MD = dilate(dilate(dilate(dilate(mask_L_MD_3D)))) - dilate(dilate(dilate(mask_L_MD_3D)));

# ╔═╡ 10adebdd-a69d-4bac-8521-38b3d70a1ca8
@bind h4 overlay_mask_bind(ring_mask_L_MD)

# ╔═╡ ef793dbc-63fc-4530-8456-2c4e35726fb7
overlay_mask_plot(arr, ring_mask_L_MD, h4, "ring mask")

# ╔═╡ d4434bfd-4c33-4cdb-84fe-2c21b1dfccb9
md"""
### Calculations
"""

# ╔═╡ 3d3ecbef-c5fb-487c-ba29-d6d835cebba5
begin
	single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
	s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
end

# ╔═╡ f135c026-17d6-4b08-83d6-fdc5c4628b61
S_Obj_MD = intensity(400)

# ╔═╡ ab1fcaa7-6296-4a05-b343-c39443967a5b
begin
	alg_L_MD = Integrated(arr[mask_L_MD_3D])
	vol_l_md = score(s_bkg_L_MD, S_Obj_MD, pixel_size, alg_L_MD)
end

# ╔═╡ aeb0e00c-5196-4a80-a299-38f9218b7d59
begin
	ρ_MD = 0.4 # mg/mm^3
	mass_l_md = score(s_bkg_L_MD, S_Obj_MD, pixel_size, ρ_MD, alg_L_MD)
end

# ╔═╡ 39e51df8-5e9d-4eea-bb48-c0939209b95e
md"""
## Low Density
"""

# ╔═╡ a4a68ce0-ad10-41ef-918b-092fc6e381f4
begin
	mask_L_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_LD_3D[:, :, z] = mask_L_LD
	end
end;

# ╔═╡ 3240bd88-513e-43c8-8068-9b0ce88f64d0
md"""
#### Dilated mask
"""

# ╔═╡ be02edba-fe11-4c79-a783-d7a22ee79043
dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D));

# ╔═╡ 288ad98d-5472-4e59-b5d3-45c9fbdc1a47
@bind i2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ b9d8d488-8b33-4941-967a-0cf9c019c91b
overlay_mask_plot(arr, dilated_mask_L_LD, i2, "dilated mask")

# ╔═╡ a6d95216-8d80-4101-af62-1a1b1f1d936c
md"""
#### Ring (background) mask
"""

# ╔═╡ f30d2b4d-23df-41c2-aaf1-28878fd28a40
ring_mask_L_LD = dilate(dilate(dilate(dilate(mask_L_LD_3D)))) - dilate(dilate(dilate(mask_L_LD_3D)));

# ╔═╡ 0129e1a3-3cc2-4e61-a325-1abf356b9718
@bind i4 overlay_mask_bind(ring_mask_L_LD)

# ╔═╡ 309ca624-fefe-4bcb-81e3-e5f4ebb9f773
overlay_mask_plot(arr, ring_mask_L_LD, i4, "ring mask")

# ╔═╡ 7c969f9e-dc47-4478-910d-a91ce541ec17
md"""
### Calculations
"""

# ╔═╡ 536a1a1f-43f2-4279-979e-a510f1e97223
begin	
	single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
	s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
end

# ╔═╡ 1b895ec8-ed17-437b-a355-81d8039c4393
S_Obj_LD = intensity(200)

# ╔═╡ 7b5234eb-e4d5-410a-9740-20c23c62baa6
begin
	alg_L_LD = Integrated(arr[mask_L_LD_3D])
	vol_l_ld = score(s_bkg_L_LD, S_Obj_LD, pixel_size, alg_L_LD)
end

# ╔═╡ 1054ff92-52c9-408e-bc8d-194585288052
begin
	ρ_LD = 0.2
	mass_l_ld = score(s_bkg_L_LD, S_Obj_LD, pixel_size, ρ_LD, alg_L_LD)
end

# ╔═╡ 57bf1a8f-6f62-4d37-991f-83843e2c7e73
md"""
# Score Medium Inserts
"""

# ╔═╡ 68f7258e-f6b1-44e3-82fa-606f65bf5c0f
md"""
## High Density
"""

# ╔═╡ 8b5b6433-cd38-4c18-8176-110e2a33965c
begin
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_HD_3D[:, :, z] = mask_M_HD
	end
end;

# ╔═╡ 28ec799f-fcfd-40ea-bddb-0c1450d14d39
md"""
#### Dilated mask
"""

# ╔═╡ 1e341397-f0c4-4a92-8b12-669f59380f9e
dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 17175b11-5ddc-4e78-b21a-b817049be2ca
@bind j2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ e9f0f552-c88e-44f5-afbc-88672cb34093
overlay_mask_plot(arr, dilated_mask_M_HD, j2, "dilated mask")

# ╔═╡ b4e30b34-ef4c-4cf4-86a2-5e4a81a41705
md"""
#### Ring (background) mask
"""

# ╔═╡ cf195f98-70b4-4ba2-bab6-d27e925070e6
ring_mask_M_HD = dilate(dilate(dilate(dilate(dilate(mask_M_HD_3D))))) - dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 0301bf1e-acd5-4d1e-bfbe-0f6813d72357
@bind j4 overlay_mask_bind(ring_mask_M_HD)

# ╔═╡ bf6c5954-bcef-4b38-bf6d-2850c81283e6
overlay_mask_plot(arr, ring_mask_M_HD, j4, "ring mask")

# ╔═╡ f6dc3d73-6bd7-4224-a94b-ea195f63ce2c
md"""
### Calculations
"""

# ╔═╡ 6ba5d457-6bfc-4f9d-8a83-1e6c6f033a49
begin
	single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
	s_bkg_M_HD = mean(single_arr[single_ring_mask_M_HD])
end

# ╔═╡ 0d467ffd-8cc8-4c3b-8729-e08bce2dcf5d
begin
	alg_M_HD = Integrated(arr[mask_M_HD_3D])
	vol_m_hd = score(s_bkg_M_HD, S_Obj_HD, pixel_size, alg_M_HD)
end

# ╔═╡ 5c1ece3e-4cea-4107-9cea-0c23aaf0bbe1
mass_m_hd = score(s_bkg_M_HD, S_Obj_HD, pixel_size, ρ_HD, alg_M_HD)

# ╔═╡ d2c408bc-6b14-4fca-a277-da342f0b68da
md"""
## Medium Density
"""

# ╔═╡ 79425b93-a29e-40f3-980f-0e5774f75a6b
begin
	mask_M_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_MD_3D[:, :, z] = mask_M_MD
	end
end;

# ╔═╡ 0478b474-d354-4c62-a56a-f868d70a042b
md"""
#### Dilated mask
"""

# ╔═╡ 93c337b6-3398-4998-8cde-c12a411e23df
dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))));

# ╔═╡ 5c7df635-9388-421d-8929-f0f2fa3f4862
@bind k2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ 8d987296-7c94-44a9-ac43-bba5d455bbc9
overlay_mask_plot(arr, dilated_mask_M_MD, k2, "dilated mask")

# ╔═╡ b7b2917f-86af-4591-b298-64941ec91665
md"""
#### Ring (background) mask
"""

# ╔═╡ 63c468be-281a-4cfd-a546-260239c89fbd
ring_mask_M_MD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))));

# ╔═╡ 7cfbb85a-49f2-4826-a814-7eb9f5696ce4
@bind k4 overlay_mask_bind(ring_mask_M_MD)

# ╔═╡ 20d9dd5d-d4f1-4c04-99c4-9d18fff34d04
overlay_mask_plot(arr, ring_mask_M_MD, k4, "ring mask")

# ╔═╡ edeea6af-ad4a-430d-97d7-5e62d7c68aab
md"""
### Calculations
"""

# ╔═╡ 8dc152a0-00d7-46e7-8368-33c2a36da2f6
begin
	single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
	s_bkg_M_MD = mean(single_arr[single_ring_mask_M_MD])
end

# ╔═╡ 5ffeba59-c088-45ee-9af8-0a593eb32dce
begin
	alg_M_MD = Integrated(arr[mask_M_MD_3D])
	vol_m_md = score(s_bkg_M_MD, S_Obj_MD, pixel_size, alg_M_MD)
end

# ╔═╡ 37126474-6677-45b0-a425-50031b592e8b
mass_m_md = score(s_bkg_M_MD, S_Obj_MD, pixel_size, ρ_MD, alg_M_MD)

# ╔═╡ f7e87875-b9fc-49cb-8b81-2a75d1be11fb
md"""
## Low Density
"""

# ╔═╡ d1ef2035-7854-44f3-9286-dc5b5b36097a
begin
	mask_M_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_LD_3D[:, :, z] = mask_M_LD
	end
end;

# ╔═╡ a345439c-83d7-4c45-8e98-77db1b36f279
md"""
#### Dilated mask
"""

# ╔═╡ 339600a4-236e-46b3-9612-f30ee7f6ba78
dilated_mask_M_LD = dilate(dilate(dilate(dilate(mask_M_LD_3D))));

# ╔═╡ fccb3dba-69ba-4782-9181-65ac303813b2
@bind l2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ a9dd7902-1871-4ea8-8f58-0e50adf82f0b
overlay_mask_plot(arr, dilated_mask_M_LD, l2, "dilated mask")

# ╔═╡ 297bf17a-1a57-4422-b15d-599eef64b125
md"""
#### Ring (background) mask
"""

# ╔═╡ c4b622c7-b037-4839-a61e-759f43dafeee
ring_mask_M_LD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ e8a8f9cf-7359-4131-b307-c8b55bf73546
@bind l4 overlay_mask_bind(ring_mask_M_LD)

# ╔═╡ 35e417b7-a71e-4827-b154-0419888952e2
overlay_mask_plot(arr, ring_mask_M_LD, l4, "ring mask")

# ╔═╡ 81c064f4-cff9-4ed7-b2a8-a620cf15ccaf
md"""
### Calculations
"""

# ╔═╡ 3e62017a-90e3-4948-8d16-effe14e5f5af
begin
	single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
	s_bkg_M_LD = mean(single_arr[single_ring_mask_M_LD])
end

# ╔═╡ dcb8d1b7-5e9e-4078-b86a-a8c572a870a9
begin
	alg_M_LD = Integrated(arr[mask_M_LD_3D])
	vol_m_ld = score(s_bkg_M_LD, S_Obj_LD, pixel_size, alg_M_LD)
end

# ╔═╡ 1f579a15-3efc-4bc0-ba4c-dcc11dcb4c63
mass_m_ld = score(s_bkg_M_LD, S_Obj_LD, pixel_size, ρ_LD, alg_M_LD)

# ╔═╡ e8ec60e9-5ee6-4085-8c1e-907aa2d4624a
md"""
# Score Small Inserts
"""

# ╔═╡ b6525a66-35ae-4a74-937d-db51fa7a3024
md"""
## High Density
"""

# ╔═╡ e3280173-b6ac-4cc1-945b-ec1c394c9260
begin
	mask_S_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_HD_3D[:, :, z] = mask_S_HD
	end
end;

# ╔═╡ 559669c9-8539-4d83-ac0b-3135b127028b
md"""
#### Dilated mask
"""

# ╔═╡ b883f972-ce50-4400-8bc5-edce46a6689d
dilated_mask_S_HD = dilate(dilate(dilate(dilate(mask_S_HD_3D))));

# ╔═╡ f0e695a3-eb78-47ba-8e4a-6af332e1b496
@bind m2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ 069d2662-24c6-49af-be33-8d269b05bc9b
overlay_mask_plot(arr, dilated_mask_S_HD, m2, "dilated mask")

# ╔═╡ bbaff240-af5c-4621-97ba-613040420f6e
md"""
#### Ring (background) mask
"""

# ╔═╡ da33cf36-6c3d-4c94-a1f1-59a71ae06acb
ring_mask_S_HD = dilate(dilate(dilate(dilate(dilate(mask_S_HD_3D))))) - dilate(dilate(dilate(dilate(mask_S_HD_3D))));

# ╔═╡ a48fc90b-f99d-4408-be54-d170977e394f
@bind m4 overlay_mask_bind(ring_mask_S_HD)

# ╔═╡ 2563b25a-a386-4f25-b88e-d5db49a1473a
overlay_mask_plot(arr, ring_mask_S_HD, m4, "ring mask")

# ╔═╡ 6fb7e314-c5cc-4c90-90c0-9b9db04cfe1b
md"""
### Calculations
"""

# ╔═╡ fe263281-fdd2-4025-a75c-11701059ea78
begin
	single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
	s_bkg_S_HD = mean(single_arr[single_ring_mask_S_HD])
end

# ╔═╡ 2ee67d3b-b2a7-4605-9bbb-06caa5c888e2
begin
	alg_S_HD = Integrated(arr[mask_S_HD_3D])
	vol_s_hd = score(s_bkg_S_HD, S_Obj_HD, pixel_size, alg_S_HD)
	if vol_s_hd < 0
		vol_s_hd = 0
	end
end

# ╔═╡ b040556f-5afe-4539-be19-ef5eaecb8448
mass_s_hd = score(s_bkg_S_HD, S_Obj_HD, pixel_size, ρ_HD, alg_S_HD)

# ╔═╡ ca315565-2e7f-435a-be3f-05d4893cfcce
md"""
## Medium Density
"""

# ╔═╡ 164dd82f-47bd-4255-9507-f0386cd94b40
begin
	mask_S_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_MD_3D[:, :, z] = mask_S_MD
	end
end;

# ╔═╡ 35fa87bb-ce0e-4505-b667-6de69e3641ed
md"""
#### Dilated mask
"""

# ╔═╡ 48ce3652-c8ef-4fcc-b284-1cb8bb720da8
dilated_mask_S_MD = dilate(dilate(dilate(dilate(mask_S_MD_3D))));

# ╔═╡ f1a0c8a1-b79d-435f-907d-64b5671e65c2
@bind n2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ 7e7d92bb-d74a-470a-a0f2-e48b62ae56e0
overlay_mask_plot(arr, dilated_mask_S_MD, n2, "dilated mask")

# ╔═╡ e3c667a9-b3ff-494f-943e-8c797f31bd91
md"""
#### Ring (background) mask
"""

# ╔═╡ 4c7c4aeb-2f16-4aa8-acc9-1e9be473b57f
ring_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D))))) - dilate(dilate(dilate(dilate(mask_S_MD_3D))));

# ╔═╡ b6c84ed9-573b-472b-be9b-fe4b91bbf40b
@bind n4 overlay_mask_bind(ring_mask_S_MD)

# ╔═╡ f1f32ea8-2a17-478c-9ca1-c5ce7afe3912
overlay_mask_plot(arr, ring_mask_S_MD, n4, "ring mask")

# ╔═╡ b315f712-7630-4f5f-af85-bdfc46d4d5f7
md"""
### Calculations
"""

# ╔═╡ 47e47506-08de-4991-90fa-35e3e6b18c65
begin
	single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
	s_bkg_S_MD = mean(single_arr[single_ring_mask_S_MD])
end

# ╔═╡ c8756384-6c90-4819-b0dc-2798eca244a8
begin
	alg_S_MD = Integrated(arr[mask_S_MD_3D])
	vol_s_md = score(s_bkg_S_MD, S_Obj_HD, pixel_size, alg_S_MD)
	if vol_s_md < 0
		vol_s_md = 0
	end
end

# ╔═╡ 55a9dfbf-6f8a-40ab-9e23-5fe5231a250f
mass_s_md = score(s_bkg_S_MD, S_Obj_HD, pixel_size, ρ_MD, alg_S_MD)

# ╔═╡ c10d0942-ba67-4c71-ae53-1e2ed5020dae
md"""
## Low Density
"""

# ╔═╡ c384cf03-98fb-4c5e-b426-487b688e1e69
begin
	mask_S_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_LD_3D[:, :, z] = mask_S_LD
	end
end;

# ╔═╡ 23189d22-61e7-4c30-b7ac-50d20692ecd5
md"""
#### Dilated mask
"""

# ╔═╡ 75b50ca5-f1f5-494a-85f1-33881db475b3
dilated_mask_S_LD = dilate(dilate(dilate(dilate(mask_S_LD_3D))));

# ╔═╡ 4233274a-f1ec-4cd9-b68f-c68f2f10dfa0
@bind o2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ 976fc996-2620-45a7-aacc-da7da24bfc0b
overlay_mask_plot(arr, dilated_mask_S_LD, o2, "dilated mask")

# ╔═╡ 76e2b877-4ebf-4550-b796-52dfbf55d361
md"""
#### Ring (background) mask
"""

# ╔═╡ 3b7ddaa4-4e99-4710-b2fe-26be6288edbc
ring_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D))))) - dilate(dilate(dilate(dilate(mask_S_LD_3D))));

# ╔═╡ eb3bbceb-2198-438a-9f25-18e074b1423e
@bind o4 overlay_mask_bind(ring_mask_S_LD)

# ╔═╡ 0fbac767-2689-435e-a123-a3ea84f0a294
overlay_mask_plot(arr, ring_mask_S_LD, o4, "ring mask")

# ╔═╡ 8b4ba108-bea3-4850-a1fc-0fad0cd7a6ea
md"""
### Calculations
"""

# ╔═╡ 4a6bd3b1-53ed-48f4-916f-e1665aab1664
begin
	single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
	s_bkg_S_LD = mean(single_arr[single_ring_mask_S_LD])
end

# ╔═╡ ac73fe19-8ad0-44a7-bb95-95af48f170d5
begin
	alg_S_LD = Integrated(arr[mask_S_LD_3D])
	vol_s_ld = score(s_bkg_S_LD, S_Obj_HD, pixel_size, alg_S_LD)
	if vol_s_ld < 0
		vol_s_ld = 0
	end
end

# ╔═╡ 14e7d218-b103-413a-9ce3-51f051243c56
mass_s_ld = score(s_bkg_S_LD, S_Obj_HD, pixel_size, ρ_LD, alg_S_LD)

# ╔═╡ 4ba64fcd-b8d3-4045-8a47-86d1705d19ff
md"""
# Results
"""

# ╔═╡ fa08de3f-27b6-4ff3-917e-73ea692296b0
md"""
### Volume
"""

# ╔═╡ dc3ae6f7-8fa7-42ca-b2f2-75e42befa5ce
inserts = [
	"Low Density",
	"Medium Density",
	"High Density"
]

# ╔═╡ 8162eea8-e609-490c-8f8f-0a058a3c2d8b
ground_truth_volume_large = [
	98.2,
	98.2,
	98.2,
] # mm^3

# ╔═╡ 75530e55-95e8-498f-82e0-3e4d7bac8ed4
calculated_volume_large = [
	vol_l_ld,
	vol_l_md,
	vol_l_hd
]

# ╔═╡ 6a094720-d715-4e81-ad6a-e23efa13d545
ground_truth_volume_medium = [
	21.2,
	21.2,
	21.2
]

# ╔═╡ 514c21a8-2e4e-47ed-b6e5-53d67c049ebb
calculated_volume_medium = [
	vol_m_ld,
	vol_m_md,
	vol_m_hd
]

# ╔═╡ bd11ad26-f577-45cb-beea-dc9231cbf0a1
ground_truth_volume_small = [
	0.8,
	0.8,
	0.8
]

# ╔═╡ 75cea2c9-46f8-42af-a92f-5124318c21c2
calculated_volume_small = [
	vol_s_ld,
	vol_s_md,
	vol_s_hd
]

# ╔═╡ 62340699-9e33-4e45-a926-5c6d84e8b4ce
df1 = DataFrame(
	inserts = inserts,
	ground_truth_volume_large = ground_truth_volume_large,
	calculated_volume_large = calculated_volume_large,
	ground_truth_volume_medium = ground_truth_volume_medium,
	calculated_volume_medium = calculated_volume_medium,
	ground_truth_volume_small = ground_truth_volume_small,
	calculated_volume_small = calculated_volume_small
)

# ╔═╡ ba56da4a-e854-481f-9088-4be7c8be425d
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

# ╔═╡ 82d5fdb6-0734-4eff-b329-ac640c40bc3e
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

# ╔═╡ cf8232e0-bab8-4e9b-9c85-2b5df5bbe5be
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

# ╔═╡ d758c9db-39e1-42cb-9f3f-e8cb962c925f
md"""
### Mass
"""

# ╔═╡ e564c3e3-058a-49e9-9a79-68375ab4482b
ground_truth_mass_large = [
	19.6,
	39.3,
	78.5
] # mg

# ╔═╡ 1cf8f266-448e-4ece-bfed-41db8fb23526
calculated_mass_large = [
	mass_l_ld,
	mass_l_md,
	mass_l_hd
]

# ╔═╡ 22ffc892-6f90-4172-94f6-06d39237715d
ground_truth_mass_medium = [
	4.2,
	8.5,
	17.0
]

# ╔═╡ 2031a3c1-db8b-4a31-b1cb-3caf2b7537c3
calculated_mass_medium = [
	mass_m_ld,
	mass_m_md,
	mass_m_hd
]

# ╔═╡ 21935abe-702f-4bd3-8cee-dddad07d23e2
ground_truth_mass_small = [
	0.2,
	0.3,
	0.6
]

# ╔═╡ fdf4804b-af2a-42ba-9d12-1d4e0a9c4496
calculated_mass_small = [
	mass_s_ld,
	mass_s_md,
	mass_s_hd
]

# ╔═╡ 338ad5b0-7c19-4b14-a513-cdabb0ef05f4
df2 = DataFrame(
	inserts = inserts,
	ground_truth_mass_large = ground_truth_mass_large,
	calculated_mass_large = calculated_mass_large,
	ground_truth_mass_medium = ground_truth_mass_medium,
	calculated_mass_medium = calculated_mass_medium,
	ground_truth_mass_small = ground_truth_mass_small,
	calculated_mass_small = calculated_mass_small
)

# ╔═╡ 12bf244c-d9f8-456e-b818-81d024f451b2
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

# ╔═╡ 925a05e7-10b8-484d-a44c-9b2e532b63aa
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

# ╔═╡ d333096e-a985-4abb-9b32-a13dc82d824b
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

# ╔═╡ 812299b1-c58a-40c6-8d93-80becfcc008e
md"""
### Save Results
"""

# ╔═╡ c25e5b1d-053f-49b7-ad1f-f4c417aa491c
df_final = leftjoin(df1, df2, on=:inserts)

# ╔═╡ 5947da89-8010-48ac-a5b3-3e9371280288
if ~isdir(string(cd(pwd, "..") , "/data/output/", VENDER, "2"))
	mkdir(string(cd(pwd, "..") , "/data/output/", VENDER, "2"))
end

# ╔═╡ 443b7def-1fb5-463c-bd88-8975c6f59501
output_path = string(cd(pwd, "..") , "/data/output/", VENDER, "2", "/", scan, ".csv")

# ╔═╡ 67bd57ff-7f71-48c7-9e3b-92e6b7e5bc30
CSV.write(output_path, df_final)

# ╔═╡ Cell order:
# ╠═d693d25a-40c8-443d-8cab-264b611fdbb8
# ╠═d04784c6-71b9-4c9b-b7a9-267557518ce4
# ╟─10f3154b-79be-4a0f-b1aa-77d751b15f77
# ╠═9fcf4eca-1593-494e-95f2-e992c9107ebe
# ╟─f03483ff-7567-4147-9ec3-722eb93ee498
# ╠═e528161c-e925-4891-9c92-66fcdcda7b3e
# ╠═8103fab5-f51b-4aed-b877-84764adb58ad
# ╠═879b7ea4-91a7-4423-8e14-14b383c4f1a6
# ╠═8ab76bd2-33dc-4834-bc82-b69f9a4bdec3
# ╠═042c7d6f-9427-4dfa-8712-f7d07f2ffd8c
# ╟─852c5afc-fdae-4b6f-9b5e-706cb79483ed
# ╟─bc91a832-2c72-4ef5-a0c3-90a9069fb1ea
# ╟─abb3fbbe-b7e1-450e-a020-7dcde6fcabd4
# ╟─132a7e0f-1748-4178-9db4-4b42f900b715
# ╟─f97db31e-7308-47be-9f78-2b745b95e5bd
# ╠═f655ef59-87a9-43af-8938-a760060ee2a9
# ╠═13a28169-e4a8-4e10-85aa-4f910f8b3e08
# ╠═5752445b-1311-4ada-a78f-78110bb3a1d9
# ╠═b27f97a9-3cbb-4db7-9c1b-e1d2f28135ef
# ╠═0789b8c8-1472-4c8e-9458-d28c48217074
# ╠═55e8902c-44b1-477f-9a5f-d6c8065b9a95
# ╟─19e46b98-c5f3-4372-a7bc-63e00e739573
# ╠═90463a54-7ae1-44e4-a55c-4406fa359206
# ╠═65edf7a9-ee57-4cf5-aa08-03a9dc2f504c
# ╠═c5bba353-1293-46ae-8bdd-f3dbb3437934
# ╟─2c257af8-8c92-4a9b-b135-fa86c57fc0a1
# ╠═8e31710e-a3b6-4a8f-a00c-8966900fca15
# ╠═c22d1c6f-47e4-4797-9e88-b89b84faf386
# ╠═3cead08a-a32b-4c19-bc71-ea4b88f59e93
# ╟─5010860c-112c-4ae5-8044-cf95e207b9c5
# ╟─4911c60b-928e-4f31-9bed-374a238cbc38
# ╠═66429a9e-86a8-42b1-b042-9aa7bb159e40
# ╠═7b36195d-0579-43bb-9305-9b62164338bf
# ╠═63f9e598-a8b1-42d2-a958-d369d2bc77ac
# ╠═a4bdaca0-ad29-4845-9879-e72462d2204d
# ╠═28c810fc-012a-423e-8e5b-613b2f8ee664
# ╠═c1afe8a1-24d9-477e-8f98-126c789b2801
# ╠═b7c82ec4-9f9e-47e6-a0ae-7d0cf526bc8f
# ╠═4e654494-bbdc-45af-82fd-d26bcb8d106b
# ╟─ebdb78c0-3b17-42b8-80cd-6eff134805b0
# ╠═e24af383-13dc-4cd4-92a2-2e35b44e0910
# ╠═57d3fb05-2516-4e9d-9acd-6263dd572f02
# ╠═0220131d-cd4c-4f4a-aa02-726685da39f0
# ╠═4f032201-37c5-4823-b314-cc0e8761e196
# ╠═a0954527-e227-40aa-a7ae-5bd2e1a44f48
# ╠═fcb1b8de-cba2-4989-b71f-5e1fd1591297
# ╠═a8590648-5378-48dc-b9bd-c46d77e27aa6
# ╠═8f360e64-11d8-420e-bbb2-a55c06650af1
# ╟─a8773915-c2be-4127-b1db-299be8781f7c
# ╠═eba83ef8-efa8-4cad-b6dc-aebdca73bee2
# ╠═dea725bd-de75-44b8-84b3-f393f66adc5c
# ╠═f0c50fef-3bfe-430e-b060-9981b8d417d7
# ╟─5710f522-d980-473b-9c1d-90c5c9ead5af
# ╠═b33e0312-63d8-44d2-9ca4-bd07192dcff5
# ╠═1b0f7d16-f1cf-4a5f-80c0-4fa3ee8e3e36
# ╟─8b95766c-4f50-4923-8384-9aea67860d24
# ╠═f44689fc-41c1-49a8-9198-bf4f46a7cee7
# ╟─17211561-e3c6-4b13-8c4d-75f8a05b284b
# ╠═871f1466-0be9-4d7f-8797-8ec54cc8cb7a
# ╠═23963e27-0957-4751-9ca2-b3817f7a485a
# ╠═d1704206-df85-4fdb-8ba0-8cf5784bb1b3
# ╟─c91ea7f7-ffd5-41a1-80f1-952daa00e2d0
# ╠═3fff8462-1e59-4692-a107-4cfeaa3d435a
# ╠═483b5682-d35f-4f82-8b45-c679fab9bba9
# ╠═9dfa799a-98e3-4c57-9317-17c21fac0f8b
# ╟─74b3220b-7b75-4a45-b76f-ede69567d942
# ╠═20253e73-f846-4094-9196-8509e7aa207f
# ╠═5e7e2d95-59fb-467b-946f-d6a0917ef971
# ╠═20d94be3-5244-4762-8f6a-9b7a397c7ef6
# ╠═c0ffcf34-0ac4-49aa-9687-4a867c3683d4
# ╠═64ecc3de-c80c-4aa3-b8c6-8f539c40936d
# ╟─32f5c372-f092-4e2c-ae9a-e1127ed36209
# ╠═efe7eadd-0be4-46c7-8e05-1cf0c78b8c39
# ╟─ff0c4b5e-201b-495b-951a-fb283fa8df02
# ╠═35d5e53a-e543-4288-a5bc-d41bb0b0dd85
# ╠═550c3434-bb91-4721-bfa4-b3f7fddab95f
# ╠═f80dfb77-0a01-42a8-9d8e-b4f5b78a10bb
# ╠═e6d277f1-7cb8-4887-829a-a991cbfa32bc
# ╠═7ae02410-b40e-414a-ad44-e1d62b48281f
# ╠═10adebdd-a69d-4bac-8521-38b3d70a1ca8
# ╠═ef793dbc-63fc-4530-8456-2c4e35726fb7
# ╟─d4434bfd-4c33-4cdb-84fe-2c21b1dfccb9
# ╠═3d3ecbef-c5fb-487c-ba29-d6d835cebba5
# ╠═f135c026-17d6-4b08-83d6-fdc5c4628b61
# ╠═ab1fcaa7-6296-4a05-b343-c39443967a5b
# ╠═aeb0e00c-5196-4a80-a299-38f9218b7d59
# ╟─39e51df8-5e9d-4eea-bb48-c0939209b95e
# ╠═a4a68ce0-ad10-41ef-918b-092fc6e381f4
# ╟─3240bd88-513e-43c8-8068-9b0ce88f64d0
# ╠═be02edba-fe11-4c79-a783-d7a22ee79043
# ╠═288ad98d-5472-4e59-b5d3-45c9fbdc1a47
# ╠═b9d8d488-8b33-4941-967a-0cf9c019c91b
# ╟─a6d95216-8d80-4101-af62-1a1b1f1d936c
# ╠═f30d2b4d-23df-41c2-aaf1-28878fd28a40
# ╠═0129e1a3-3cc2-4e61-a325-1abf356b9718
# ╠═309ca624-fefe-4bcb-81e3-e5f4ebb9f773
# ╟─7c969f9e-dc47-4478-910d-a91ce541ec17
# ╠═536a1a1f-43f2-4279-979e-a510f1e97223
# ╠═1b895ec8-ed17-437b-a355-81d8039c4393
# ╠═7b5234eb-e4d5-410a-9740-20c23c62baa6
# ╠═1054ff92-52c9-408e-bc8d-194585288052
# ╟─57bf1a8f-6f62-4d37-991f-83843e2c7e73
# ╟─68f7258e-f6b1-44e3-82fa-606f65bf5c0f
# ╠═8b5b6433-cd38-4c18-8176-110e2a33965c
# ╟─28ec799f-fcfd-40ea-bddb-0c1450d14d39
# ╠═1e341397-f0c4-4a92-8b12-669f59380f9e
# ╠═17175b11-5ddc-4e78-b21a-b817049be2ca
# ╠═e9f0f552-c88e-44f5-afbc-88672cb34093
# ╟─b4e30b34-ef4c-4cf4-86a2-5e4a81a41705
# ╠═cf195f98-70b4-4ba2-bab6-d27e925070e6
# ╠═0301bf1e-acd5-4d1e-bfbe-0f6813d72357
# ╠═bf6c5954-bcef-4b38-bf6d-2850c81283e6
# ╟─f6dc3d73-6bd7-4224-a94b-ea195f63ce2c
# ╠═6ba5d457-6bfc-4f9d-8a83-1e6c6f033a49
# ╠═0d467ffd-8cc8-4c3b-8729-e08bce2dcf5d
# ╠═5c1ece3e-4cea-4107-9cea-0c23aaf0bbe1
# ╟─d2c408bc-6b14-4fca-a277-da342f0b68da
# ╠═79425b93-a29e-40f3-980f-0e5774f75a6b
# ╟─0478b474-d354-4c62-a56a-f868d70a042b
# ╠═93c337b6-3398-4998-8cde-c12a411e23df
# ╠═5c7df635-9388-421d-8929-f0f2fa3f4862
# ╠═8d987296-7c94-44a9-ac43-bba5d455bbc9
# ╠═b7b2917f-86af-4591-b298-64941ec91665
# ╠═63c468be-281a-4cfd-a546-260239c89fbd
# ╠═7cfbb85a-49f2-4826-a814-7eb9f5696ce4
# ╠═20d9dd5d-d4f1-4c04-99c4-9d18fff34d04
# ╠═edeea6af-ad4a-430d-97d7-5e62d7c68aab
# ╠═8dc152a0-00d7-46e7-8368-33c2a36da2f6
# ╠═5ffeba59-c088-45ee-9af8-0a593eb32dce
# ╠═37126474-6677-45b0-a425-50031b592e8b
# ╠═f7e87875-b9fc-49cb-8b81-2a75d1be11fb
# ╠═d1ef2035-7854-44f3-9286-dc5b5b36097a
# ╠═a345439c-83d7-4c45-8e98-77db1b36f279
# ╠═339600a4-236e-46b3-9612-f30ee7f6ba78
# ╠═fccb3dba-69ba-4782-9181-65ac303813b2
# ╠═a9dd7902-1871-4ea8-8f58-0e50adf82f0b
# ╠═297bf17a-1a57-4422-b15d-599eef64b125
# ╠═c4b622c7-b037-4839-a61e-759f43dafeee
# ╠═e8a8f9cf-7359-4131-b307-c8b55bf73546
# ╠═35e417b7-a71e-4827-b154-0419888952e2
# ╠═81c064f4-cff9-4ed7-b2a8-a620cf15ccaf
# ╠═3e62017a-90e3-4948-8d16-effe14e5f5af
# ╠═dcb8d1b7-5e9e-4078-b86a-a8c572a870a9
# ╠═1f579a15-3efc-4bc0-ba4c-dcc11dcb4c63
# ╠═e8ec60e9-5ee6-4085-8c1e-907aa2d4624a
# ╠═b6525a66-35ae-4a74-937d-db51fa7a3024
# ╠═e3280173-b6ac-4cc1-945b-ec1c394c9260
# ╠═559669c9-8539-4d83-ac0b-3135b127028b
# ╠═b883f972-ce50-4400-8bc5-edce46a6689d
# ╠═f0e695a3-eb78-47ba-8e4a-6af332e1b496
# ╠═069d2662-24c6-49af-be33-8d269b05bc9b
# ╠═bbaff240-af5c-4621-97ba-613040420f6e
# ╠═da33cf36-6c3d-4c94-a1f1-59a71ae06acb
# ╠═a48fc90b-f99d-4408-be54-d170977e394f
# ╠═2563b25a-a386-4f25-b88e-d5db49a1473a
# ╠═6fb7e314-c5cc-4c90-90c0-9b9db04cfe1b
# ╠═fe263281-fdd2-4025-a75c-11701059ea78
# ╠═2ee67d3b-b2a7-4605-9bbb-06caa5c888e2
# ╠═b040556f-5afe-4539-be19-ef5eaecb8448
# ╠═ca315565-2e7f-435a-be3f-05d4893cfcce
# ╠═164dd82f-47bd-4255-9507-f0386cd94b40
# ╠═35fa87bb-ce0e-4505-b667-6de69e3641ed
# ╠═48ce3652-c8ef-4fcc-b284-1cb8bb720da8
# ╠═f1a0c8a1-b79d-435f-907d-64b5671e65c2
# ╠═7e7d92bb-d74a-470a-a0f2-e48b62ae56e0
# ╠═e3c667a9-b3ff-494f-943e-8c797f31bd91
# ╠═4c7c4aeb-2f16-4aa8-acc9-1e9be473b57f
# ╠═b6c84ed9-573b-472b-be9b-fe4b91bbf40b
# ╠═f1f32ea8-2a17-478c-9ca1-c5ce7afe3912
# ╠═b315f712-7630-4f5f-af85-bdfc46d4d5f7
# ╠═47e47506-08de-4991-90fa-35e3e6b18c65
# ╠═c8756384-6c90-4819-b0dc-2798eca244a8
# ╠═55a9dfbf-6f8a-40ab-9e23-5fe5231a250f
# ╠═c10d0942-ba67-4c71-ae53-1e2ed5020dae
# ╠═c384cf03-98fb-4c5e-b426-487b688e1e69
# ╠═23189d22-61e7-4c30-b7ac-50d20692ecd5
# ╠═75b50ca5-f1f5-494a-85f1-33881db475b3
# ╠═4233274a-f1ec-4cd9-b68f-c68f2f10dfa0
# ╠═976fc996-2620-45a7-aacc-da7da24bfc0b
# ╠═76e2b877-4ebf-4550-b796-52dfbf55d361
# ╠═3b7ddaa4-4e99-4710-b2fe-26be6288edbc
# ╠═eb3bbceb-2198-438a-9f25-18e074b1423e
# ╠═0fbac767-2689-435e-a123-a3ea84f0a294
# ╠═8b4ba108-bea3-4850-a1fc-0fad0cd7a6ea
# ╠═4a6bd3b1-53ed-48f4-916f-e1665aab1664
# ╠═ac73fe19-8ad0-44a7-bb95-95af48f170d5
# ╠═14e7d218-b103-413a-9ce3-51f051243c56
# ╠═4ba64fcd-b8d3-4045-8a47-86d1705d19ff
# ╠═fa08de3f-27b6-4ff3-917e-73ea692296b0
# ╠═dc3ae6f7-8fa7-42ca-b2f2-75e42befa5ce
# ╠═8162eea8-e609-490c-8f8f-0a058a3c2d8b
# ╠═75530e55-95e8-498f-82e0-3e4d7bac8ed4
# ╠═6a094720-d715-4e81-ad6a-e23efa13d545
# ╠═514c21a8-2e4e-47ed-b6e5-53d67c049ebb
# ╠═bd11ad26-f577-45cb-beea-dc9231cbf0a1
# ╠═75cea2c9-46f8-42af-a92f-5124318c21c2
# ╠═62340699-9e33-4e45-a926-5c6d84e8b4ce
# ╠═ba56da4a-e854-481f-9088-4be7c8be425d
# ╠═82d5fdb6-0734-4eff-b329-ac640c40bc3e
# ╠═cf8232e0-bab8-4e9b-9c85-2b5df5bbe5be
# ╟─d758c9db-39e1-42cb-9f3f-e8cb962c925f
# ╠═e564c3e3-058a-49e9-9a79-68375ab4482b
# ╠═1cf8f266-448e-4ece-bfed-41db8fb23526
# ╠═22ffc892-6f90-4172-94f6-06d39237715d
# ╠═2031a3c1-db8b-4a31-b1cb-3caf2b7537c3
# ╠═21935abe-702f-4bd3-8cee-dddad07d23e2
# ╠═fdf4804b-af2a-42ba-9d12-1d4e0a9c4496
# ╠═338ad5b0-7c19-4b14-a513-cdabb0ef05f4
# ╠═12bf244c-d9f8-456e-b818-81d024f451b2
# ╠═925a05e7-10b8-484d-a44c-9b2e532b63aa
# ╠═d333096e-a985-4abb-9b32-a13dc82d824b
# ╟─812299b1-c58a-40c6-8d93-80becfcc008e
# ╠═c25e5b1d-053f-49b7-ad1f-f4c417aa491c
# ╠═5947da89-8010-48ac-a5b3-3e9371280288
# ╠═443b7def-1fb5-463c-bd88-8975c6f59501
# ╠═67bd57ff-7f71-48c7-9e3b-92e6b7e5bc30

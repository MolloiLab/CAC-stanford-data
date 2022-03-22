### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ d5e88f20-a3de-11ec-2c6b-350b7ea5c841
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
		Pkg.add(url="https://github.com/Dale-Black/PhantomSegmentation.jl")
		Pkg.add(url="https://github.com/Dale-Black/CalciumScoring.jl")
		Pkg.add("CairoMakie")
		Pkg.add("HypothesisTests")
		Pkg.add("Colors")
	end
	
	using PlutoUI
	using Statistics
	using StatsBase: quantile!, rmsd
	using ImageMorphology
	using ImageFiltering
	using CSV
	using DataFrames
	using GLM
	using DICOM
	using DICOMUtils
	using PhantomSegmentation
	using CalciumScoring
	using CairoMakie
	using HypothesisTests
	using Colors
end

# ╔═╡ ada1a311-b8cc-4376-9811-104737119f81
TableOfContents()

# ╔═╡ 0e222b32-34f5-4355-b6f2-cc22b984496e
md"""
# Integrated
"""

# ╔═╡ afd300fc-e1e4-4a81-b80d-e56c080cc309
md"""
## Canon
"""

# ╔═╡ fec52734-5ec4-42e6-9257-494ebe0ae3a7
canon_path = "/Users/daleblack/Google Drive/dev/MolloiLab/CAC-stanford-data/data/output/Canon_Aquilion_One_Vision_integrated_score_script";

# ╔═╡ 3b0b0d7d-7551-4116-b570-92840e4747cb
df_i_canon = CSV.read(string(canon_path, "/output.csv"), DataFrame);

# ╔═╡ 33a88630-0059-496f-ab43-6f5a62dd672f
inserts = repeat([200, 400, 800], 10)

# ╔═╡ 8bdc88e6-a8da-408f-ba72-76ccb0dea5ea
md"""
## GE
"""

# ╔═╡ 55940006-c781-4fa1-b2e9-5d60dd2a27cb
ge_path = "/Users/daleblack/Google Drive/dev/MolloiLab/CAC-stanford-data/data/output/GE_Revolution_integrated_score_script";

# ╔═╡ 9877afc3-c563-4b52-a472-c64abd201848
df_i_ge = CSV.read(string(ge_path, "/output.csv"), DataFrame);

# ╔═╡ 667bf5e4-b73f-41e3-a429-6e72a5a913eb
md"""
## Philips
"""

# ╔═╡ 3cbc7db2-c0bd-43e3-9e14-a1582f14f260
p_path = "/Users/daleblack/Google Drive/dev/MolloiLab/CAC-stanford-data/data/output/Philips_Brilliance_iCT_integrated_score_script";

# ╔═╡ d65feef0-e76d-4904-a1dd-fe17ce690823
df_i_p = CSV.read(string(p_path, "/output.csv"), DataFrame);

# ╔═╡ 28bcc00c-ee82-4140-93d1-fee53b30d9ff
md"""
# Agatston/Mass
"""

# ╔═╡ fbf02079-6286-42c7-a4b4-010acce061d5
md"""
## Canon
"""

# ╔═╡ 50e90e4d-e584-4628-9577-a29eb83acacc
canon_path_agat = "/Users/daleblack/Google Drive/dev/MolloiLab/CAC-stanford-data/data/output/Canon_Aquilion_One_Vision_agatston_mass_score_script";

# ╔═╡ 24c33f8e-adb7-4f35-b2a1-5dc0572b1ffd
df_a_canon = CSV.read(string(canon_path_agat, "/output.csv"), DataFrame);

# ╔═╡ 23da2c11-27db-4860-8c7a-fd6406a31856
md"""
## GE
"""

# ╔═╡ 7e0af03c-0ffc-4762-ac6e-624791c5f286
ge_path_agat = "/Users/daleblack/Google Drive/dev/MolloiLab/CAC-stanford-data/data/output/GE_Revolution_agatston_mass_score_script";

# ╔═╡ fc2a9ae1-95a6-49c1-91d2-58b4664af146
df_a_ge = CSV.read(string(ge_path_agat, "/output.csv"), DataFrame);

# ╔═╡ 1cb3defe-9b63-4aa5-92ac-4f4deb2e3c20
md"""
## Philips
"""

# ╔═╡ 526bf7c0-bf2f-43fc-aca5-c46c99eff6de
p_path_agat = "/Users/daleblack/Google Drive/dev/MolloiLab/CAC-stanford-data/data/output/Philips_Brilliance_iCT_agatston_mass_score_script";

# ╔═╡ 63e9dbbb-8add-4de1-8d30-d871eaa311c0
df_a_p = CSV.read(string(p_path_agat, "/output.csv"), DataFrame);

# ╔═╡ b7e2164b-063c-49e2-b3b8-ac561df22eb0
md"""
# RMSD
"""

# ╔═╡ b7a2b1bb-1bcd-4ee8-9577-a8634b980944
md"""
## Large Inserts
"""

# ╔═╡ 18c64198-b1f4-47b4-a08d-552a1b99e715
int_rms_large_c = rmsd(df_i_canon[!, :calculated_mass_large], df_i_canon[!, :ground_truth_mass_large])

# ╔═╡ aca94cc3-73b5-4808-babb-fdb97a160fe0
int_rms_large_g = rmsd(df_i_ge[!, :calculated_mass_large], df_i_ge[!, :ground_truth_mass_large])

# ╔═╡ fa492c39-948c-4e88-b6cd-376565bbcb4f
int_rms_large_p = rmsd(df_i_p[!, :calculated_mass_large], df_i_p[!, :ground_truth_mass_large])

# ╔═╡ 90b0a22d-d30a-4699-bea9-c74edddda55f
agat_rms_large_c = rmsd(df_a_canon[!, :calculated_mass_large], df_a_canon[!, :ground_truth_mass_large])

# ╔═╡ 79fcc2d0-60d6-4664-98ca-3284ea90e399
agat_rms_large_g = rmsd(df_a_ge[!, :calculated_mass_large], df_a_ge[!, :ground_truth_mass_large])

# ╔═╡ df5b3e6c-e5fb-4a24-ab9d-8c2687dc5fba
agat_rms_large_p = rmsd(df_a_p[!, :calculated_mass_large], df_a_p[!, :ground_truth_mass_large])

# ╔═╡ 3c23f3e0-de5b-4d5c-bc9a-3a80f7d62489
md"""
## Medium Inserts
"""

# ╔═╡ bece3d60-e5a1-4aee-a719-3f9a6d4f7464
int_rms_medium_c = rmsd(df_i_canon[!, :calculated_mass_medium], df_i_canon[!, :ground_truth_mass_medium])

# ╔═╡ 6b583951-ee79-4de2-8d45-a1eed776a867
int_rms_medium_g = rmsd(df_i_ge[!, :calculated_mass_medium], df_i_ge[!, :ground_truth_mass_medium])

# ╔═╡ 97b723f6-95ef-4cc2-9329-dc17de6d00f5
int_rms_medium_p = rmsd(df_i_p[!, :calculated_mass_medium], df_i_p[!, :ground_truth_mass_medium])

# ╔═╡ 25eb7c73-6e0a-4459-8f22-b8787017abba
agat_rms_medium_c = rmsd(df_a_canon[!, :calculated_mass_medium], df_a_canon[!, :ground_truth_mass_medium])

# ╔═╡ 68632cf8-d048-4b18-b76c-86ba515f25f7
agat_rms_medium_g = rmsd(df_a_ge[!, :calculated_mass_medium], df_a_ge[!, :ground_truth_mass_medium])

# ╔═╡ e227d63b-a9ac-4144-8885-d39f5131a527
agat_rms_medium_p = rmsd(df_a_p[!, :calculated_mass_medium], df_a_p[!, :ground_truth_mass_medium])

# ╔═╡ 6a69d11d-f495-4289-b543-7ea78fcfd41c
md"""
## Small Inserts
"""

# ╔═╡ 502d16ab-85da-47a8-95c3-1a28dde03360
int_rms_small_c = rmsd(df_i_canon[!, :calculated_mass_small], df_i_canon[!, :ground_truth_mass_small])

# ╔═╡ a57316d8-a5a8-4347-9888-e656a8d5019d
int_rms_small_g = rmsd(df_i_ge[!, :calculated_mass_small], df_i_ge[!, :ground_truth_mass_small])

# ╔═╡ dfc49d5e-e18a-42e3-bbbe-8f162f1c9eda
int_rms_small_p = rmsd(df_i_p[!, :calculated_mass_small], df_i_p[!, :ground_truth_mass_small])

# ╔═╡ 04a32b4d-8c5d-4475-b2b5-8b878ba01661
agat_rms_small_c = rmsd(df_a_canon[!, :calculated_mass_small], df_a_canon[!, :ground_truth_mass_small])

# ╔═╡ c6b59690-3051-43c5-b678-0d9f437d4577
agat_rms_small_g = rmsd(df_a_ge[!, :calculated_mass_small], df_a_ge[!, :ground_truth_mass_small])

# ╔═╡ 82a305a9-fc59-4ecc-bbf8-bd07bca52b39
agat_rms_small_p = rmsd(df_a_p[!, :calculated_mass_small], df_a_p[!, :ground_truth_mass_small])

# ╔═╡ 509e0e2e-0cee-49dc-9c26-bd44d22153c8
md"""
## Average RMSD
"""

# ╔═╡ b714cc63-d714-44dc-9070-4be6107e7c22
md"""
### Large
"""

# ╔═╡ 8ea31ce1-0edd-4b61-9398-b6ff75a43f71
avg_l_i = mean([int_rms_large_c, int_rms_large_g, int_rms_large_p])

# ╔═╡ 176f21d0-247f-4492-8d7a-1ef97bc6d8b5
avg_l_a = mean([agat_rms_large_c, agat_rms_large_g, agat_rms_large_p])

# ╔═╡ 60986505-c141-46ce-873e-29924b93b967
begin
	large_ttest = OneSampleTTest(df_i_canon[!, :calculated_mass_large], df_a_canon[!, :calculated_mass_large])
	large_pvalue = pvalue(large_ttest)
	if large_pvalue < 0.05
		large_pvalue = "≤ 0.05"
	else
		large_pvalue = string(large_pvalue)
	end
end

# ╔═╡ fa0dfd59-2e4c-467a-8068-1df9f8923c3e
large_std_i = std(df_i_canon[!, :calculated_mass_large])

# ╔═╡ 60ecda95-c33a-4bed-b9ad-4d9742a07c98
large_std_a = std(df_a_canon[!, :calculated_mass_large])

# ╔═╡ 56446bb1-43a9-4cbf-a9b3-bdb4599f248a
md"""
### Medium
"""

# ╔═╡ a6d41ade-7d36-487c-b731-272b3eadc1ba
avg_m_i = mean([int_rms_medium_c, int_rms_medium_g, int_rms_medium_p])

# ╔═╡ 87c5314e-027e-4d4d-9234-0af49919a1f5
avg_m_a = mean([agat_rms_medium_c, agat_rms_medium_g, agat_rms_medium_p])

# ╔═╡ fa07ebba-5719-403e-90ac-cecace8a4695
begin
	medium_ttest = OneSampleTTest(df_i_canon[!, :calculated_mass_medium], df_a_canon[!, :calculated_mass_medium])
	medium_pvalue = pvalue(medium_ttest)
	if medium_pvalue < 0.05
		medium_pvalue = "≤ 0.05"
	else
		medium_pvalue = string(medium_pvalue)
	end
end

# ╔═╡ 94a6fa3a-4f04-4f52-9640-c504fca27525
medium_std_i = std(df_i_canon[!, :calculated_mass_medium])

# ╔═╡ 93a296c8-6765-43b9-904d-eadc34479867
medium_std_a = std(df_a_canon[!, :calculated_mass_medium])

# ╔═╡ 28981589-b284-40c9-9916-5843a0bb8118
md"""
### Small
"""

# ╔═╡ 0731371d-6511-4f5d-b318-6a4122fa30c5
avg_s_i = mean([int_rms_small_c, int_rms_small_g, int_rms_small_p])

# ╔═╡ 2e7a5c57-c91f-4780-ba8d-bc59a530da11
avg_s_a = mean([agat_rms_small_c, agat_rms_small_g, agat_rms_small_p])

# ╔═╡ 6ac3ee1b-5365-4a3e-a2e2-1e1dfa7762f7
begin
	small_ttest = OneSampleTTest(df_i_canon[!, :calculated_mass_small], df_a_canon[!, :calculated_mass_small])
	small_pvalue = pvalue(small_ttest)
	if small_pvalue < 0.05
		small_pvalue = "≤ 0.05"
	else
		small_pvalue = string(round(small_pvalue, digits=2))
	end
end

# ╔═╡ 771f119b-eb79-4638-b016-ad5252fabb66
small_std_i = std(df_i_canon[!, :calculated_mass_small])

# ╔═╡ f94b1309-a68f-48f1-88d2-30112224450c
small_std_a = std(df_a_canon[!, :calculated_mass_small])

# ╔═╡ 174cc8bc-1a74-43eb-aee7-8a3a629e1eb9
md"""
## Total
"""

# ╔═╡ 25f30d11-2e53-43f0-a565-70a33ed357d4
rmsd_tot_i = mean([int_rms_large_c, int_rms_medium_c, int_rms_small_c])

# ╔═╡ 1181f896-0491-43d0-bce1-44a9314faa65
rmsd_tot_a = mean([agat_rms_large_c, agat_rms_medium_c, agat_rms_small_c])

# ╔═╡ 119a83d5-c016-4282-941d-faffc1c57f90
md"""
# Canon Only
"""

# ╔═╡ 79a62129-c9ea-430a-862a-d6bc86658bb3
inserts2 = ["Small Inserts", "Medium Inserts", "Large Inserts"]

# ╔═╡ 54d0e363-f838-4712-b46f-423291bf58a8
df = DataFrame(
	"Inserts" => inserts2,
	"iCAC RMSD" => [int_rms_small_c, int_rms_medium_c, int_rms_large_c],
	"aCAC RMSD" => [agat_rms_small_c, agat_rms_medium_c, agat_rms_large_c],
	"P Values (iCAc vs aCAC)" => [small_pvalue, medium_pvalue, large_pvalue],
	"Variance (iCAC)" => [small_std_i, medium_std_i, large_std_i],
	"Variance (aCAC)" => [small_std_a, medium_std_a, large_std_a]
)

# ╔═╡ d5f64fac-e9c7-4137-a43a-a1588e001f08
begin
	f1a = Figure()
	ax1a = Axis(f1a[1, 1])

	scatter!(inserts, df_i_canon[!, :calculated_mass_large], label="Calculated Mass")
	scatter!(inserts, df_i_canon[!, :ground_truth_mass_large], label="Known Mass")
	
	ax1a.title = "Large Inserts"
	ax1a.ylabel = "iCAC Mass Score (mg)"
	# ax1a.xlabel = "Density (mg/cc)"
	ax1a.xticks = [200, 400, 800]

	ylims!(ax1a, 0, 100)

	ax2a = Axis(f1a[1, 2])

	scatter!(inserts, df_i_canon[!, :calculated_mass_medium], label="calculated_mass_medium")
	scatter!(inserts, df_i_canon[!, :ground_truth_mass_medium], label="ground_truth_mass_medium")
	
	ax2a.title = "Medium Inserts"
	# ax2a.xlabel = "Density (mg/cc)"
	ax2a.xticks = [200, 400, 800]

	ylims!(ax2a, 0, 25)
	
	ax3a = Axis(f1a[1, 3])

	scatter!(inserts, df_i_canon[!, :calculated_mass_small], label="calculated_mass_small")
	scatter!(inserts, df_i_canon[!, :ground_truth_mass_small], label="ground_truth_mass_small")
	
	ax3a.title = "Small Inserts"
	# ax3a.xlabel = "Density (mg/cc)"
	ax3a.xticks = [200, 400, 800]

	ylims!(ax3a, 0, 2)
	
	# f1a[1, 3] = Legend(f1a, ax3a, framevisiblae = false)
	ax1b = Axis(f1a[2, 1])

	scatter!(inserts, df_a_canon[!, :calculated_mass_large], label="Calculated Mass")
	scatter!(inserts, df_a_canon[!, :ground_truth_mass_large], label="Known Mass")
	
	# ax1b.title = "Large Inserts"
	ax1b.ylabel = "aCAC Mass Score (mg)"
	# ax1b.xlabel = "Density (mg/cc)"
	ax1b.xticks = [200, 400, 800]

	ylims!(ax1b, 0, 100)

	ax2b = Axis(f1a[2, 2])

	scatter!(inserts, df_a_canon[!, :calculated_mass_medium], label="calculated_mass_medium")
	scatter!(inserts, df_a_canon[!, :ground_truth_mass_medium], label="ground_truth_mass_medium")
	
	# ax2b.title = "Medium Inserts"
	ax2b.xlabel = "Density (mg/cc)"
	ax2b.xticks = [200, 400, 800]

	ylims!(ax2b, 0, 25)
	
	ax3b = Axis(f1a[2, 3])

	scatter!(inserts, df_a_canon[!, :calculated_mass_small], label="calculated_mass_small")
	scatter!(inserts, df_a_canon[!, :ground_truth_mass_small], label="ground_truth_mass_small")
	
	# ax3b.title = "Small Inserts"
	# ax3b.xlabel = "Density (mg/cc)"
	ax3b.xticks = [200, 400, 800]

	ylims!(ax3b, 0, 2)
	
	f1a
end

# ╔═╡ eedbec30-34e2-4f1a-8e81-733fdef72302
# save("larger.pdf", f1a, pt_per_unit = 2)

# ╔═╡ b818b581-e36f-413a-a16c-1783f88c1fb9
md"""
## Bar chart
"""

# ╔═╡ 35d577a7-be59-4524-b7a1-630ad8c9515d
md"""
### Small Only
"""

# ╔═╡ 11cfb2cf-4159-4a1c-b651-21f432de0117
num_nonzero_i = length(findall(x -> x != 0, df_i_canon[!, :calculated_mass_small]))

# ╔═╡ c73cb2c4-f640-4822-9f92-ca0332170b75
num_zero_i = length(findall(x -> x == 0, df_i_canon[!, :calculated_mass_small]))

# ╔═╡ 059fd708-3aca-4d90-a831-1a48c6a79f9c
num_nonzero_a = length(findall(x -> x != 0, df_a_canon[!, :calculated_mass_small]))

# ╔═╡ 819272ec-5300-47cb-bec8-46617d4b5e17
num_zero_a = length(findall(x -> x == 0, df_a_canon[!, :calculated_mass_small]))

# ╔═╡ 15e5bf9b-ae71-4523-a02e-554a312707fe
begin
	x = [1, 1, 2, 2]
	height = [num_zero_i, num_zero_a, num_nonzero_i, num_nonzero_a]
	grp = [1, 2, 1, 2]
end

# ╔═╡ b5374ee9-607d-4b03-82f6-a9bef97c272d
colors = cgrad(:tab10)

# ╔═╡ ab1445fd-9462-40ee-b105-1e0db7748f06
md"""
### Total
"""

# ╔═╡ e68a4a0d-b533-4484-8f50-fb752268e907
tot_i = hcat(df_i_canon[!, :calculated_mass_small], df_i_canon[!, :calculated_mass_medium], df_i_canon[!, :calculated_mass_large]);

# ╔═╡ 6c9d6fe2-9e35-46b6-acd9-61b5feb76d5e
tot_a = hcat(df_a_canon[!, :calculated_mass_small], df_a_canon[!, :calculated_mass_medium], df_a_canon[!, :calculated_mass_large]);

# ╔═╡ 5d1257f3-854a-4222-8065-ad9256f1365a
num_nonzero_i_tot = length(findall(x -> x != 0, tot_i))

# ╔═╡ c1f584cc-26a6-42af-898e-d5d9d6c1e111
num_zero_i_tot = length(findall(x -> x == 0, tot_i))

# ╔═╡ 36cd1c4b-69c0-4bb1-9d42-7cb895fac23d
num_nonzero_a_tot = length(findall(x -> x != 0, tot_a))

# ╔═╡ c3838058-b3ab-4975-b5d0-d5e7233d4d85
num_zero_a_tot = length(findall(x -> x == 0, tot_a))

# ╔═╡ f2834ca0-967e-4896-ae1d-1e5bf908089b
height2 = [num_zero_i_tot, num_zero_a_tot, num_nonzero_i_tot, num_nonzero_a_tot]

# ╔═╡ 54c5a455-1521-40ef-acca-28a9ef276bac


# ╔═╡ 383214b7-8b7e-47e1-8bca-da6a5d737c2e
md"""
## IHS/AS vs Known Mass
"""

# ╔═╡ 0093205e-f8cf-4689-b756-823edd57a454
begin
	fInt = Figure()
	axInt = Axis(fInt[1, 1])

	scatter!(df_i_canon[!, :ground_truth_mass_large], df_i_canon[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df_i_canon[!, :ground_truth_mass_medium], df_i_canon[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df_i_canon[!, :ground_truth_mass_small], df_i_canon[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity")
	
	axInt.title = "IHU Mass vs Known Mass"
	axInt.ylabel = "Calculated Mass (mg)"
	axInt.xlabel = "Known Mass (mg)"
	axInt.xticks = [0, 25, 50, 75, 100]
	axInt.yticks = [0, 25, 50, 75, 100]

	xlims!(axInt, 0, 100)
	ylims!(axInt, 0, 100)
	
	fInt[1, 3] = Legend(fInt, axInt, framevisible = false)

	axAgat = Axis(fInt[1, 2])

	scatter!(df_a_canon[!, :ground_truth_mass_large], df_a_canon[!, :calculated_mass_large], label="Mass: Large Inserts")
	scatter!(df_a_canon[!, :ground_truth_mass_medium], df_a_canon[!, :calculated_mass_medium], label="Mass: Medium Inserts")
	scatter!(df_a_canon[!, :ground_truth_mass_small], df_a_canon[!, :calculated_mass_small], label="Mass: Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity Line")
	
	axAgat.title = "AS Mass vs Known Mass"
	axAgat.xlabel = "Known Mass (mg)"
	axAgat.xticks = [0, 25, 50, 75, 100]
	axAgat.yticks = [0, 25, 50, 75, 100]

	xlims!(axAgat, 0, 100)
	ylims!(axAgat, 0, 100)

	axbar2 = Axis(fInt[2, 1])

	barplot!(axbar2, x, height2, dodge = grp, color =  [colors[1], colors[2], colors[1], colors[2]])

	axbar2.xticks = (1:2, ["CAC = 0", "CAC > 0"])
	axbar2.ylabel = "Number of Inserts"
	axbar2.title = "Calcium Scores (All Inserts)"

	labels2 = ["Integrated Score", "Agatston Score"]
	elements2 = [PolyElement(polycolor = colors[i]) for i in 1:length(labels2)]
	title2 = "Scoring Techniques"

	axbar1 = Axis(fInt[2, 2])

	barplot!(axbar1, x, height, dodge = grp, color =  [colors[1], colors[2], colors[1], colors[2]])

	axbar1.xticks = (1:2, ["CAC = 0", "CAC > 0"])
	axbar1.title = "Calcium Scores (Small Inserts)"

	labels = ["Integrated HU Score", "Agatston Score"]
	elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
	title = "Scoring Techniques"
	
	Legend(fInt[2, 3], elements, labels, title)
	
	fInt
end

# ╔═╡ 7f4b8ddc-ff5c-48db-a159-3f21f8d52939
# save("scct2.png", fInt, pt_per_unit = 3)

# ╔═╡ b1990d31-40d4-4acc-8bcf-00b6744ac8dc
md"""
# CANON - Calibration Types
"""

# ╔═╡ 45982b4c-4b55-40b2-b290-301cea893e4d
md"""
## Point Calibration
"""

# ╔═╡ 28818e58-6568-43e8-bfa9-3eba9ef9dd4b
canon_path2 = "/Users/daleblack/Google Drive/dev/MolloiLab/CAC-stanford-data/data/output/Canon_Aquilion_One_Vision4";

# ╔═╡ 29ba89eb-339e-4329-a34c-1e9737660b8b
df_i_canon2 = CSV.read(string(canon_path2, "/full.csv"), DataFrame);

# ╔═╡ 27359b10-27df-4fe6-bbac-23501958bd3b
begin
	fInt2 = Figure()
	axInt2 = Axis(fInt2[1, 1])

	scatter!(df_i_canon2[!, :ground_truth_mass_large], df_i_canon2[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df_i_canon2[!, :ground_truth_mass_medium], df_i_canon2[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df_i_canon2[!, :ground_truth_mass_small], df_i_canon2[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity")
	
	axInt2.title = "IHU Mass vs Known Mass"
	axInt2.ylabel = "Calculated Mass (mg)"
	axInt2.xlabel = "Known Mass (mg)"
	axInt2.xticks = [0, 25, 50, 75, 100]
	axInt2.yticks = [0, 25, 50, 75, 100]

	xlims!(axInt2, 0, 100)
	ylims!(axInt2, 0, 100)
	
	fInt2[1, 2] = Legend(fInt2, axInt2, framevisible = false)
	fInt2
end

# ╔═╡ 82212a13-576b-4299-92f8-434fcca599c1
md"""
## Line Calibration
"""

# ╔═╡ 2550f774-423a-44fa-860c-638290caf2cf
canon_path3 = "/Users/daleblack/Google Drive/dev/MolloiLab/CAC-stanford-data/data/output/Canon_Aquilion_One_Vision3";

# ╔═╡ 9a165cf4-9574-43ba-9e49-34153a473b49
df_i_canon3 = CSV.read(string(canon_path3, "/full.csv"), DataFrame);

# ╔═╡ e6ee6ed7-3f45-44d2-8aa5-2d96b8b257e8
begin
	fInt3 = Figure()
	axInt3 = Axis(fInt3[1, 1])

	scatter!(df_i_canon3[!, :ground_truth_mass_large], df_i_canon3[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df_i_canon3[!, :ground_truth_mass_medium], df_i_canon3[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df_i_canon3[!, :ground_truth_mass_small], df_i_canon3[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity")
	
	axInt3.title = "IHU Mass vs Known Mass"
	axInt3.ylabel = "Calculated Mass (mg)"
	axInt3.xlabel = "Known Mass (mg)"
	axInt3.xticks = [0, 25, 50, 75, 100]
	axInt3.yticks = [0, 25, 50, 75, 100]

	xlims!(axInt3, 0, 100)
	ylims!(axInt3, 0, 100)
	
	fInt3[1, 2] = Legend(fInt3, axInt3, framevisible = false)
	fInt3
end

# ╔═╡ Cell order:
# ╠═d5e88f20-a3de-11ec-2c6b-350b7ea5c841
# ╠═ada1a311-b8cc-4376-9811-104737119f81
# ╟─0e222b32-34f5-4355-b6f2-cc22b984496e
# ╟─afd300fc-e1e4-4a81-b80d-e56c080cc309
# ╠═fec52734-5ec4-42e6-9257-494ebe0ae3a7
# ╠═3b0b0d7d-7551-4116-b570-92840e4747cb
# ╠═33a88630-0059-496f-ab43-6f5a62dd672f
# ╟─8bdc88e6-a8da-408f-ba72-76ccb0dea5ea
# ╠═55940006-c781-4fa1-b2e9-5d60dd2a27cb
# ╠═9877afc3-c563-4b52-a472-c64abd201848
# ╟─667bf5e4-b73f-41e3-a429-6e72a5a913eb
# ╠═3cbc7db2-c0bd-43e3-9e14-a1582f14f260
# ╠═d65feef0-e76d-4904-a1dd-fe17ce690823
# ╟─28bcc00c-ee82-4140-93d1-fee53b30d9ff
# ╟─fbf02079-6286-42c7-a4b4-010acce061d5
# ╠═50e90e4d-e584-4628-9577-a29eb83acacc
# ╠═24c33f8e-adb7-4f35-b2a1-5dc0572b1ffd
# ╟─23da2c11-27db-4860-8c7a-fd6406a31856
# ╠═7e0af03c-0ffc-4762-ac6e-624791c5f286
# ╠═fc2a9ae1-95a6-49c1-91d2-58b4664af146
# ╟─1cb3defe-9b63-4aa5-92ac-4f4deb2e3c20
# ╠═526bf7c0-bf2f-43fc-aca5-c46c99eff6de
# ╠═63e9dbbb-8add-4de1-8d30-d871eaa311c0
# ╟─b7e2164b-063c-49e2-b3b8-ac561df22eb0
# ╟─b7a2b1bb-1bcd-4ee8-9577-a8634b980944
# ╠═18c64198-b1f4-47b4-a08d-552a1b99e715
# ╠═aca94cc3-73b5-4808-babb-fdb97a160fe0
# ╠═fa492c39-948c-4e88-b6cd-376565bbcb4f
# ╠═90b0a22d-d30a-4699-bea9-c74edddda55f
# ╠═79fcc2d0-60d6-4664-98ca-3284ea90e399
# ╠═df5b3e6c-e5fb-4a24-ab9d-8c2687dc5fba
# ╟─3c23f3e0-de5b-4d5c-bc9a-3a80f7d62489
# ╠═bece3d60-e5a1-4aee-a719-3f9a6d4f7464
# ╠═6b583951-ee79-4de2-8d45-a1eed776a867
# ╠═97b723f6-95ef-4cc2-9329-dc17de6d00f5
# ╠═25eb7c73-6e0a-4459-8f22-b8787017abba
# ╠═68632cf8-d048-4b18-b76c-86ba515f25f7
# ╠═e227d63b-a9ac-4144-8885-d39f5131a527
# ╟─6a69d11d-f495-4289-b543-7ea78fcfd41c
# ╠═502d16ab-85da-47a8-95c3-1a28dde03360
# ╠═a57316d8-a5a8-4347-9888-e656a8d5019d
# ╠═dfc49d5e-e18a-42e3-bbbe-8f162f1c9eda
# ╠═04a32b4d-8c5d-4475-b2b5-8b878ba01661
# ╠═c6b59690-3051-43c5-b678-0d9f437d4577
# ╠═82a305a9-fc59-4ecc-bbf8-bd07bca52b39
# ╟─509e0e2e-0cee-49dc-9c26-bd44d22153c8
# ╟─b714cc63-d714-44dc-9070-4be6107e7c22
# ╠═8ea31ce1-0edd-4b61-9398-b6ff75a43f71
# ╠═176f21d0-247f-4492-8d7a-1ef97bc6d8b5
# ╠═60986505-c141-46ce-873e-29924b93b967
# ╠═fa0dfd59-2e4c-467a-8068-1df9f8923c3e
# ╠═60ecda95-c33a-4bed-b9ad-4d9742a07c98
# ╟─56446bb1-43a9-4cbf-a9b3-bdb4599f248a
# ╠═a6d41ade-7d36-487c-b731-272b3eadc1ba
# ╠═87c5314e-027e-4d4d-9234-0af49919a1f5
# ╠═fa07ebba-5719-403e-90ac-cecace8a4695
# ╠═94a6fa3a-4f04-4f52-9640-c504fca27525
# ╠═93a296c8-6765-43b9-904d-eadc34479867
# ╟─28981589-b284-40c9-9916-5843a0bb8118
# ╠═0731371d-6511-4f5d-b318-6a4122fa30c5
# ╠═2e7a5c57-c91f-4780-ba8d-bc59a530da11
# ╠═6ac3ee1b-5365-4a3e-a2e2-1e1dfa7762f7
# ╠═771f119b-eb79-4638-b016-ad5252fabb66
# ╠═f94b1309-a68f-48f1-88d2-30112224450c
# ╟─174cc8bc-1a74-43eb-aee7-8a3a629e1eb9
# ╠═25f30d11-2e53-43f0-a565-70a33ed357d4
# ╠═1181f896-0491-43d0-bce1-44a9314faa65
# ╟─119a83d5-c016-4282-941d-faffc1c57f90
# ╠═79a62129-c9ea-430a-862a-d6bc86658bb3
# ╟─54d0e363-f838-4712-b46f-423291bf58a8
# ╟─d5f64fac-e9c7-4137-a43a-a1588e001f08
# ╠═eedbec30-34e2-4f1a-8e81-733fdef72302
# ╟─b818b581-e36f-413a-a16c-1783f88c1fb9
# ╟─35d577a7-be59-4524-b7a1-630ad8c9515d
# ╠═11cfb2cf-4159-4a1c-b651-21f432de0117
# ╠═c73cb2c4-f640-4822-9f92-ca0332170b75
# ╠═059fd708-3aca-4d90-a831-1a48c6a79f9c
# ╠═819272ec-5300-47cb-bec8-46617d4b5e17
# ╠═15e5bf9b-ae71-4523-a02e-554a312707fe
# ╠═b5374ee9-607d-4b03-82f6-a9bef97c272d
# ╟─ab1445fd-9462-40ee-b105-1e0db7748f06
# ╠═e68a4a0d-b533-4484-8f50-fb752268e907
# ╠═6c9d6fe2-9e35-46b6-acd9-61b5feb76d5e
# ╠═5d1257f3-854a-4222-8065-ad9256f1365a
# ╠═c1f584cc-26a6-42af-898e-d5d9d6c1e111
# ╠═36cd1c4b-69c0-4bb1-9d42-7cb895fac23d
# ╠═c3838058-b3ab-4975-b5d0-d5e7233d4d85
# ╠═f2834ca0-967e-4896-ae1d-1e5bf908089b
# ╠═54c5a455-1521-40ef-acca-28a9ef276bac
# ╟─383214b7-8b7e-47e1-8bca-da6a5d737c2e
# ╟─0093205e-f8cf-4689-b756-823edd57a454
# ╠═7f4b8ddc-ff5c-48db-a159-3f21f8d52939
# ╟─b1990d31-40d4-4acc-8bcf-00b6744ac8dc
# ╟─45982b4c-4b55-40b2-b290-301cea893e4d
# ╠═28818e58-6568-43e8-bfa9-3eba9ef9dd4b
# ╠═29ba89eb-339e-4329-a34c-1e9737660b8b
# ╟─27359b10-27df-4fe6-bbac-23501958bd3b
# ╟─82212a13-576b-4299-92f8-434fcca599c1
# ╠═2550f774-423a-44fa-860c-638290caf2cf
# ╠═9a165cf4-9574-43ba-9e49-34153a473b49
# ╟─e6ee6ed7-3f45-44d2-8aa5-2d96b8b257e8

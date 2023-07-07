using Base: Float64
include("utils.jl")
include("settings.jl")
include("SolverDLRA.jl")

using PyPlot
using NPZ
using DelimitedFiles


s = Settings(); # create settings class with 501 x 501 spatial cells and a maximal rank of 100


############## read data ###################

u = npzread("data/SolutionUFullNx$(s.NCellsX)Ny$(s.NCellsY)N$(s.nPN)tEnd$(s.tEnd).npy");
B = npzread("data/SolutionBFullNx$(s.NCellsX)Ny$(s.NCellsY)N$(s.nPN)tEnd$(s.tEnd).npy");
rankInTime = npzread("data/SolutionRankDLRANx$(s.NCellsX)Ny$(s.NCellsY)N$(s.nPN)eps$(s.epsAdapt)tEnd$(s.tEnd).npy");
uDLRA = npzread("data/SolutionUDLRANx$(s.NCellsX)Ny$(s.NCellsY)N$(s.nPN)eps$(s.epsAdapt)tEnd$(s.tEnd).npy");
BDLRA = npzread("data/SolutionBDLRANx$(s.NCellsX)Ny$(s.NCellsY)N$(s.nPN)eps$(s.epsAdapt)tEnd$(s.tEnd).npy");
mass = npzread("data/SolutionMassFullNx$(s.NCellsX)Ny$(s.NCellsY)N$(s.nPN)tEnd$(s.tEnd).npy");
massDLRA = npzread("data/SolutionMassDLRANx$(s.NCellsX)Ny$(s.NCellsY)N$(s.nPN)eps$(s.epsAdapt)tEnd$(s.tEnd).npy");


################################################################
######################### execute code #########################
################################################################

################### run explicit Euler ###################
u1 = Vec2Mat(s.NCellsX,s.NCellsY,u);
B1 = Vec2Mat(s.NCellsX,s.NCellsY,B);
idx = findall(B1 .< 0); B1[idx] .= NaN;
u2 = Vec2Mat(s.NCellsX,s.NCellsY,uDLRA);
B2 = Vec2Mat(s.NCellsX,s.NCellsY,BDLRA);
idx = findall(B2 .< 0); B2[idx] .= NaN;

############################################################
######################### plotting #########################
############################################################

################### plot scalar fluxes full ###################

## explicit Euler
fig = figure("scalar_flux",figsize=(10,10),dpi=100)
ax = gca()
im1 = pcolormesh(s.xMid,s.yMid, 4.0*pi*sqrt(2)*u1[:,:,1],vmin=0.0)
color_bar = fig.colorbar(im1, ax=ax, pad=0.03, shrink = 0.71)
ax.tick_params("both",labelsize=20) 
ax.set(adjustable = "box", aspect = "equal")
plt.xlabel(L"x_1", fontsize=20)
plt.ylabel(L"x_2", fontsize=20)
plt.title(L"$\Phi$, full", fontsize=25)
show()
tight_layout();
savefig("full_scalar_flux_exp_Euler_nx$(s.NCellsX)_tEnd$(s.tEnd).png")

## explicit Euler
fig = figure("T",figsize=(10,10),dpi=100)
ax = gca()
im2 = pcolormesh(s.xMid,s.yMid, 4.0*pi*sqrt(2)*B1.^(1/4),vmin=0.0)
color_bar = fig.colorbar(im2, ax=ax, pad=0.03, shrink = 0.71)
ax.tick_params("both",labelsize=20)
ax.set(adjustable = "box", aspect = "equal")
plt.xlabel(L"x_1", fontsize=20)
plt.ylabel(L"x_2", fontsize=20)
plt.title(L"$T$, full", fontsize=25)
show()
tight_layout();
savefig("full_temperature_exp_Euler_nx$(s.NCellsX)_tEnd$(s.tEnd).png")


################### plot scalar fluxes DLRA ###################

## explicit Euler
fig = figure("scalar_flux_DLRA",figsize=(10,10),dpi=100)
ax = gca()
im1DLRA = pcolormesh(s.xMid,s.yMid, 4.0*pi*sqrt(2)*u2[:,:,1],vmin=0.0)
color_bar = fig.colorbar(im1DLRA, ax=ax, pad=0.03, shrink = 0.71)
ax.tick_params("both",labelsize=20) 
ax.set(adjustable = "box", aspect = "equal")
plt.xlabel(L"x_1", fontsize=20)
plt.ylabel(L"x_2", fontsize=20)
plt.title(L"$\Phi$, DLRA", fontsize=25)
show()
tight_layout();
savefig("DLRA_scalar_flux_exp_Euler_nx$(s.NCellsX)_tEnd$(s.tEnd).png")

## explicit Euler
fig = figure("T_DLRA",figsize=(10,10),dpi=100)
ax = gca()
im2DLRA = pcolormesh(s.xMid,s.yMid, 4.0*pi*sqrt(2)*B2.^(1/4),vmin=0.0)
color_bar = fig.colorbar(im2DLRA, ax=ax, pad=0.03, shrink = 0.71)
ax.tick_params("both",labelsize=20)
ax.set(adjustable = "box", aspect = "equal")
plt.xlabel(L"x_1", fontsize=20)
plt.ylabel(L"x_2", fontsize=20)
plt.title(L"$T$, DLRA", fontsize=25)
show()
tight_layout();
savefig("DLRA_temperature_exp_Euler_nx$(s.NCellsX)_tEnd$(s.tEnd).png")


################### plot rank for DLRA ###################

# plot rank in time
t = rankInTime[1,:];
fig, ax = subplots()
ax.plot(t,rankInTime[2,:], color="blue", linewidth=2, label="DLRA", alpha=0.9)
ax.set_xlim([0.0,t[end]])
ax.set_ylim([0.0,maximum(rankInTime[2,2:end])+0.1])
ax.tick_params("both",labelsize=15)
if s.tEnd == 1.5
    ax.set_xticks(range(0, 1.5, step=0.25))
end
xlabel(L"t",fontsize=15)
ylabel(L"$r$",fontsize=15)
tight_layout()
show()
savefig("rank_exp_Euler_nx$(s.NCellsX)_tEnd$(s.tEnd).png")


################### plot mass error ###################

t = mass[1,:];
massInit = massDLRA[2,1];
fig, ax = subplots()
ax.plot(t,abs.(mass[2,:] .- massInit)/massInit, color="red", linestyle = "--", linewidth=2, label="full", alpha=0.9)
ax.plot(t,abs.(massDLRA[2,:] .- massInit)/massInit, color="blue", linestyle = ":", linewidth=2, label="DLRA", alpha=0.9)
ax.legend(loc="upper left")
ax.set_xlim([0.0,t[end]])
ax.tick_params("both",labelsize=15)
if s.tEnd == 1.5
    ax.set_xticks(range(0, 1.5, step=0.25))
end
ax.set_yscale("log")
xlabel(L"t",fontsize=15)
ylabel(L"$\frac{|m^0-m^n|}{\Vert m^0\Vert}$",fontsize=15)
tight_layout()
show()
savefig("masserror_exp_Euler_nx$(s.NCellsX)_tEnd$(s.tEnd).png")


############# plot everything in two plots ################

fig, axs = subplots(2, 2, figsize=(12, 12), constrained_layout=true)

im1 = axs[1].pcolormesh(s.xMid,s.yMid, 4.0*pi*sqrt(2)*u1[:,:,1],vmin=0.0)
color_bar = fig.colorbar(im1, ax=axs[1], pad=-0.015, shrink = 0.73)
axs[1].tick_params("both",labelsize=15) 
axs[1].set(adjustable = "box", aspect = "equal")
axs[1].set_xticks(range(-1, 1, step=0.5))
axs[1].set_yticks(range(-1, 1, step=0.5))
axs[1].set_xlabel(L"x_1", fontsize=15)
axs[1].set_ylabel(L"x_2", fontsize=15)
axs[1].set_title(L"$\Phi$, full", fontsize=15)

im3 = axs[3].pcolormesh(s.xMid,s.yMid, 4.0*pi*sqrt(2)*u2[:,:,1],vmin=0.0)
color_bar = fig.colorbar(im3, ax=axs[3], pad=-0.015, shrink = 0.73)
axs[3].tick_params("both",labelsize=15) 
axs[3].set(adjustable = "box", aspect = "equal")
axs[3].set_xticks(range(-1, 1, step=0.5))
axs[3].set_yticks(range(-1, 1, step=0.5))
axs[3].set_xlabel(L"x_1", fontsize=15)
axs[3].set_ylabel(L"x_2", fontsize=15)
axs[3].set_title(L"$\Phi$, DLRA", fontsize=15)

im2 = axs[2].pcolormesh(s.xMid,s.yMid, 4.0*pi*sqrt(2)*B1.^(1/4),vmin=0.0)
color_bar = fig.colorbar(im2, ax=axs[2], pad=-0.001, shrink = 0.73)
axs[2].tick_params("both",labelsize=15)
axs[2].set(adjustable = "box", aspect = "equal")
axs[2].set_xticks(range(-1, 1, step=0.5))
axs[2].set_yticks(range(-1, 1, step=0.5))
axs[2].set_xlabel(L"x_1", fontsize=15)
axs[2].set_ylabel(L"x_2", fontsize=15)
axs[2].set_title(L"$T$, full", fontsize=15)

im4 = axs[4].pcolormesh(s.xMid,s.yMid, 4.0*pi*sqrt(2)*B2.^(1/4),vmin=0.0)
color_bar = fig.colorbar(im4, ax=axs[4], pad=-0.001, shrink = 0.73)
axs[4].tick_params("both",labelsize=15)
axs[4].set(adjustable = "box", aspect = "equal")
axs[4].set_xticks(range(-1, 1, step=0.5))
axs[4].set_yticks(range(-1, 1, step=0.5))
axs[4].set_xlabel(L"x_1", fontsize=15)
axs[4].set_ylabel(L"x_2", fontsize=15)
axs[4].set_title(L"$T$, DLRA", fontsize=15)

savefig("allplots_exp_Euler_nx$(s.NCellsX)_tEnd$(s.tEnd).png")
show()

########### second plot #############

fig, axs = subplots(1, 2, figsize=(12,4), constrained_layout=true)

### plot rank in time ###
t = rankInTime[1,:];
axs[1].plot(t,rankInTime[2,:], color="blue", linewidth=2, label="DLRA", alpha=0.9)
axs[1].set_xlim([0.0,t[end]])
axs[1].tick_params("both",labelsize=15)
if s.tEnd == 1.5
    axs[1].set_xticks(range(0, 1.5, step=0.25))
end
axs[1].set_xlabel(L"t",fontsize=15)
axs[1].set_ylabel(L"$r$",fontsize=15)

### plot mass error in time ###

t = mass[1,:];
massInit = massDLRA[2,1];

axs[2].plot(t,abs.(mass[2,:] .- massInit)/massInit, color="red", linestyle = "--", linewidth=2, label="full", alpha=0.9)
axs[2].plot(t,abs.(massDLRA[2,:] .- massInit)/massInit, color="blue", linestyle = ":", linewidth=2, label="DLRA", alpha=0.9)
axs[2].legend(loc="upper left")
axs[2].set_xlim([0.0,t[end]])
axs[2].tick_params("both",labelsize=15)
if s.tEnd == 1.5
    axs[2].set_xticks(range(0, 1.5, step=0.25))
end
axs[2].set_yscale("log")
axs[2].set_xlabel(L"t",fontsize=15)
axs[2].set_ylabel(L"$\frac{|m^0-m^n|}{\Vert m^0\Vert}$",fontsize=15)

savefig("allplots2_exp_Euler_nx$(s.NCellsX)_tEnd$(s.tEnd).png")
show()

println("main finished")

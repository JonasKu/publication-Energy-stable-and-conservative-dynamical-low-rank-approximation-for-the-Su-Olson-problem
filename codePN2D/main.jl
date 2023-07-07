using Base: Float64
include("utils.jl")
include("settings.jl")
include("SolverDLRA.jl")

using PyPlot
using NPZ
using DelimitedFiles

close("all")

s = Settings(); # create settings class with 501 x 501 spatial cells and a maximal rank of 100

################################################################
######################### execute code #########################
################################################################

################### run explicit Euler ###################
solver = SolverDLRA(s);
@time tEnd, uDLRA, BDLRA, rankInTime, massDLRA = SolveDLRA(solver);
@time tEnd, u, B, mass = Solve(solver);
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
savefig("DLRA_scalar_flux_exp_Euler_nx$(s.NCellsX)_eps$(s.epsAdapt)_tEnd$(s.tEnd).png")

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
savefig("DLRA_temperature_exp_Euler_nx$(s.NCellsX)_eps$(s.epsAdapt)_tEnd$(s.tEnd).png")


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
savefig("rank_exp_Euler_nx$(s.NCellsX)_eps$(s.epsAdapt)_tEnd$(s.tEnd).png")


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
savefig("masserror_exp_Euler_nx$(s.NCellsX)_eps$(s.epsAdapt)_tEnd$(s.tEnd).png")



################### save data #####################

#npzwrite("data/SolutionUFullNx$(s.NCellsX)Ny$(s.NCellsY)N$(s.nPN)tEnd$(s.tEnd).npy", u);
#npzwrite("data/SolutionBFullNx$(s.NCellsX)Ny$(s.NCellsY)N$(s.nPN)tEnd$(s.tEnd).npy", B);
#npzwrite("data/SolutionMassFullNx$(s.NCellsX)Ny$(s.NCellsY)N$(s.nPN)tEnd$(s.tEnd).npy", mass);
npzwrite("data/SolutionUDLRANx$(s.NCellsX)Ny$(s.NCellsY)N$(s.nPN)eps$(s.epsAdapt)tEnd$(s.tEnd).npy", uDLRA);
npzwrite("data/SolutionBDLRANx$(s.NCellsX)Ny$(s.NCellsY)N$(s.nPN)eps$(s.epsAdapt)tEnd$(s.tEnd).npy", BDLRA);
npzwrite("data/SolutionRankDLRANx$(s.NCellsX)Ny$(s.NCellsY)N$(s.nPN)eps$(s.epsAdapt)tEnd$(s.tEnd).npy", rankInTime);
npzwrite("data/SolutionMassDLRANx$(s.NCellsX)Ny$(s.NCellsY)N$(s.nPN)eps$(s.epsAdapt)tEnd$(s.tEnd).npy", massDLRA);

println("main finished")

include("settings.jl")
include("Solver.jl")

using DelimitedFiles
using NPZ
using PyPlot

# axes_grid1 = pyimport("mpl_toolkits.axes_grid1")

s = Settings();

### read data ###

if s.problem == "SuOlsonTestcase"
    f = npzread("data/SuOlsonSolutionfFullNx$(s.Nx)N$(s.nPN)tEnd$(s.tEnd).npy");
    u = npzread("data/SuOlsonSolutionUFullNx$(s.Nx)N$(s.nPN)tEnd$(s.tEnd).npy");
    B = npzread("data/SuOlsonSolutionBFullNx$(s.Nx)N$(s.nPN)tEnd$(s.tEnd).npy");
    mass = npzread("data/SuOlsonSolutionMassFullNx$(s.Nx)N$(s.nPN)tEnd$(s.tEnd).npy");
    fDLRA = npzread("data/SuOlsonSolutionfDLRANx$(s.Nx)N$(s.nPN)tEnd$(s.tEnd).npy");
    uDLRA = npzread("data/SuOlsonSolutionUDLRANx$(s.Nx)N$(s.nPN)tEnd$(s.tEnd).npy");
    BDLRA = npzread("data/SuOlsonSolutionBDLRANx$(s.Nx)N$(s.nPN)tEnd$(s.tEnd).npy");
    rankInTime = npzread("data/SuOlsonSolutionRankDLRANx$(s.Nx)N$(s.nPN)tEnd$(s.tEnd).npy");
    massDLRA = npzread("data/SuOlsonSolutionMassDLRANx$(s.Nx)N$(s.nPN)tEnd$(s.tEnd).npy");
else 
    f = npzread("data/SuOlsonSolutionfFullNx$(s.Nx)N$(s.nPN)tEnd$(s.tEnd).npy");
    u = npzread("data/PlaneSourceSolutionUFullNx$(s.Nx)N$(s.nPN)tEnd$(s.tEnd).npy");
    B = npzread("data/PlaneSourceSolutionBFullNx$(s.Nx)N$(s.nPN)tEnd$(s.tEnd).npy");
    mass = npzread("data/PlaneSourceSolutionMassFullNx$(s.Nx)N$(s.nPN)tEnd$(s.tEnd).npy");
    fDLRA = npzread("data/SuOlsonSolutionfDLRANx$(s.Nx)N$(s.nPN)tEnd$(s.tEnd).npy");
    uDLRA = npzread("data/PlaneSourceSolutionUDLRANx$(s.Nx)N$(s.nPN)tEnd$(s.tEnd).npy");
    BDLRA = npzread("data/PlaneSourceSolutionBDLRANx$(s.Nx)N$(s.nPN)tEnd$(s.tEnd).npy");
    rankInTime = npzread("data/PlaneSourceSolutionRankDLRANx$(s.Nx)N$(s.nPN)tEnd$(s.tEnd).npy");
    massDLRA = npzread("data/PlaneSourceSolutionMassDLRANx$(s.Nx)N$(s.nPN)tEnd$(s.tEnd).npy");
end


### plotting ###

### plot scalar flux ###
fig, ax = subplots()
ax.plot(s.xMid,u[:,1], color="red", linestyle = "--", linewidth=2, label="full", alpha=0.9)
ax.plot(s.xMid,uDLRA[:,1], color="blue", linestyle = ":", linewidth=2, label="DLRA", alpha=0.9)
ax.legend(loc="upper right")
if s.problem == "SuOlsonTestcase"
    ax.set_xlim([0.0,3.0]);
    ax.set_xticks(range(0, 3, step=0.5))
else
    ax.set_xlim([-10.0,10.0]);
    ax.set_xticks(range(-10, 10, step=5.0))
end
ax.tick_params("both",labelsize=15) 
xlabel(L"x",fontsize=15)
ylabel(L"$\Phi$",fontsize=15)
tight_layout()
show()
if s.problem == "SuOlsonTestcase"
    savefig("SuOlson_scalar_flux_nx$(s.Nx).png")
else 
    savefig("PlaneSource_scalar_flux_nx$(s.Nx).png")
end

### plot temperature ###
fig, ax = subplots()
ax.plot(s.xMid,B.^(1/4), color="red", linestyle = "--", linewidth=2, label="full", alpha=0.9)
ax.plot(s.xMid,BDLRA.^(1/4), color="blue", linestyle = ":", linewidth=2, label="DLRA", alpha=0.9)
ax.legend(loc="upper right")
if s.problem == "SuOlsonTestcase"
    ax.set_xlim([0.0,3.0]);
    ax.set_xticks(range(0, 3, step=0.5))
else
    ax.set_xlim([-10.0,10.0]);
    ax.set_xticks(range(-10, 10, step=5.0))
end
ax.tick_params("both",labelsize=15) 
xlabel(L"x",fontsize=15)
ylabel(L"$T$",fontsize=15)
tight_layout()
show()
if s.problem == "SuOlsonTestcase"
    savefig("SuOlson_temperature_nx$(s.Nx).png")
else 
    savefig("PlaneSource_temperature_nx$(s.Nx).png")
end

### plot function f ###
fig, ax = subplots()
im = pcolormesh(s.xMid, s.yMid, f')
color_bar = fig.colorbar(im, ax=ax, pad=0.03)
ax.tick_params("both",labelsize=15) 
plt.xlabel("x", fontsize=15)
plt.ylabel(L"\mu", fontsize=15)
plt.title("f full", fontsize=15)
show()
tight_layout();
if s.problem == "SuOlsonTestcase"
    savefig("SuOlson_f_nx$(s.Nx)_ny$(s.Ny).png")
else 
    savefig("PlaneSource_f_nx$(s.Nx)_ny$(s.Ny).png")
end

fig, ax = subplots()
imDLRA = pcolormesh(s.xMid, s.yMid, fDLRA')
color_bar = fig.colorbar(imDLRA, ax=ax, pad=0.03)
ax.tick_params("both",labelsize=15) 
plt.xlabel("x", fontsize=15)
plt.ylabel(L"\mu", fontsize=15)
plt.title("f DLRA", fontsize=15)
show()
tight_layout();
if s.problem == "SuOlsonTestcase"
    savefig("SuOlson_fDLRA_nx$(s.Nx)_ny$(s.Ny).png")
else 
    savefig("PlaneSource_fDLRA_nx$(s.Nx)_ny$(s.Ny).png")
end

### plot rank in time ###
t = rankInTime[1,:];
fig, ax = subplots()
ax.plot(t,rankInTime[2,:], color="blue", linewidth=2, label="DLRA", alpha=0.9)
ax.set_xlim([0.0,t[end]])
ax.tick_params("both",labelsize=15)
if s.problem == "SuOlsonTestcase"
    ax.set_yticks(range(0, 22, step=5.0))
else
    ax.set_yticks(range(0, 25, step=5.0))
end
if s.problem == "SuOlsonTestcase"
    ax.set_xticks(range(0, 3, step=0.5))
else
    ax.set_xticks(range(0,8, step=1.0))
end
xlabel(L"t",fontsize=15)
ylabel(L"$r$",fontsize=15)
tight_layout()
show()
if s.problem == "SuOlsonTestcase"
    savefig("SuOlson_rank_nx$(s.Nx).png")
else 
    savefig("PlaneSource_rank_nx$(s.Nx).png")
end

### plot mass error in time ###
t = mass[1,:];
massInit = mass[2,1];
fig, ax = subplots()
ax.plot(t,abs.(mass[2,:] .- massInit)/massInit, color="red", linestyle = "--", linewidth=2, label="full", alpha=0.9)
ax.plot(t,abs.(massDLRA[2,:] .- massInit)/massInit, color="blue", linestyle = ":", linewidth=2, label="DLRA", alpha=0.9)
ax.legend(loc="lower right")
ax.set_xlim([0.0,t[end]])
ax.tick_params("both",labelsize=15)
if s.problem == "SuOlsonTestcase"
    ax.set_xticks(range(0, 3, step=0.5))
else
    ax.set_xticks(range(0,8, step=1.0))
end
ax.set_yscale("log")
xlabel(L"t",fontsize=15)
ylabel(L"$\frac{|m^0-m^n|}{\Vert m^0\Vert}$",fontsize=15)
tight_layout()
show()
if s.problem == "SuOlsonTestcase"
    savefig("SuOlson_masserror_nx$(s.Nx).png")
else 
    savefig("PlaneSource_masserror_nx$(s.Nx).png")
end


### plot everything in one plot ###

if s.problem == "SuOlsonTestcase"
    fig, axs = subplots(3, 2, figsize=(12, 12), constrained_layout=true)
    # plt.subplots_adjust(hspace=0.4, wspace=0.3)

    im1 = axs[1].pcolormesh(s.xMid, s.yMid, f')
    color_bar = fig.colorbar(im1, ax=axs[1], pad=-0.05)
    axs[1].tick_params("both",labelsize=15) 
    axs[1].set_yticks(range(-1, 1, step=1.0), fontsize=15)
    axs[1].set_xticks(range(-10, 10, step=5.0), fontsize=15)
    axs[1].set_xlabel("x", fontsize=15)
    axs[1].set_ylabel(L"\mu", fontsize=15)
    axs[1].set_title("f full", fontsize=15)

    im4 = axs[4].pcolormesh(s.xMid, s.yMid, fDLRA')
    color_bar = fig.colorbar(im4, ax=axs[4], pad=-0.05)
    axs[4].tick_params("both",labelsize=15)
    axs[4].set_yticks(range(-1, 1, step=1.0), fontsize =15)
    axs[4].set_xticks(range(-10, 10, step=5.0), fontsize=15)
    axs[4].set_xlabel("x", fontsize=15)
    axs[4].set_ylabel(L"\mu", fontsize=15)
    axs[4].set_title("f DLRA", fontsize=15)

    axs[2].plot(s.xMid,u[:,1], color="red", linestyle = "--", linewidth=2, label="full", alpha=0.9)
    axs[2].plot(s.xMid,uDLRA[:,1], color="blue", linestyle = ":", linewidth=2, label="DLRA", alpha=0.9)
    axs[2].legend(loc="upper right")
    axs[2].set_xlim([0.0,3.0])
    axs[2].set_xticks(range(0, 3, step=0.5))
    axs[2].tick_params("both",labelsize=15)
    axs[2].set_xlabel(L"x",fontsize=15)
    axs[2].set_ylabel(L"$\Phi$",fontsize=15)

    axs[5].plot(s.xMid,B.^(1/4), color="red", linestyle = "--", linewidth=2, label="full", alpha=0.9)
    axs[5].plot(s.xMid,BDLRA.^(1/4), color="blue", linestyle = ":", linewidth=2, label="DLRA", alpha=0.9)
    axs[5].legend(loc="upper right")
    axs[5].set_xlim([0.0,3.0]);
    axs[5].set_xticks(range(0, 3, step=0.5))
    axs[5].tick_params("both",labelsize=15) 
    axs[5].set_xlabel(L"x",fontsize=15)
    axs[5].set_ylabel(L"$T$",fontsize=15)

    # plot rank in time
    t = rankInTime[1,:];
    axs[3].plot(t,rankInTime[2,:], color="blue", linewidth=2, label="DLRA", alpha=0.9)
    axs[3].set_xlim([0.0,t[end]])
    axs[3].tick_params("both",labelsize=15)
    axs[3].set_yticks(range(0, 22, step=5.0))
    axs[3].set_xticks(range(0,3, step=1.0))
    axs[3].set_xlabel(L"t",fontsize=15)
    axs[3].set_ylabel(L"$r$",fontsize=15)

    fig.delaxes(axs[3,2]) 

    savefig("SuOlson_allplots_nx$(s.Nx).png")
    
else

    fig, axs = subplots(3, 2, figsize=(12, 12), constrained_layout=true)

    im1 = axs[1].pcolormesh(s.xMid, s.yMid, f')
    color_bar = fig.colorbar(im1, ax=axs[1], pad=-0.001)
    axs[1].tick_params("both",labelsize=15) 
    axs[1].set_yticks(range(-1, 1, step=1.0), fontsize=15)
    axs[1].set_xticks(range(-10, 10, step=5.0),fontsize=15)
    axs[1].set_xlabel("x", fontsize=15)
    axs[1].set_ylabel(L"\mu", fontsize=15)
    axs[1].set_title("f full", fontsize=15)

    im4 = axs[4].pcolormesh(s.xMid, s.yMid, fDLRA')
    color_bar = fig.colorbar(im4, ax=axs[4], pad=-0.001)
    axs[4].tick_params("both",labelsize=15) 
    axs[4].tick_params("both",labelsize=15)
    axs[4].set_yticks(range(-1, 1, step=1.0), fontsize=15)
    axs[4].set_xticks(range(-10, 10, step=5.0),fontsize=15)
    axs[4].set_xlabel("x", fontsize=15)
    axs[4].set_ylabel(L"\mu", fontsize=15)
    axs[4].set_title("f DLRA", fontsize=15)

    axs[2].plot(s.xMid,u[:,1], color="red", linestyle = "--", linewidth=2, label="full", alpha=0.9)
    axs[2].plot(s.xMid,uDLRA[:,1], color="blue", linestyle = ":", linewidth=2, label="DLRA", alpha=0.9)
    axs[2].legend(loc="upper right")
    axs[2].set_xlim([-10.0,10.0]);
    axs[2].set_xticks(range(-10, 10, step=5.0))
    axs[2].tick_params("both",labelsize=15)
    axs[2].set_xlabel(L"x",fontsize=15)
    axs[2].set_ylabel(L"$\Phi$",fontsize=15)

    axs[5].plot(s.xMid,B.^(1/4), color="red", linestyle = "--", linewidth=2, label="full", alpha=0.9)
    axs[5].plot(s.xMid,BDLRA.^(1/4), color="blue", linestyle = ":", linewidth=2, label="DLRA", alpha=0.9)
    axs[5].legend(loc="upper right")
    axs[5].set_xlim([-10.0,10.0])
    axs[5].set_xticks(range(-10, 10, step=5.0))
    axs[5].tick_params("both",labelsize=15) 
    axs[5].set_xlabel(L"x",fontsize=15)
    axs[5].set_ylabel(L"$T$",fontsize=15)

    # plot rank in time
    t = rankInTime[1,:];
    axs[3].plot(t,rankInTime[2,:], color="blue", linewidth=2, label="DLRA", alpha=0.9)
    axs[3].set_xlim([0.0,t[end]])
    axs[3].tick_params("both",labelsize=15)
    axs[3].set_yticks(range(0, 25, step=5.0))
    axs[3].set_xticks(range(0,8, step=1.0))
    axs[3].set_xlabel(L"t",fontsize=15)
    axs[3].set_ylabel(L"$r$",fontsize=15)

    # plot mass error in time
    t = mass[1,:];
    massInit = mass[2,1];
    axs[6].plot(t,abs.(mass[2,:] .- massInit)/massInit, color="red", linestyle = "--", linewidth=2, label="full", alpha=0.9)
    axs[6].plot(t,abs.(massDLRA[2,:] .- massInit)/massInit, color="blue", linestyle = ":", linewidth=2, label="DLRA", alpha=0.9)
    axs[6].set_xlim([0.0,t[end]])
    axs[6].legend(loc="lower right")
    axs[6].tick_params("both",labelsize=15)
    axs[6].set_xticks(range(0,8, step=1.0))
    axs[6].set_yscale("log")
    axs[6].set_xlabel(L"t",fontsize=15)
    axs[6].set_ylabel(L"$\frac{|m^0-m^n|}{\Vert m^0\Vert}$",fontsize=15)

    savefig("PlaneSource_allplots_nx$(s.Nx).png")

end

# subplots_adjust(wspace=0.2, hspace=10.0)

show()

println("main finished")

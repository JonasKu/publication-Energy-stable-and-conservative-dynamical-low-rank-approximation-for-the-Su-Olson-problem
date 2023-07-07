include("settings.jl")
include("Solver.jl")

using DelimitedFiles
using NPZ
using PyPlot

close("all")

s = Settings();

############################
cfl = 1.0
solver = Solver(s);
@time tEnd, u,B, mass,f = Solve(solver);
@time tEnd, uDLRA, BDLRA, rankInTime, massDLRA, fDLRA = SolveDLRA(solver);

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

println("main finished")

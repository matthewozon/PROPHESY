## load the packages used in the estimation
# plotting
using PyPlot
using PyCall
rc("text", usetex=true)
rc("figure",max_open_warning=50)
color_array = ["tab:blue"; "tab:orange"; "tab:green"; "tab:red"; "tab:purple"; "tab:brown"; "tab:pink"; "tab:gray"; "tab:olive"; "tab:cyan"; "magenta"; "yellow"; "hotpink"; "darkmagenta"; "chartreuse"; "deepskyblue"; "navy"; "darkcyan"; "crimson"; "firebrick"]; 

# data manipulation (loading, writing, etc)
using Printf
using XPSfile
using DataFrames
using Query

# scientific package from the official Julia repositories
using LinearAlgebra
using StatsBase


# modeling XPS
using XPSpack # experiment model (geometry factor and cross section estimation)
using MINOTAUR

# inversion package
using XPSinv
using XPSsampling


SAVE_DATA = true
if SAVE_DATA
    using XLSX
end
PLOT_FIG  = true
PLOT_DATA = false

# folders where the data files are
data_folder = "../../../data/TK2/";
data_folderC1s = string(data_folder,"C1s/")

FLAG_CC = false
FLAG_COSO3 = !FLAG_CC

# flags
FLAG_PLOT = true;
FLAG_SAVE_PLOT = false;

# geometry setup
λe0 = 2.0e-3;         # reference penetration depth in μm
δr = 2.5e-3           # transition to vacuum layer thickness 
k0 = 10;              # compute the measurement model over a distance of k0*λe0
Nr = 101; # 201; #    # number of discretization points in the radial dimension
Nθ = 256;             # number of discretization points in the polar angle dimension
Ny = 256;             # number of discretization points in the cylinder axis dimension
L  = 50.0;            # height of the irradiated sample (the vertical extent of the beam is more like 20μm instead of 100μm)
μ0 = 10.0;            # radius of the cylinder
x0 = sqrt(3.0)*5000.0 # (x0,y0,z0) are the coordinates of the analyzer's apperture
y0 = 0.0;
z0 = 5000.0;
# smooth edge parameters
Δr_water = 0.25e-3
κ_transition = 5.0




# make sure that δr>dboundary;
# bulk depth and boundary distance
d0 = 2.0*1.92e-3-Δr_water;    # bulk
d0 = 10.0*1.92e-3-Δr_water;    # bulk
dboundary = 1.92e-3+Δr_water; # boundary 

if (δr<=dboundary)
    δr = dboundary
end



# spacial discretization 
r = collect(range(μ0-k0*λe0,μ0+δr,length=Nr));
θ0 = atan(x0,z0)
θ = collect(range(θ0-π/2.0,θ0+π/2.0,Nθ));
y = collect(range(-L/2.0,L/2.0,length=Ny));

# smooth edge cut 
N0 = findlast(r.-μ0.<=-d0); # bulk
NB = findlast(r.-μ0.<=dboundary); # boundary
N_trunc = NB-N0+1;

# standard deviation of the known values
σB     = 0.1;        

# standard deviation of the smoothness a priori
if FLAG_CC
    σd     = 0.5 # 1.0 # 0.1
else
    σd     = 0.5 # 0.1 # 0.01
end
cor_len_lowres = 2.5;

# amplitude of the communication mecanism for sampling the a posteriori model
if FLAG_CC
    σw     = 0.5e-2
else
    σw     = 1.0e-3
end
cor_len_sampling = N_trunc/5;
Ns      = 1000000;
Ns_burn =  100000;

# optimization 
τ0 = 1.0e1 # 
N_max_iter = 200000# 0#00; 
r_n_tol=0.000001; # 0.001
r_y_tol=0.000001; # 0.001
W_stop = collect(LinRange(1.0,10.0,N_trunc));
x00 = collect(LinRange(1.0,0.0,N_trunc))


# unit conversion constant (some of the quantities are in μm, some in L and some in Mbarn)
NA = 6.022e23;
κ_simple_units = 1.0e-37*NA; # simplified model
κ_units        = 1.0e-25*NA; # original model



# column names fo the fitted results
photon_sym = Symbol("Photon energy");
ph_flu_sym = Symbol("Photon flux");
bind_sym   = Symbol("Binding energy");
shift_sym  = Symbol("Peak shift");
gauss_sym  = Symbol("FWHM(G)");
loren_sym  = Symbol("FWHM(L)");
area_sym   = Symbol("Area");
σtot_sym   = Symbol("Sigma");
α_al_sym   = Symbol("Alignment");

# regular expression to navigate the files
regData    = r"Eph ?= ?"                  # pattern found in the data sheet's name that does not appear in other sheets
regFit     = r"[Ff]itt?ing(_| )?results?" # pattern found in fit sheet's name
regEph     = r"=[0-9]*eV"                 # pattern for reading the photon energy
sortBy     = Symbol("Photon energy")      # column name for the photon energy  (sort the data by increasing order of photon energy)
thenSortBy = Symbol("Binding energy")     # colmun name for the binding energy (for the data sharing the same photon energy, sort them by binding energy)


# alignement from C1s data
dt_load = @elapsed include("loadDataAndAlignment.jl")
println("Loading and computing the measurement models: ",dt_load," s")

# now, you can run the inversion with the data

# y_data_1 # SNR OK
# y_data_2 # SNR OK ish
# y_data_3 # SNR... well, there's no need to try this one
# σ_noise  # standard deviation of the normalized noise


##
## run the inversion using only y_data_1
##
include("profileReconstruction.jl")

##
## plot the results
##
if PLOT_FIG
    figure(figsize=[10, 5])
    ax1 = subplot(121)
    l_ρ, = plot(1.0e3*(r.-μ0),ρC1s_bulk*ρ_cp,color=color_array[1])
    σρ_est = [σB*ones(Cdouble,N0-1); stdρ_HI; σB*ones(Cdouble,Nr-NB)]
    μρ_est = ρ_cp[:];
    l_ρ_std =  fill_between(1.0e3*(r.-μ0),ρC1s_bulk*(μρ_est-σρ_est),ρC1s_bulk*(μρ_est+σρ_est),alpha=0.5,color=color_array[1])
    l_λ = Array{PyObject,1}(undef,Ndata)
    for i in 1:Ndata
        l_λ[i], = plot(-[λ_all[i]; λ_all[i]],[0.0; 1.0], color=color_array[i+1])
    end
    l_stoi, = plot(1.0e3*(r.-μ0),ρC1s_bulk*MINOTAUR.f_logistic.(r,μ0,Δr_water;A=1.0), color=color_array[6])
    # l_s, = plot(-[0.0; 0.0],[0.0; 1.0], color=color_array[Ndata+2]) # ;"sharp edge surface"
    # xlim(1.0e3*(r[1].-μ0),2.5e3*(r[end].-μ0))
    # xlim(-d0*1.0e3,dboundary*1.0e3)
    # xlim(-1.1maximum([d0*1.0e3; λ_all]),1.0e3δr)
    xlim(-1.1maximum([d0*1.0e3; λ_all; 5.0]),1.0e3δr)
    ylim(-0.01ρC1s_bulk,max(1.6ρC1s_bulk,1.1maximum(ρC1s_bulk*(ρ_est+stdρ_HI))))
    xlabel("distance [nm]",fontsize=14); 
    ylabel("C1s molar concentration [M]",fontsize=14) 
    xticks(fontsize=14); yticks(fontsize=14); 
    # ax1 = gca()
    ax1.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
    ax1.yaxis.offsetText.set_size(14)
    ax1.xaxis.offsetText.set_size(14)
    # legend(fontsize=12)
    legend([(l_ρ,l_ρ_std);l_λ;l_stoi],["estimates\$\\pm\\sigma\$ [M]"; string.("\$\\lambda_e\$ = ",floor.(100λ_all)/100," [nm]");"stoichiometric profile"],fontsize=12,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)

    ax2 = subplot(122)
    if false
        scatter(λ_all,ρC1s_bulk*y_data_1)
        ylim(0.9ρC1s_bulk*minimum(y_data_1),1.1ρC1s_bulk*maximum(y_data_1))
        xlabel("IMFP [nm]",fontsize=14); 
        ylabel("peak area [mol]",fontsize=14)
        xticks(fontsize=14); yticks(fontsize=14); 
        ax2.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
        ax2.yaxis.offsetText.set_size(14)
        ax2.xaxis.offsetText.set_size(14)
        legend(["data"],fontsize=12,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)
    else
        ρ_solvent = MINOTAUR.f_logistic.(r,μ0,Δr_water;A=1.0)
        # l_excess, = plot(1.0e3*(r[1:NB].-μ0),ρ_cp[1:NB]./ρ_solvent[1:NB],color=color_array[1])
        l_excess, = plot(1.0e3*(r.-μ0),ρ_cp./ρ_solvent,color=color_array[1])
        l_excess_std =  fill_between(1.0e3*(r.-μ0),(μρ_est-σρ_est)./ρ_solvent,(μρ_est+σρ_est)./ρ_solvent,alpha=0.5,color=color_array[1])
        xlim(-1.1maximum([d0*1.0e3; λ_all; 5.0]),1.0e3δr)
        # ylim(-1.0,400.0)
        ylim(-1.0,800.0)
        xlabel("distance [nm]",fontsize=14); 
        ylabel("C1s relative excess",fontsize=14) 
        xticks(fontsize=14); yticks(fontsize=14); 
        # ax1 = gca()
        ax1.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
        ax1.yaxis.offsetText.set_size(14)
        ax1.xaxis.offsetText.set_size(14)
        legend([(l_excess,l_excess_std)],["estimated excess\$\\pm\\sigma\$"],fontsize=12,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)
    end

    tight_layout(pad=0.2, w_pad=0.2, h_pad=0.2)
    ax1.text(-0.1, 0.97, "a)", transform=ax1.transAxes,fontsize=16)
    ax2.text(0.05, 0.1+0.75, replace(data_filesC1s[idx_file][1:end-5],"_"=>" "), transform=ax2.transAxes,fontsize=16)
    ax2.text(-0.1, 0.97, "b)", transform=ax2.transAxes,fontsize=16)

    if FLAG_CC
        savefig(string(data_filesC1s[idx_file][1:end-5],"_reconstruction_and_data_smooth_edge_CC.png"))
        savefig(string(data_filesC1s[idx_file][1:end-5],"_reconstruction_and_data_smooth_edge_CC.pdf"))
    else
        savefig(string(data_filesC1s[idx_file][1:end-5],"_reconstruction_and_data_smooth_edge_COSO3.png"))
        savefig(string(data_filesC1s[idx_file][1:end-5],"_reconstruction_and_data_smooth_edge_COSO3.pdf"))
    end

    figure()
    plot(deltaU)
end

if SAVE_DATA
    figure()
    # scatter(1.0e8α_al_noise,1.0e8(α_al_noise_est-α_al_noise))
    scatter(1.0e8α_al_noise,100*(α_al_noise_est-α_al_noise)./α_al_noise)
    xlabel("\$\\alpha_{\\mathrm{stoichiometric}}\$ [cm\$^{-2}\$]")
    ylabel("relative variation [\\%]")

    include("save_data.jl")
end
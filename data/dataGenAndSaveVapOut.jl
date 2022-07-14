## load the packages used in the estimation
# plotting
using PyPlot
rc("text", usetex=true)
rc("figure",max_open_warning=50)
using myPlot

# data manipulation (loading, writing, etc)
using Printf
using CSV
using DataFrames

# scientific package from the official Julia repositories
using LinearAlgebra
using StatsBase
using Interpolations

# implemented scientific packages
using utilsFun  # for the softMax functions

# modeling XPS
using XPSpack # experiment model (geometry factor and cross section estimation)
using ATTIRE  # kinetic energy analyzer



# tags
LOW_RES   = true               # set to true for computing the low resolution measurement models
MARG_UN   = (false & LOW_RES)    # set to true for computing the mean measurement operator as well as the covarainces
SHORT_RANGE = false              # select either wide range of attenuation lengths (false) or a restricted range more similar to experimental setup (true)
SIMULATE_DATA = true            # set true for simulating some data (several nose level)
SAVE_DATA = (false & SIMULATE_DATA)                # set to true for saving the generated variables
SAVE_MODEL = false               # set to true to save the models

MODEL_5   = false               # select the number of attenuation lengths probed
MODEL_10  = false
MODEL_20  = true

SAMPLE_MODEL = (false & LOW_RES) # if set to true, generate plenty of models by drawing attenuation lengths (only in low resolution, and does not compute the marginalization)
N_model_sample = 100;           # 100 should be more than enough for 5 attenuation length, but maybe not for 20

FLAG_0001 = (true & SIMULATE_DATA)               # selection of the profile (one must be true and the others false)
FLAG_0002 = (false & SIMULATE_DATA)
FLAG_0003 = (false & SIMULATE_DATA)
FLAG_0004 = (false & SIMULATE_DATA)

save_folder = "./";

# geometry setup
λe0 = 2.0e-3;        # reference penetration depth in μm
δr = 2.0e-3          # transition to vacuum layer thickness (let's set about 1 nm)
k0 = 10;             # compute the measurement model over a distance of k0*λe0
Nr = 201;            # number of discretization points in the radial dimension
Nr_lowres = 101 ;    # for the low resolution model
Nθ = 256;            # number of discretization points in the polar angle dimension
Ny = 256;            # number of discretization points in the cylinder axis dimension
μ0 = 10.0; # 20.0;   # radius of the cylinder
L = 100.0;           # height of the irradiated sample
x0 = sqrt(2.0)*100.0 # (x0,y0,z0) are the coordinates of the analyzer's apperture
y0 = 0.0;
z0 = 100.0;

# soon add the option of sphere or plane geometry?
save_folder = string(save_folder,"cylinder_radius_",μ0,"/")


# spacial discretization 
r = collect(range(μ0-k0*λe0,μ0+δr,length=Nr));
r_lowres = collect(range(μ0-k0*λe0,μ0+δr,length=Nr_lowres));
θ0 = atan(x0,z0)
θ = collect(range(θ0-π/2.0,θ0+π/2.0,Nθ));
y = collect(range(-L/2.0,L/2.0,length=Ny));

# concentration profiles (4 different cases)
ρ0 = 1.0
ρ_vac = 0.0
r_th  = 2.346; # 2.0


if FLAG_0001
    # ρA_1 = logistic.(1000.0reverse(μ0.-r).-r_th,ρ_vac,ρ0,2.0); # r_th=2.0
    ρA_1 = logistic.(1000.0reverse(δr.+μ0.-r).-r_th,ρ_vac,ρ0,2.0);
    exp_tag     = "0001"
end
if FLAG_0002
    # ρA_1 = logistic.(1000.0reverse(μ0.-r).-r_th,ρ_vac,ρ0,2.0) .+ 2.0exp.(-(1000.0reverse(μ0.-r).-1.0).^2. /(2.0*0.25^2)); # r_th=2.0
    ρA_1 = logistic.(1000.0reverse(δr.+μ0.-r).-r_th,ρ_vac,ρ0,2.0) .+ 2.0exp.(-(1000.0reverse(μ0.-r_surf).-0.0).^2. /(2.0*0.25^2));
    exp_tag     = "0002"
end
if FLAG_0003
    # ρA_1 = logistic.(1000.0reverse(μ0.-r).-r_th,ρ_vac,ρ0,2.0) .+ exp.(-(1000.0reverse(μ0.-r).-1.5).^2. /(2.0*0.5^2)); # r_th=2.0
    ρA_1 = logistic.(1000.0reverse(δr.+μ0.-r).-r_th,ρ_vac,ρ0,2.0) .+ exp.(-(1000.0reverse(μ0.-r_surf).-0.5).^2. /(2.0*0.5^2));
    exp_tag     = "0003"
end
if FLAG_0004
    # ρA_1 = exp.(-(1000.0reverse(μ0.-r).-2.5).^2. /(2.0*0.5^2));
    ρA_1 = exp.(-(1000.0(reverse(δr.+μ0.-r).-δr)).^2. /(2.0*0.5^2));
    exp_tag     = "0004"
end

# figure()
# plot(1000.0reverse(μ0.-r),ρA_1)

# measurement operator (only the geometrical term since the other comes as a multiplicative scalar estimated from the data)
if MODEL_5                                   # number of measurement (penetration depth)
    Ndata = 5;
end
if MODEL_10
    Ndata = 10;
end
if MODEL_20
    Ndata = 20;
end

if SHORT_RANGE                                          # bounds of the attenuation lengths
    λe1 = 1.3
    λe2 = 2.5
    save_folder = string(save_folder,"eal_",Ndata,"_restricted_range/")
else
    λe1 = 0.5;
    λe2 = 5.5; 
    save_folder = string(save_folder,"eal_",Ndata,"/")
end
λe = 1.0e-3collect(range(λe1,λe2,Ndata));              # attenuation length range

##
## compute model
##

H_highres = zeros(Cdouble,Ndata,Nr);                   # high resolution operator used for the data simulation (the attenuation length used for the computation of this operator are the reference ones)
H_lowres  = zeros(Cdouble,Ndata,Nr_lowres);
for i in 1:Ndata
    # high resolution measurement operator (using true λe values)
    H_highres[i,:],_,_,_,_ = cylinder_gain_H(r,θ,y,x0,y0,z0,μ0,λe[i]);
    H_lowres[i,:],_,_,_,_  = cylinder_gain_H(r_lowres,θ,y,x0,y0,z0,μ0,λe[i]);
end
H_highres = reverse(H_highres,dims=2);
H_lowres = reverse(H_lowres,dims=2);
if SAVE_MODEL
    mkpath(save_folder)
    CSV.write(string(save_folder,"radial_discretization.csv"),DataFrame(reverse(r)',:auto);header=true)
    CSV.write(string(save_folder,"attenuation_length.csv"),DataFrame(λe',:auto);header=false)
    CSV.write(string(save_folder,"H_highres.csv"),DataFrame(H_highres,:auto);header=true)
end


##
## geometry and concentration factor 
##

M_HCj = H_highres*ρA_1;

## 
## cross section density
## 

dKe = 0.05;
Be = collect(286.0:dKe:298.0);
Nspectrum = length(Be);
μBe = [290.2; 292.0; 293.0] # assume that the peaks are at the same location for every photon energy

σ_be = [0.45; 0.25; 0.6];

p_peak = zeros(Cdouble,Ndata,3);
p_peak[:,1] = 0.85 .+ (0.77-0.85)*(λe.-λe[1])./(λe[end]-λe[1]);
p_peak[:,2] = 0.125 .+ (0.12-0.125)*(λe.-λe[1])./(λe[end]-λe[1]);
p_peak[:,3] = 1.0 .- (p_peak[:,1]+p_peak[:,2]);
sum(p_peak,dims=2)

figure()
for i in 1:Ndata
    σ_peak_1 = (1.0/sqrt(2.0π*σ_be[1]^2))*exp.(-0.5*(Be.-μBe[1]).^2/(2.0σ_be[1]^2));
    σ_peak_2 = (1.0/sqrt(2.0π*σ_be[2]^2))*exp.(-0.5*(Be.-μBe[2]).^2/(2.0σ_be[2]^2));
    σ_peak_3 = (1.0/sqrt(2.0π*σ_be[3]^2))*exp.(-0.5*(Be.-μBe[3]).^2/(2.0σ_be[3]^2));
    plot(Be,p_peak[i,1]*σ_peak_1+p_peak[i,2]*σ_peak_2+p_peak[i,3]*σ_peak_3)
end


function σ_cs(hν::Cdouble,Ke::Cdouble,μKe::Cdouble;μKe0::Cdouble=50.0,μKe1::Cdouble=1200.0)
    Be = hν-reverse(Ke);
    # partial cross section (one for each chemical state)
    σ_peak_1 = (1.0/sqrt(2.0π*σ_be[1]^2))*exp.(-0.5*(Be-μBe[1])^2/(2.0σ_be[1]^2));
    σ_peak_2 = (1.0/sqrt(2.0π*σ_be[2]^2))*exp.(-0.5*(Be-μBe[2])^2/(2.0σ_be[2]^2));
    σ_peak_3 = (1.0/sqrt(2.0π*σ_be[3]^2))*exp.(-0.5*(Be-μBe[3])^2/(2.0σ_be[3]^2));
    # quantity of chemical states
    p1 = 0.85 .+ (0.77-0.85)*(μKe.-μKe0)./(μKe1-μKe0);
    p2 = 0.125 .+ (0.12-0.125)*(μKe.-μKe0)./(μKe1-μKe0);
    p3 = 1.0-(p1+p2);

    # cross section value (only for hν ∈ [295,1500.0])
    XPSpack.σ_C1s_interp[hν]*(p1*σ_peak_1+p2*σ_peak_2+p3*σ_peak_3)
end
#TODO: save cross section density and total cross section (check Yeah 1985 values)

##
## simple analyzer model
##

# for a given photon energy νj, measure a spectrum
hν = collect(LinRange(365.0,1500.0,Ndata));
dhν = hν.*((1.0/25000.0)*(hν.<500.0) + (1.0/15000.0)*(hν.>=500.0));
Fνj = 1.0e3*ones(Cdouble,Ndata);
j = 1 # Ndata
j = Ndata
Tj   = collect(LinRange(5.0,10.0,Ndata));       # one transmission factor per photon energy
# μKe  = collect(LinRange(50.0,1200.0,Ndata));    # one central kinetic energy per photon energy
σ_ke = 2.0*dKe*collect(LinRange(1.0,2.0,Ndata)); # one spread per pass energy
Keij = reverse(hν[j] .- Be) ;
μKe = 0.5*(Keij[1]+Keij[end]);


# φi  = Array{Function,1}(undef,Nspectrum);
# for i in 1:Nspectrum
#     φi[i] = (x::Cdouble->(Tj[j]/sqrt(2π*σ_ke[j]^2))*exp(-(x-Keij[i])^2/(2σ_ke[j]^2)))
# end

φi = Φi(Keij,σ_ke[j],Tj[j])

# TODO: plot some φi
figure()
for i in 10:20
    plot(Keij,φi[i].(Keij))
end

# TODO: compute a spectral spread of the beamline for a given photon energy (does not really matter at low energy, but becomes a problem at high energy)
# function F_density(hν::Cdouble,hνj::Cdouble,Δν::Cdouble)
#     (Fνj[j]/sqrt(2π*Δν^2))*exp(-(hν-hνj)^2/(2Δν^2))
# end

hνplot = collect(LinRange(hν[j]-3.0hν[j]/20000.0,hν[j]+3.0hν[j]/20000.0,11))
# F_density.(hνplot,hν[j],hν[j]/20000.0)


figure(); plot(hνplot,sourceSpread.(hνplot,hν[j],dhν[j],hν[j]))
# figure(); plot(hνplot,F_density.(hνplot,hν[j],hν[j]/20000.0))


# discretization of the spectral integrations: source and analyzer
Gm = dKe*ones(Cdouble,Nspectrum);
Gm[1] = 0.5dKe;
Gm[end] = 0.5dKe
hνjl = collect(hν[j]-5*dKe:dKe:hν[j]+5*dKe)
NdensF = length(hνjl)
Fl = dKe*ones(Cdouble,NdensF);
Fl[1] = 0.5dKe;
Fl[end] = 0.5dKe

# compute the signal measured
Aij = zeros(Cdouble,Nspectrum,NdensF); # discretization of the spread of the source and analyzer
Sj = Array{Cdouble,1}(undef,Nspectrum);
Ssignal = Array{Cdouble,2}(undef,Nspectrum,Ndata);
for i in 1:Ndata
    for m in 1:Nspectrum
        for l in 1:NdensF
            Aij[m,l] = φi[i](Keij[m])*sourceSpread(hνjl[l],hν[j],dhν[j],Fνj[j])*σ_cs(hνjl[l],Keij[m],μKe);
        end
    end
    Sj[i] = Gm'*Aij*Fl;
    Ssignal[i,j] = Sj[i]*M_HCj[j];
end

# plot spread of the signal due to the source and the measurement device
figure(); imshow(Aij); colorbar()

# compare the signals, without background
# figure(); plot(Keij,Sj)
# figure(); plot(Keij,σ_cs.(hν[j],Keij,μKe))
figure(); plot(Keij,Sj/Fνj[j]); plot(Keij,Tj[j]*σ_cs.(hν[j],Keij,μKe))


##
## add background and noise
##

# parameters used for the simulation of the inelastic background 
# Here, the background is made up of electrons that undergo inelastic collisions,
# at least one so that they don't appear in the sharp peak, but at lower kinetic energy
# the model is fairly simple and arbitrary, but it's better than no the background at all
BeC1s = mean(Be);
ΔBeC1s = (Be[end]-Be[1])/5;
SbgC1s = 0.5Keij./(1.0.+exp.((Keij.-(hν[j]-BeC1s))./ΔBeC1s)); # the inelastic collision background

# add the noise on the signal that hit the sensor, i.e. signal with background
SC1snoise = countElectrons(SbgC1s+reverse(Ssignal[:,j]))
SC1s = SC1snoise - SbgC1s; # noisy signal without background -> negative values appears

# plot signals w.r.t. the kinetic energy
figure(); plot(Keij,reverse(Ssignal[:,j])); scatter(Keij,SC1s)
figure(); plot(Keij,SbgC1s); plot(Keij,SbgC1s+reverse(Ssignal[:,j])); scatter(Keij,rand.(Poisson.(SbgC1s+reverse(Ssignal[:,j]))))
# plot the signal w.r.t. the binding energy
figure(); plot(Be,Ssignal[:,j]); scatter(Be,reverse(SC1snoise-SbgC1s))

# Kedumdum = collect(0.0:1.5maximum(Keij));
# Sdumdum = 900.0exp.(-0.5*(Kedumdum.-(hν[j]-BeC1s)).^2);
# figure(); plot(Kedumdum,0.2Kedumdum./(1.0.+exp.((Kedumdum.-(hν[j]-BeC1s))./ΔBeC1s))+Sdumdum)

"""
    (Fνj::Cdouble,hνj::Cdouble,Δνj::Cdouble): photon beam's parameters
    (Ki::Array{Cdouble,1},ΔKi::Cdouble,T::Cdouble): analyzer's parameters
    (r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0,y0,z0,μ0,λe): geometry factor (parameters)
    (hν::Array{Cdouble,1},Ki::Array{Cdouble,1},Be0::Array{Cdouble,1},σ_cs_0::Array{Cdouble,1}): cross section
    (ρ::Array{Cdouble,1}): concentration profile
"""
function simulateSpectrum(Fνj::Cdouble,hνj::Cdouble,Δνj::Cdouble,
    Ki::Array{Cdouble,1},ΔKi::Cdouble,T::Cdouble,
    Be0::Array{Cdouble,1},σ_cs_0::Array{Cdouble,1},
    r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,λe::Cdouble,
    ρ::Array{Cdouble,1})
    # efficiency functions of the analyzer
    Ndata = length(Ki); # number of readings in a spectrum
    φi = Φi(Ki,ΔKi,T);

    # photon energy discretization space
    nν = 5 # this should depend on the relative value ΔKi and Δνj
    Δhν = ΔKi;
    hνjl = collect(hνj-nν*Δhν:Δhν:hνj+nν*Δhν); # discretization of the photon energy space

    # number of discretization point in the kinetic energy space and in the photon energy space
    Nspectrum = length(Ki); # NOTE: the discretization of the kinetic enegy space does not have to coincide with the center of the channel, it can be whatever subdivision of that space
    NdensF = length(hνjl);
    # discretization of the integral: piecewise linear basis function
    Gm = dKe*ones(Cdouble,Nspectrum); Gm[1] = 0.5dKe; Gm[end] = 0.5dKe;
    Fl = Δhν*ones(Cdouble,NdensF); Fl[1] = 0.5Δhν; Fl[end] = 0.5Δhν;

    # define an interpolating tool whose nodes are Be0::Array{Cdouble,1},σ_cs_0::Array{Cdouble,1}
    σ_cs_interp = extrapolate(interpolate((Be0,), σ_cs_0, Gridded(Linear())),Line())
    # for each hν, compute the binding energies for Be = hν-reverse(Ki);
    Be = hνjl[1] .- reverse(Ki);
    σ_cs_interp[Be]

    # compute the "convolution" of the corss section by the spread of the light source and the spread of the kinetic energy analyzer
    Aij = zeros(Cdouble,Nspectrum,NdensF); # discretization of the spread of the source and analyzer
    Sj = Array{Cdouble,1}(undef,Nspectrum);
    # F_dens = sourceSpread.(hνjl,hνj,Δνj,Fνj)
    # σ_tot = σ_C1s_interp[hνjl] 
    # σ_val = Array{Cdouble,NdensF,Nspectrum}
    # for l in 1:NdensF
    #     # interpolator for the current photon energy
    #     Be = hνjl[l] .- reverse(Ki);
    #     σ_val[l,:] = σ_tot[l]*σ_cs_interp[Be] # no need to repeat that computation at each iteration of the i loop, once at first is enough
    # end
    # TODO: there's a lot to optimize overhere! the efficiency functions can be evaluated once at the beginning, the corss section density can be evaluated once in the beginning too (see lines above)
    for i in 1:Ndata 
        # φi_val = φi[i].(Ki)
        for m in 1:Nspectrum
            for l in 1:NdensF
                # total cross section for the given photon energy 
                σ_tot =  σ_C1s_interp[hνjl[l]] 
                # interpolator for the current photon energy
                Be = hνjl[l] .- reverse(Ki);
                σ_val = σ_tot*σ_cs_interp[Be[m]] # no need to repeat that computation at each iteration of the i loop, once at first is enough
                Aij[m,l] = φi[i](Ki[m])*sourceSpread(hνjl[l],hνj,Δνj,Fνj)*σ_val # σ_cs(hνjl[l],Keij[m],μKe);
            end
        end
        Sj[i] = Gm'*Aij*Fl;
    end

    # geometry factor
    H_highres,_,_,_,_ = cylinder_gain_H(r,θ,y,x0,y0,z0,μ0,λe);
    
    # compute geometry signal
    S_geom = H_highres*ρ;

    # compute total signal
    S_geom*Sj
end


# TODO: save cross section and image of the cross section
# TODO: add noise

if false
    save_folder_data = string(save_folder,exp_tag,"/")
    # generate some data (data point and covariance)
    Nnoise = 200;
    SNRmodulation = 1.0 .+ (λe.-λe[1])./(λe[end]-λe[1]);
    SNR_level = [5.0; 10.0; 50.0; 100.0; 500.0; 1000.0];
    Signal_level = abs.(H_highres*ρA_1);


    for k in 1:length(SNR_level)
        # figure()
        y_data = zeros(Cdouble,Nnoise,Ndata);
        ΓIsqrt = diagm(Signal_level./(SNR_level[k]*SNRmodulation));
        y_data[1,:] = H_highres*ρA_1;
        for i in 2:Nnoise
            y_data[i,:] = y_data[1,:] + ΓIsqrt*randn(Cdouble,Ndata);
            # scatter(collect(1:Ndata),y_data[i,:])
        end
        y_data[y_data.<0.0] = -y_data[y_data.<0.0]; # just making sure the integrated peak area are positive values
        if SAVE_DATA
            mkpath(string(save_folder_data,"/noise_level_",σ_level[k],"/"))
            CSV.write(string(save_folder_data,"/noise_level_",σ_level[k],"/","repeated_data.csv"),DataFrame([λe';y_data],:auto);header=true)
        end
    end
    if SAVE_DATA
        CSV.write(string(save_folder_data,"concentration_profile.csv"),DataFrame(ρA_1',:auto);header=true)
    end

end


# find the data files
folder_content = readdir(data_folderS2p);
list_match = match.(r"^S2p.*xlsx$",folder_content)
data_filesS2p = folder_content[list_match.!=nothing]

# S2p
# idx_file = 1; # BS: not enough data (3 data points, do not use for profile reconstruction)
# idx_file = 2; # OK
# idx_file = 3; # not so good, probably problem in the data... or maybe not, idk 
# idx_file = 4; # OK
idx_file = 5; # OK


fileName = data_filesS2p[idx_file][1:end-5];

ρS2p_out  = 0.0                                                                  # [M] concentration of C 1s outside of the sharp edge volume (at the outter boundary)
NS2p_per_SDS = 1.0
ρS2p_bulk = NS2p_per_SDS*1.0e-3parse(Cdouble,match(r"SDS_[0-9]*mM",fileName).match[5:end-2]); # [M] bulk concentration of C 1s

# dictAllData,df_fit,Ndata = dataAndFit_xlsx2df(string(data_folderS2p,fileName,".xlsx");         regData=regData, regFit=regFit, regEph=regEph, sortBy=sortBy, thenSortBy=thenSortBy);
dictAllData,df_fit,Ndata = dataAndFit_xlsx2df_missing(string(data_folderS2p,fileName,".xlsx"); regData=regData, regFit=regFit, regEph=regEph, sortBy=sortBy, thenSortBy=thenSortBy);

# list all the photon energies
df_Eph = select(df_fit,photon_sym) |> @unique() |> @orderby(_[photon_sym]) |> DataFrame
# put everything in nice boxes
dictPeak = Dict();
for i in df_Eph[!,photon_sym]
    local df = df_fit |> @filter(_[photon_sym]==i) |> DataFrame
    dictPeak[Symbol(string("hν_",i))] = df # (df,names(df))
end


# alignement
α_al_noise       = zeros(Cdouble,Ndata);
α_al_noise_sharp = zeros(Cdouble,Ndata);
α_al_ratio       = zeros(Cdouble,Ndata);
T                = zeros(Cdouble,Ndata);
Fν               = zeros(Cdouble,Ndata);
σ_S2p            = zeros(Cdouble,Ndata);
H_geom           = zeros(Cdouble,Ndata,Nr);
H_geom_sharp     = zeros(Cdouble,Ndata,Nr);
σ_data           = zeros(Cdouble,Ndata);
y_peak_1         = zeros(Cdouble,Ndata);
y_peak_2         = zeros(Cdouble,Ndata);
λ_all            = zeros(Cdouble,Ndata);

if FLAG_PLOT
    figure(figsize=[12, 10])
end
for i in 1:Ndata
    # get fitting data into arrays
    local plot_sym = Symbol(string("hν_",df_Eph[!,photon_sym][i]));
    # local Be       = dictAllData[plot_sym].Wavelength;
    local Be       = collect(skipmissing(dictAllData[plot_sym].Wavelength));
    local dKe      = median(abs.(Be[2:end]-Be[1:end-1]))
    local Npeaks, BePeak, σePeak, σePeakL, AePeak, σ_peak = curveFromFit(dictPeak[plot_sym],Be,true;bind_sym=bind_sym, shift_sym=shift_sym, gauss_sym=gauss_sym, area_sym=area_sym, loren_sym=loren_sym);

    # get peak areas
    y_peak_1[i] = dKe*sum(collect(skipmissing(dictAllData[plot_sym].Curve1))); # dKe*AePeak[1] -> C-C first peak 
    y_peak_2[i] = dKe*sum(collect(skipmissing(dictAllData[plot_sym].Curve2))); # dKe*AePeak[2] -> C-C second peak

    # plot data and fits
    if FLAG_PLOT
        local ax = subplot(2,2,i)
        plot(Be,collect(skipmissing(dictAllData[plot_sym].Raw_spectrum)),label="Data")
        # plot(Be,collect(skipmissing(dictAllData[plot_sym].Background)).+dKe*dropdims(AePeak'*σ_peak,dims=1),label="fitted spectra peaks")
        plot(Be,collect(skipmissing(dictAllData[plot_sym].Background.+dictAllData[plot_sym].Curve1.+dictAllData[plot_sym].Curve2)),label="fitted spectra curves") # S2p
        ylim(0.0)
        xlabel("binding energy [eV]",fontsize=14); 
        ylabel("spectrum [count]",fontsize=14) 
        xticks(fontsize=14); yticks(fontsize=14); 
        legend(fontsize=14)
        ax.invert_xaxis()
        ax.text(0.1, 0.5, string("Eph = ",df_Eph[!,photon_sym][i]," [eV]"), transform=ax.transAxes,fontsize=14)
    end
    
    # measurement model
    local λe = mean(1.0e-3dictPeak[plot_sym][!,:IMFP][1:2])
    local σ_all = collect(skipmissing(dictAllData[plot_sym].Curve1.+dictAllData[plot_sym].Curve2)); # S2p
    σ_all = σ_all/(dKe*sum(σ_all));
    # compute the geomtry factor
    H_geom_sharp[i,:],_,_,_,_ = cylinder_gain_H(r,θ,y,x0,y0,z0,μ0,λe);
    H_geom[i,:],_,_,_,_       = cylinder_gain_smooth_H(r,θ,y,x0,y0,z0,μ0,Δr_water,λe;κ=κ_transition,Nτ=20)
    # estimate the alignment parameter
    local Sbg                = collect(skipmissing(dictAllData[plot_sym].Background));
    local S_noisy            = collect(skipmissing(dictAllData[plot_sym].Raw_spectrum));
    local ρ_sharp            = ρS2p_bulk*ones(Cdouble,Nr)
    local ρ                  = ρS2p_bulk*MINOTAUR.f_logistic.(r,μ0,Δr_water;A=1.0); # for the estimation of the alignement parameter, assume a concentration profile similar to that of water 
    α_al_noise_sharp[i],_    = noiseAndParameterEstimation(σ_all,H_geom_sharp[i,:],S_noisy,Sbg,ρ_sharp)
    α_al_noise_sharp[i]      = α_al_noise_sharp[i]/(κ_units*σ_S2p_exp(convert(Cdouble,df_Eph[!,photon_sym][i]))*dictPeak[plot_sym][!,ph_flu_sym][1])
    α_al_noise[i],noise_data = noiseAndParameterEstimation(σ_all,H_geom[i,:],S_noisy,Sbg,ρ)
    α_al_noise[i]            = α_al_noise[i]/(κ_units*σ_S2p_exp(convert(Cdouble,df_Eph[!,photon_sym][i]))*dictPeak[plot_sym][!,ph_flu_sym][1])
    α_al_ratio[i]            = dictPeak[plot_sym][!,α_al_sym][1]
    σ_data[i]                = sqrt.(var(noise_data))
    Fν[i]                    = dictPeak[plot_sym][!,ph_flu_sym][1];
    T[i]                     = 1.0;
    σ_S2p[i]                 = σ_S2p_exp(convert(Cdouble,df_Eph[!,photon_sym][i]));
    λ_all[i]                 = 1.0e3λe # dictPeak[plot_sym][!,:IMFP][2];
end

y_data_1 = y_peak_1./(α_al_noise.*T.*Fν.*σ_S2p.*ρS2p_bulk*κ_units) 
y_data_2 = y_peak_2./(α_al_noise.*T.*Fν.*σ_S2p.*ρS2p_bulk*κ_units)
σ_noise  = σ_data./(α_al_noise.*T.*Fν.*σ_S2p.*ρS2p_bulk*κ_units)

y_S    = y_data_1 + y_data_2;

P1 = y_data_1./(y_data_1+y_data_2);
P2 = y_data_2./(y_data_1+y_data_2);

P_S    = P1+P2;

if FLAG_PLOT
    tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)
    if FLAG_SAVE_PLOT
        savefig(string(fileName,"_plot.pdf"))
        savefig(string(fileName,"_plot.png"))
    end
end

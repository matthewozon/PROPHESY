# find the data files
folder_content = readdir(data_folderC1s);
list_match = match.(r"^C1s.+xlsx$",folder_content)
data_filesC1s = folder_content[list_match.!=nothing]

# go through the list of file and load the data and compute the alignment parameter
α_noiseC1s = Dict();
α_ratioC1s = Dict();
for idx_file in 1:eachindex(data_filesC1s)
    fileName = data_filesC1s[idx_file][1:end-5];
    ρC1s_bulk = 1.0e-3parse(Cdouble,match(r"SDS_[0-9]*mM",fileName).match[5:end-2]); # extract concentration in mM
    dictAllData,df_fit,Ndata = dataAndFit_xlsx2df_missing(string(data_folderC1s,fileName,".xlsx"); regData=regData, regFit=regFit, regEph=regEph, sortBy=sortBy, thenSortBy=thenSortBy);

    # list all the photon energies
    df_Eph = select(df_fit,photon_sym) |> @unique() |> @orderby(_[photon_sym]) |> DataFrame
    # put everything in nice boxes
    dictPeak = Dict();
    for i in df_Eph[!,photon_sym]
        local df = df_fit |> @filter(_[photon_sym]==i) |> DataFrame
        dictPeak[Symbol(string("hν_",i))] = df # (df,names(df))
    end

    # alignement
    α_al_noise = zeros(Cdouble,Ndata);
    α_al_ratio = zeros(Cdouble,Ndata);
    global idx = 1

    if FLAG_PLOT
        figure(figsize=[12, 10])
    end
    for i in 1:Ndata
        # get fitting data into arrays
        local plot_sym  = Symbol(string("hν_",df_Eph[!,photon_sym][i]));
        local Be        = dictAllData[plot_sym].Wavelength;
        local Be_nomiss = collect(skipmissing(Be));
        # local Ke_mean   = convert(Float64,df_Eph[!,photon_sym][i]) + mean(Be_nomiss);
        # local T_mean    = T_r4000.(Ke_mean,E_pass);
        local dKe       = median(abs.(Be_nomiss[2:end]-Be_nomiss[1:end-1])) # this is an estimate of the channel bandwidth (FWHM)
        local N_avg     = dictAllData[plot_sym].averaged_from[1];
        local Δt        = (dictAllData[plot_sym].time_per_bin[1]*N_avg)/1000.0 # total integration time in s
        # local Npeaks, BePeak, σePeak, σePeakL, AePeak, σ_peak = curveFromFit(dictPeak[plot_sym],Be_nomiss,true;bind_sym=bind_sym, shift_sym=shift_sym, gauss_sym=gauss_sym, area_sym=area_sym, loren_sym=loren_sym);

        local raw_tot    = N_avg*collect(skipmissing(dictAllData[plot_sym].Raw_spectrum))/10.0
        local Curve1_tot = N_avg*collect(skipmissing(dictAllData[plot_sym].Curve1))/10.0
        local Curve2_tot = N_avg*collect(skipmissing(dictAllData[plot_sym].Curve2))/10.0
        local Curve3_tot = N_avg*collect(skipmissing(dictAllData[plot_sym].Curve3))/10.0
        local Bg_tot     = N_avg*collect(skipmissing(dictAllData[plot_sym].Background))/10.0
        
        # plot data and fits
        if FLAG_PLOT
            local ax = subplot(2,2,i)
            plot(Be_nomiss,raw_tot,label="Data")
            # plot(Be_nomiss,Bg_tot.+dKe*dropdims(AePeak'*σ_peak,dims=1),label="fitted spectra peaks")
            plot(Be_nomiss,Bg_tot.+Curve1_tot.+Curve2_tot.+Curve3_tot,label="fitted spectra curves")
            ylim(0.0)
            xlabel("binding energy [eV]",fontsize=14); 
            ylabel("spectrum [count]",fontsize=14) 
            xticks(fontsize=14); yticks(fontsize=14); 
            legend(fontsize=14)
            ax.invert_xaxis()
            ax.text(0.1, 0.5, string("Eph = ",df_Eph[!,photon_sym][i]," [eV]"), transform=ax.transAxes,fontsize=14)
        end
        
        # measurement model
        local λe = 1.0e-3dictPeak[plot_sym][!,:IMFP][2]
        local σ_all = Curve1_tot .+ Curve2_tot .+ Curve3_tot
        σ_all = σ_all/(dKe*sum(σ_all));
        # compute the geomtry factor
        local H_geom,H_rθy,Arn,Aθj,Ayk = cylinder_gain_H(r,θ,y,x0,y0,z0,μ0,λe);
        # estimate the alignment parameter
        local Sbg = Bg_tot # collect(skipmissing(dictAllData[plot_sym].Background))/10.0;
        local S_noisy = raw_tot # collect(skipmissing(dictAllData[plot_sym].Raw_spectrum))/10.0;
        local ρ = ρC1s_bulk*ones(Cdouble,Nr)
        # compute the mean value of the multiplicative factor τ_k
        α_al_noise[idx],_,_    = noiseAndParameterEstimation(σ_all,H_geom,S_noisy,Sbg,ρ)

        # divide by all known parameters to get the mean value αT_k in the right units
        α_al_noise[idx]      = α_al_noise[idx]/(κ_units*σ_C1s_exp(convert(Cdouble,df_Eph[!,photon_sym][i]))*dictPeak[plot_sym][!,ph_flu_sym][1]*dKe*Δt)
        α_al_ratio[idx]      = dictPeak[plot_sym][!,α_al_sym][1]
        idx = idx + 1;
    end
    if FLAG_PLOT
        tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)
        if FLAG_SAVE_PLOT
            savefig(string(fileName,"_plot.pdf"))
            savefig(string(fileName,"_plot.png"))
        end
    end
    
    # push data in dictionaries
    α_noiseC1s[Symbol(fileName)] = α_al_noise
    α_ratioC1s[Symbol(fileName)] = α_al_ratio
end


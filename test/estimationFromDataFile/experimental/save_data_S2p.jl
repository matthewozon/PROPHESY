dfSampleParam       = DataFrame("radius"=>μ0, "water_transition_length"=>Δr_water, "height"=>L,
                        "boundary_distance"=>dboundary, "bulk_distance"=>d0)
dfDiscParam         = DataFrame("distance_min"=>r[1],"distance_max"=>r[end],"distance_Ndisc"=>Nr,
                        "height_min"=>y[1],"height_max"=>y[end],"height_Ndisc"=>Ny,
                        "theta0_angle"=>θ0, "angle_Ndisc"=>Nθ)
dfTruncation        = DataFrame("first_index"=>N0,"first_depth"=>r[N0].-μ0,
                                "last_index"=>NB,"last_depth"=>r[NB].-μ0,
                                "Ntrunc"=>N_trunc)
dfAnalyzerPosition  = DataFrame("x0"=>x0,"y0"=>y0,"z0"=>z0)
dfReconstructionSetup = DataFrame("bulk_relative_conc_std"=>σB,"boundary_relative_conc_std"=>σB,
                                "bulk_conc"=>ρS2p_bulk,"out_of_sample_conc"=>ρS2p_out,
                                "regularization_amplitude"=>σd,"regularization_index_correlation_length"=>cor_len_lowres,
                                "uncertainty_sampling_amplitude"=>σw,"uncertainty_sampling_correlation_length"=>cor_len_sampling,
                                "uncertainty_sampling_N"=>Ns,"uncertainty_sampling_burni_in"=>Ns_burn)

dfOptimization      = DataFrame("algorithm"=>"Pock and Chambolle ALG2","τ0"=>τ0,"max_iter"=>N_max_iter,
                        "rel_tol_state"=>r_n_tol,"rel_tol_data"=>r_y_tol,"weight_stop"=>W_stop,"intial_state"=>x00,"distance"=>r[N0:NB])


Ke = [mean(dictPeak[Symbol(string("hν_",toto))][!,Symbol("Kinetic energy")][1:2]) for toto in df_Eph[!,photon_sym]]
peak_proba = P_S;
normalized_data = y_S;
saveFileName = string(fileName,"_S.xlsx")

dfDataAndMeta       = DataFrame("photon_energy"=>convert.(Float64,df_Eph[!,1]),"attenuation_length"=>λ_all,"kinetic_energy"=>Ke,
                        "total_cross_section"=>σ_S2p,"photon_flux"=>Fν,"transmition_value"=>T,"peak_proba"=>peak_proba,
                        "alignement_ratio"=>α_al_ratio, "alignement_sharp"=>α_al_noise_sharp,
                        "alignement_smooth"=>α_al_noise,"alignement_estimate"=>α_al_noise_est,"normalized_data"=>normalized_data)

dfFileName          = DataFrame("filename"=>fileName)
dfProfile           = DataFrame("distance"=>r,"profile"=>ρS2p_bulk*ρ_cp)
dfCovMat            = DataFrame(["distance"=>r[N0:NB];[string("distance_",string(i,pad=3)) => Γρ_HI[i-N0+1,:] for i in N0:NB]])
dfHsharp            = DataFrame([Symbol("conc_profile_APE")=>ρS2p_bulk*ones(Cdouble,Nr);[Symbol(string("hν_",df_Eph[!,photon_sym][i])) => H_geom_sharp[i,:] for i in 1:Ndata]])
dfHsmooth           = DataFrame([Symbol("conc_profile_APE")=>ρS2p_bulk*MINOTAUR.f_logistic.(r,μ0,Δr_water;A=1.0);[Symbol(string("hν_",df_Eph[!,photon_sym][i])) => H_geom[i,:] for i in 1:Ndata]])

XLSX.writetable(saveFileName, "profile reconstruction" => dfProfile, "profile covmat"=>dfCovMat,
    "sharp model"=>dfHsharp, "smooth model"=>dfHsmooth, "file name"=>dfFileName,"data and meta data"=>dfDataAndMeta,
    "sample parameters"=>dfSampleParam, "discretization parameters"=>dfDiscParam,"truncation parameters"=>dfTruncation,
    "analyzer coordinates"=>dfAnalyzerPosition, "reconstruction setup"=>dfReconstructionSetup, "optimization"=>dfOptimization)


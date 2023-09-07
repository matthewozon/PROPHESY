##
## inversion model
##

function profile_reconstruction(y_data::Array{Cdouble,1},σ_noise::Array{Cdouble,1},H_geom::Array{Cdouble,2},ρS2p_bulk::Cdouble,ρS2p_out::Cdouble,σB::Cdouble,N0::Int64,NB::Int64,N_trunc::Int64,Nr::Int64,τ0::Cdouble,N_max_iter::Int64,r_n_tol::Cdouble,r_y_tol::Cdouble,W_stop::Array{Cdouble,1},x00::Array{Cdouble,1})
    # # slice the model (3 terms: boundary, surface and bulk)
    H0 = H_geom[:,NB+1:end]; # boundary
    H_tilde = H_geom[:,N0:NB];
    Hb = H_geom[:,1:N0-1]; # bulk
    Hnot = [H0 Hb];
    Hnot1 = sum(Hnot;dims=2);


    # data correction
    Δy = dropdims(sum(Hb,dims=2),dims=2)*ρS2p_bulk/ρS2p_bulk; # in the bulk (deeper than 5 nm) if the data are normalized by the bulk concentration, then the normalized concentration in the bulk is 1, otherwise ρB
    δy = dropdims(sum(H0,dims=2),dims=2)*ρS2p_out/ρS2p_bulk;                          # outside the sample the concentration is nearly 0
    y_tilde = y_data-(Δy+δy);
    


    # regularization (smoothness: applied as sparsity in the second order difference)
    DN = D2nd(N_trunc+4);
    Db = DN[:,1:2];
    D_tilde = DN[:,3:end-2];
    D0 = DN[:,end-1:end];


    # correction of the regularization "data"
    Δyd = -dropdims(sum(Db,dims=2),dims=2)*ρS2p_bulk/ρS2p_bulk;
    δyd = -dropdims(sum(D0,dims=2),dims=2)*ρS2p_out/ρS2p_bulk;
    yd  = Δyd+δyd;


    # smoosh together the several part of the model into augmented operators
    Htrunc = [H_tilde; D_tilde];                                     # conditional to data and measurement model


    #
    # covariances: the crafed Bayesian models assumes that some covariance are known, i.e. measurement noise, smoothness, known values and measurement operator (the last only in the marginalized case)
    #

    # measurement noise covariance
    ΓI = diagm(σ_noise.^2);
    ΓItrunc = ΓI + σB^2*Hnot1*Hnot1';
    ΓIinv = inv(ΓItrunc);

    # covariance matrix for the a priori distribution (second order difference)
    Γprior = zeros(Cdouble,Nr,Nr);
    for i in 1:Nr
        Γprior[i,i] =  1.0
        for j in i+1:Nr
            Γprior[i,j] = Γprior[i,i]*exp(-(i-j)^2/(0.5*cor_len_lowres^2));
            Γprior[j,i] = Γprior[i,j];
        end
    end

    Γd = (N_trunc/Ndata)*(σd^2)*Γprior[N0-1:NB+1,N0-1:NB+1]
    Γd_inv = inv(Γd);


    ##
    ## reconstruction
    ##
    dt_reconstruction = @elapsed ρ_est,_,N_last = alg2_cp_quad_LM(x00,y_tilde,yd,Htrunc,ΓItrunc,Γd,W_stop;τ0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);

    println("profile reconstruction, dt = ", dt_reconstruction)

    ##
    ## sampling the a posteriori model
    ##
    w = σw*ones(Cdouble,N_trunc); 
    Γsqrt = real(sqrt(corrCovariance(w;cor_len=cor_len_sampling))); 
    deltaU = zeros(Cdouble,Ns);

    # dt_sampling = @elapsed global μρ_HI,Γρ_HI,deltaU = samplePosteriorMeanAndCov(ρ_cp[N0:NB],Γsqrt,y_tilde,yd,ΓIinv,Γd_inv,H_tilde,D_tilde;Ns=Ns,Nburn=Ns_burn);
    dt_sampling = @elapsed global μρ_HI,Γρ_HI,deltaU = samplePosteriorMeanAndCov(ρ_est,Γsqrt,y_tilde,yd,ΓIinv,Γd_inv,H_tilde,D_tilde;Ns=Ns,Nburn=Ns_burn);

    println("sampling a posteriori model, dt = ", dt_sampling)


    ρ_est,N_last,Γρ_HI # ,Γd_inv,y_tilde,yd,ΓIinv
end


y_data = y_S;
ρ_est,N_last,Γρ_HI = profile_reconstruction(y_data,σ_noise,H_geom,ρS2p_bulk,ρS2p_out,σB,N0,NB,N_trunc,Nr,τ0,N_max_iter,r_n_tol,r_y_tol,W_stop,x00); 

ρ_cp = [(ρS2p_bulk/ρS2p_bulk)*ones(Cdouble,N0-1); ρ_est; (ρS2p_out/ρS2p_bulk)*ones(Cdouble,Nr-NB)];
stdρ_HI = sqrt.(diag(Γρ_HI));




# estimate the alignment parameter using the estimated concentration profile 
α_al_noise_est = zeros(Cdouble,Ndata);
for i in 1:Ndata
    local plot_sym = Symbol(string("hν_",df_Eph[!,photon_sym][i]));
    local Be       = collect(skipmissing(dictAllData[plot_sym].Wavelength));
    local dKe      = median(abs.(Be[2:end]-Be[1:end-1]))
    local σ_all = collect(skipmissing(dictAllData[plot_sym].Curve1.+dictAllData[plot_sym].Curve2)); # S2p
    σ_all = σ_all/(dKe*sum(σ_all));
    local Sbg                = collect(skipmissing(dictAllData[plot_sym].Background));
    local S_noisy            = collect(skipmissing(dictAllData[plot_sym].Raw_spectrum));
    α_al_noise_est[i],_    = noiseAndParameterEstimation(σ_all,H_geom[i,:],S_noisy,Sbg,ρS2p_bulk*ρ_cp)
    α_al_noise_est[i]      = α_al_noise_est[i]/(κ_units*σ_S2p_exp(convert(Cdouble,df_Eph[!,photon_sym][i]))*dictPeak[plot_sym][!,ph_flu_sym][1])
end


# if the first estimate of the alignment parameter is too bad, loop until the alignement parameter becomes somewhat acceptable
if (false & (sqrt(sum((α_al_noise-α_al_noise_est).^2)) /sqrt(sum(α_al_noise.^2))>0.5))

    for _ in 1:5
        # update data and noise 
        y_tmp = y_data.*α_al_noise./α_al_noise_est
        σ_tmp = σ_noise.*α_al_noise./α_al_noise_est

        # reconstruction using the updated data and noise 
        global ρ_est,N_last,Γρ_HI = profile_reconstruction(y_tmp,σ_tmp,H_geom,ρS2p_bulk,ρS2p_out,σB,N0,NB,N_trunc,Nr,τ0,N_max_iter,r_n_tol,r_y_tol,W_stop,x00); 
        
        # update profile 
        global ρ_cp = [(ρS2p_bulk/ρS2p_bulk)*ones(Cdouble,N0-1); ρ_est; (ρS2p_out/ρS2p_bulk)*ones(Cdouble,Nr-NB)];
        
        # update alignement parameter
        # global α_al_noise_est = zeros(Cdouble,Ndata);
        for i in 1:Ndata
            local plot_sym = Symbol(string("hν_",df_Eph[!,photon_sym][i]));
            local Be       = collect(skipmissing(dictAllData[plot_sym].Wavelength));
            local dKe      = median(abs.(Be[2:end]-Be[1:end-1]))
            local σ_all = collect(skipmissing(dictAllData[plot_sym].Curve1.+dictAllData[plot_sym].Curve2)); # S2p
            σ_all = σ_all/(dKe*sum(σ_all));
            local Sbg                = collect(skipmissing(dictAllData[plot_sym].Background));
            local S_noisy            = collect(skipmissing(dictAllData[plot_sym].Raw_spectrum));
            α_al_noise_est[i],_    = noiseAndParameterEstimation(σ_all,H_geom[i,:],S_noisy,Sbg,ρS2p_bulk*ρ_cp)
            α_al_noise_est[i]      = α_al_noise_est[i]/(κ_units*σ_S2p_exp(convert(Cdouble,df_Eph[!,photon_sym][i]))*dictPeak[plot_sym][!,ph_flu_sym][1])
        end
    end

    stdρ_HI = sqrt.(diag(Γρ_HI));
end 




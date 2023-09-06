##
## inversion model
##


# # slice the model (3 terms: boundary, surface and bulk)
H0 = H_geom[:,NB+1:end]; # boundary
H_tilde = H_geom[:,N0:NB];
Hb = H_geom[:,1:N0-1]; # bulk
Hnot = [H0 Hb];
Hnot1 = sum(Hnot;dims=2);


# data correction
Δy = dropdims(sum(Hb,dims=2),dims=2)*ρC1s_bulk/ρC1s_bulk; # in the bulk (deeper than 5 nm) if the data are normalized by the bulk concentration, then the normalized concentration in the bulk is 1, otherwise ρB
δy = dropdims(sum(H0,dims=2),dims=2)*ρC1s_out/ρC1s_bulk;                          # outside the sample the concentration is nearly 0
if FLAG_CC
    y_tilde = y_CC-(Δy+δy);
else
    y_tilde = y_COSO3-(Δy+δy);
end



# regularization (smoothness: applied as sparsity in the second order difference)
DN = D2nd(N_trunc+4);
Db = DN[:,1:2];
D_tilde = DN[:,3:end-2];
D0 = DN[:,end-1:end];


# correction of the regularization "data"
Δyd = -dropdims(sum(Db,dims=2),dims=2)*ρC1s_bulk/ρC1s_bulk;
δyd = -dropdims(sum(D0,dims=2),dims=2)*ρC1s_out/ρC1s_bulk;
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
N = N_trunc


##
## reconstruction
##
W_stop = collect(LinRange(1.0,10.0,N));
τ0 = 1.0e1 # 
x00 = 0.5ones(Cdouble,N); # since the concentration is normalized by the bulk concentration, the initial state is taken as uniform with value 1/2
N_max_iter = 200000# 0#00; 
r_n_tol=0.000001; # 0.001
r_y_tol=0.000001; # 0.001
ρ_cp    = zeros(Cdouble,Nr);
μρ_HI = zeros(Cdouble,N);
Γρ_HI = zeros(Cdouble,N,N);

ρ_est,_,N_last = alg2_cp_quad_LM(x00,y_tilde,yd,Htrunc,ΓItrunc,Γd,W_stop;τ0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);

# ρ_cp = [(ρC1s_bulk/ρC1s_bulk)*ones(Cdouble,N0-1); ρ_est; ρC1s_out/ρC1s_bulk];
dt_reconstruction = @elapsed ρ_cp = [(ρC1s_bulk/ρC1s_bulk)*ones(Cdouble,N0-1); ρ_est; (ρC1s_out/ρC1s_bulk)*ones(Cdouble,Nr-NB)];

println("profile reconstruction, dt = ", dt_reconstruction)

##
## sampling the a posteriori model
##
if SAMPLING
    w = σw*ones(Cdouble,N); 
    Γsqrt = real(sqrt(corrCovariance(w;cor_len=20.0))); # 5.0
    Ns      = 1000000;
    Ns_burn =  100000;
    deltaU = zeros(Cdouble,Ns);

    dt_sampling = @elapsed global μρ_HI,Γρ_HI,deltaU = samplePosteriorMeanAndCov(ρ_cp[N0:NB],Γsqrt,y_tilde,yd,ΓIinv,Γd_inv,H_tilde,D_tilde;Ns=Ns,Nburn=Ns_burn);

    println("sampling a posteriori model, dt = ", dt_sampling)

    global stdρ_HI = sqrt.(diag(Γρ_HI));
end



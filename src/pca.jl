using MultivariateStats
using Interpolations

function interp_one_obs_order_for_pca_1d(data::EOST, log_lambda_out::AbstractRange{T1}; order::Integer, obs::Integer,
      output::AbstractArray{T2,1} = zeros(eltype(data.flux),length(log_lambda_out)) ,
      z_bc::Real = zero(T1),
      blaze::Poly{T4} = Poly([one(T1)])
      ) where { EOST <: AbstractEchelleObservationSet, T1<:Real, T2<:Real, T4<:Real}
    @assert length(log_lambda_out) == length(output)
    @assert 1<=order<=num_orders(data)
    @assert 1<=obs<=num_obs(data)
    pixel_range = 1:size(data.flux,1)
    interpolator = extrapolate(
            interpolate( (log.(lambda_ssb_from_obs(data.lambda[:,order,obs],z_bc)), ),
            data.flux[:,order,obs]./polyval(blaze,pixel_range), Gridded(Linear()) ),
                  Flat() )
    output .= interpolator(log_lambda_out)
    return output
end

function interp_obs_order_for_pca_1d(data::EOST, log_lambda_out::AbstractRange{T1}; order::Integer,
      output::AbstractArray{T2,2} = zeros(eltype(data.flux),(length(log_lambda_out),num_obs(data))),
      z_bc::AbstractVector{T3} = zeros(T1,num_obs(data)),
      blaze::Poly{T4} = Poly([one(T1)])
      ) where { EOST <: AbstractEchelleObservationSet, T1<:Real, T2<:Real, T3<:Real, T4<:Real}
   obs_range = 1:num_obs(data)
   map(i->interp_one_obs_order_for_pca_1d(data,log_lambda_out,order=order,obs=i,
      output=view(output,:,i), z_bc=z_bc[i], blaze=blaze), obs_range)
   return output
end


function fit_pca_obs_order_interpolated(data::EOST, log_lambda_out::AbstractRange{T1}; order::Integer,
         z_bc::AbstractVector{T3} = zeros(T1,num_obs(data)),
         maxoutdim::Integer=10, pratio::Real=1-1e-6
         ) where { EOST <: AbstractEchelleObservationSet, T1<:Real, T3<:Real }
   obs_range = 1:num_obs(data)
   obs_interp = interp_obs_order_for_pca_1d(data, log_lambda_out, order=order, z_bc=z_bc)
   pca = fit(PCA,obs_interp, maxoutdim=10,pratio=1-1e-6)
   return pca
end

#=
data_for_pca = interp_obs_order_for_pca_1d(data,res2d[1].log_lambda_range, order=1)
pca1 = fit(PCA,data_for_pca, maxoutdim=10,pratio=1-1e-6)
principalvars(pca1)./tvar(pca1)
=#

#=
function interp_obs_for_pca_1d(data::EOST, log_lambda_out::AbstractRange{T1}; order::Integer, obs::Integer,
    output::AbstractArray{T2,1} = zeros(eltype(data.flux),(length(log_lambda_out),num_obs(data) ) ) ,
    z_bc::AbstractVector{T3} = zeros(T1,num_obs(data)),
    pixel_range::AbstractRange{T4} = 1:size(data.lambda,1),
    blaze::Poly{T5} = Poly([one(T2)])
    ) where { EOST <: AbstractEchelleObservationSet, T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Real}
=#

using Statistics     # For computing sum/mean/median over pixels

# Size of data
function num_pixels_per_order(data::EOT) where {EOT <: AbstractEchelleObservation}
   size(data.flux,1)
end
function num_orders(data::EOT) where {EOT <: AbstractEchelleObservation}
   size(data.flux,2)
end

function num_pixels_per_order(data::EOST) where {EOST <: AbstractEchelleObservationSet}
   size(data.flux,1)
end
function num_orders(data::EOST) where {EOST <: AbstractEchelleObservationSet}
   size(data.flux,2)
end

function num_obs(data::EOST) where {EOST <: AbstractEchelleObservationSet}
   size(data.flux,3)
end

# summary statistics for data
function flux_mean_over_pixels(data::AbstractEchelleObservationSet)
   reshape( sum(data.flux./data.var,dims=3)./sum(one(eltype(data.var))./data.var,dims=3),
            (size(data.flux,1),size(data.flux,2)) )
end

function flux_median_over_pixels(data::AbstractEchelleObservationSet)
   reshape(median(data.flux,dims=3),(size(data.flux,1),size(data.flux,2)))
end

function lambda_mean_over_pixels(data::AbstractEchelleObservationSet)
   reshape(mean(data.lambda,dims=3),(size(data.lambda,1),size(data.lambda,2)))
end

function lambda_median_over_pixels(data::AbstractEchelleObservationSet)
   reshape(median(data.lambda,dims=3),(size(data.lambda,1),size(data.lambda,2)))
end

# Default templates for star and tellurics
""" default_telluric_template(n; num_basis_vectors = 0)
   Returns array with (num_basis_vectors, n) zeros. """
function default_telluric_template(n::Integer; num_basis_vectors::Integer = 0)
   @assert 0 <= min_num_telluric_basis_vectors <= num_basis_vectors <= max_num_telluric_basis_vectors
   @assert 1 <= n <= max_template_length
   zeros(num_basis_vectors,n)
end

""" default_stellar_template(n; num_basis_vectors = 0)
   Returns array with one row of n 1's followed by num_basis_vectors rows of zeros. """
function default_stellar_template(n::Integer; num_basis_vectors::Integer = 0)
   @assert 0 <= min_num_star_basis_vectors-1 <= num_basis_vectors <= max_num_star_basis_vectors
   @assert 1 <= n <= max_template_length
   result = vcat(ones(1,n),zeros(num_basis_vectors,n))
   return result
end

# Doppler shift wavelengths
"""Calculate lambda in solar system barycenter frame given observed lambda and barycentric correction)."""
function lambda_ssb_from_obs(lambda_obs, z_bc)
   lambda_obs * (1+z_bc)
end

"""Calculate log lambda in solar system barycenter frame given observed log lambda and barycentric correction)."""
function log_lambda_ssb_from_obs(log_lambda_obs, z_bc)
   log_lambda_obs + log1p(z_bc)
end

"""Calculate lambda in source frame given observed lambda, barycentric correction and source z relative to solar system barycenter)."""
function lambda_source(lambda_obs, z_bc; z_source = zero(z_bc) )
   lambda_obs * (1+z_bc) / (1+z_source) # TODO: Check EQNs
end

"""Calculate lambda in source frame given observed lambda, barycentric correction and source z relative to solar system barycenter)."""
function log_lambda_source(log_lambda_obs, z_bc; z_source = zero(z_bc) )
   log_lambda_obs + log1p(z_bc) - log1p(z_source) # TODO: Check EQNs
end

function make_model_wavelength_range(data::EOT, order::OrdersType;
            oversample_factor::Real = 1.0,
            num_knots::Integer=round(Int64,oversample_factor*num_pixels_per_order(data))
            ) where { EOT <: AbstractEchelleObservation }
      range(log(minimum(data.lambda[1,order])), stop=log(maximum(data.lambda[end,order])), length=num_knots )
end

function make_model_wavelength_range(data::EOST, order::OrdersType;
            oversample_factor::Real = 1.0,
            num_knots::Integer=round(Int64,oversample_factor*num_pixels_per_order(data))
            ) where { EOST <: AbstractEchelleObservationSet }
      range(log(minimum(data.lambda[1,order,:])), stop=log(maximum(data.lambda[end,order,:])), length=num_knots )
end

function make_model_wavelength_range(data::EOST;
            oversample_factor::Real = 1.0,
            num_knots::Integer=round(Int64,oversample_factor*num_orders(data)*num_pixels_per_order(data))
            ) where { EOST <: AbstractEchelleObservationSet }
      range( log(minimum(data.lambda[1,1,:])),
             stop = log(maximum(data.lambda[end,end,:])),
             length = num_knots )
end

# Structure to define file format (works for HPF, see if need to generalize for EXPRES/HARPSN)
struct ManifestFormatEntry <: AbstractManifestFormatEntry
    name::Symbol
    type::Type
    cols::UnitRange{Int64}
end

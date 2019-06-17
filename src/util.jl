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
export lambda_ssb_from_obs, log_lambda_ssb_from_obs, lambda_source, log_lambda_source

#----

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
export num_pixels_per_order, num_orders, num_obs

#---

# Utility functions, so can work with either an EcheelObservation or an EcheelObservationSet
function get_min_log_lambda(data::T, pixel_range::AbstractRange = 1:num_pixels_per_order(data)
      ) where T<:AbstractEchelleObservationSet
   pixel_lo = max(1,pixel_range[1] - 1)
   log_lambda_lo = log(minimum(data.lambda[pixel_lo,order,:]))
end
function get_max_log_lambda(data::T, pixel_range::AbstractRange = 1:num_pixels_per_order(data)
      ) where T<:AbstractEchelleObservationSet
   pixel_hi = min(num_pixels_per_order(data), pixel_range[end]+1)
   log_lambda_hi = log(maximum(data.lambda[pixel_hi,order,:]))
end

function get_min_log_lambda(data::T, pixel_range::AbstractRange = 1:num_pixels_per_order(data)
      ) where T<:AbstractEchelleObservation
   pixel_lo = max(1,pixel_range[1] - 1)
   log_lambda_lo = log(minimum(data.lambda[pixel_lo,order]))
end
function get_max_log_lambda(data::T, pixel_range::AbstractRange = 1:num_pixels_per_order(data)
      ) where T<:AbstractEchelleObservation
   pixel_hi = min(num_pixels_per_order(data), pixel_range[end]+1)
   log_lambda_hi = log(maximum(data.lambda[pixel_hi,order]))
end
export get_min_log_lambda, get_max_log_lambda

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
export flux_mean_over_pixels, lambda_mean_over_pixels
export flux_median_over_pixels, lambda_median_over_pixels

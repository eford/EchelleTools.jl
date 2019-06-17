@reexport module Blaze
#using Reexport

import ..AbstractEchelleObservation, ..AbstractEchelleObservationSet

abstract type AbstractBlazeType end
abstract type AbstractBlaze1DType <: AbstractBlazeType end
abstract type AbstractBlaze2DType <: AbstractBlazeType end
export AbstractBlazeType, AbstractBlaze1DType, AbstractBlaze2DType

# Default functions for AbstractBlazeTypes
min_log_lambda(blaze::AbstractBlazeType) = blaze.log_lambda_lo
max_log_lambda(blaze::AbstractBlazeType) = blaze.log_lambda_hi
order(blaze::AbstractBlazeType) = blaze.order
pixel_range(blaze::AbstractBlazeType) = blaze.pixel_range
export order, min_log_lambda, max_log_lambda, pixel_range

# Constants from calling module
import ..default_blaze_degree# , ..max_blaze_degree
import ..max_num_pixels_per_order
import ..default_pixel_buffer_lo, ..default_pixel_buffer_hi

# Util functions on data from EchelleTools
import ..num_pixels_per_order, ..num_orders, ..num_obs

include("common.jl")
include("blaze_map.jl")
include("blaze_polynomial.jl")

end # module Blaze

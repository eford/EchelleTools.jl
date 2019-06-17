#@reexport module BlazeMap

#export fit_blaze_1d, make_flat_blaze, make_blaze, calc_blaze, eval_blaze
#export calc_normalizations, renormalize

#import ..AbstractBlazeType, ..AbstractBlaze1DType, ..AbstractBlaze2DType
import ..AbstractEchelleObservation, ..AbstractEchelleObservationSet
import ..default_pixel_buffer_lo, ..default_pixel_buffer_hi
import ..max_num_pixels_per_order
import ..num_pixels_per_order, ..num_orders, ..num_obs
#import ..calc_blaze
#using ..BlazeCommon

struct BlazeMap1D{T1,T2,T3} <: AbstractBlaze1DType where {T1<:Integer, T2<:Real, T3<:Real}
   order::T1
   log_lambda_lo::T2
   log_lambda_hi::T2
   pixel_range::UnitRange{T1}
   model::AbstractArray{T3,1}
end
export BlazeMap1D

function make_flat_blaze(;num_pixels_per_order::Integer = max_num_pixels_per_order,
      pixel_range::AbstractRange = 1:num_pixels_per_order,
      order::Integer = 0, log_lambda_lo::Real = 0.0, log_lambda_hi::Real = 0.0)
      return BlazeMap1D(order, log_lambda_lo, log_lambda_hi, pixel_range, ones(length(pixel_range)) )
end
make_flat_blaze_map_1d = make_flat_blaze
export make_flat_blaze_map_1d

# TODO: Change calc_blaze to eval_blaze
function calc_blaze(b::BlazeMap1D{T1,T2,T3}, pixel::T4) where {T1<:Integer, T2<:Real, T3<:Real, T4<:Real}
      #println("# pix_rng.lo=", b.pixel_range[1], " pix=", pixel, " pix_rng.hi=", b.pixel_range[end], "len(pix_rng)=", length(b.pixel_range))
      #@assert (b.pixel_range[1]-1) <= pixel <= (b.pixel_range[end] +1)
      @assert 1 <= pixel <= length(b.model)
      pixel_lo = max(1, floor(Int64,pixel))
      pixel_hi = min(pixel_lo+1, b.pixel_range[end])
      w_hi = (pixel-pixel_lo)
      w_lo = 1-w_hi
      w_lo*b.model[pixel_lo] + w_hi*b.model[pixel_hi]
end

function calc_blaze(b::BlazeMap1D{T1,T2,T3}, pixel::T4) where {T1<:Integer, T2<:Real, T3<:Real, T4<:Integer}
      @assert b.pixel_range[1] <= pixel <= b.pixel_range[end]
      b.model[pixel]
end
export calc_blaze

function calc_blaze_at_log_lambda(b::BT, log_lambda::Real) where {BT<:AbstractBlaze1DType }
      pixel = (log_lambda-b.log_lambda_lo)/(b.log_lambda_hi-b.log_lambda_lo)*(length(b.pixel_range)-1) + b.pixel_range[1]
      calc_blaze(b,pixel)
end
function calc_blaze_at_log_lambda(b::BT, log_lambda::Union{AbstractRange{T},AbstractArray{T,1}}) where {BT<:AbstractBlaze1DType, T<:Real }
      map(ll->calc_blaze_at_log_lambda(b,ll),log_lambda)
end
export calc_blaze_at_log_lambda



#=
function calc_blaze(b::BT, pixel_range::AbstractRange{T2}) where {BT<:AbstractBlaze1DType, T2<:Integer}
      @assert b.pixel_range[1]-0.5 <= pixel_range[1] <= b.pixel_range[end]+0.5
      @assert b.pixel_range[1]-0.5 <= pixel_range[end] <= b.pixel_range[end]+0.5
      b.model[pixel_range]
end

=#

struct BlazeMap2D{T1,T2,T3} <: AbstractBlaze2DType where {T1<:Integer, T2<:Real, T3<:Real}
   blaze_list::Vector{BlazeMap1D{T1,T2,T3}}
end
export BlazeMap2D

function make_flat_blaze_map_2d(;num_pixels_per_order::Integer = max_num_pixels_per_order,
      num_orders::T1 = max_num_orders,
      pixel_range::AbstractArray{AbstractRange} = fill(1:num_pixels_per_order, num_orders),
      orders::AbstractArray{T1,1} = collect(1:num_orders),
      log_lambda_lo::AbstractArray{T2,1} = zeros(num_orders), log_lambda_hi::AbstractArray{T2,1} = zeros(num_orders)
      ) where { T1<:Integer, T2<:Real }
      @assert length(pixel_range) == length(orders) == length(log_lambda_lo) == length(log_lambda_hi) >= 1
      blaze_map = ones(length(pixel_range[1]))
      return BlazeMap2D( [
            BlazeMap1D(orders[i], log_lambda_lo[i], log_lambda_hi[i], pixel_range[i], blaze_map )
                  for i in 1:length(orders) ] )
end
export make_flat_blaze_map_2d

function covert_blaze_map_into_blaze_polynomial(b::BMT; degree = default_blaze_degree
   ) where { BMT <: BlazeMap1D }
   poly_model = fit_blaze_1d(b.pixel_range, b.model[b.pixel_range], degree=degree)
   BlazePolynomial(b.order,b.lambda_lo,b.lambda_hi,b.pixel_range,poly_model)
end

function covert_blaze_map_into_blaze_polynomial(blaze_list::BMT; degree = default_blaze_degree
      ) where { BMT <: BlazeMap2D }
   BlazeMap2D( [ covert_blaze_map_into_blaze_polynomial(blaze_list[i],degree=degree) for i in 1:length(blaze_list) ] )
end
export covert_blaze_map_into_blaze_polynomial

#include("common.jl")

#end # module BlazeMap

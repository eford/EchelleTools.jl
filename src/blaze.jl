@reexport module Blaze

#export BlazePolynomial, BlazePolynomialList, BlazeMap1D, BlazeMap2D


# Default functions for AbstractBlazeTypes
min_log_lambda(blaze::AbstractBlazeType) = blaze.log_lambda_lo
max_log_lambda(blaze::AbstractBlazeType) = blaze.log_lambda_hi
order(blaze::AbstractBlazeType) = blaze.order
pixel_range(blaze::AbstractBlazeType) = blaze.pixel_range
export order, min_log_lambda, max_log_lambda, pixel_range

#= immutable
function update_log_lambda_min_max!(blaze::AbstractBlazeType, data::DataT) where DataT <: AbstractEchelleObservationOrSet
   r = EchelleTools.make_model_wavelength_range(data)
   blaze.log_lambda_lo = r[1]
   blaze.log_lambda_hi = r[end]
end
=#

#@reexport module BlazePolynomial

using Polynomials     # For blaze function


""" Type for storing a polynomial model for one order's blaze function. """
struct BlazePolynomial{T1,T2} <: AbstractBlaze1DType where {T1<:Integer, T2<:Real}
   order::T1
   log_lambda_lo::T2
   log_lambda_hi::T2
   pixel_range::UnitRange{T1}
   model::Poly{T2}
end
export BlazePolynomial

struct BlazePolynomialList{T1,T2} <: AbstractBlaze2DType where {T1<:Integer, T2<:Real}
   blaze_list::Vector{BlazePolynomial{T1,T2}}
end
export BlazePolynomialList

function fit_blaze_1d(x::AbstractRange{T1}, flux::AbstractArray{T2,1}; degree::Int = default_blaze_degree
      ) where {T1<:Real, T2<:Real}
   polyfit(x,flux,degree) # replaced x.-mean(x)
end
export fit_blaze_1d

function make_flat_blaze_polynomial(;num_pixels_per_order::Integer = max_num_pixels_per_order,
      pixel_range::AbstractRange = 1:num_pixels_per_order,
      order::Integer = 0, degree::Integer = 0, log_lambda_lo::Real = 0.0, log_lambda_hi::Real = 0.0)
      @assert 0 <= degree <= max_degree_blaze
      coeff = vcat([1.0],zeros(degree))
      return BlazePolynomial(order, log_lambda_lo, log_lambda_hi, pixel_range, Poly(coeff) )
end


function make_flat_blaze_polynomial(data::EOST, order::Integer;
      pixel_range::AbstractRange = 1:num_pixels_per_order(data), degree::Integer = 0
      ) where { EOST <: Union{AbstractEchelleObservationSet,AbstractEchelleObservation} }
   @assert 1 <= order <= max_num_orders
   @assert 0 <= degree <= max_degree_blaze
   pixel_lo = max(1,pixel_range[1] - 1)
   pixel_hi = min(num_pixels_per_order(data), pixel_range[end]+1)
   log_lambda_lo = log(minimum(data.lambda[pixel_lo,order,:]))
   log_lambda_hi = log(maximum(data.lambda[pixel_hi,order,:]))
   make_flat_blaze_polynomial(order=order, log_lambda_lo=log_lambda_lo, log_lambda_hi=log_lambda_hi,
                  pixel_range=pixel_range, degree=degree)
end
export make_flat_blaze_polynomial

function make_flat_blaze_polynomial_2d(data::EOST;
      pixel_range::AbstractRange = 1:num_pixels_per_order(data), degree::Integer = 0,
      ) where { EOST <: Union{AbstractEchelleObservationSet,AbstractEchelleObservation} }
   map(order-> make_flat_blaze_polynomial(pixel_range=pixel_range, order=order, degree=degree), 1:num_orders(data) )
end
export make_flat_blaze_polynomial_2d

function make_blaze(data::EOST, order::Integer;
      pixel_range::AbstractRange = 1:size(data.lambda,1),
      pixel_buffer_lo::Integer = default_pixel_buffer_lo,
      pixel_buffer_hi::Integer = default_pixel_buffer_hi,
      degree::Integer = default_blaze_degree,
      verbose::Integer = 0
      ) where { EOST <: AbstractEchelleObservationSet }
      @assert 0 <= degree <= max_degree_blaze
      @assert 1 <= order <= num_orders(data)
      max_itterations_fit_blaze = 20
      # TODO: Insert code to truncate pixel_range to enforce positive blaze
      # blaze = fit_blaze_1d(pixel_range, vec(mean(data.flux[pixel_range,order,:],dims=2)), degree = degree)
      #println("# Order = ", order, " pixel_range = ", pixel_range)
      new_pixel_range = pixel_range
      last_pixel_range = 0:0
      local blaze
      i = 0
      while !((2*abs(new_pixel_range[1] - last_pixel_range[1]) <= pixel_buffer_lo) && (2*abs(new_pixel_range[end] - last_pixel_range[end]) <= pixel_buffer_hi )) &&
                  (i<max_itterations_fit_blaze)
            last_pixel_range = new_pixel_range
            blaze = fit_blaze_1d(new_pixel_range, vec(mean(data.flux[new_pixel_range,order,:],dims=2)), degree = degree)
            pixel_lo = max(pixel_range[1],  findfirst(x->x>0,polyval(blaze,pixel_range)) ) +pixel_buffer_lo
            pixel_hi = min(pixel_range[end],findlast( x->x>0,polyval(blaze,pixel_range)) ) -pixel_buffer_hi
            if verbose > 1
                  println("# order: ", order, " it=", i," range: ", pixel_lo, ":", pixel_hi, " blaze = ", polyval(blaze,max(pixel_lo-1,1)), " ", polyval(blaze,pixel_lo), " ... ", polyval(blaze, pixel_hi) )
            end
            new_pixel_range = pixel_lo:pixel_hi
            i += 1
      end
      if verbose > -1 && !((2*abs(new_pixel_range[1] - last_pixel_range[1]) <= pixel_buffer_lo) && (2*abs(new_pixel_range[end] - last_pixel_range[end]) <= pixel_buffer_hi ))
            println("# Warning: Reached maximum number of itterations fitting blaze for order ", order)
            println("# last_pixel_range = ", last_pixel_range)
            println("# new_pixel_range = ", new_pixel_range)
      end
      # blaze = blaze / mean(polyval(blaze,new_pixel_range))
      pixel_lo = max(1,new_pixel_range[1] - 1)
      pixel_hi = min(num_pixels_per_order(data), new_pixel_range[end]+1)
      log_lambda_lo = log(minimum(data.lambda[pixel_lo,order,:]))
      log_lambda_hi = log(maximum(data.lambda[pixel_hi,order,:]))
      return BlazePolynomial(order, log_lambda_lo, log_lambda_hi, new_pixel_range, blaze )
end
export make_blaze

function make_blaze_polynomial(data::EOST;
      pixel_range::AbstractRange = 1:size(data.lambda,1), degree::Integer = default_blaze_degree
      ) where { EOST <: AbstractEchelleObservationSet}
      map(order-> make_blaze_polynomial(data, order, pixel_range=pixel_range, degree=degree), 1:num_orders(data) )
end
export make_blaze_polynomial

# TODO: Change calc_blaze to eval_blaze
function calc_blaze(b::BlazePolynomial, pixel::Real)
      @assert b.pixel_range[1] <= pixel <= b.pixel_range[end]
      polyval(b.model, pixel)
end
export calc_blaze

function calc_blaze(b::BlazePolynomial, pixel_range::AbstractRange{T2}) where {T2<:Real}
      @assert b.pixel_range[1]-0.5 <= pixel_range[1] <= b.pixel_range[end]+0.5
      @assert b.pixel_range[1]-0.5 <= pixel_range[end] <= b.pixel_range[end]+0.5
      blaze = polyval(b.model, pixel_range)
end

function calc_blaze(b::BlazePolynomial) where {T2<:Real}
      blaze = polyval(b.model, b.pixel_range)
end

function calc_blaze(blaze_list::BlazeLT, data::EOT;
      output::AbstractArray{T2,2} = zeros(eltype(data.flux),size(data.lambda))
      ) where { T2<:Real,
      BlazeLT <: AbstractBlaze2DType,
      #BlazeT <: AbstractBlazeType, BlazeLT <: AbstractArray{BlazeT,1},
      EOT <: AbstractEchelleObservation  }
      @assert length(blaze_list) == num_orders(data)
      for blaze in blaze_list
            output[blaze.pixel_range,blaze.order] .= calc_blaze(blaze)
      end
      return output
end

function calc_blaze(blaze_list::BlazeLT, data::EOST;
      output::AbstractArray{T2,3} = zeros(eltype(data.flux),size(data.lambda))
      ) where { T2<:Real,
      BlazeLT <: AbstractBlaze2DType,
      #BlazeT <: AbstractBlazeType, BlazeLT <: AbstractArray{BlazeT,1},
      EOST <: AbstractEchelleObservationSet  }
      @assert length(blaze_list) == num_orders(data)
      for blaze in blaze_list
            output[blaze.pixel_range,blaze.order,:] .= calc_blaze(blaze)
      end
      return output
end

function calc_normalizations(blaze_list::BlazeLT, data::EOST;
      output::OutputT = zeros(eltype(data.flux),num_obs(data))
      ) where { T2<:Real,
      BlazeLT <: AbstractBlaze2DType,
      # BlazeT <: AbstractBlazeType, BlazeLT <: AbstractArray{BlazeT,1},
      EOST <: AbstractEchelleObservationSet, OutputT <: AbstractArray{T2,1}  }
   @assert length(blaze_list) == num_orders(data)
   @assert length(output) == num_obs(data)
   output .= map(obsid->
      mean(map(order->
            mean(data.flux[blaze_list[order].pixel_range,order,obsid] ./
                  calc_blaze(blaze_list[order]) )
            , 1:num_orders(data) ))
      , 1:num_obs(data))
end
export calc_normalizations

function renormalize(blaze::BlazePolynomial, norm::Real)
      return BlazePolynomial(blaze.order, blaze.log_lambda_lo, blaze.log_lambda_hi,
                  blaze.pixel_range, blaze.model*norm )
end
export renormalize

#end # module BlazePolynomial

#@reexport module BlazeMap

struct BlazeMap1D{T1,T2,T3} <: AbstractBlaze1DType where {T1<:Integer, T2<:Real, T3<:Real}
   order::T1
   log_lambda_lo::T2
   log_lambda_hi::T2
   pixel_range::UnitRange{T1}
   model::AbstractArray{T3,1}
end
export BlazeMap1D

function make_flat_blaze_map_1d(;num_pixels_per_order::Integer = max_num_pixels_per_order,
      pixel_range::AbstractRange = 1:num_pixels_per_order,
      order::Integer = 0, log_lambda_lo::Real = 0.0, log_lambda_hi::Real = 0.0)
      return BlazeMap1D(order, log_lambda_lo, log_lambda_hi, pixel_range, ones(length(pixel_range)) )
end
export make_flat_blaze_map_1d

# TODO: Change calc_blaze to eval_blaze
function calc_blaze(b::BT, pixel::Real) where {BT<:AbstractBlaze1DType }
      @assert b.pixel_range[1]-1 <= pixel <= b.pixel_range[end] +1
      pixel_lo = maximum(1, floor(Int64,pixel))
      pixel_hi = minimum(pixel_lo+1, b.pixel_range[end])
      w_hi = (pixel-pixel_lo)
      w_lo = 1-w_hi
      w_lo*b.model[pixel_lo] + w_hi*b.model[pixel_hi]
end

function calc_blaze(b::BT, pixel::Integer) where {BT<:AbstractBlaze1DType }
      @assert b.pixel_range[1] <= pixel <= b.pixel_range[end]
      b.model[pixel]
end

function calc_blaze(b::BT, pixel_range::AbstractRange{T2}) where {BT<:AbstractBlaze1DType, T2<:Real}
      map(p->calc_blaze(b,p), pixel_range)
end

function calc_blaze(b::BT, pixel_range::AbstractRange{T2}) where {BT<:AbstractBlaze1DType, T2<:Integer}
      @assert b.pixel_range[1]-0.5 <= pixel_range[1] <= b.pixel_range[end]+0.5
      @assert b.pixel_range[1]-0.5 <= pixel_range[end] <= b.pixel_range[end]+0.5
      b.model[pixel_range]
end

function calc_blaze(b::BT) where {BT<:AbstractBlaze1DType, T2<:Real}
      calc_blaze(b, b.pixel_range)
end
export calc_blaze


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
export covert_blaze_map_into_blaze_polynomial

function covert_blaze_map_into_blaze_polynomial(blaze_list::BMT; degree = default_blaze_degree
      ) where { BMT <: BlazeMap2D }
   BlazeMap2D( [ covert_blaze_map_into_blaze_polynomial(blaze_list[i],degree=degree) for i in 1:length(blaze_list) ] )
end
export covert_blaze_map_into_blaze_polynomial

#end # module BlazeMap

end # module Blaze

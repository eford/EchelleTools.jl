using Polynomials     # For blaze function

""" Type for storing info about model of blaze function. """
struct Blaze{T1,T2} <: AbstractBlazeType where {T1<:Integer, T2<:Real}
   order::T1
   pixel_range::UnitRange{T1}
   model::Poly{T2}
end
export Blaze

function fit_blaze_1d(x::AbstractRange{T1}, flux::AbstractArray{T2,1}; degree::Int = default_blaze_degree
      ) where {T1<:Real, T2<:Real}
   polyfit(x,flux,degree) # replaced x.-mean(x)
end

function make_flat_blaze(;num_pixels_per_order::Integer = max_num_pixels_per_order, pixel_range::AbstractRange = 1:num_pixels_per_order, order::Integer = 0, degree::Integer = 0)
      @assert 0 <= degree <= max_degree_blaze
      coeff = vcat([1.0],zeros(degree))
      return Blaze(order, pixel_range, Poly(coeff) )
end

function make_flat_blaze(data::EOST, order::Integer;
      pixel_range::AbstractRange = 1:num_pixels_per_order(data), degree::Integer = 0
      ) where { EOST <: Union{AbstractEchelleObservationSet,AbstractEchelleObservation} }
   @assert 1 <= order <= max_num_orders
   @assert 0 <= degree <= max_degree_blaze
   make_flat_blaze(pixel_range=pixel_range, order=order, degree=degree)
end

function make_flat_blaze(data::EOST;
      pixel_range::AbstractRange = 1:num_pixels_per_order(data), degree::Integer = 0
      ) where { EOST <: Union{AbstractEchelleObservationSet,AbstractEchelleObservation} }
   map(order-> make_flat_blaze(pixel_range=pixel_range, order=order, degree=degree), 1:num_orders(data) )
end

function make_blaze(data::EOST, order::Integer;
      pixel_range::AbstractRange = 1:size(data.lambda,1), degree::Integer = default_blaze_degree
      ) where { EOST <: AbstractEchelleObservationSet }
      @assert 0 <= degree <= max_degree_blaze
      @assert 1 <= order <= num_orders(data)
      # TODO: Insert code to truncate pixel_range to enforce positive blaze
      # blaze = fit_blaze_1d(pixel_range, vec(mean(data.flux[pixel_range,order,:],dims=2)), degree = degree)
      #println("# Order = ", order, " pixel_range = ", pixel_range)
      new_pixel_range = pixel_range
      last_pixel_range = 0:0
      local blaze
      while new_pixel_range != last_pixel_range
            last_pixel_range = new_pixel_range
            blaze = fit_blaze_1d(new_pixel_range, vec(mean(data.flux[new_pixel_range,order,:],dims=2)), degree = degree)
            pixel_lo = max(pixel_range[1],  pixel_range[findfirst(x->x>0,polyval(blaze,pixel_range))] )
            pixel_hi = min(pixel_range[end],pixel_range[findlast( x->x>0,polyval(blaze,pixel_range))] )
            #println(pixel_lo, ":", pixel_hi, " blaze = ", polyval(blaze,max(pixel_lo-1,1)), " ", polyval(blaze,pixel_lo), " ... ", polyval(blaze, pixel_hi) )
            new_pixel_range = pixel_lo:pixel_hi
      end
      # blaze = blaze / mean(polyval(blaze,new_pixel_range))
      return Blaze(order, new_pixel_range, blaze )
end

function make_blaze(data::EOST;
      pixel_range::AbstractRange = 1:size(data.lambda,1), degree::Integer = default_blaze_degree
      ) where { EOST <: AbstractEchelleObservationSet}
      map(order-> make_blaze(data, order, pixel_range=pixel_range, degree=degree), 1:num_orders(data) )
end

# TODO: Change calc_blaze to eval_blaze
function calc_blaze(b::Blaze, pixel::Real)
      @assert b.pixel_range[1] <= pixel <= b.pixel_range[end]
      polyval(b.model, pixel)
end

function calc_blaze(b::Blaze, pixel_range::AbstractRange{T2}) where {T2<:Real}
      #@assert b.pixel_range[1] <= pixel_range[1] <= b.pixel_range[end]
      #@assert b.pixel_range[1] <= pixel_range[end] <= b.pixel_range[end]
      blaze = polyval(b.model, pixel_range)
end

function calc_blaze(b::Blaze) where {T2<:Real}
      blaze = polyval(b.model, b.pixel_range)
end

function calc_blaze(blaze_list::BlazeLT, data::EOST;
      output::AbstractArray{T2,3} = zeros(eltype(data.flux),size(data.lambda))
      ) where { T2<:Real,
      BlazeT <: AbstractBlazeType, BlazeLT <: AbstractArray{BlazeT,1},
      EOST <: AbstractEchelleObservationSet  }
      @assert length(blaze_list) == num_orders(data)
      for blaze in blaze_list
            output[:,blaze.order,:] .= calc_blaze(blaze,1:num_pixels_per_order(data))
      end
      return output
end

function calc_normalizations(blaze_list::BlazeLT, data::EOST;
      output::OutputT = zeros(eltype(data.flux),num_obs(data))
      ) where { T2<:Real,
      BlazeT <: AbstractBlazeType, BlazeLT <: AbstractArray{BlazeT,1},
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

struct PSFDeltaFunction <: AbstractPSF
end
export PSFDeltaFunction

#=
@reexport module BlazeCommon

import ..AbstractEchelleObservation, ..AbstractEchelleObservationSet
import ..AbstractBlazeType, ..AbstractBlaze1DType, ..AbstractBlaze2DType
import ..max_num_pixels_per_order
import ..default_pixel_buffer_lo, ..default_pixel_buffer_hi

# Util functions on data from EchelleTools
import ..num_pixels_per_order, ..num_orders, ..num_obs
=#

using Statistics

function calc_blaze(b::BT, pixel_range::AbstractRange{T2}) where {BT<:AbstractBlaze1DType, T2<:Real}
      map(p->calc_blaze(b,p), pixel_range)
end

function calc_blaze(b::BT) where {BT<:AbstractBlaze1DType}
      calc_blaze(b, b.pixel_range)
end

function calc_blaze(blaze_list::BlazeLT, data::EOT;
      output::AbstractArray{T2,2} = zeros(eltype(data.flux),size(data.lambda))
      ) where { T2<:Real,
      BlazeLT <: AbstractBlaze2DType,
      #BlazeT <: AbstractBlazeType, BlazeLT <: AbstractArray{BlazeT,1},
      EOT <: AbstractEchelleObservation  }
      @assert length(blaze_list.blaze_list) == num_orders(data)
      for blaze in blaze_list.blaze_list
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
      @assert length(blaze_list.blaze_list) == num_orders(data)
      for blaze in blaze_list.blaze_list
            output[blaze.pixel_range,blaze.order,:] .= calc_blaze(blaze)
      end
      return output
end
export calc_blaze


function calc_normalizations(blaze_list::BlazeLT, data::EOST;
      output::OutputT = zeros(eltype(data.flux),num_obs(data))
      ) where { T2<:Real,
      BlazeLT <: AbstractBlaze2DType,
      # BlazeT <: AbstractBlazeType, BlazeLT <: AbstractArray{BlazeT,1},
      EOST <: AbstractEchelleObservationSet, OutputT <: AbstractArray{T2,1}  }
   @assert length(blaze_list.blaze_list) == num_orders(data)
   @assert length(output) == num_obs(data)
   output .= map(obsid->
      mean(map(order->
            mean(data.flux[blaze_list.blaze_list[order].pixel_range,order,obsid] ./
                  calc_blaze(blaze_list.blaze_list[order]) )
            , 1:num_orders(data) ))
      , 1:num_obs(data))
end
export calc_normalizations

#end # module BlazeComm

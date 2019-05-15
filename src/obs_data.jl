#include("abstract_types.jl")

# Structures for reading in data
""" EchelleObservation contains observations for one echelle observation.
Arrays organized into (pixels,orders).
z_bc = Barrycentric correction to convert from observer frame to solar system barrycenter frame.
"""
struct EchelleObservation{T1,T2,T3,T4,T5} <: AbstractEchelleObservation where { T1<:Number, T2<:Number, T3<:Number, T4<:Number, T5<:AbstractInstrument }
   lambda::AbstractArray{T1,2}
   flux::AbstractArray{T2,2}
   var::AbstractArray{T3,2}
   z_bc::T4
   instrument::T5
   metadata::Dict{Symbol,Any}

   function EchelleObservation{T1,T2,T3,T4,T5}(;lambda::LambdaT,flux::FluxT,var::VarT,
       z_bc::T4, instrument::T5 = T5(), metadata::AbstractDict = Dict{Symbol,Any}()
       )  where { T1<:Number, T2<:Number, T3<:Number, T4<:Number, T5<:AbstractInstrument,
       LambdaT <: AbstractArray{T1,2}, FluxT <: AbstractArray{T2,2},
       VarT <: AbstractArray{T3,2}
       }
      @assert size(lambda)==size(flux)==size(var)
      @assert size(lambda,1) <= max_num_pixels_per_order
      @assert size(lambda,2) <= max_num_orders
      new(lambda,flux,var,z_bc,instrument,metadata)
   end
end
export EchelleObservation

""" EchelleObservationSet contains observations for multiple echelle observations.
Arrays organized into (pixels,orders,observations).
z_bc = Barrycentric correction to convert from observer frame to solar system barrycenter frame.
"""
struct EchelleObservationSet{T1,T2,T3, #= T4, =# T5} <: AbstractEchelleObservationSet where {T1<:Number, T2<:Number, T3<:Number, #= T4<:Number, =# T5<:AbstractInstrument}
   lambda::Array{T1,3}
   flux::Array{T2,3}
   var::Array{T3,3}
   #= z_bc::Array{T4,1} =#
   instrument::T5
   metadata::Dict{Symbol,Any}
   metadata_list::Array{Dict{Symbol,Any},1}

   function EchelleObservationSet{T1,T2,T3, #= T4, =# T5}(;lambda::AbstractArray{T1,3},
         flux::AbstractArray{T2,3},var::AbstractArray{T3,3}, instrument::T5 = T5(),
         metadata::AbstractDict = Dict{Symbol,Any}(), metadata_list::AbstractArray{Dict{Symbol,Any},1} = Dict{Symbol,Any}[]
          )  where {T1<:Number, T2<:Number, T3<:Number, #= T4<:Number, =# T5<:AbstractInstrument}
      @assert size(lambda)==size(flux)==size(var)
      @assert size(lambda,1) <= max_num_pixels_per_order
      @assert size(lambda,2) <= max_num_orders
      #@assert length(z_bc)==size(lambda,3)
      new(lambda,flux,var,#= z_bc, =# instrument,metadata,metadata_list)
   end

end
export EchelleObservationSet

""" Preallocate arrays for EchelleObservationSet.  Takes obs_size (pixels_per_order,num_orders) and num_obs.  Optionally, instrument and metadata."""
function EchelleObservationSet{T1,T2,T3, #= T4, =# T5}(obs_size::Tuple{Int64,Int64},num_obs::Int = 0;
      instrument::T5 = T5(), metadata::AbstractDict = Dict{Symbol,Any}()
      )  where {T1<:Number, T2<:Number, T3<:Number, #=T4<:Number, =# T5<:AbstractInstrument}
   @assert 0 <= num_obs <= max_num_obs
   lambda = Array{T1,3}(undef,obs_size[1],obs_size[2],num_obs)
   flux   = Array{T2,3}(undef,obs_size[1],obs_size[2],num_obs)
   var    = Array{T3,3}(undef,obs_size[1],obs_size[2],num_obs)
   #=z_bc   = zeros(T4,num_obs) =#
   metadata_list = Array{Dict{Symbol,Any},1}(undef,num_obs)
   EchelleObservationSet{eltype(lambda),eltype(flux),eltype(var),#= eltype(z_bc), =# typeof(instrument) }(lambda=lambda,flux=flux,var=var,#= z_bc=z_bc, =# instrument=instrument,metadata=metadata,metadata_list=metadata_list)
end

function z_bc( eos::EchelleObservationSet{T1,T2,T3, #= T4, =# T5}, obsid::Integer
   ) where {T1<:Number, T2<:Number, T3<:Number, #= T4<:Number, =# T5<:AbstractInstrument}
   @assert 1 <= obsid <= num_obs(eos)
   return metadata_list[obsid][:z_bc]
end

""" Extract one EchelleObservation from an EchelleObservationSet and obs number."""
function EchelleObservation(eos::EchelleObservationSet{T1,T2,T3, #= T4, =# T5}, obs::Integer) where {T1<:Number, T2<:Number, T3<:Number, #= T4<:Number, =# T5<:AbstractInstrument}
   # md = merge(eos.metadata,eos.metadata_list[obs]) # TODO: Should we merge metadata from obs set?
   EchelleObservation(lambda=view(eos.lambda,:,:,obs), flux=view(eos.flux,:,:,obs), var=view(eos.var,:,:,obs), z_bc=z_bc(eos,obs), instrument=eos.instrument, metadata=eos.metadata_list[obs] )
end

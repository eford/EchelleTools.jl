# TODO: Need to think about how to refactor so general for NEID, EXPRES, HARPS-N,... 

# Define file format
abstract type AbstractManifestFormatEntry end

struct ManifestFormatEntry <: AbstractManifestFormatEntry
    name::Symbol
    type::Type
    cols::UnitRange{Int64}
end

OrdersType = Union{Int64,UnitRange{Int64},AbstractVector{Int64}}

# Abstract Types  
abstract type AbstractEchelleObservation end
abstract type AbstractEchelleObservationSet end

struct EchelleObservation{T1,T2,T3} <: AbstractEchelleObservation where {T1<:Number, T2<:Number, T3<:Number}
   lambda::Array{T1,2}
   flux::Array{T2,2}
   var::Array{T3,2}
   metadata::Dict{Symbol,Any}

   function EchelleObservation{T1,T2,T3}(;lambda::AbstractArray{T1,2},flux::AbstractArray{T2,2},var::AbstractArray{T3,2}, metadata::AbstractDict = Dict{Symbol,Any}() )  where {T1<:Number, T2<:Number, T3<:Number}
      @assert size(lambda)==size(flux)==size(var)
      new(lambda,flux,var,metadata)
   end
end

struct EchelleObservationSet{T1,T2,T3} <: AbstractEchelleObservationSet where {T1<:Number, T2<:Number, T3<:Number}
   lambda::Array{T1,3}
   flux::Array{T2,3}
   var::Array{T3,3}
   metadata::Dict{Symbol,Any}
   metadata_list::Array{Dict{Symbol,Any},1}

   function EchelleObservationSet{T1,T2,T3}(;lambda::AbstractArray{T1,3},flux::AbstractArray{T2,3},var::AbstractArray{T3,3}, metadata::AbstractDict = Dict{Symbol,Any}(), metadata_list::AbstractArray{Dict{Symbol,Any},1} = Dict{Symbol,Any}[] )  where {T1<:Number, T2<:Number, T3<:Number}
      @assert size(lambda)==size(flux)==size(var)
      new(lambda,flux,var,metadata,metadata_list)
   end

end

function EchelleObservationSet{T1,T2,T3}(obs_size::Tuple{Int64,Int64},num_obs::Int = 0; metadata::AbstractDict = Dict{Symbol,Any}() )  where {T1<:Number, T2<:Number, T3<:Number}
      lambda = Array{T1,3}(undef,obs_size[1],obs_size[2],num_obs)
      flux   = Array{T2,3}(undef,obs_size[1],obs_size[2],num_obs)
      var    = Array{T3,3}(undef,obs_size[1],obs_size[2],num_obs)
      metadata_list = Array{Dict{Symbol,Any},1}(undef,num_obs)
      EchelleObservationSet{T1,T2,T3}(lambda=lambda,flux=flux,var=var,metadata=metadata,metadata_list=metadata_list)
end


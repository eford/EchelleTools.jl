# Abstract Types
#abstract type AbstractInstrument end  # in constants

abstract type AbstractEchelleObservationOrSet end
abstract type AbstractEchelleObservation <: AbstractEchelleObservationOrSet end
abstract type AbstractEchelleObservationSet <: AbstractEchelleObservationOrSet end
abstract type AbstractBlazeType end
abstract type AbstractPSF end
abstract type AbstractSpectralOrderModel end
abstract type AbstractSpectralModel end
abstract type AbstractObservationModel end
abstract type AbstractObsSeriesModel end
abstract type AbstractManifestFormatEntry end
export  AbstractEchelleObservationOrSet, AbstractEchelleObservation, AbstractEchelleObservationSet
export AbstractBlazeType, AbstractPSF
export AbstractSpectralOrderModel, AbstractSpectralModel
export AbstractObservationModel, AbstractManifestFormatEntry

# Type Aliases
#const AbstractEchelleObservationOrSet = Union{AbstractEchelleObservation, AbstractEchelleObservationSet}
const OrdersType = Union{Int64,UnitRange{Int64},AbstractVector{Int64}}
const Wavelength1DType{T1} = Union{AbstractArray{T1,1}, AbstractRange{T1} } where T1<:Real
const Wavelength2DType{T1} = Union{AbstractArray{T1,2}, AbstractArray{AbstractRange{T1},1} } where T1<:Real
export OrdersType, Wavelength1DType, Wavelength2DType

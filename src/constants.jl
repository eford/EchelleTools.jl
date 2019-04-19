# abstract type AbstractInstrument end  # in types.jl
struct HPFType    <: AbstractInstrument end
struct NEIDType   <: AbstractInstrument end
struct EXPRESType <: AbstractInstrument end
struct HARPSType  <: AbstractInstrument end
struct HARPSNType <: AbstractInstrument end

const instrument_symbol_list = [:HPF]  # Eventually add:  :HARPS, :HARPSN, :EXPRES, :NEID


# abstract type AbstractInstrument end  # in abstract_types.jl
struct HPFType    <: AbstractInstrument end
struct NEIDType   <: AbstractInstrument end
struct EXPRESType <: AbstractInstrument end
struct HARPSType  <: AbstractInstrument end
#struct HARPSNType <: AbstractInstrument end

export AbstractInstrument, HPFType, NEIDType, HARPSType
#export HARPSNType

#instrument_symbol_list = Symbol[]
instrument_symbol_list = Symbol[:HPF, :NEID, :EXPRES, :HARPSN, :HARPS ]
export instrument_symbol_list

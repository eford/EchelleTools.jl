# abstract type AbstractInstrument end  # in types.jl
abstract type AbstractInstrument end
struct HPFType    <: AbstractInstrument end
struct NEIDType   <: AbstractInstrument end
struct EXPRESType <: AbstractInstrument end
struct HARPSType  <: AbstractInstrument end
struct HARPSNType <: AbstractInstrument end

const instrument_symbol_list = [:HPF, :NEID, :EXPRES, :HARPSN, :HARPS ]

export AbstractInstrument, HPFType, NEIDType, HARPSType, HARPSNType, instrument_symbol_list

const max_num_pixels_per_order = 8096
const max_num_orders = 128
const max_template_length = max_num_orders * max_num_pixels_per_order
const max_num_obs = 1024
const max_degree_blaze = 5
const max_num_basis_vectors = 16
const min_num_star_basis_vectors = 1
const max_num_star_basis_vectors = max_num_basis_vectors
const min_num_telluric_basis_vectors = 0
const max_num_telluric_basis_vectors = max_num_basis_vectors
export max_num_pixels_per_order, max_num_orders, max_template_length, max_num_obs, max_degree_blaze
export max_num_basis_vectors, min_num_star_basis_vectors, max_num_star_basis_vectors, min_num_telluric_basis_vectors, max_num_telluric_basis_vectors

const default_blaze_degree = 3
export default_blaze_degree

const speed_of_light_mps = 299792458.0 # m/s
export speed_of_light_mps

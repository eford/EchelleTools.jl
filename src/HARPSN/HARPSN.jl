@reexport module HARPSN

# Import types from EchelleTools
#import ..ManifestFormatEntry
import ..EchelleObservation, ..EchelleObservationSet, ..OrdersType
import ..AbstractInstrument #, ..HARPSNType
import ..instrument_symbol_list

using Reexport
@reexport using ..Blaze
#using ..Blaze.BlazeMap

struct HARPSNType <: AbstractInstrument end
push!(instrument_symbol_list,:HARPSN)

include("constants.jl")

# Set paths # TODO: Move to config file
base_path_default = "SET ME IN customizations.jl"
manifest_subdir_default = "."
data_subdir_default = "e2ds"
include("customizations.jl")
base_path = base_path_default

include("manifests.jl")
get_harpsn_base_path = Manifests.get_base_path
set_harpsn_base_path! = Manifests.set_base_path!
export get_harpsn_base_path, set_harpsn_base_path!
read_harpsn_manifest_file = Manifests.read
read_harpsn_manifest_files = Manifests.read
get_harpsn_data_filenames = Manifests.get_data_filenames
get_harpsn_data_filenames_and_bsrv = Manifests.get_data_filenames_and_bsrv
export read_harpsn_manifest_file, read_harpsn_manifest_files, get_harpsn_data_filenames
export get_harpsn_data_filenames_and_bsrv

include("readdata.jl")
read_harpsn_data = read_data
read_harpsn_blaze_list = read_blaze_list
export read_harpsn_data
export read_harpsn_blaze_list

end # module HARPSN

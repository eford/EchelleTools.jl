module HARPSN

# Import types from EchelleTools
#import ..ManifestFormatEntry
import ..EchelleObservation, ..EchelleObservationSet, ..OrdersType
import ..AbstractInstrument, ..HARPSNType

include("constants.jl")

# Set paths # TODO: Move to config file
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
export read_harpsn_data

end # module HARPSN

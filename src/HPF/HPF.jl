module HPF

# Import types from EchelleTools
import ..ManifestFormatEntry
import ..EchelleObservation, ..EchelleObservationSet, ..OrdersType
import ..AbstractInstrument, ..HPFType

include("constants.jl")

# Set paths # TODO: Move to config file
include("customizations.jl")
base_path = base_path_default

include("manifests.jl")
get_hpf_base_path = Manifests.get_base_path
set_hpf_base_path! = Manifests.set_base_path!
export get_hpf_base_path, set_hpf_base_path!
read_hpf_manifest_file = Manifests.read
read_hpf_manifest_files = Manifests.read
get_hpf_data_filenames = Manifests.get_data_filenames
export read_hpf_manifest_file, read_hpf_manifest_files, get_hpf_data_filenames

include("readdata.jl")
read_hpf_data = read_data
export read_hpf_data

end # module HPF


module EchelleTools

# Specify paths to manifest and data files
include("customizations.jl")
base_path = base_path_default
export get_base_path, set_base_path!

include("types.jl")
export EchelleObservation, EchelleObservationSet

include("constants.jl")

include("hpf/manifests.jl")
read_hpf_manifest_file = HPFManifests.read
read_hpf_manifest_files = HPFManifests.read
#read_hpf_manifest_file = HPFManifests.read_manifest_file 
#read_hpf_manifest_files = HPFManifests.read_manifest_files
get_hpf_data_filenames = HPFManifests.get_data_filenames
export read_hpf_manifest_file, read_hpf_manifest_files, get_hpf_data_filenames

include("hpf/readdata.jl")
read_hpf_data = HPFData.read_data
export read_hpf_data

include("analyze.jl")
export mean_spectrum, median_spectrum

end # module

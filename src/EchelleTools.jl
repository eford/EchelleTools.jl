module EchelleTools

include("types.jl")
export EchelleObservation, EchelleObservationSet

include("constants.jl")

include("HPF/HPF.jl")
read_hpf_manifest_file = HPF.read_hpf_manifest_file
read_hpf_manifest_files = HPF.read_hpf_manifest_files
get_hpf_data_filenames = HPF.get_hpf_data_filenames
read_hpf_data = HPF.read_hpf_data
export read_hpf_manifest_file, read_hpf_manifest_files, get_hpf_data_filenames, read_hpf_data

include("analyze.jl")
export mean_spectrum, median_spectrum

end # module EchelleTools

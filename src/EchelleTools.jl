module EchelleTools
using Reexport

include("constants.jl")

include("abstract_types.jl")

include("util.jl")
export num_pixels_per_order, num_orders, num_obs
export flux_mean_over_pixels, lambda_mean_over_pixels
export flux_median_over_pixels, lambda_median_over_pixels
export default_telluric_template, default_stellar_template
export make_model_wavelength_range
export lambda_ssb_from_obs, log_lambda_ssb_from_obs, lambda_source, log_lambda_source
export ManifestFormatEntry

include("blaze.jl")
export Blaze
export fit_blaze_1d, make_flat_blaze, make_blaze, calc_blaze, eval_blaze
export calc_normalizations

include("obs_data.jl")
export EchelleObservation, EchelleObservationSet

include("model.jl")
export SpectralOrderModel, SpectralModel
export get_order

export ObservationModel, ObsSeriesModel
export eval_model, eval_model_baseline

#=
include("analyze.jl")
export mean_spectrum, median_spectrum
export lambda_ssb_from_obs, log_lambda_ssb_from_obs
export make_blaze_and_template_1d, make_blaze_and_template_2d
export make_stellar_basis, default_telluric_model
export make_template_1d
export make_spectral_order_model, make_spectral_model
export make_observations_model
# export model_spectrum_1d, model_spectrum_2d

include("pca.jl")
export interp_one_obs_order_for_pca_1d, fit_pca_obs_order_interpolated
=#

#=
include("toy_data.jl")
export spectral_line, make_toy_data
=#

# Instrument specific files
include("HARPSN/HARPSN.jl")
read_harpsn_manifest_file = HARPSN.read_harpsn_manifest_file
read_harpsn_manifest_files = HARPSN.read_harpsn_manifest_files
get_harpsn_data_filenames = HARPSN.get_harpsn_data_filenames
get_harpsn_data_filenames = HARPSN.get_harpsn_data_filenames
read_harpsn_data = HARPSN.read_harpsn_data
export read_harpsn_manifest_file, read_harpsn_manifest_files, get_harpsn_data_filenames, read_harpsn_data
export get_harpsn_data_filenames_and_bsrv

#=
include("HPF/HPF.jl")
read_hpf_manifest_file = HPF.read_hpf_manifest_file
read_hpf_manifest_files = HPF.read_hpf_manifest_files
get_hpf_data_filenames = HPF.get_hpf_data_filenames
read_hpf_data = HPF.read_hpf_data
export read_hpf_manifest_file, read_hpf_manifest_files, get_hpf_data_filenames, read_hpf_data
=#

end # module EchelleTools

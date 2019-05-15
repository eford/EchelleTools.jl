using Dates
using EchelleTools
using Statistics

(fns,zbsrv) = EchelleTools.HARPSN.get_harpsn_data_filenames_and_bsrv(date_start=Date(2015,1,1),date_stop=Date(2015,12,1))
data = EchelleTools.HARPSN.read_data(fns)
zlist = map(d->d[:z_bc],data.metadata_list)
flat_blaze_list = make_flat_blaze(data)
blaze_list = make_blaze(data)
sm = SpectralModel(data)
sm = SpectralModel(data, blaze_list=blaze_list)
normalizations = calc_normalizations(blaze_list,data)
obsm = ObsSeriesModel(num_obs(data), z_bc=zlist, normalization=normalizations)
som = SpectralOrderModel(sm,1)
om = ObservationModel(obsm,1)
m11 = eval_model(log.(data.lambda[:,1,1]), som, om)[som.blaze.pixel_range]
m1 = eval_model(log.(data.lambda[:,:,1]), sm, om)
m1 = eval_model(log.(data.lambda), sm, obsm)


order = 1
obsid = 1
obs_i = EchelleTools.ObservationModel( star_score=view(obsm.star_score,:,obsid), telluric_score=view(obsm.telluric_score,:,obsid), z_bc=obsm.z_bc[obsid], normalization=obsm.normalization[obsid])
llr = make_model_wavelength_range(EchelleObservation(data,obsid), order)
som = EchelleTools.make_spectral_order_model(data, llr, order)

llr = make_model_wavelength_range(data)
sm = make_spectral_model(data, z_bc=zlist)

evalm = eval_model(llr,sm[1],obs_i)


obsm = make_observations_model(data, sm, z_bc=zlist)
evalm = eval_model(llr,sm[1],1)

#=
obsid = 1
obs_i = EchelleTools.ObservationModel( star_score=view(obsm.star_score,:,obsid), telluric_score=view(obsm.telluric_score,:,obsid), z_bc=obsm.z_bc[obsid], normalization=obsm.normalization[obsid])
=#

#=
using Plots
pyplot()
plot(exp.(llr),evalm)
scatter!(data.lambda[:,1,1],data.flux[:,1,1],markersize=1)
=#



#=
import EchelleTools.AbstractBlazeType
import EchelleTools.AbstractEchelleObservationSet
import EchelleTools.SpectralOrderModel
import EchelleTools.num_obs
import EchelleTools.num_orders
using Interpolations
import EchelleTools.calc_blaze
import EchelleTools.make_stellar_basis
import EchelleTools.default_telluric_template
import EchelleTools.make_flat_blaze

=#

#=
som = EchelleTools.SpectralOrderModel(llr)
obsm = EchelleTools.ObservationsModel(length(zlist))
(blaze1, templ1) = make_blaze_and_template_1d(harpsn_data,order=1,z_bc=zlist)

templ = EchelleTools.make_template_1d(harpsn_data,llr, 1)
starb = EchelleTools.default_stellar_template(length(llr))
telluricb = EchelleTools.default_telluric_template(length(llr))
pixrng = 1:length(llr)
blaz = make_flat_blaze(pixel_range=pixrng)
m = EchelleTools.SpectralOrderModel(llr, star_basis=starb, telluric_basis=telluricb, blaze=blaz )
=#

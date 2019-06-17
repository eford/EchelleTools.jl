using EchelleTools
using Dates
(fns,zbsrv) = EchelleTools.HARPSN.get_harpsn_data_filenames_and_bsrv(date_start=Date(2015,1,1),date_stop=Date(2015,9,1))
data = EchelleTools.HARPSN.read_data(fns)
zlist = map(d->d[:z_bc],data.metadata_list)
blaze_list = read_blaze_list(data) # Need to inject lambda limits from somewhere into blaze_list

sm = SpectralModel(data, blaze_list=blaze_list)
normalizations =  EchelleTools.Blaze.calc_normalizations(blaze_list,data)
obsm = ObsSeriesModel(num_obs(data), z_bc=zlist, normalization=normalizations)
order = 30
obsid = 1
om = ObservationModel(obsm,obsid)
som = SpectralOrderModel(sm,order)

m11 = eval_model(som, om)[som.blaze.pixel_range]
m1 = eval_model(sm, om, log.(data.lambda[:,:,obsid]), )
m = eval_model(sm, obsm, log.(data.lambda))
#chisq =  eval_model_fit(data, sm, obsm) # WARNING: slow
num_pixels_summed = mapreduce(i->length(SpectralOrderModel(sm,i).blaze.pixel_range),+,1:num_orders(data))

using Plots
pyplot()
plot(data.lambda[:,order,obsid],m2[:,order,obsid]./obsm.normalization[obsid],markersize=1)
plot!(data.lambda[:,order,obsid+1],m2[:,order,obsid+1]./obsm.normalization[obsid+1],markersize=1)
plot!(data.lambda[:,order,obsid+2],m2[:,order,obsid+2]./obsm.normalization[obsid+2],markersize=1)
scatter!(data.lambda[:,order,obsid],data.flux[:,order,obsid]./obsm.normalization[obsid],markersize=1)
scatter!(data.lambda[:,order,obsid+1],data.flux[:,order,obsid+1]./obsm.normalization[obsid+1],markersize=1)
scatter!(data.lambda[:,order,obsid+2],data.flux[:,order,obsid+2]./obsm.normalization[obsid+2],markersize=1)


#=
flat_blaze_list = make_flat_blaze_polynomial_2d(data)
blaze_list = make_blaze_polynomial_2d(data)
sm = SpectralModel(data)
normalizations = calc_normalizations(blaze_list,data)
ave_normalization = mean(normalizations)
blaze_list = map(i->renormalize(blaze_list[i],ave_normalization), 1:length(blaze_list))
blaze_list = read_blaze_list()
normalizations ./= ave_normalization
=#


#=
using Statistics

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
=#

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

using Interpolations # For making initial guess of template
using Polynomials

function make_blaze_and_template_1d(data::EOST; order::Integer, oversample_factor::Real = 1.0, num_knots::Integer=round(Int64,oversample_factor*num_pixels_per_order(data)),
      z_bc::AbstractVector{T3} = zeros(eltype(data.lambda),num_obs(data))
      ) where { EOST <: AbstractEchelleObservationSet, T3<:Real}
      #=
      pixel_range = 1:size(data.lambda,1)
      blaze = fit_blaze_1d(pixel_range, vec(mean(data.flux[:,order,:],dims=2)), degree = 3)
      blaze = blaze / mean(polyval(blaze,pixel_range))
      =#
      log_lambda_range = make_model_wavelength_range(data, order=order)
      blaze = make_blaze_1d(data, order=order)
      # log_lambda_range = range(log(minimum(data.lambda[1,order,:])), stop=log(maximum(data.lambda[end,order,:])), )
      template = make_template_1d(data, log_lambda_range, order, z_bc=z_bc, blaze=blaze) #pixel_range=blaze.pixel_range, blaze=blaze.model)
      #template = make_template_1d(data, log_lambda_range, order, z_bc=z_bc, pixel_range=blaze.pixel_range, blaze=blaze.model).template
      #return (pixel_range=blaze.pixel_range, blaze=blaze.model, log_lambda_range=log_lambda_range, template=template )
      return blaze, template
end

function make_blaze_and_template_2d(data::EOST;
      orders::Union{AbstractVector{T1},AbstractRange{T1}} = 1:size(data.lambda,2),
      oversample_factor::Real = 1.0, num_knots::Integer=round(Int64,oversample_factor*num_pixels_per_order(data)),
      z_bc::AbstractVector{T3} = zeros(eltype(data.lambda),num_obs(data))
      ) where { EOST <: AbstractEchelleObservationSet, T1<:Integer, T3<:Real}
   map(order->make_blaze_and_template_1d(data,order=order,z_bc=z_bc, num_knots=num_knots), orders) # TODO: TEST
end


""" make_stellar_basis(mean_spectrum; num_basis_vectors = 0)
   Returns array with one row containing the mean_spectrum followed by num_basis_vectors rows of zeros. """
function make_stellar_basis(mean_spectrum::AbstractArray{T,1}; num_basis_vectors::Integer = 0) where {T<:Real}
   @assert 1 <= length(mean_spectrum) <= max_template_length
   @assert 0 <= num_basis_vectors <= max_num_basis_vectors
   output = zeros(T,1+num_basis_vectors,length(mean_spectrum))
   output[1,:] .= mean_spectrum
   return output
end

function make_template_1d(data::EOST, log_lambda_out::AbstractRange{T1}, order::Integer;
    template::AbstractArray{T2,1} = zeros(T1,length(log_lambda_out)),
    z_bc::AbstractVector{T3} = zeros(T1,num_obs(data)),
    pixel_range::AbstractRange{T4} = 1:size(data.lambda,1),
    blaze::AbstractBlazeType = make_flat_blaze(data, pixel_range=pixel_range, order=order),
    num_star_basis_vectors::Integer = 0, num_telluric_basis_vectors::Integer = 0
    ) where { EOST <: AbstractEchelleObservationSet, T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Real}
    @assert num_obs(data) == length(z_bc)
    @assert length(log_lambda_out) == length(template)
    @assert 1<=order<=num_orders(data)
    interpolators = [ interpolate( (log.(lambda_ssb_from_obs(data.lambda[:,order,j],z_bc[j])), ),
         data.flux[:,order,j], Gridded(Linear()) ) for j in 1:num_obs(data) ]

    k_start = findfirst(k->(log_lambda_out[k]) >= log(minimum(data.lambda[1,order,:])),   1:length(log_lambda_out)  )
    if isnothing(k_start)
      k_start = 1
    end
    k_stop  = findlast( k->(log_lambda_out[k]) <= log(maximum(data.lambda[end,order,:])), 1:length(log_lambda_out)  )
    if isnothing(k_stop)
      k_stop = length(log_lambda_out)
    end
    weight = vec(sum(data.flux[:,order,:],dims=1))/size(data.flux,1) # average counts per pixel in each observation
    for k in k_start:k_stop
       weight_sum = 0.0
       for j in 1:num_obs(data)
         if interpolators[j].knots[1][1] <= log_lambda_out[k] <=  interpolators[j].knots[end][end]
            template[k] += interpolators[j](log_lambda_out[k]) # normalized flux times weight
            weight_sum  += weight[j]
         end
       end # for j over observations
       if weight_sum > 0
          # TODO: WARNING: Blaze calculation assumes length(log_lambda_out) == num_pixels_per_order
          pixel = blaze.pixel_range.start+length(blaze.pixel_range)*(k-1)/length(log_lambda_out)
          # println("k=",k," pixel=",pixel)
          b = calc_blaze(blaze,pixel) # polyval(blaze,k) #-mean(pixel_range))
          template[k] /= b * weight_sum
       end
   end # for k over log_lambda_out
   return template
end

function make_spectral_order_model(data::EOST, log_lambda_out::AbstractRange{T1}, order::Integer;
   pixel_range::AbstractRange{T4} = 1:size(data.lambda,1),
   blaze::AbstractBlazeType = make_blaze_1d(data, order=order),
   z_bc::AbstractVector{T3} = zeros(T1,num_obs(data)),
    template::AbstractArray{T2,1} = make_template_1d(data,log_lambda_out,order, z_bc=z_bc, pixel_range=pixel_range, blaze=blaze ), # zeros(T1,length(log_lambda_out)),
    num_star_basis_vectors::Integer = 0, num_telluric_basis_vectors::Integer = 0
    ) where { EOST <: AbstractEchelleObservationSet, T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Real}
    @assert num_obs(data) == length(z_bc)
    @assert length(log_lambda_out) == length(template)
    @assert 1<=order<=num_orders(data)
    # template = make_template_1d(data,log_lambda_out,order, z_bc=z_bc, pixel_range=pixel_range, blaze=blaze )
    model = SpectralOrderModel(log_lambda_out,
      star_basis = make_stellar_basis(template,num_basis_vectors=num_star_basis_vectors),
      telluric_basis = default_telluric_model(length(log_lambda_out),num_basis_vectors=num_telluric_basis_vectors),
      pixel_range=blaze.pixel_range, blaze=blaze )
   return model
#   return (pixel_range=pixel_range, blaze=blaze, log_lambda_range=log_lambda_out,template=template)

end

function make_spectral_model(data::EOST; #, log_lambda_out::AbstractRange{T1};
    pixel_range = 1:size(data.lambda,1),
    z_bc::AbstractVector{T3} = zeros(T1,num_obs(data)),
    num_star_basis_vectors::Integer = 0, num_telluric_basis_vectors::Integer = 0
    ) where { EOST <: AbstractEchelleObservationSet, T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Real}
    @assert num_obs(data) == length(z_bc)

    model = map( order->make_spectral_order_model(data,make_model_wavelength_range(data,order),order,
               pixel_range=pixel_range, z_bc=z_bc,
               num_star_basis_vectors=num_star_basis_vectors,
               num_telluric_basis_vectors=num_telluric_basis_vectors), 1:size(data.lambda,2) )
end

#=
import EchelleTools.AbstractEchelleObservationSet, EchelleTools.SpectralOrderModelAbstract
import EchelleTools.min_num_star_basis_vectors
import EchelleTools.max_num_star_basis_vectors
import EchelleTools.min_num_telluric_basis_vectors
import EchelleTools.max_num_telluric_basis_vectors
=#

function make_observations_model(data::DataT, som::AbstractArray{SOMT,1};
   orders::AbstractRange = 1:size(data.lambda,2),
   num_star_basis_vectors::Integer = 1,
   num_telluric_basis_vector::Integer = 0,
   z_bc::AbstractArray{T1,1} = zeros(num_obs),
   )  where {T1<:Real,
      DataT <: AbstractEchelleObservationSet, SOMT <: SpectralOrderModelAbstract
      }
   @assert min_num_star_basis_vectors <= num_star_basis_vectors <= max_num_star_basis_vectors
   @assert min_num_telluric_basis_vectors <= num_telluric_basis_vector <= max_num_telluric_basis_vectors
   @assert length(z_bc) == size(data.lambda,3)
   num_pixels = size(data.lambda,1)
   num_orders = size(data.lambda,2)
   num_obs = size(data.lambda,3)
   star_score = default_star(num_pixels,num_basis_vectors=num_star_basis_vectors)
   telluric_score = default_telluric_model(num_pixels,num_basis_vectors=num_telluric_basis_vector)
   normalization = zeros(num_obs)
   for i in 1:num_obs
      # TODO: specify pixel_range

      normalization[i] = mapreduce(
         order->mean( view(data.flux,:,order,i) ./
                        eval_model_baseline(view(data.lambda,:,order,i), som[order]) ),
                        +, orders) / length(orders)
   end
   ObservationsModel{eltype(star_score),eltype(telluric_score),eltype(z_bc),eltype(normalization)}(star_score,telluric_score,z_bc,normalization)
end

# Rest is no longer used?
function model_spectrum_1d(log_lambdas_template::AbstractRange{T1}, template::AbstractArray{T2,1},
      log_lambdas_out::Union{AbstractRange{T3},AbstractVector{T3}}, z_bc::Real;
      pixel_range::AbstractRange{T4} = 1:length(log_lambdas_out),
      blaze::Poly{T5} = Poly([one(T2)]),
      output::AbstractVector{T6} = zeros(T2,length(log_lambdas_out)) ) where {T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Real, T6<:Real}
   @assert length(log_lambdas_template) == length(template)
   @assert length(log_lambdas_out) == length(pixel_range) == length(output)

   fb = polyval(blaze, pixel_range )
   interpolator = scale(interpolate( template, BSpline(Linear())), log_lambdas_template)
   ll = log_lambda_ssb_from_obs.(log_lambdas_out, -z_bc)
   fm = interpolator.(ll)
   output .= fb.*fm
end


function model_spectrum_2d(log_lambdas_template::AbstractVector{ART1},
      template::AbstractArray{AVT2,1}, log_lambdas_out::AbstractArray{T3,2}, z_bc::Real;
      orders::Union{AbstractVector{T4},AbstractRange{T4}} = 1:size(log_lambdas_out,2),
      pixel_range::AbstractRange{T5} = 1:size(log_lambdas_out,1), # TODO: WARNING Assumes all orders have same pixel_range
      blaze::AbstractVector{Poly{T6}} = fill(Poly([one(eltype(eltype(AVT2)))]),size(log_lambdas_out,1)),
      output::AbstractArray{T7,2} = zeros(eltype(eltype(AVT2)),size(log_lambdas_out))
      ) where {T1<:Real, ART1<:AbstractRange{T1}, T2<:Real, AVT2<:AbstractVector{T2}, T3<:Real, T4<:Integer, T5<:Real, T6<:Real, T7<:Real}

   map(order->model_spectrum_1d(log_lambdas_template[order], template[order],
         view(log_lambdas_out,:,order), z_bc, pixel_range=pixel_range, blaze=blaze[order],
         output=view(output,:,order)), orders) # TODO: TEST
   return output
end


#=
ENV["PYTHON"] = "/usr/bin/python3"
using Pkg
Pkg.activate(".")
using Revise
using EchelleTools
using DataFrames, FileIO # JLD2 #
#JLD2.@load("../../Data/GJ_436.jld2")
data = load("../../Data/GJ_436.jld2")["data"]
res2d = make_blaze_and_template_2d(data)
orders = 1:12
llt = map(i->res2d[i].log_lambda_range,orders)
t = map(i->res2d[i].template,orders)
llo = log.(data.lambda[:,orders,1])
b = map(i->res2d[i].blaze,orders)
out = model_spectrum_2d(llt,t,llo,0.0,orders=orders,blaze=b )

using Polynomials, Interpolations
import EchelleTools.num_obs, EchelleTools.num_orders, EchelleTools.num_pixels_per_order
import EchelleTools.model_spectrum_2d
include("src/analyze.jl")

=#

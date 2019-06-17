using Interpolations

""" Model for one order of an echelle spectrum
Arrays organized into (basis_component,pixel).
"""
struct SpectralOrderModel{T1<:Real, T2<:Union{Real,Base.TwicePrecision{T1}}, T3<:Union{Real,Base.TwicePrecision{T1}},
      T4<:Real, T5<:Real, BlazeT<:AbstractBlaze1DType} <: AbstractSpectralOrderModel
   log_lambda::StepRangeLen{T1,T2,T3}
   star_basis::AbstractArray{T4,2}
   telluric_basis::AbstractArray{T5,2}
   blaze::BlazeT
   psf::AbstractPSF
   #=
   function SpectralOrderModel{T1,T2,T3,T4,T5}(log_lambda::LambdaT,
      star_basis::StarBasisT,
      telluric_basis::TelluricBasisT,
      blaze::BlazeT{T4,T5}
      )  where {    T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Real,
         LambdaT <: AbstractRange{T1}, StarBasisT <: AbstractArray{T2,2},
         TelluricBasisT <: AbstractArray{T3,2},
         #PixelRangeT <: AbstractRange{T4},
         BlazeT <: AbstractBlazeType
         }
      @assert length(log_lambda)==size(star_basis,2)==size(telluric_basis,2)
      @assert 1 <= size(star_basis,1) <= max_num_basis_vectors
      @assert 0 <= size(telluric_basis,1) <= max_num_basis_vectors
      new(log_lambda,star_basis,telluric_basis,blaze) # pixel_range,blaze)
   end
   =#
end

function SpectralOrderModel(data::ObsT, order::Integer;
   log_lambda::LambdaT, # TODO FIX = make_model_wavelength_range(data,order),
   star_basis::StarBasisT = default_stellar_template(length(log_lambda)), #num_pixels_per_order(data)),
   telluric_basis::TelluricBasisT = default_telluric_template(length(log_lambda)), #num_pixels_per_order(data)),
   pixel_range::PixelRangeT = 1:num_pixels_per_order(data),
   degree::Integer = default_blaze_degree,
   blaze::BlazeT = make_flat_blaze_polynomial(data,order,degree=degree,pixel_range=pixel_range),
   psf::PSFT = PSFDeltaFunction()
   )  where {    T1<:Real, T2<:Real, T3<:Real, T4<:Real,
      ObsT <: AbstractEchelleObservationOrSet,
      LambdaT <: AbstractRange{T1},
      StarBasisT <: AbstractArray{T2,2}, TelluricBasisT <: AbstractArray{T3,2},
      PixelRangeT <: AbstractRange{T4}, BlazeT <: AbstractBlazeType,
      PSFT <: AbstractPSF
      }
   # TODO: move asserts to inner constructor?
   @assert 1 <= order <= num_orders(data)
   @assert 1 <= length(log_lambda) <= max_template_length
   @assert size(star_basis,2) == size(telluric_basis,2)
   @assert 1 <= size(star_basis,2)
   @assert min_num_star_basis_vectors <= size(star_basis,1) <= max_num_star_basis_vectors
   @assert min_num_telluric_basis_vectors <= size(telluric_basis,1) <= max_num_telluric_basis_vectors
   @assert blaze.order == order

   SpectralOrderModel(log_lambda,star_basis,telluric_basis,blaze,psf)
end

function get_order(som::SOMT) where { SOMT <: AbstractSpectralOrderModel }
   som.blaze.order
end

""" Model for one order of an echelle spectrum
Arrays organized into (basis_component,pixel).
"""
struct SpectralModel{T1<:Real, T2<:Real, T3<:Real, BlazeT<: AbstractBlaze2DType } <: AbstractSpectralModel
   log_lambda::StepRangeLen{T1,Base.TwicePrecision{T1},Base.TwicePrecision{T1} }
   star_basis::AbstractArray{T2,2}
   telluric_basis::AbstractArray{T3,2}
   blaze_list::BlazeT
   psf::AbstractPSF
end

function SpectralModel(data::ObsT;
   log_lambda::LambdaT = make_model_wavelength_range(data),
   star_basis::StarBasisT = default_stellar_template(length(log_lambda)), # num_orders(data)*num_pixels_per_order(data)),
   telluric_basis::TelluricBasisT = default_telluric_template(length(log_lambda)), #num_orders(data)*num_pixels_per_order(data)),
   pixel_range::PixelRangeT = 1:num_pixels_per_order(data),
   degree::Integer = default_blaze_degree,
   blaze_list::BlazeT = make_flat_blaze_polynomial_2d(data,pixel_rage=pixel_range,degree=degree),
      #map(order->make_flat_blaze(data,order,degree=degree,pixel_range=pixel_range),1:num_orders(data)),
   psf::PSFT = PSFDeltaFunction()
   )  where {    T1<:Real, T2<:Real, T3<:Real, T4<:Real,
      ObsT <: AbstractEchelleObservationOrSet,
      LambdaT <: AbstractRange{T1},
      StarBasisT <: AbstractArray{T2,2}, TelluricBasisT <: AbstractArray{T3,2},
      PixelRangeT <: AbstractRange{T4}, BlazeT <: AbstractBlaze2DType, # BT<: AbstractBlazeType, BlazeT <: AbstractArray{BT,1},
      PSFT <: AbstractPSF
      }
   # TODO: Move asserts to inner constructor?
   @assert 1 <= length(log_lambda) <= max_template_length
   @assert length(log_lambda)==size(star_basis,2) == size(telluric_basis,2)
   @assert 1 <= size(star_basis,1)
   @assert min_num_star_basis_vectors <= size(star_basis,1) <= max_num_star_basis_vectors
   @assert min_num_telluric_basis_vectors <= size(telluric_basis,1) <= max_num_telluric_basis_vectors
   @assert all( 1 .<= map(b->b.order,blaze_list.blaze_list) .<= num_orders(data) )
   @assert length(blaze_list.blaze_list) == num_orders(data)  # TODO: Could relax, but would need to remove assumption index=blaze.order

   SpectralModel(log_lambda,star_basis,telluric_basis,blaze_list,psf)
end

function SpectralOrderModel(model::SMT, order::Integer; max_abs_z::Real = 1e-4 ) where { SMT <: AbstractSpectralModel}
   @assert 1 <= order <= length(model.blaze_list.blaze_list)
   log_lambda_lo = min_log_lambda(model.blaze_list.blaze_list[order]) - log1p(max_abs_z)
   log_lambda_hi = max_log_lambda(model.blaze_list.blaze_list[order]) + log1p(max_abs_z)
   # TODO OPT: Above use model.log_lambda's step to figure out how much to move idx without taking log's
   idx_lo = log_lambda_lo > model.log_lambda[1] ? findfirst(l->l>=log_lambda_lo, model.log_lambda) : 1
   idx_hi = log_lambda_hi < model.log_lambda[end] ? findlast(l->l<=log_lambda_hi, model.log_lambda) : length(model.log_lambda)
   len = idx_hi-idx_lo+1
   #log_lambda_narrow = StepRangeLen{Float64,Float64,Base.TwicePrecision{Float64}}(model.log_lambda[1],model.log_lambda.step,len,idx_lo-1)
   log_lambda_narrow = StepRangeLen{Float64,Float64,Base.TwicePrecision{Float64}}(model.log_lambda[idx_lo],model.log_lambda.step,len)
   #log_lambda_narrow = range(model.log_lambda[idx_lo],step=model.log_lambda.step,length=len)
   #log_lambda_narrow = range(model.log_lambda[idx_lo],model.log_lambda.step,len)
   #log_lambda_narrow = StepRangeLen(model.log_lambda[1],model.log_lambda.step,len,idx_lo-1)
   SpectralOrderModel(log_lambda_narrow,view(model.star_basis,:,idx_lo:idx_hi),view(model.telluric_basis,:,idx_lo:idx_hi),model.blaze_list.blaze_list[order],model.psf)
end

# Default templates for star and tellurics
""" default_telluric_template(n; num_basis_vectors = 0)
   Returns array with (num_basis_vectors, n) zeros. """
function default_telluric_template(n::Integer; num_basis_vectors::Integer = 0)
   @assert 0 <= min_num_telluric_basis_vectors <= num_basis_vectors <= max_num_telluric_basis_vectors
   @assert 1 <= n <= max_template_length
   zeros(num_basis_vectors,n)
end

""" default_stellar_template(n; num_basis_vectors = 0)
   Returns array with one row of n 1's followed by num_basis_vectors rows of zeros. """
function default_stellar_template(n::Integer; num_basis_vectors::Integer = 0)
   @assert 0 <= min_num_star_basis_vectors-1 <= num_basis_vectors <= max_num_star_basis_vectors
   @assert 1 <= n <= max_template_length
   result = vcat(ones(1,n),zeros(num_basis_vectors,n))
   return result
end

function make_model_wavelength_range(data::EOT, order::OrdersType;
            oversample_factor::Real = 1.0, max_abs_z::Real = 1e-4,
            num_knots::Integer=round(Int64,oversample_factor*(1+max_abs_z)^2*num_pixels_per_order(data))
            ) where { EOT <: AbstractEchelleObservation }
      min_log_lambda = log(minimum(data.lambda[1,order])  / (1+max_abs_z) )
      max_log_lambda = log(maximum(data.lambda[end,order]) * (1+max_abs_z) )
      range(min_log_lambda, stop=max_log_lambda, length=num_knots )
end

function make_model_wavelength_range(data::EOST, order::OrdersType;
            oversample_factor::Real = 1.0, max_abs_z::Real = 1e-4,
            num_knots::Integer=round(Int64,oversample_factor*(1+max_abs_z)^2*num_pixels_per_order(data))
            ) where { EOST <: AbstractEchelleObservationSet }
      min_log_lambda = log(minimum(data.lambda[1,order])  / (1+max_abs_z) )
      max_log_lambda = log(maximum(data.lambda[end,order]) * (1+max_abs_z) )
      range(min_log_lambda, stop=max_log_lambda, length=num_knots )
end

function make_model_wavelength_range(data::EOST;
            oversample_factor::Real = 1.0, max_abs_z::Real = 1e-4,
            num_knots::Integer=round(Int64,oversample_factor*(1+max_abs_z)^2*num_orders(data)*num_pixels_per_order(data))
            ) where { EOST <: AbstractEchelleObservationSet }
      min_log_lambda = log(minimum(data.lambda[1,1,:])      / (1+max_abs_z) )
      max_log_lambda = log(maximum(data.lambda[end,end,:])  * (1+max_abs_z) )
      range(min_log_lambda, stop = max_log_lambda, length = num_knots )
end


""" Model for one echelle observation (potentially multiple orders)
Arrays over (basis_component).
"""
struct ObservationModel{T1<:Real, T2<:Real, T3<:Real, T4<:Real } <: AbstractObservationModel
   star_score::Array{T1,1}
   telluric_score::Array{T2,1}
   z_bc::T3
   normalization::T4

   function ObservationModel{T1,T2,T3,T4}( star_score::AbstractArray{T1,1},
      telluric_score::AbstractArray{T2,1}, z_bc::T3, normalization::T4
      )  where {T1<:Real, T2<:Real, T3<:Real, T4<:Real }
      @assert min_num_star_basis_vectors <= size(star_score,1) <= max_num_star_basis_vectors
      @assert min_num_telluric_basis_vectors <= size(telluric_score,1) <= max_num_telluric_basis_vectors
      new(star_score,telluric_score,z_bc,normalization)
   end

end
export ObservationModel

function ObservationModel(;
   num_star_basis_vectors::Integer = 1,
   num_telluric_basis_vector::Integer = 0,
   star_score::AbstractArray{T1,1} = zeros(num_star_basis_vectors),
   telluric_score::AbstractArray{T2,1} = zeros(num_telluric_basis_vector),
   z_bc::T3 = 0.0,
   normalization::T4 = 1.0
   )  where {T1<:Real, T2<:Real, T3<:Real, T4<:Real }
   ObservationModel{eltype(star_score),eltype(telluric_score),typeof(z_bc),typeof(normalization)}(star_score,telluric_score,z_bc,normalization)
end

""" Model for a time series of echelle observations (potentially multiple orders)
Arrays of scores over (basis_component,observation).
Arrays of z and normalization over (observation).
"""
struct ObsSeriesModel{T1<:Real, T2<:Real, T3<:Real, T4<:Real } <: AbstractObsSeriesModel
   star_score::Array{T1,2}
   telluric_score::Array{T2,2}
   z_bc::Array{T3,1}
   normalization::Array{T4,1}
   function ObsSeriesModel{T1,T2,T3,T4}( star_score::AbstractArray{T1,2},
      telluric_score::AbstractArray{T2,2}, z_bc::AbstractArray{T3,1}, normalization::AbstractArray{T4,1}
      )  where {T1<:Real, T2<:Real, T3<:Real, T4<:Real }
      @assert min_num_star_basis_vectors <= size(star_score,1) <= max_num_star_basis_vectors
      @assert min_num_telluric_basis_vectors <= size(telluric_score,1) <= max_num_telluric_basis_vectors
      @assert size(star_score,2) == size(telluric_score,2) == length(z_bc) == length(normalization)
      @assert length(z_bc) <= max_num_obs
      new(star_score,telluric_score,z_bc,normalization)
   end
end
export ObsSeriesModel

function ObsSeriesModel(num_obs::Integer;
   num_star_basis_vectors::Integer = 1,
   num_telluric_basis_vector::Integer = 0,
   star_score::AbstractArray{T1,2} = vcat(ones(1,num_obs),zeros(num_star_basis_vectors-1,num_obs)),
   telluric_score::AbstractArray{T2,2} = zeros(num_telluric_basis_vector,num_obs),
   z_bc::AbstractArray{T3,1} = zeros(num_obs),
   normalization::AbstractArray{T4,1} = ones(num_obs)
   )  where {T1<:Real, T2<:Real, T3<:Real, T4<:Real }
   @assert min_num_star_basis_vectors <= size(star_score,1) <= max_num_star_basis_vectors
   @assert min_num_telluric_basis_vectors <= size(telluric_score,1) <= max_num_telluric_basis_vectors
   @assert size(star_score,2) == size(telluric_score,2) == length(z_bc) == length(normalization)

   ObsSeriesModel{eltype(star_score),eltype(telluric_score),eltype(z_bc),eltype(normalization)}(star_score,telluric_score,z_bc,normalization)
end

function ObservationModel(obs_models::ObsSeriesModel, i::Integer) where {ObsSeriesModel <: AbstractObsSeriesModel}
   @assert 1 <= i <= length(obs_models.normalization)
   ObservationModel( star_score=view(obs_models.star_score,:,i), telluric_score=view(obs_models.telluric_score,:,i), z_bc=obs_models.z_bc[i], normalization=obs_models.normalization[i] )
end

function num_obs(osm::OSMT) where { OSMT <: AbstractObsSeriesModel }
  return length(osm.normalization)
end

#=
function eval_model_baseline(log_lambda_obs::LambdaT,
      spec_model::SOMT;
      pixel_range_obs::AbstractRange = 1:length(log_lambda_obs),
      output::OutputT = zeros(eltype(spec_model.star_basis),length(log_lambda_obs))
      ) where {  T1 <: Real, T2 <: Real,
      LambdaT <: AbstractArray{T1,1} #= Wavelength1DType{T1} =#, SOMT <: AbstractSpectralOrderModel,
      OMT <: AbstractObservationModel, OutputT <: AbstractVector{T2} }
   @assert length(log_lambda_obs) == length(pixel_range_obs)
   @assert length(log_lambda_obs) == length(output)
   # Compute each term on a grid
   blaze = calc_blaze(spec_model.blaze,pixel_range_obs) # TODO: Use sub-pixels once integrating over PSF
   star_on_grid = spec_model.star_basis[1,:] # * obs.star_score
   interp_star = scale(interpolate( star_on_grid, BSpline(Cubic(Line(OnGrid())))), spec_model.log_lambda)
   # Interpolate each term to observed wavelength
   flux_star = interp_star.(log_lambda_obs)
   telluric_factor = 1.0
   # TODO: Replace product of interpolations with integration over PSF
   output .= obs.normalization .* blaze .* telluric_factor .* flux_star
end
=#

function eval_model_star(spec_model::SOMT, obs_model::OMT;
      log_lambda_obs::LambdaT = spec_model.log_lambda, # view(spec_model.log_lambda,2:(length(spec_model.log_lambda)-1)),
      output::OutputT = zeros(eltype(spec_model.star_basis),length(log_lambda_obs))
      ) where {  T1 <: Real, T2 <: Real,
      SOMT <: AbstractSpectralOrderModel, OMT <: AbstractObservationModel,
      LambdaT <: AbstractArray{T1,1} #=Wavelength1DType{T1}=#, #= PixelRangeT <: AbstractRange, =#
      OutputT <: AbstractVector{T2} }
   @assert length(log_lambda_obs) == length(output)
   @assert size(spec_model.star_basis,1) == size(obs_model.star_score,1)
   @assert size(spec_model.telluric_basis,1) == size(obs_model.telluric_score,1)
   # Compute each term on a grid
   star_on_grid = spec_model.star_basis' * obs_model.star_score
   # TODO: Replace with my interpolator to avoid having to allocate memory to store values at knots
   interp_star = scale(extrapolate(interpolate( star_on_grid, BSpline(Cubic(Line(OnGrid())))),Line()), spec_model.log_lambda)
   # Interpolate each term to observed wavelength
   log_lambda_source = log_lambda_ssb_from_obs.(log_lambda_obs, -obs_model.z_bc)
   flux_star = interp_star.(log_lambda_source)
end

function eval_model_telluric(spec_model::SOMT, obs_model::OMT;
      log_lambda_obs::LambdaT = spec_model.log_lambda,
      output::OutputT = zeros(eltype(spec_model.star_basis),length(log_lambda_obs))
      ) where {  T1 <: Real, T2 <: Real,
      SOMT <: AbstractSpectralOrderModel, OMT <: AbstractObservationModel,
      LambdaT <: AbstractArray{T1,1} #=Wavelength1DType{T1}=#, #= PixelRangeT <: AbstractRange, =#
      OutputT <: AbstractVector{T2} }
   @assert length(log_lambda_obs) == length(output)
   @assert size(spec_model.star_basis,1) == size(obs_model.star_score,1)
   # Compute each term on a grid
   telluric_on_grid = 1.0 .- spec_model.telluric_basis' * obs_model.telluric_score
   # TODO: Replace with my interpolator to avoid having to allocate memory to store values at knots
   interp_telluric = scale(interpolate( telluric_on_grid, BSpline(Cubic(Line(OnGrid())))), spec_model.log_lambda)
   # Interpolate each term to observed wavelength
   telluric_factor = interp_telluric.(log_lambda_obs)
end

function eval_model_star_telluric(spec_model::SOMT, obs_model::OMT;
      log_lambda_obs::LambdaT = spec_model.log_lambda,
      output::OutputT = zeros(eltype(spec_model.star_basis),length(log_lambda_obs))
      ) where {  T1 <: Real, T2 <: Real,
      SOMT <: AbstractSpectralOrderModel, OMT <: AbstractObservationModel,
      LambdaT <: AbstractArray{T1,1} #=Wavelength1DType{T1}=#, OutputT <: AbstractVector{T2} }
   @assert length(log_lambda_obs) == length(output)
   @assert size(spec_model.telluric_basis,1) == size(obs_model.telluric_score,1)
   output .= eval_model_star(spec_model,obs_model,log_lambda_obs=log_lambda_obs) .*
      eval_model_telluric(spec_model,obs_model,log_lambda_obs=log_lambda_obs) .*
      obs_model.normalization
end

function eval_model_blaze(spec_model::SOMT #=, obs_model::OMT=#;
      pixel_range_obs::PixelRangeT = spec_model.blaze.pixel_range,
      output::OutputT = zeros(eltype(spec_model.star_basis),length(pixel_range_obs))
      ) where {  #= T1 <: Real, =# T2 <: Real,
      SOMT <: AbstractSpectralOrderModel, OMT <: AbstractObservationModel,
      PixelRangeT <: AbstractRange, OutputT <: AbstractVector{T2} }
   @assert length(pixel_range_obs) == length(output)
   output .= calc_blaze(spec_model.blaze,pixel_range_obs)
end

function eval_model(spec_model::SOMT, obs_model::OMT;
      #pixel_range_obs::PixelRangeT = spec_model.blaze.pixel_range,
      log_lambda_obs::LambdaT = spec_model.log_lambda,
      output::OutputT = zeros(eltype(spec_model.star_basis),length(log_lambda_obs))
      ) where { T1 <: Real, T2 <: Real,
      SOMT <: AbstractSpectralOrderModel, OMT <: AbstractObservationModel,
      LambdaT <: AbstractArray{T1,1} #=Wavelength1DType{T1}=# ,  #=PixelRangeT <: AbstractRange, =#
      OutputT <: AbstractVector{T2} }
   @assert length(output) == length(log_lambda_obs)
   # @assert length(log_lambda_obs) == length(spec_model.blaze.pixel_range) # ERROR TODO: UNDERSTAND WHY?
   #blaze = eval_model_blaze(spec_model)
   blaze = calc_blaze_at_log_lambda(spec_model.blaze, log_lambda_obs)
   incoming_flux = eval_model_star_telluric(spec_model,obs_model,log_lambda_obs=log_lambda_obs)
   output .= blaze .* incoming_flux
end

function eval_model(spec_model::SMT, obs_model::OMT, log_lambda_obs::LambdaT;
      #= pixel_range_obs::PixelRangeTA = fill(1:size(log_lambda_obs,1),size(log_lambda_obs,2)), =#
      output::OutputT = zeros(eltype(spec_model.star_basis),size(log_lambda_obs))
      ) where {  T1 <: Real, #= T2 <: Real, =# T3 <: Real,
      LambdaT <: AbstractArray{T1,2}, SMT <: AbstractSpectralModel,
      OMT <: AbstractObservationModel, #= PixelRangeT <: AbstractRange{T2},
      PixelRangeTA <: AbstractArray{PixelRangeT,1}, =#
      OutputT <: AbstractArray{T3,2} }
   @assert size(output,2) == size(log_lambda_obs,2)
   map(i->eval_model(SpectralOrderModel(spec_model,i),obs_model,
               log_lambda_obs=view(log_lambda_obs,:,i),
               #= pixel_range_obs=pixel_range_obs[i], =#
               output=view(output,:,i)), 1:size(log_lambda_obs,2) )
   return output
end

function eval_model(spec_model::SMT, obs_series_model::OSMT, log_lambda_obs::LambdaT;
      pixel_range_obs::PixelRangeTA = fill(1:size(log_lambda_obs,1),size(log_lambda_obs,2)),
      output::OutputT = zeros(eltype(spec_model.star_basis),size(log_lambda_obs))
      ) where {  T1 <: Real, T2 <: Real, T3 <: Real,
      LambdaT <: AbstractArray{T1,3}, SMT <: AbstractSpectralModel,
      OSMT <: AbstractObsSeriesModel, PixelRangeT <: AbstractRange{T2},
      PixelRangeTA <: AbstractArray{PixelRangeT,1},
      OutputT <: AbstractArray{T3,3} }
   @assert size(output,3) == size(log_lambda_obs,3)
   map(obsid -> eval_model(spec_model,ObservationModel(obs_series_model,obsid),
                           view(log_lambda_obs,:,:,obsid),
                           #pixel_range_obs=pixel_range_obs,
                           output=view(output,:,:,obsid)),
               1:num_obs(obs_series_model) )

   return output
end

#=
function eval_model(log_lambda_obs::LambdaT,
      spec_model::AbstractArray{SOMT,1}, obs_model::OMT;
      pixel_range_obs::AbstractArray{AbstractRange,1} = fill(1:size(log_lambda_obs,1),size(log_lambda_obs,2)),
      output::OutputT = zeros(eltype(model.star_basis),size(log_lambda_obs,1),size(log_lambda_obs,2))
      ) where {  T1 <: Real, T2 <: Real,
      LambdaT <: AbstractArray{AbstractRange{T1},1}, SOMT <: AbstractSpectralOrderModel,
      OMT <: AbstractObservationModel, OutputT <: AbstractArray{T2,2} }
   map(i->eval_model(log_lambda_obs[i],spec_model[i],obs_model,
         pixel_range=pixel_range_obs[i],output=view(output,:,i)), 1:size(log_lambda_obs,2) )
   return output
end
=#

function make_stellar_template(data::EOST, sm::SMT, obsm::OSMT
   ) where { EOST <: AbstractEchelleObservationSet, SMT <: AbstractSpectralModel, OSMT <: AbstractObsSeriesModel }
   model = eval_model(log.(data.lambda), sm, obsm)
   template = reshape(mean(data.flux ./ model, dims=3), num_pixels_per_order(data), num_orders(data) )
end

function eval_model_fit(data::EOST, spec_model::SMT, obs_model::OSMT
      #;output::OutputT = zeros(eltype(spec_model.star_basis),size(data.lambda))
      ) where {  T2<:Real, T3 <: Real,
      EOST <: AbstractEchelleObservationSet, SMT <: AbstractSpectralModel,
      OSMT <: AbstractObsSeriesModel, PixelRangeT <: AbstractRange{T2},
      PixelRangeTA <: AbstractArray{PixelRangeT,1},
      OutputT <: AbstractArray{T3,2} }
   sum_sq = 0.
   for obsid in 1:num_obs(data)
      obsm = ObservationModel(obs_model,obsid)
      for i in 1:num_orders(data)
         pixel_range = spec_model.blaze_list.blaze_list[i].pixel_range
         sum_sq += sum( ( data.flux[pixel_range] .-
                  eval_model(log.(view(data.lambda,pixel_range,i,obsid)),
                     SpectralOrderModel(spec_model,i),obsm,pixel_range=pixel_range ) ).^2 ./
                        data.var[pixel_range,i,obsid] )
      end
   end
   return sum_sq
end

#=
function eval_model_fit(data::EOST, sm::SMT, obsm::OSMT
   ) where { EOST <: AbstractEchelleObservationSet, SMT <: AbstractSpectralModel, OSMT <: AbstractObsSeriesModel }
   model = eval_model(log.(data.lambda), sm, obsm)
   data.flux .- model
   template = reshape(mean(data.flux ./ model, dims=3), num_pixels_per_order(data), num_orders(data) )
end
=#

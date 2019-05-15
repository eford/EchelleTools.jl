#module HARPSNData
#export read_data

import ..EchelleObservation, ..EchelleObservationSet, ..OrdersType
import ..AbstractInstrument, ..HARPSNType
#import ..harpsn_flux_hdu,  ..harpsn_all_orders
using FITSIO

"""
read_data(filename; orders)
Example:
   data = read_data(filename_list[1])
"""
function read_data(filename::String; orders::OrdersType = harpsn_all_orders )
   local fits
   try
      fits = FITS(filename)
   catch
      @error("# Can't read " * filename)
   end
   if typeof(orders) <: Integer
      orders = orders:orders
   end
   hdr = read_header(fits[1])
   # Get header info for estimating flux uncertainties
   gain = hdr["TNG DRS CCD CONAD"]
   ron = hdr["TNG DRS CCD SIGDET"]
   # Get barycentric Earth RV correction (should be positive in morning when observer is approaching sun)
   berv = hdr["TNG DRS BERV"]
   berv /= speed_of_light_mps
   flux = read(fits[harpsn_flux_hdu])[:,orders]
   flux *= gain # in ADU
   #flux .*= gain
   var = copy(flux)
   flux[flux .< zero(eltype(flux))] .= zero(eltype(flux))
   var .+= ron^2
   npix,nord = size(flux)
   # Extract pixel->wavelegnth polynomials (thank ACC)
   coeff = zeros(4,nord)
   for ord in 0:(nord-1)
      for poly in 0:3
         n = ord*4+poly
         coeff[1+poly,1+ord] = hdr["TNG DRS CAL TH COEFF LL" * string(n)]
         #push!(coeff,hdr["TNG DRS CAL TH COEFF LL" * string(n)])
         #n += 1
      end
   end
   xpix = 0:(npix-1)
   xarray =  vcat(ones(length(xpix))', xpix', (xpix.^2)', (xpix.^3)' )
   lambda = zeros(npix,length(orders))
   for ord in orders
      lambda[:,ord] .= xarray' * coeff[:,ord]
   end
   metadata = Dict{Symbol,Any}(:Filename=>filename,:Instrument=>:HARPSN,:z_bc=>berv)
   return EchelleObservation{eltype(lambda),eltype(flux),eltype(var),typeof(berv),HARPSNType}(flux=flux,var=var,lambda=lambda,z_bc=berv,metadata=metadata)
end

"""
   read_data!(data_out, obsid, filename; orders)
Writes to data_out[:,:,i]'s flux, var and lambda.
Creates metadata_list including :Filename and :z_bc (redshift of observer relative to solar system barcenter)

Example:
   data = read_data!(data_out, obsid, filename, orders=10:15)
"""
function read_data!(data::EchelleObservationSet{T1,T2,T3,#= T4, =# HARPSNType}, i::Int, filename::String;
      orders::OrdersType = harpsn_all_orders
      ) where {T1<:Number, T2<:Number, T3<:Number, #=T4<:Number,  HARPSNType <:AbstractInstrument =#}
   @assert 1 <= i <= size(data.flux,3)
   local fits
   try
      fits = FITS(filename)
   catch
      @error("# Can't read " * filename)
   end
   if typeof(orders) <: Integer
      orders = orders:orders
   end
   obs = read_data(filename, orders=orders)
   data.flux[:,:,i] .= obs.flux
   data.var[:,:,i] .= obs.var
   data.lambda[:,:,i] .= obs.lambda
   data.metadata_list[i] = Dict{Symbol,Any}(:Filename=>filename,:z_bc=>obs.metadata[:z_bc])
   return data
end


"""
   read_data(filenames; orders, z_bsrv)
Returns struct with flux, var, lambda, and metadata_list including :Filename and :z_bc
If pass z_bsrv, then adds that to the barrycentric correction stored in metadata_list.

Example:
   data = read_data(filename_list, orders=10:15)
"""
function read_data(filename_list::AbstractArray{String,1};
      orders::OrdersType = harpsn_all_orders,
      z_bsrv::AbstractArray = zeros(length(filename_list)) )
   @assert length(filename_list) >= 1
   data_first = read_data(filename_list[1], orders=orders)
   data = EchelleObservationSet{Float64,Float64,Float64,HARPSNType}(size(data_first.flux), length(filename_list), metadata=Dict(:Instrument=>:HARPSN) )
   data.lambda[:,:,1] = data_first.lambda
   data.flux[:,:,1]   = data_first.flux
   data.var[:,:,1]    = data_first.var
   for (i,filename) in enumerate(filename_list)
          read_data!(data, i, filename, orders=orders)
          data.metadata_list[i][:z_bc] += z_bsrv[i]
   end
   return data
end

#end # module

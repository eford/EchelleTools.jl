module HPFData
export read_data

import ..EchelleObservation, ..EchelleObservationSet, ..OrdersType
import ..hpf_flux_hdu, ..hpf_var_hdu, ..hpf_lambda_hdu, ..hpf_all_orders
using FITSIO

#const OrdersType = Union{Int64,UnitRange{Int64},AbstractVector{Int64}}

"""
read_data(filename; orders)
Example:
   data = read_data(filename_list[1])
"""
function read_data(filename::String; orders::OrdersType = hpf_all_orders )
   fits = FITS(filename)
   if typeof(orders) <: Integer
      orders = orders:orders
   end
   flux = read(fits[hpf_flux_hdu])[:,orders]
   var = read(fits[hpf_var_hdu])[:,orders]
   lambda = read(fits[hpf_lambda_hdu])[:,orders]
   metadata = Dict{Symbol,Any}(:Filename=>filename,:Instrument=>:HPF)
   return EchelleObservation{eltype(lambda),eltype(flux),eltype(var)}(flux=flux,var=var,lambda=lambda,metadata=metadata)
end

function read_data!(data::EchelleObservationSet{T1,T2,T3}, i::Int, filename::String; orders::OrdersType = hpf_all_orders ) where {T1<:Number, T2<:Number, T3<:Number}
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
   data.flux[:,:,i] = read(fits[hpf_flux_hdu])[:,orders]
   data.var[:,:,i] = read(fits[hpf_var_hdu])[:,orders]
   data.lambda[:,:,i] = read(fits[hpf_lambda_hdu])[:,orders]
   data.metadata_list[i] = Dict{Symbol,Any}(:Filename=>filename)
   return data
end

"""
read_data(filenames; orders)
Example:
   data = read_data(filename_list, orders=10:15)
"""
function read_data(filename_list::AbstractArray{String,1}; orders::OrdersType = hpf_all_orders )
   @assert length(filename_list) >= 1
   data_first = read_data(filename_list[1], orders=orders)
   data = EchelleObservationSet{Float64,Float32,Float32}(size(data_first.flux), length(filename_list), metadata=Dict(:Instrument=>:HPF) )
   data.lambda[:,:,1] = data_first.lambda
   data.flux[:,:,1]   = data_first.flux
   data.var[:,:,1]    = data_first.var
   for (i,filename) in enumerate(filename_list)
          read_data!(data, i, filename, orders=orders)
   end
   return data
end

#=
read_hpf_data(filename::String; orders::OrdersType = hpf_all_orders ) = read_data(filename, orders=orders)
read_hpf_data(filename_list::AbstractArray{String,1}; orders::OrdersType = hpf_all_orders ) = read_data(filename_list, orders=orders)
=#

end # module



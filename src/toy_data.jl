#function spectral_line(x::AbstractArray; center::Real=0.0, width::Real=0.0, depth::Real=1.0, max_dist::Real=5.0 )
#   return spectral_line.(x,center=center,width=width,depth=depth, max_dist=max_dist)
#end


function spectral_line(x::Real; center::Real=0.0, width::Real=0.0, depth::Real=1.0, max_dist::Real=5.0 )
   return (x-center)/width > max_dist ? 1.0 : 1.0-depth*exp(-0.5*((x-center)/width)^2)
end

function spectral_line(x::Real; center::AbstractArray, width::AbstractArray, depth::AbstractArray, max_dist::Real=5.0 )
   @assert length(center) == length(width) == length(depth)
   output = one(typeof(x))
   for i in 1:length(center)
      output *= spectral_line(x,center=center[i],width=width[i],depth=depth[i], max_dist=max_dist)
   end
   return output
end

function make_toy_data(obs_size::Tuple{Int64,Int64},num_obs::Int = 0)
   instrument = :HPF
   data = EchelleObservationSet(obs_size, num_obs, instrument=instrument)
   lambda_min = 4000.0
   lambda_max = 8000.0
   num_lines = 1000
   line_widths = fill(5.0,num_lines)
   line_depths = fill(0.5,num_lines)
   line_centers = lambda_min.+(lambda_max-lambda_min).*rand(num_lines)

   snr = 100
   (num_obs, num_orders) = obs_size
   for i in 1:num_orders
      data.lambda[:,i,:] .= range(lambda_min+(lambda_max-lambda_min)*(i-1)/num_orders,
         stop=lambda_min+(lambda_max-lambda_min)*(i)/num_orders,length=obs_size[1])
   end
   data.flux .= spectral_line.(data.lambda, center=line_centers, width=line_widths, depth=line_depths, max_dist=5.0)
   data.var .= 1.0/snr
   data.flux .*= 1.0 .+ randn(size(data.flux))./snr
   return data
end

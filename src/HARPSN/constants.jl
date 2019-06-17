import ..speed_of_light_mps

const harpsn_all_orders = 1:69
const harpsn_flux_hdu = 1

const harpsn_rv_orders = vcat(collect(1:54), collect(56:60), collect(62:63), collect(66:68) )

default_blaze_degree = 4
max_blaze_degree = 6
default_pixel_buffer_lo = 64
default_pixel_buffer_hi = 8

#const harpsn_var_hdu = 5   # computed from flux and ron
#const harpsn_lambda_hdu = 8 # computed from polynomial coefficients in heard

#= Since in CSV w/ header, don't need to specify columns
const harpsn_manifest_format = [
		    ManifestFormatEntry(:UTTimestamp,String, 4:26),
                    ManifestFormatEntry(:UTDate,String,     30:37),
                    ManifestFormatEntry(:ObsNum,Int,        43:46),
                    ManifestFormatEntry(:Frame,String,      50:89),
                    ManifestFormatEntry(:ObsType,String,    97:99),
                    ManifestFormatEntry(:iTime,Float64,    100:117),
                    ManifestFormatEntry(:QProg,String,     118:145),
                    ManifestFormatEntry(:Object,String,    146:190)
                    ]
=#

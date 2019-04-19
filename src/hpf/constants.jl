const hpf_all_orders = 1:28
const hpf_flux_hdu = 2
const hpf_var_hdu = 5
const hpf_lambda_hdu = 8

const hpf_manifest_format = [ 
		    ManifestFormatEntry(:UTTimestamp,String, 4:26),
                    ManifestFormatEntry(:UTDate,String,     30:37),
                    ManifestFormatEntry(:ObsNum,Int,        43:46),
                    ManifestFormatEntry(:Frame,String,      50:89),
                    ManifestFormatEntry(:ObsType,String,    97:99),
                    ManifestFormatEntry(:iTime,Float64,    100:117),
                    ManifestFormatEntry(:QProg,String,     118:145),
                    ManifestFormatEntry(:Object,String,    146:190)
                    ]


# EchelleTools.jl

## Install
Using julia >v1.0
```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/eford/EchelleTools.jl"))
```

## Tools Avaliable
- Read Manifest Files
```julia
df = read_hpf_manifest_files() # Read them all!
df = read_hpf_manifest_files(filter_object=r"GJ.699",filter_obstype=r"Sci")   # Filter manifests
```
- Manipulate Manifest in Dataframe form
```julia
df_gj699 = df[findall(s->occursin(r"GJ.*699",s),df[:Object]),:]
```
- Translate manifest into filenames
```julia
filename_list = get_hpf_data_filenames(df)
```

- Extract data
```julia
data = read_hpf_data(filename_list[1])
data_3d = read_hpf_data(filename_list,orders=5:20)
```

- Analyze data
```julia
mean_spec = mean_spectrum(data)
```


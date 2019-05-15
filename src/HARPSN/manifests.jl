module Manifests

import ..base_path, ..base_path_default, ..manifest_subdir_default, ..data_subdir_default
import ..speed_of_light_mps

export read
export get_data_filenames, get_data_filenames_and_bsrv
export get_base_path, set_base_path!, get_manifest_path, set_manifest_path!, get_data_path, set_data_path!

using CSV
using FITSIO
using DataFrames
using Dates
using Printf

# Specify paths to manifest and data files
manifest_path_default = joinpath(base_path_default,manifest_subdir_default)
data_path_default = joinpath(base_path_default,data_subdir_default)

get_manifest_path() = manifest_path_default
set_manifest_path!(path::String) = manifest_path_default = path

get_data_path() = data_path_default
set_data_path!(path::String) = data_path_default = path

get_base_path() = base_path
function set_base_path!(path::String)
   base_path = path
   set_manifest_path!( joinpath(base_path,manifest_subdir_default) )
   set_data_path!( joinpath(base_path,data_subdir_default) )
end

# Initialize base_path
set_base_path!(base_path_default)

# Read in manifest files
function read(filename::String; dvtol=0.0001, qthresh=0.99, date_start::Date, date_stop::Date = date_start, verbose = false)
   df = CSV.read(filename)
   df = df[.!ismissing.(df[:JD]),:]  # clean missing values
   vcorr = df[:RVfinal] - df[:RVhel]
   cond1 = abs.(vcorr) .< dvtol
   cond2 = df[:qualflag] .> qthresh
   date_fmt = DateFormat("d/m/y")
   df[:Date] = Date.(df[:Date],date_fmt)
   cond3 = date_start .<= df[:Date] .<= date_stop
   df = df[cond1 .& cond2 .& cond3,:]
   num_spectra = size(df,1)
   if verbose
      println("# Found ", num_spectra, " HARPS-N spectra matching requirements in manifest.")
   end
   #bsrv = df[:RVhel] .- df[:RVbary]
   return df
end

"""Get data filenames from a manifest DataFrame.  Optioanlly specify path to data."""
function get_data_filenames(df::DataFrame; path::String = get_data_path(), suffix::String = "_e2ds_A.fits", verbose=false )
   #num_files = size(df,1)
   #files = Vector{String}()
   #df[:File] = map(df->joinpath(path,,@sprintf("%4d-%02d",year(df[1]),month(df[1])), @sprintf("%4d-%02d-%02d",year(df[1]),month(df[1]),day(df[1])), df[2] * "_e2ds_A.fits"),zip(df[:,:Date],df[:,:File]))
   df[:File] = map(df->joinpath(path,df[2] * suffix), zip(df[:,:Date],df[:,:File]) )
   return df[:,:File]
   #=
    for i in 1:num_files
       obsnum = Printf.@sprintf("%04d",df[i,:ObsNum])
       subdir = joinpath(df[i,:UTDate],obsnum)
       dirname = joinpath(path,subdir)
       file_pattern = Regex("(" * df[i,:Frame][1:end-5] * ".*fits)")
       local files_here
       try
          files_here = readdir(dirname)
       catch
          if verbose
             println("# Skipping " * string(file_pattern) * " not found in " * subdir )
          end
          continue
       end
       files_matched = map(f->match(file_pattern,f),files_here)
       idx_keep = findall(x->x!=nothing,files_matched)
       append!(files,map(m->joinpath(dirname,m.match),files_matched[idx_keep]))
    end
    =#
end

"""Get data filenames and corresponding barycentric solar radial velocities from searching all manifest entries.  Optioanlly specify path to data, constraints on which data files to read."""
function get_data_filenames_and_bsrv(; manifest_path::String = get_manifest_path(), manifest_filename::String = "Sun_harpsn_index.csv",
      data_path::String = get_data_path(), verbose=false,
      dvtol=0.0001, qthresh=0.99, date_start=Date(2015,01,01), date_stop::Date = date_start)
   df = read(joinpath(manifest_path,manifest_filename), dvtol=dvtol, qthresh=qthresh, date_start=date_start, date_stop=date_stop)
   bsrv = (df[:RVhel] .- df[:RVbary])  ./ speed_of_light_mps
   filenames = get_data_filenames(df, path=data_path, verbose=verbose)
   found_on_disk = isfile.(filenames)
   return (filenames=filenames[found_on_disk], z_bsrv=bsrv[found_on_disk])
end

"""Get data filenames from searching all manifest entries.  Optioanlly specify path to data, constraints on which data files to read."""
function get_data_filenames(; manifest_path::String = get_manifest_path(), manifest_filename::String = "Sun_harpsn_index.csv",
      data_path::String = get_data_path(), verbose=false,
      dvtol=0.0001, qthresh=0.99, date_start=Date(2015,01,01), date_stop::Date = date_start)
   output = get_data_filenames_and_bsrv(manifest_path=manifest_path, manifest_filename=manifest_filename,data_path=data_path,verbose=verbose,
               dvtol=dvtol, qthresh=qthresh, date_start=date_start, date_stop=date_stop)
   return output.filenames
end

end

module HPFManifests

import ..ManifestFormatEntry
 import ..hpf_manifest_format
import ..base_path, ..base_path_default, ..manifest_subdir_default, ..data_subdir_default

export read
export get_data_filenames
export get_base_path, set_base_path!, get_manifest_path, set_manifest_path!, get_data_path, set_data_path!
#export get_format, set_format!

using DataFrames
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

#manifest_format = hpf_manifest_format
#get_format() = manifest_format
#set_format!(format::Vector{ManifestFormatEntry}) = manifest_formation = format

# Read in manifest files
function read(filename::String; filter_object::Regex = r"", filter_obstype::Regex = r"", verbose=false, manifest_format = hpf_manifest_format )
    df_types = map(t->Union{t,Missing},getproperty.(manifest_format,:type) )
    manifest = DataFrame(df_types,getproperty.(manifest_format,:name),0)

    manifest_lines = readlines(filename)
    for line in manifest_lines[2:end]
        line_parsed = Dict()
        for entry in manifest_format
            name = entry.name
            type = entry.type
            if entry.type == String
                line_parsed[entry.name] = strip(line[entry.cols])
            elseif entry.type <: Number
                val = tryparse(entry.type,line[entry.cols])
                if val == nothing  
                   val = missing
                end
                line_parsed[entry.name] = val
            else
                @error("What am I supposed to do with an entry of type" * string(entry.type) * "?")
            end
        end # for entry
        if occursin(filter_object,line_parsed[:Object]) && occursin(filter_obstype,line_parsed[:ObsType])
            push!(manifest,line_parsed)
        end
    end # for line
   return manifest
end

function read(; path::String = get_manifest_path(), filter_files::Regex = r"", filter_object::Regex = r"", filter_obstype::Regex = r"", manifest_format = hpf_manifest_format, verbose=false )
    df_types = map(t->Union{t,Missing},getproperty.(manifest_format,:type) )
    manifest = DataFrame(df_types,getproperty.(manifest_format,:name),0)
    manifest_files = readdir(path)
    for manifest_filename in manifest_files
        if ! occursin(filter_files,manifest_filename )
            continue
        end
        manifest_this = read(joinpath(path,manifest_filename), filter_object=filter_object, filter_obstype=filter_obstype, verbose=verbose)
        append!(manifest,manifest_this)
    end
    return manifest
end

"""Get data filenames from a manifest DataFrame.  Optioanlly specify path to data."""
function get_data_filenames(df::DataFrame; path::String = get_data_path(), verbose=false )
    num_files = size(df,1)
    files = Vector{String}()
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
    return files
end

"""Get data filenames from searching all manifest entries.  Optioanlly specify path to data, filter for manifest filenames, and/or filter for object names."""
function get_data_filenames(; manifest_path::String = get_manifest_path(), data_path::String = get_data_path(), filter_files::Regex = r"", filter_object::Regex = r"", filter_obstype::Regex = r"", manifest_format = hpf_manifest_format, verbose = false  )
     df = read(path=manifest_path, filter_files=filter_files, filter_object=filter_object, filter_obstype=filter_obstype, manifest_format=manifest_format, verbose=verbose)
     get_data_filenames(df, path=data_path, verbose=verbose)
end

end


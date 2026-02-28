using CSV
using DataFrames
using GZip

"""
    read_vcf(vcf_file::String) -> Tuple{DataFrame, Vector{String}}

Read `vcf_file` and return content.

`vcf_file` may be gzipped.

## Returns
* `body::DataFrame` : The body of `vcf_file`.
* `header::Vector{String}` : The header of `vcf_file`, where each element is a header row.
"""
function read_vcf(filename)
    # Structure of vcf files:
    # * A couple of lines beginning with `##` (header)
    # * Single line beginning with `#CHROM` (tab-separated column names)
    # * All further lines are content (tab-separated)
    !(endswith(filename, ".vcf.gz") || endswith(filename, ".vcf")) && error("invalid file format")
    header = String[]
    columns = nothing
    content = []
    f = ifelse(endswith(filename, ".gz"), GZip.open(filename), open(filename))
    for line in readlines(f)
        if startswith(line, "##")
            push!(header, line)
        elseif startswith(line, "#CHROM")
            columns = split(strip(line), "\t")
        else
            push!(content, split(strip(line), "\t"))
        end
    end
    close(f)
    content = hcat(content...)
    return DataFrame([columns[j] => content[j, :] for j in eachindex(columns)]), header
end

#--------------------------------------------------------------------------------------------------

"""
    readcounts(vcf_file::String; <keyword arguments>) -> DataFrame

Import `vcf_file` and extract read counts.
"""
function readcounts(
    vcf_file;
    caller,
    tumorID = nothing,
    minvaf = 0.01,
    mindepth = 30,
    filter_value = ["PASS", "."],
    ignore_XY = false,
    remove_indels = false,
    max_autosome_number::Integer = 22,
    refgenome = "hg19"
    )
    
    caller = lowercase(caller)
    caller ∉ ["strelka", "mutect", "dkfz", "sentieon"] && error("invalid caller")
    
    vcf, _ = read_vcf(vcf_file)
    for column in ["#CHROM", "POS", "REF", "ALT", "INFO", "FORMAT"]
        column ∉ names(vcf) && error("$(column) not found in vcf file")
    end
    
    if isnothing(tumorID)
        tumorID = _tumorID(vcf, caller)
        @info "inferred tumorID = $(tumorID)"
        flush(stderr)
    end
    tumorID ∉ names(vcf) && error("tumorID not found in vcf file")
    
    println("Extracting read counts for tumorID = $(tumorID)")
    println("Total number of variants       : ", nrow(vcf))
    nd = ndigits(nrow(vcf)) # used for padding info printing
    
    if filter_value !== nothing
        if !hasproperty(vcf, :FILTER)
            @warn "no FILTER column detected -> filter step is skipped"
            flush(stderr)
        else
            if isa(filter_value, AbstractString)
                filter_value = [filter_value]
            end
            vcf = filter(:FILTER => f -> f ∈ filter_value, vcf)
            nrow(vcf) == 0 && error("no variants passed filtering")
            println("Variants passing filter        : $(lpad(nrow(vcf), nd))")
        end
    end
    
    # remove multi-allelic variants
    vcf = filter(:ALT => a -> !occursin(",", a), vcf)
    nrow(vcf) == 0 && error("no bi-allelic variants found")
    println("Bi-allelic variants            : $(lpad(nrow(vcf), nd))")
    
    # remove indels
    if remove_indels
        vcf = filter([:REF, :ALT] => (r, a) -> length(r) == length(a) == 1, vcf)
        nrow(vcf) == 0 && error("no single nucelotide variants found")
        println("Single nucleotide variants     : $(lpad(nrow(vcf), nd))")
    end
    
    df = vcf[!, ["#CHROM", "POS", "REF", "ALT"]]
    rename!(df, [:chrom, :pos, :ref, :alt])
    # add read counts and vaf to df
    df = hcat(df, _readcount_df(vcf, caller, tumorID))
    
    # select only valid chromosomes
    valid_chroms = push!(string.(1:max_autosome_number), "X", "Y")
    append!(valid_chroms, ["chr" * c for c in valid_chroms]) 
    df = filter(:chrom => c -> c ∈ valid_chroms, df)
    println("Primary contig variants        : $(lpad(nrow(df), nd))")

    if ignore_XY
        df = filter(:chrom => c -> c ∉ ["X", "Y"], df)
        println("Autosomal variants             : $(lpad(nrow(df), nd))")
    end
    
    df = filter(:vaf => v -> (v ≥ minvaf), df)
    println("Variants with vaf ≥ minvaf     : $(lpad(nrow(df), nd))")
    
    df = filter(:depth => d -> d ≥ mindepth, df)
    println("Variants with depth ≥ mindepth : $(lpad(nrow(df), nd))")

    if eltype(df.pos) <: AbstractString
        df.pos .= parse.(Int, df.pos)
    end

    df[!, :chrom] .= map(s -> startswith(s, "chr") ? s : "chr" * s, df[!, :chrom])
    
    # Fixed assignment for new columns in DataFrames
    df[!, :sampleID] .= tumorID
    df[!, :refgenome] .= refgenome
    return df
end

#--------------------------------------------------------------------------------------------------
function _tumorID(vcf, caller)
    stdcols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    sns = setdiff(names(vcf), stdcols)
    if caller == "dkfz"
        # tumorID is the column name right after the standard columns
        return sns[1]
    else
        # In case matched-normal is used, second entry will always be the tumor sample
        # If not, first sample is assumed to be tumor
        return ifelse(length(sns) > 1, sns[2], sns[1])
        #return sns[1]
    end
end


#--------------------------------------------------------------------------------------------------

function _readcount_df(vcf, caller, tumorID)
    ref_count = Int[]
    alt_count = Int[]
    depth = Int[]
    vaf = Float64[]
    
    if caller == "strelka"
        fmt = split(vcf[1, :FORMAT], ":")
        counts = hcat(split.(vcf[!, tumorID], ":")...)
        counts = map(x -> parse(Int, split(x, ",")[1]), counts)
        counts_df = DataFrame(Dict(zip(fmt, [counts[j,:] for j in eachindex(fmt)])))
        rename!(counts_df, :AU => :A, :CU => :C, :GU => :G, :TU => :T)
        for (j, row) in enumerate(eachrow(vcf))
            push!(ref_count, counts_df[j, row.REF])
            push!(alt_count, counts_df[j, row.ALT])
            push!(depth, counts_df[j, row.REF] + counts_df[j, row.ALT])
            push!(vaf, counts_df[j, row.ALT] / (counts_df[j, row.REF] + counts_df[j, row.ALT]))
        end
        
    elseif caller == "dkfz"
        for row in eachrow(vcf)
            counts = parse.(Int, split(match(r"DP4=([^;]*);",row.INFO).captures[1], ","))
            push!(ref_count, sum(counts[1:2]))
            push!(alt_count, sum(counts[3:4]))
            push!(depth, sum(counts))
            push!(vaf, sum(counts[3:4]) / sum(counts))
        end
        
    elseif caller == "mutect"
        for (j, row) in enumerate(eachrow(vcf))
            counts = parse.(Int, 
                split(
                    Dict(zip(split(row.FORMAT, ":"), split(vcf[j, tumorID], ":")))["AD"],
                    ",")
            )
            push!(ref_count, counts[1])
            push!(alt_count, sum(counts[2]))
            push!(depth, sum(counts))
            push!(vaf, sum(counts[2]) / sum(counts))
        end
        
    elseif caller == "sentieon"
        for (j, row) in enumerate(eachrow(vcf))
            values = Dict(zip(split(row.FORMAT, ":"), split(vcf[j, tumorID], ":")))
            dp = parse(Int, values["DP"])
            af = parse(Float64, values["AF"])
            push!(ref_count, round(Int, dp * (1 - af)))
            push!(alt_count, round(Int, dp * af))
            push!(depth, dp)
            push!(vaf, af)
        end
        
    else
        error("unknown caller")
    end
    
    return DataFrame(
        :ref_count    => ref_count,
        :alt_count    => alt_count,
        :depth        => depth,
        :vaf          => vaf
        )
end

###################################################################################################

function read_cnv(file; kwargs...)
    kwargs = merge(Dict(:missingstring => ["NA", "missing", "", "sub"]), kwargs)
    CSV.read(file, DataFrame; kwargs...)
end

#--------------------------------------------------------------------------------------------------

function copynumberinfo(
    fname;
    chr_col  ::Union{String, Integer, Nothing} = nothing,
    start_col::Union{String, Integer, Nothing} = nothing,
    end_col  ::Union{String, Integer, Nothing} = nothing,
    A_col    ::Union{String, Integer, Nothing} = nothing,
    B_col    ::Union{String, Integer, Nothing} = nothing,
    tcn_col  ::Union{String, Integer, Nothing} = nothing,
    merge_tolerance = 1e5,
    ignore_XY = true,
    max_cn = 4,
    tumorID = nothing,
    max_dev = 0.2,
    csv_args...)

    cn_df = read_cnv(fname; csv_args...)
    cn_df = _adjust_dfcols(cn_df, chr_col, start_col, end_col, A_col, B_col, tcn_col)
    (nrow(cn_df) == 0) && error("no copy number information in file $(fname)")
    println("Extracted $(nrow(cn_df)) segments on $(length(unique(cn_df.chrom))) chromosomes")
    nd = ndigits(nrow(cn_df))

    cn_df = dropmissing(cn_df)
    println("Segments with copy number information : $(lpad(nrow(cn_df), nd))")
    
    cn_df.tcn .= Float64.(cn_df.tcn)
    cn_df = filter(:tcn => tcn -> (tcn % 1 ≤ max_dev) || tcn % 1 ≥ 1 - max_dev, cn_df)
    cn_df.tcn .= round.(Int, cn_df.tcn)
    cn_df = filter(:A => A -> (A % 1 ≤ max_dev/2) || (A % 1 ≥ 1 - max_dev/2), cn_df)
    cn_df = filter(:B => B -> (B % 1 ≤ max_dev/2) || (B % 1 ≥ 1 - max_dev/2), cn_df)
    cn_df.A .= round.(Int, cn_df.A)
    cn_df.B .= round.(Int, cn_df.B)
    println("Segments with clonal copy number      : $(lpad(nrow(cn_df), nd))")

    cn_df = filter(:tcn => tcn -> 0 < tcn ≤ max_cn, cn_df)
    println("Segments with copy numbers $(min_cn) ≤ cn ≤ $(max_cn) : $(lpad(nrow(cn_df), nd))")
    (nrow(cn_df) == 0) && error("No segments with copy numbers $(min_cn) ≤ cn ≤ $(max_cn)")

    cn_df.chrom .= map(s -> replace(s, "chr" => ""), cn_df.chrom)

    if (ignore_XY)
        cn_df = filter(:chrom => c -> c ∉ ["X", "Y"], cn_df)
        println("Autosomal segments                    : $(lpad(nrow(cn_df), nd))")
    end

    cn_df = _merge_adjacent_segs(cn_df, merge_tolerance)
    cn_df = filter([:start, :end] => (s,e) -> s ≤ e, cn_df)
    println("Remaining segments after merging      : $(lpad(nrow(cn_df), nd))")
    
    if sum(cn_df.end .- cn_df.start) < 3e8
        @warn "Less than 10% of the genome with valid copy number information."
        flush(stderr)
    end
    
    cn_df[!, :sample] .= ifelse(isnothing(tumorID), "sample", tumorID)

    return cn_df
end

#--------------------------------------------------------------------------------------------------

function _adjust_dfcols(cn_df, chr_col, start_col, end_col, A_col, B_col, tcn_col)
    colnames = lowercase.(names(cn_df))
    errmsg(name, p="") = "failed to infer $(name) column$(p), please provide column information"
    
    function infer_colname(col, func)
        if isa(col, Integer)
            return names(cn_df)[col]
        elseif isa(col, AbstractString)
            col ∉ names(cn_df) && error("column $(col) not found")
            return col
        elseif isnothing(col)
            cands = findall(func, colnames)
            if length(cands) == 1
                return names(cn_df)[first(cands)]
            else
                return nothing
            end
        else
            error("invalid type for parameter col")
        end
    end
    infer_colname(func) = infer_colname(nothing, func)
    
    chr_col = infer_colname(chr_col, contains("chr"))
    for id in ["chr", "chrom", "chromosome"]
        chr_col !== nothing && break
        chr_col = infer_colname(chr_col, isequal(id))
    end
    isnothing(chr_col) && error(errmsg("chromosome"))
    
    for id in ["start", "chromstart"]
        start_col = infer_colname(start_col, startswith(id))
        start_col !== nothing && break
    end
    isnothing(start_col) && error(errmsg("start"))
    
    for id in ["end", "chromend"]
        end_col = infer_colname(end_col, startswith(id))
        end_col !== nothing && break
    end
    isnothing(end_col) && error(errmsg("end"))
    
    for id in ["tcn", "cnt", "copynumber", "copy number"]
        tcn_col = infer_colname(tcn_col, isequal(id))
        tcn_col !== nothing && break
    end
    
    A_col = infer_colname(A_col, contains("major"))
    if isnothing(A_col)
        A_col = infer_colname(isequal("a"))
    end
    
    B_col = infer_colname(B_col, contains("minor"))
    if isnothing(B_col)
        B_col = infer_colname(isequal("b"))
    end
    
    df = DataFrame(
        :chrom => cn_df[!, chr_col],
        :start => cn_df[!, start_col],
        :end => cn_df[!, end_col],
    )
    
    if isnothing(tcn_col) && (isnothing(A_col) || isnothing(B_col))
        error(errmsg("tcn and A/B", "s"))
    end
    
    if isnothing(tcn_col)
        df[!, :tcn] .= cn_df[!, A_col] .+ cn_df[!, B_col]
        df[!, :A] .= cn_df[!, A_col]
        df[!, :B] .= cn_df[!, B_col]
    else
        if (A_col !== nothing) && (B_col !== nothing)
            df[!, :tcn] .= cn_df[!, tcn_col]
            df[!, :A] .= cn_df[!, A_col]
            df[!, :B] .= cn_df[!, B_col]
        end
        if (A_col === nothing) && (B_col !== nothing)
            df[!, :tcn] .= cn_df[!, tcn_col]
            df[!, :A] .= cn_df[!, tcn_col] .- cn_df[!, B_col]
            df[!, :B] .= cn_df[!, B_col]
        end
        if (A_col !== nothing) && (B_col === nothing)
            df[!, :tcn] .= cn_df[!, tcn_col]
            df[!, :A] .= cn_df[!, A_col]
            df[!, :B] .= cn_df[!, tcn_col] .- cn_df[!, A_col]
        end
        if (A_col === nothing) && (B_col === nothing)
            df[!, :tcn] .= cn_df[!, tcn_col]
            df.A .= ceil(Int, df.tcn ./ 2)
            df.B .= floor(Int, df.tcn ./ 2)
        end
    end

    k = ismissing.(df.tcn)
    df[k, :tcn] .= df[k, :A] .+ df[k, :B]
    k = ismissing.(df.A)
    df[k, :A] .= df[k, :tcn] .- df[k, :B]
    k = ismissing.(df.B)
    df[k, :B] .= df[k, :tcn] .- df[k, :A]
    
    return df
end

#--------------------------------------------------------------------------------------------------

function _merge_adjacent_segs(cn_df, merge_tolerance)
    distinct_seqs(df) = prepend!(
        1 .+ findall(@. !((df.start[2:end] - df.end[1:end-1] < merge_tolerance) &&
                          (df.A[1:end-1] == df.A[2:end]) && 
                          (df.B[1:end-1] == df.B[2:end])))
        , 1)
    
    dfs = DataFrame[]
    for subdf in groupby(cn_df, :chrom)
        keep = distinct_seqs(subdf)
        df = subdf[keep, :]
        df.end .= subdf.end[push!(keep[2:end] .- 1, nrow(subdf))]
        push!(dfs, df)
    end
    return vcat(dfs...)
end

###################################################################################################

function add_copynumberinfo(vcf::DataFrame, cn::DataFrame)
    sdfs = DataFrame[]
    for chrom in unique(vcf.chrom)
        svcf = filter(:chrom => c -> c == chrom, vcf)
        scn  = filter(:chrom => c -> c == chrom, cn)
        
        svcf.tcn = fill!(zeros(Union{Int, Missing}, nrow(svcf)), missing)
        svcf.A   = fill!(zeros(Union{Int, Missing}, nrow(svcf)), missing)
        svcf.B   = fill!(zeros(Union{Int, Missing}, nrow(svcf)), missing)
        isempty(scn) && continue
        for (k, row) in enumerate(eachrow(svcf))
            for row2 in eachrow(scn)
                if row2.start ≤ row.pos ≤ row2.end
                    svcf[k, [:tcn, :A, :B]] .= row2.tcn, row2.A, row2.B
                elseif row2.end > row.pos
                    continue
                end
            end
        end
        push!(sdfs, svcf)
    end
    return dropmissing(vcat(sdfs...))
end

###################################################################################################
#                    --- BATCH PROCESSING SCRIPT ---


input_dir = raw"C:\Users\damme\Documents\uni\4_Semester\Hoefer_lab_rotation\Wilms_Tumor\vcf_files_dkfz_pipeline"
output_dir = raw"C:\Users\damme\Documents\uni\4_Semester\Hoefer_lab_rotation\Wilms_Tumor\Julia_SVF\mutect_vcf"

mkpath(output_dir)

input_files = filter(f -> endswith(f, ".gz"), readdir(input_dir, join=true))

for filepath in input_files
    filename = basename(filepath)
    
    clean_name = replace(filename, ".vcf.gz" => "")
    clean_name = replace(clean_name, "_hg19_intersect" => "") 
    clean_name = replace(clean_name, "-" => "_") 
    
    output_filename = "$(clean_name)_julia.csv"
    output_path = joinpath(output_dir, output_filename)
    
    println("\n========================================================")
    println("--- Processing: $filename ---")
    println("========================================================")
    

    df = readcounts(
        filepath, 
        caller="mutect", 
        minvaf=0.0, 
        mindepth=0, 
        ignore_XY=false
    )
    
    CSV.write(output_path, df)
    println("Saved successfully to: $output_path")
end

println("\nAll files have been successfully processed!")
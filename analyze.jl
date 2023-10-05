using Pkg
Pkg.activate("AnnotationEnv")
Pkg.instantiate()

using RNASeqTools, CSV, DataFrames, GenomicFeatures, Combinatorics

function run_tss()
    coverage_files = CoverageFiles(joinpath(@__DIR__, "coverage/"))
    coverages = [Coverage(f1, f2) for (f1, f2) in coverage_files]
    samplenames = ["tex_lcd_kf_01", "tex_hcd_kf_10", "tex_lcd_kf_02", "tex_hcd_kf_09"]

    features = Features()

    for (sname, coverage) in zip(samplenames, coverages)

        features *= maxdiffpositions(coverage; type="TSS", samplename=sname, min_diff=10)

    end

    write(joinpath(@__DIR__, "tss.gff"), features)

    df = DataFrame([(ref=refname(f), position=leftposition(f), strand=strand(f), diff=f["diff"], height=f["height"], from=f["from"]) for f in features if type(f) == "TSS"])
    CSV.write(joinpath(@__DIR__, "tss.csv"), df)
end

function make_annotation()
    df_srna = DataFrame(CSV.File(joinpath(@__DIR__, "2023-05-10_allsRNAsMGH.csv"); stringtype=String))

    features = Features(joinpath(@__DIR__, "MGH78578_plasmids_clean.gff3"); name_keys=["locus_tag", "ID"])
    for row in eachrow(df_srna)
        f = Interval(row.refname, row.left, row.right, row.strand == "positive" ? '+' : '-',
            Annotation("ncRNA", row.name, Dict("Name"=>row.name, "alternative_name"=>ismissing(row.alternative_name) ? row.name : row.alternative_name)))
        push!(features, f)
    end

    df_tss = DataFrame(CSV.File(joinpath(@__DIR__, "TSS-MGH-finished-v02.csv"); stringtype=String))

    max_scan_distance = 250
    refs = unique(df_tss.ref)
    tss_forward = Dict(r=>falses(maximum(df_tss.position[(df_tss.ref .== r) .& (df_tss.strand .== "+")])+max_scan_distance) for r in refs)
    tss_reverse = Dict(r=>falses(maximum(df_tss.position[(df_tss.ref .== r) .& (df_tss.strand .== "-")])+max_scan_distance) for r in refs)
    df_tss.decide[ismissing.(df_tss.decide)] .= "-"
    for r in refs
        forward_index = (df_tss.ref .== r) .& (df_tss.strand .== "+")
        tss_forward[r][df_tss.position[forward_index][df_tss.decide[forward_index] .== "y"]] .= true
        reverse_index = (df_tss.ref .== r) .& (df_tss.strand .== "-")
        tss_reverse[r][df_tss.position[reverse_index][df_tss.decide[reverse_index] .== "y"]] .= true
    end
    df_tss_trans = Dict((row.ref, row.position, String(row.strand))=>i for (i, row) in enumerate(eachrow(df_tss)) if row.decide == "y")

    cds_features = Features(joinpath(@__DIR__, "MGH78578_plasmids_clean.gff3"), "CDS"; name_keys=["locus_tag", "ID"])

    df_tss[:, :annotation] = repeat([""], nrow(df_tss))
    df_tss[:, :gene] = repeat([""], nrow(df_tss))
    categories = Dict(r=>falses(5, max(length(tss_forward[r]), length(tss_reverse[r]))) for r in refs)
    for cds in cds_features
        f = if strand(cds) == STRAND_POS
            tssi = tss_forward[refname(cds)]
            tss = leftposition(cds)-max_scan_distance < 1 ? Int64[] : findall(view(tssi, leftposition(cds)-max_scan_distance:leftposition(cds)-1))
            if isempty(tss)
                Interval(refname(cds), leftposition(cds)-150, leftposition(cds)-1, '+', Annotation("5UTR", name(cds), merge(params(cds), Dict("source"=>"manual"))))
            else
                offset = leftposition(cds)-max_scan_distance-1
                for i in 1:length(tss)
                    tss_position = offset+tss[i]
                    df_index = df_tss_trans[(refname(cds), tss_position, "+")]
                    df_tss.gene[df_index] = name(cds)

                    if i == length(tss)
                        categories[refname(cds)][1, tss_position] = true
                        df_tss.annotation[df_index] *= "proximal"
                    else
                        categories[refname(cds)][2, tss_position] = true
                        df_tss.annotation[df_index] *= "distal"
                    end

                    for f in eachoverlap(Interval(refname(cds), offset+first(tss), offset+first(tss), '+'), cds_features)
                        categories[refname(cds)][3, tss_position] = true
                        df_tss.annotation[df_index] *= df_tss.annotation[df_index] == "" ? "internal" : "-internal"
                    end

                    for f in eachoverlap(Interval(refname(cds), offset+first(tss), offset+first(tss), '-'), cds_features)
                        categories[refname(cds)][4, tss_position] = true
                        df_tss.annotation[df_index] *= df_tss.annotation[df_index] == "" ? "antisense" : "-antisense"
                    end
                end
                Interval(refname(cds), offset+first(tss), leftposition(cds)-1, '+', Annotation("5UTR", name(cds), merge(params(cds), Dict("source"=>"TSS"))))
            end
        else
            tssi = tss_reverse[refname(cds)]
            tss = rightposition(cds)+max_scan_distance > length(tssi) ? Int64[] : findall(view(tssi, rightposition(cds)+1:rightposition(cds)+max_scan_distance))
            if isempty(tss)
                Interval(refname(cds), rightposition(cds)+1, rightposition(cds)+150, '-', Annotation("5UTR", name(cds), merge(params(cds), Dict("source"=>"manual"))))
            else
                offset = rightposition(cds)
                for i in 1:length(tss)
                    tss_position = offset+tss[i]
                    df_index = df_tss_trans[(refname(cds), tss_position, "-")]
                    df_tss.gene[df_index] = name(cds)

                    if i == 1
                        categories[refname(cds)][1, tss_position] = true
                        df_tss.annotation[df_index] *= "proximal"
                    else
                        categories[refname(cds)][2, tss_position] = true
                        df_tss.annotation[df_index] *= "distal"
                    end

                    for f in eachoverlap(Interval(refname(cds), offset+last(tss), offset+last(tss), '-'), cds_features)
                        categories[refname(cds)][3, tss_position] = true
                        df_tss.annotation[df_index] *= df_tss.annotation[df_index] == "" ? "internal" : "-internal"
                    end

                    for f in eachoverlap(Interval(refname(cds), offset+last(tss), offset+last(tss), '+'), cds_features)
                        categories[refname(cds)][4, tss_position] = true
                        df_tss.annotation[df_index] *= df_tss.annotation[df_index] == "" ? "antisense" : "-antisense"
                    end
                end
                Interval(refname(cds), offset+1, offset+last(tss), '-', Annotation("5UTR", name(cds), merge(params(cds), Dict("source"=>"TSS"))))
            end
        end
        push!(features, f)
    end



    for r in refs
        for tss in findall(tss_forward[r])
            if !any(categories[r][:, tss])
                df_index = df_tss_trans[(r, tss, "+")]
                df_tss.annotation[df_index] *= "orphan"
                categories[r][5, tss] = true
                for f in eachoverlap(Interval(r, tss, tss, '+'), cds_features)
                    categories[r][3, tss] = true
                    df_tss.annotation[df_index] *= df_tss.annotation[df_index] == "" ? "internal" : "-internal"
                end
                for f in eachoverlap(Interval(r, tss, tss, '-'), cds_features)
                    categories[r][4, tss] = true
                    df_tss.annotation[df_index] *= df_tss.annotation[df_index] == "" ? "antisense" : "-antisense"
                end
            end
        end
        for tss in findall(tss_reverse[r])
            if !any(categories[r][:, tss])
                df_index = df_tss_trans[(r, tss, "-")]
                df_tss.annotation[df_index] *= "orphan"
                categories[r][5, tss] = true
                for f in eachoverlap(Interval(r, tss, tss, '-'), cds_features)
                    categories[r][3, tss] = true
                    df_tss.annotation[df_index] *= df_tss.annotation[df_index] == "" ? "internal" : "-internal"
                end
                for f in eachoverlap(Interval(r, tss, tss, '+'), cds_features)
                    categories[r][4, tss] = true
                    df_tss.annotation[df_index] *= df_tss.annotation[df_index] == "" ? "antisense" : "-antisense"
                end
            end
        end
    end

    add3utrs!(cds_features; utr_length=150, min_utr_length=1)
    for f in cds_features
        type(f) == "3UTR" && push!(features, f)
    end

    write(joinpath(@__DIR__, "MGH78578_plasmids_clean_full.gff3"), features)
    CSV.write(joinpath(@__DIR__, "TSS-MGH-finished-v02_annotated.csv"), df_tss)

    test_features = Features(joinpath(@__DIR__, "MGH78578_plasmids_clean_full.gff3"))
    overlap = zeros(5, 5)
    for r in refs
        for (i, row) in enumerate(eachrow(categories[r]))
            overlap[i, i] += sum(row)
        end
        for ((i1, c1), (i2, c2)) in combinations(collect(enumerate(eachrow(categories[r]))), 2)
            overlap[i1, i2] += sum(c1 .& c2)
            overlap[i2, i1] += sum(c1 .& c2)
        end
    end

    out_frame = DataFrame(
        category = ["proximal", "distal", "orphan"],
        total = [overlap[1, 1], overlap[2, 2], overlap[5, 5]],
        internal = [overlap[1, 3], overlap[2, 3], overlap[5, 3]],
        antisense = [overlap[1, 4], overlap[2, 4], overlap[5, 4]]
    )
    CSV.write(joinpath(@__DIR__, "categories.csv"), out_frame)
end

run_tss()
make_annotation()
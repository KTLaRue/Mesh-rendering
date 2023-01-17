# Usage:
# julia --project meshgen.jl -g <cube|sphere|cylinder> [-n <divisionsU>] [-m <divisionsV>] <outfile.obj>
# julia --project meshgen.jl -i <infile.obj> <outfile.obj>

using ArgParse

push!(LOAD_PATH, pwd())
using WWUMeshes

""" parse_cmdline
Parse the command line args for the meshgen tool. """
function parse_cmdline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "-g"
            help="geometry (cube, cylider, or sphere)"
        "-n"
            help="divisionsU"
            arg_type=Int
            default=32
        "-m"
            help="divisionsV"
            arg_type=Int
            default=16
        "-i"
            help="input filename"
        "outfile"
            help="output filename"
    end

    return parse_args(s)
end

""" main
Execute the appropriate behavior from WWUMeshes based on the provided arguments. """
function main()
    args = parse_cmdline()

    if args["i"] != nothing # normal estimation usage
        est_normals(args["outfile"], args["i"])
    else # mesh generation usage
        geom = args["g"]
        divisionsU = args["n"]
        divisionsV = args["m"]
        gen_mesh(args["outfile"], geom, divisionsU, divisionsV)
    end
end

main()


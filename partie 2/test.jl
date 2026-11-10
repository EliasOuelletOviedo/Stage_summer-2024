using BenchmarkTools
using PyCall
using FFTW
using Statistics
using LinearAlgebra
using Flux
using BSON
using Plots

sp_o = pyimport("scipy.optimize")

global sdt = pyimport("sdtfile") # use python package sdtfile

include("vec_to_lifetime.jl")
# include("GUI_complet.jl")

function reshape_to_vec(file, histogram_resolution)
    return sum(reshape(file, (histogram_resolution, :)), dims = 2)[:,1]
end

function open_SDT_file(filepath)
    open(filepath,"r") do io
        seek(io, 14)
        header = read(io, 12)

        seek(io, header[12]*16^2 + header[11] + 82)
        infos = read(io, 2)
        histogram_resolution = infos[2]*16^2 + infos[1]

        if header[8]*16^2+header[7] == 512 || header[8]*16^2+header[7] == 128
            seek(io, header[2]*16^2+header[1]+22)

            vector = read(io)
            new_vector = UInt16[]

            for i in 1:2:length(vector)
                push!(new_vector, vector[i+1]*16^2 + vector[i])
            end

            return new_vector, histogram_resolution
        else
            seek(io, header[2]*16^2+header[1]+2)

            file_info = read(io, 20)

            shift_1 = file_info[8]*16^6 + file_info[7]*16^4 + file_info[6]*16^2 + file_info[5] - file_info[2]*16^2 - file_info[1]
            shift_2 = file_info[20]*16^6 + file_info[19]*16^4 + file_info[18]*16^2 + file_info[17]

            buffer = IOBuffer(read(io, shift_1))
            zip_file = ZipFile.Reader(buffer).files[1]

            file_data = Vector{UInt8}(undef, shift_2)
            read!(zip_file, file_data)

            vector = reshape_to_vec(reinterpret(UInt16, file_data), histogram_resolution)
    
            return vector, histogram_resolution
        end
    end
end

function open_SDT_file_2(filepath)
    open(filepath, "r") do io
        seek(io, 14)
        header = read!(io, Vector{UInt8}(undef, 12))

        seek(io, header[12]*0x100 + header[11] + 82)
        infos = read!(io, Vector{UInt8}(undef, 2))
        histogram_resolution = infos[2]*0x100 + infos[1]

        seek(io, header[2]*0x100 + header[1] + 22)
        vector = read(io)
        n = div(length(vector), 2)
        new_vector = reinterpret(UInt16, vector[1:2*n])
        return new_vector, histogram_resolution
    end
end

function open_SDT_file_3(filepath)
    file = sdt.SdtFile(filepath)
    data = file.data[1]
    return data, length(data)
end

function get_irf(; channel=1)
    """
    Loads an .sdt IRF file. If irf_filepath.txt exists, load the IRF from this filepath, otherwise, open a file dialog create irf_filepath.txt
    """

    # println("Getting IRF file...")

    if isfile("irf_filepath.txt")
        filepath = open(f->read(f, String), "irf_filepath.txt")
    else
        filepath = pick_file()

        open("irf_filepath.txt", "w") do io
            write(io, filepath)
        end
    end
    sdt = pyimport("sdtfile") # use python package sdtfile
    file = sdt.SdtFile(filepath)
    data = convert.(Int64, file.data[channel])
    println(open_SDT_file_2(filepath))
    median_irf = round(median(data)) # get median
    data = data.-round(Int, 1.25*median_irf) # remove noise from the IRF
    data[data .<= 0] .= 0 # set negative values to 0
    times = convert.(Float64, vec(file.times[1]).*1e9) # convert from seconds to nanoseconds
    data = [times, vec(data)]


    # println("Done")
    # println()

    return data
end

# function get_irf_2()
#     if isfile("irf_filepath.txt")
#         filepath = open(f->read(f, String), "irf_filepath.txt")
#     else
#         filepath = pick_file()
    
#         open("irf_filepath.txt", "w") do io
#             write(io, filepath)
#         end
#     end
    
#     data = open_SDT_file_2(filepath)[1]
# end

filepath = "/Users/eliasouellet-oviedo/Documents/Stage2/Codes/test/puff_full_c01000.sdt"

# sdt = pyimport("sdtfile") # use python package sdtfile
# file = sdt.SdtFile(filepath)
# println(file.data[1])

@benchmark open_SDT_file_3(filepath)
# @benchmark open_SDT_file_2(filepath)

# file = get_irf()
# file = open_SDT_file_3(filepath)
# println(file)



# if isfile("irf_filepath.txt")
#     filepath = open(f->read(f, String), "irf_filepath.txt")
# else
#     filepath = pick_file()

#     open("irf_filepath.txt", "w") do io
#         write(io, filepath)
#     end
# end

# test = open_SDT_file_2(filepath)[1]

# plot(file[1], test)
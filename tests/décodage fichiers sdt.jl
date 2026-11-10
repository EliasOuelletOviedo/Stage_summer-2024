
using PyCall
using StructTypes
using Base64
using DelimitedFiles
using ZipFile
using Base.Filesystem
# using JuliaDB
# using Plots
# using GR

# gr()

struct FileHeader
    revision::Int16
    info_offs::Int32
    info_length::Int16
    setup_offs::Int32
    setup_length::UInt16
    data_block_offs::Int32
    no_of_data_blocks::Int16
    data_block_length::UInt32
    meas_desc_block_offs::Int32
    no_of_meas_desc_blocks::Int16
    meas_desc_block_length::Int16
    header_valid::UInt16
    reserved1::UInt32
    reserved2::UInt16
    chksum::UInt16
end

struct BlockHeader
    data_offs_ext::UInt8
    next_block_offs_ext::UInt8
    data_offs::UInt32
    next_block_offs::UInt32
    block_type::UInt16
    meas_desc_block_no::Int16
    lblock_no::UInt32
    block_length::UInt32
end

filepath = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-07-17/Single_256/Single_256_c01.sdt"
# filepath = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-07-17/Single_256/Single_256_c01.sdt"
# filepath = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-07-17/Image_64/Image_64_c10.sdt"
# filepath = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-07-17/Single_64/Single_64_c01.sdt"
# filepath2 = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-06-10/single/single_c001.sdt"
# filepath2 = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/Fichiers_sdt/1.sdt"
# folderpath = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/Fichiers_sdt"
# folderpath2 = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-06-10/single"
# folderpath2 = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/Fichiers_sdt"

# filepath = "/Users/eliasouellet-oviedo/Documents/Stage2/Partie 2/Test_2_7/test_c00001.sdt"
filepath = "/Users/eliasouellet-oviedo/Documents/Stage2/Partie 2/Puffs (5ms et 20ms à 40 Hz)/test6_c0001.sdt"

open(filepath, "r") do io
    # read(io)
    seek(io, 14)
    header = read!(io, Vector{UInt8}(undef, 12))
    
    seek(io, header[12]*0x100 + header[11] + 82)
    infos = read!(io, Vector{UInt8}(undef, 2))
    histogram_resolution = infos[2]*0x100 + infos[1]
    
    seek(io, header[12]*0x100 + header[11] + 215)
    infos_2 = read!(io, Vector{UInt8}(undef, 4))
    time = reinterpret(Float32, infos_2)

    seek(io, header[2]*0x100 + header[1] + 22)
    vector = read(io)
    
    n = div(length(vector), 2)
    new_vector = reinterpret(UInt16, vector[1:2*n])
    println(length(new_vector), "\n", histogram_resolution)



    # seek(io, 14)
    # header = read(io, 12)

    # seek(io, header[12]*16^2 + header[11])
    # time = read(io, 18)
    # println(time)

    # seek(io, header[12]*16^2 + header[11] + 82)
    # infos = read(io, 2)
    # histogram_resolution = infos[2]*16^2 + infos[1]

    # seek(io, header[2]*16^2 + header[1] + 2)
    # file_info = read(io, 20)

    # shift_1 = file_info[8]*16^6 + file_info[7]*16^4 + file_info[6]*16^2 + file_info[5] - file_info[2]*16^2 - file_info[1]
    # shift_2 = file_info[20]*16^6 + file_info[19]*16^4 + file_info[18]*16^2 + file_info[17]

    # buffer = IOBuffer(read(io, shift_1))
    # zip_file = ZipFile.Reader(buffer).files[1]

    # # Pre-allocate storage for file_data
    # file_data = Vector{UInt8}(undef, shift_2)
    # read!(zip_file, file_data)

    # println()
end

# println("\nIMAGE")

# for file in readdir(folderpath)[2:end]
#     open(joinpath(folderpath, file),"r") do io
#         # seek(io, 14)
#         header = read(io, 42)
#         # seek(io, 1)
#         # test = read(io, 2)
#         # println(header[41:42])
#         # println(test[2]*16^2+test[1])
#         println(header[42]*16^2+header[41])
#     end
# end

# println("\nSINGLE")

# for file in readdir(folderpath2)[2:end]
#     open(joinpath(folderpath2, file),"r") do io
#         # seek(io, 14)
#         header = read(io, 42)
#         # println(header[41:42])
#         println(header[42]*16^2+header[41])
#     end
# end


# open(filepath,"r") do io
#     # println(read(io, 42))

#     seek(io, 14)
#     header = read(io, 8)
#     # println(header)

#     seek(io, 24)
#     test = read(io, 2)
#     println(test)

#     seek(io, test[2]*16^2 + test[1]+82)
#     infos = read(io, 2)
#     println(infos[2]*16^2 + infos[1])
# end

# open(filepath2,"r") do io
#     seek(io, 14)
#     header = read(io, 2)
#     println(header)
# end

# open(filepath,"r") do io
#     seek(io, 14)
#     header = read(io, 2)
#     # seek(io, header[2]*16^2+header[1] + 22)

#     # println(read(io, 22))

#     # println(header)

#     # file_info = read(io, 22)

#     # println(file_info)

#     # vector = read(io)
#     # new_vector = UInt16[]

#     # for i in 1:2:length(vector)
#     #     push!(new_vector, vector[i+1] * 16^2 + vector[i])
#     # end

#     # println((new_vector))

#     # zip = IOBuffer(read(io, file_info[10]*16^6 + file_info[9]*16^4 + file_info[8]*16^2 + file_info[7] - (file_info[6]*10^6 + file_info[5]*16^4 + file_info[4]*16^2 + file_info[3])))

#     # array = reinterpret(UInt16, read(ZipFile.Reader(zip).files[1], file_info[22]*16^6 + file_info[21]*16^4 + file_info[20]*16^2 + file_info[19]))

#     # println(size(array))
# end




# file = open(filepath, "r")
# seek(file, 14)
# header = read(file, 2)
# seek(file, header[2]*16^2+header[1])

# file_info = read(file, 22)

# println(file_info[10]*16^6 + file_info[9]*16^4 + file_info[8]*16^2 + file_info[7] - (file_info[6]*10^6 + file_info[5]*16^4 + file_info[4]*16^2 + file_info[3]))

# zip = IOBuffer(read(file, file_info[10]*16^6 + file_info[9]*16^4 + file_info[8]*16^2 + file_info[7] - (file_info[6]*10^6 + file_info[5]*16^4 + file_info[4]*16^2 + file_info[3])))

# @time data = reinterpret(UInt16, read(ZipFile.Reader(zip).files[1]))
# println(size(data))

# data_array = reshape(data, (1024, 1024, 256))

# close(file)

# open(filepath,"r") do io
#     seek(io, 14)
#     header = read(io, 2)
#     seek(io, header[2]*16^2+header[1])

#     file_info = read(io, 22)

#     println(file_info[10]*16^6 + file_info[9]*16^4 + file_info[8]*16^2 + file_info[7] - (file_info[6]*10^6 + file_info[5]*16^4 + file_info[4]*16^2 + file_info[3]))

#     # @time zip = IOBuffer(read(io, file_info[10]*16^6 + file_info[9]*16^4 + file_info[8]*16^2 + file_info[7] - (file_info[6]*10^6 + file_info[5]*16^4 + file_info[4]*16^2 + file_info[3])))

#     # @time unzip = ZipFile.Reader(zip).files[1]
    
#     # @time read_unzip = read(unzip, 536_870_912)

#     # @time data = reinterpret(UInt16, read_unzip)

#     # @time data_array = permutedims(reshape(data, (256, 1024, 1024))[:,:,begin:2:end], [3, 2, 1])

#     # println(size(data))
#     # println(size(data_array))

#     # heatmap(sum(data_array, dims = 3)[:,:,1], c=cgrad([:black, :white]), yflip = true, aspect_ratio = :equal)

#     # writedlm("SDT.csv", data, ",")
# end

# open(filepath,"r") do io
#     seek(io, 14)
#     header = read(io, 2)
#     seek(io, header[2]*16^2+header[1])

#     file_info = read(io, 22)

#     println(file_info[10]*16^6 + file_info[9]*16^4 + file_info[8]*16^2 + file_info[7] - (file_info[6]*10^6 + file_info[5]*16^4 + file_info[4]*16^2 + file_info[3]))

#     @time zip = IOBuffer(read(io, file_info[10]*16^6 + file_info[9]*16^4 + file_info[8]*16^2 + file_info[7] - (file_info[6]*10^6 + file_info[5]*16^4 + file_info[4]*16^2 + file_info[3])))

#     @time unzip = ZipFile.Reader(zip).files[1]
    
#     @time read_unzip = read(unzip, 536_870_912)

#     @time data = reinterpret(UInt16, read_unzip)

#     println(size(data))
# end




# a do zf
#     print(ZipFile.Reader(test))
# end

# seek(file, header.data_block_offs)
# bh_buffer = read(file, sizeof(BlockHeader))
# bh = reinterpret(FileHeader, bh_buffer)[1]


# println(header_buffer[16]*256+header_buffer[15])
# println(header_buffer[15])

# println(header_buffer)
# println(next[10]*16_777_216 + next[9]*65_536 + next[8]*256 + next[7] - (next[6]*16_777_216 + next[5]*65_536 + next[4]*256 + next[3]))
# println(next)
# print(peutetre)


# 

# file = open(filepath, "r")
# # println(size(read(file)))

# x = read(file)
# # println(size(x))
# println(x[1:100])

# println(x[1] + x[2] + x[3] + x[4])
# println(string.(x[1:1000000], base=16))
# println(hex2num.(x[1:1000]))
# println(read(file)[1:1000])

# print(parse(Int, "0xab"[3:4], base=16))

# println(chop.(read(file)[1:1000], head=2))
# println(parse.(Int, string.(x[1:100], base=16), base=16))


# plot(parse.(Int, string.(x[10000:11000], base=16), base=16))
# println(parse.(Int, string.(read(file), base=16), base=16))


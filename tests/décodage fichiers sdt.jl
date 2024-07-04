
using PyCall
using StructTypes
using Base64
using DelimitedFiles
using ZipFile
# using JuliaDB
using Plots
# using GR

gr()

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

filepath = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/Fichier_test/single_c003.sdt"


open(filepath,"r") do io
    seek(io, 14)
    header = read(io, 2)
    seek(io, header[2]*16^2+header[1] + 22)

    # println(header)

    # file_info = read(io, 22)

    # println(file_info)

    vector = read(io)
    new_vector = UInt16[]

    for i in 1:2:length(vector)
        push!(new_vector, vector[i+1] * 16^2 + vector[i])
    end

    println((new_vector))

    # zip = IOBuffer(read(io, file_info[10]*16^6 + file_info[9]*16^4 + file_info[8]*16^2 + file_info[7] - (file_info[6]*10^6 + file_info[5]*16^4 + file_info[4]*16^2 + file_info[3])))

    # array = reinterpret(UInt16, read(ZipFile.Reader(zip).files[1], file_info[22]*16^6 + file_info[21]*16^4 + file_info[20]*16^2 + file_info[19]))

    # println(size(array))
end




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


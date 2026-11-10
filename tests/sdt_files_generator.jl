include("SdtFile.jl")
using .SdtFile
using FFTW
using Statistics
using ProgressBars

function poisson_rand(lambda::Float64)::Int
    lambda <= 0 && return 0
    if lambda < 30
        L = exp(-lambda)
        k = 0
        p = 1.0
        while true
            k += 1
            p *= rand()
            p <= L && return k - 1
        end
    end
    return max(0, round(Int, lambda + sqrt(lambda) * randn()))
end

function load_irf_counts(irf_path::String)
    s = SdtFile.read_sdt(irf_path)
    counts = Float64.(vec(s.blocks[1].data))
    counts .-= round(median(counts))
    counts[counts .<= 0] .= 0
    return counts
end

function decay_curve(t_ns::Vector{Float64}, lifetimes::Vector{Float64}, amplitudes::Vector{Float64})
    length(lifetimes) == length(amplitudes) || error("lifetimes and amplitudes must have the same length")
    y = zeros(Float64, length(t_ns))
    for (tau, A) in zip(lifetimes, amplitudes)
        y .+= A .* exp.(-t_ns ./ tau)
    end
    return y
end

function convolve_circular(irf::Vector{Float64}, decay::Vector{Float64})
    y = real.(ifft(fft(irf) .* fft(decay)))
    return y ./ sum(y)
end

function write_u16!(buf::Vector{UInt8}, pos::Int, v::Integer)
    buf[pos]   = UInt8(v & 0xff)
    buf[pos+1] = UInt8((v >> 8) & 0xff)
end

function write_u32!(buf::Vector{UInt8}, pos::Int, v::Integer)
    buf[pos]   = UInt8(v & 0xff)
    buf[pos+1] = UInt8((v >> 8) & 0xff)
    buf[pos+2] = UInt8((v >> 16) & 0xff)
    buf[pos+3] = UInt8((v >> 24) & 0xff)
end

function encode_block_header(data_offs::Int, next_block_offs::Int, block_type::Integer, meas_desc_block_no::Integer, lblock_no::Integer, block_length::Integer)
    buf = zeros(UInt8, 22)
    write_u32!(buf, 3, data_offs)
    write_u32!(buf, 7, next_block_offs)
    write_u16!(buf, 11, block_type)
    write_u16!(buf, 13, meas_desc_block_no)
    write_u32!(buf, 15, lblock_no)
    write_u32!(buf, 19, block_length)
    return buf
end

"""
Splices a second MEASURE_INFO entry (cloned from the first) and a second DECAY_BLOCK after the
template's own data block, mirroring how real dual-detector bh hardware lays out multi-channel
files (as in puffs_20ms_2_10X.sdt: one FILE_HEADER, N contiguous MEASURE_INFO blocks, N chained
BLOCK_HEADER+DATA blocks). Returns the patched skeleton bytes (data payloads still zeroed) and the
1-based byte offset of each channel's data region within it.
"""
function two_channel_skeleton(template_bytes::Vector{UInt8}, hdr::SdtFile.FileHeader, n_bins::Int)
    meas_off = Int(hdr.meas_desc_block_offset)
    meas_len = Int(hdr.meas_desc_block_length)
    meas_end = meas_off + meas_len

    prefix = template_bytes[1:meas_end]
    meas_desc_copy = template_bytes[meas_off+1:meas_end]

    n_bytes = n_bins * sizeof(UInt16)
    header1_start = meas_end + meas_len
    data1_start = header1_start + 22
    header2_start = data1_start + n_bytes
    data2_start = header2_start + 22

    header1 = encode_block_header(data1_start, header2_start, 0x0003, 0, 1, n_bytes)
    header2 = encode_block_header(data2_start, 0, 0x0003, 1, 2, n_bytes)
    placeholder = zeros(UInt8, n_bytes)

    bytes = vcat(prefix, meas_desc_copy, header1, placeholder, header2, placeholder)

    write_u32!(bytes, 15, header1_start)  # data_block_offset
    write_u16!(bytes, 19, 2)              # no_of_data_blocks
    write_u16!(bytes, 29, 2)              # no_of_meas_desc_blocks
    write_u32!(bytes, 35, 2)              # reserved1

    return bytes, [data1_start + 1, data2_start + 1]
end

"""
Generates `n` synthetic .sdt files in `output_dir`, named "\$(str)_0001.sdt", "\$(str)_0002.sdt", ...
Each file clones the binary layout of `template_path` (same header/info/setup blocks) and simulates
one 256-bin DECAY_BLOCK per channel: a multi-exponential decay convolved with an IRF, scaled to a
target total photon count and drawn with Poisson shot noise. `photons` (and `photons2`) may be a
single value (shared by all files) or a vector of length `n` (one per file).

By default only channel 1 is written (byte-identical layout to `template_path`). Pass `photons2` to
also write a second channel with its own `lifetimes2`/`amplitudes2`/`photons2` (and optionally its
own `irf_path2`) — this appends a second MEASURE_INFO + DECAY_BLOCK, matching how real dual-detector
bh acquisitions (e.g. puffs_20ms_2_10X.sdt) store multiple channels in a single file. Note that
`open_SDT_file`/`vec_to_lifetime` in GUI/functions.jl currently only ever read channel 1; reading
channel 2 back out needs a small extension there.
"""
function generate_sdt_files(output_dir::String, str::String, n::Int, photons;
        template_path::String="/Users/eliasouellet-oviedo/Documents/Stages/Stage2/Partie2/Données/2026-03-12/cell1/60_1/single_c1.sdt",
        irf_path::String="/Users/eliasouellet-oviedo/Documents/Maîtrise/IRF/2025-11-19_2.sdt",
        lifetimes::Vector{Float64}=[3.0],
        amplitudes::Vector{Float64}=[1.0],
        photons2=nothing,
        lifetimes2::Vector{Float64}=[3.0],
        amplitudes2::Vector{Float64}=[1.0],
        irf_path2::Union{String,Nothing}=nothing)
    mkpath(output_dir)

    template_bytes = read(template_path)
    template = SdtFile.read_sdt(template_bytes, basename(template_path))
    blk = template.blocks[1]
    blk.block_type_info.compressed && error("Template data block is compressed; byte-patching requires an uncompressed DECAY_BLOCK")

    n_bins = length(blk.data)
    t_ns = template.times[1] .* 1e9

    irf = load_irf_counts(irf_path)
    length(irf) == n_bins || error("IRF has $(length(irf)) bins, template has $n_bins")
    shape = convolve_circular(irf, decay_curve(t_ns, lifetimes, amplitudes))

    if photons2 === nothing
        base_bytes = template_bytes
        offsets = [blk.data_offset]
        shapes = [shape]
        photon_specs = Any[photons]
    else
        irf2 = load_irf_counts(something(irf_path2, irf_path))
        length(irf2) == n_bins || error("Second IRF has $(length(irf2)) bins, template has $n_bins")
        shape2 = convolve_circular(irf2, decay_curve(t_ns, lifetimes2, amplitudes2))

        base_bytes, offsets = two_channel_skeleton(template_bytes, template.header, n_bins)
        shapes = [shape, shape2]
        photon_specs = Any[photons, photons2]
    end

    pad = max(4, ndigits(n))

    for i in ProgressBar(1:n)
        out_bytes = copy(base_bytes)
        for (offset, ch_shape, pspec) in zip(offsets, shapes, photon_specs)
            p = pspec isa AbstractVector ? pspec[i] : pspec
            counts = UInt16.(clamp.(poisson_rand.(p .* ch_shape), 0, Int(typemax(UInt16))))
            copyto!(out_bytes, offset, reinterpret(UInt8, counts))
        end

        filename = string(str, "_", lpad(i, pad, '0'), ".sdt")
        write(joinpath(output_dir, filename), out_bytes)
    end

    return nothing
end

dir = "/Users/eliasouellet-oviedo/Desktop"

# generate_sdt_files(joinpath(dir, "single_channel"), "file", 1000, 5000; lifetimes=[4.0], amplitudes=[1.0])
generate_sdt_files(joinpath(dir, "double_channel"), "file", 9999, 1000; lifetimes=[3.5], amplitudes=[1.0], photons2=4000, lifetimes2=[4.5], amplitudes2=[1.0])

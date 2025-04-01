using Printf
using NativeFileDialog

function get_folder_path(str)
    if isfile(str*"_folder_path.txt")
        folderpath = open(f->read(f, String), str*"_folder_path.txt")

        if isempty(folderpath)
            folderpath = pick_folder()

            open(str*"_folder_path.txt", "w") do io
                write(io, folderpath)
            end
        end
    else
        folderpath = pick_folder()

        open(str*"_folder_path.txt", "w") do io
            write(io, folderpath)
        end
    end

    if isempty(folderpath)
        println("Folder does not exist")

        folderpath = get_folder_path(str)
    end

    return folderpath
end

function verify_folder(folderpath)
    if isdir(joinpath(folderpath, "C1")) || isdir(joinpath(folderpath, "C2")) || isdir(joinpath(folderpath, "C3"))
        folderpath_correct = true
    elseif isdir(joinpath(folderpath, "Bliq VMS", "C1")) || isdir(joinpath(folderpath, "Bliq VMS", "C2")) || isdir(joinpath(folderpath, "Bliq VMS", "C3"))
        folderpath = joinpath(folderpath, "Bliq VMS")

        folderpath_correct = true
    else
        println("Folder error")
        folderpath_correct = false
    end

    return folderpath, folderpath_correct
end

function create_new_folder(base_name, destination_folderpath, first, last)
    for i in first:last
        folder_path = joinpath(destination_folderpath, @sprintf("%s-%03d", base_name, i))

        if !ispath(folder_path)
            mkdir(folder_path)
            mkdir(joinpath(folder_path, "Bliq VMS"))
            mkdir(joinpath(folder_path, "Bliq VMS", "C1"))
            mkdir(joinpath(folder_path, "Bliq VMS", "C2"))
            mkdir(joinpath(folder_path, "Bliq VMS", "C3"))
            # println("Created folder: $(@sprintf("%s-%03d", base_name, i))")

            return joinpath(folder_path, "Bliq VMS")
        end
    end

    error("All folder names from %03d to %03d already exist.", first, last)
end

function get_file_number(filename)
    # try
    #     println(fieldnames(typeof(match(r"T(.*?).tif", filename))))
    #     println(match(r"T(.*?).tif", filename).captures[1])
    #     println(parse(Int64, match(r"T(.*?).tif", filename).captures[1]))
    # catch e
    #     println(e)
    # end

    return match(r"T(.*?).tif", filename) !== nothing ? parse(Int64, match(r"T(.*?).tif", filename).captures[1]) : nothing
end

source_folderpath = get_folder_path("source")
destination_folderpath = get_folder_path("destination")

source, folder_correct = verify_folder(source_folderpath)

if folder_correct
    print("Folder correct")

    destination = create_new_folder("manual_VMS", destination_folderpath, 1, 999)

    source_channels = [joinpath(source, "C1"), joinpath(source, "C2"), joinpath(source, "C3")]
    destination_channels = [joinpath(destination, "C1"), joinpath(destination, "C2"), joinpath(destination, "C3")]
    available_channels = [true, true, true]
    source_files = []

    global default_channel = 1

    for n in 1:3
        try
            source_file = readdir(source_channels[n])

            if source_file[1] == ".DS_Store"
                deleteat!(source_file, 1)
            end

            sort!(source_file, by=get_file_number)

            push!(source_files, source_file)
        catch e
            push!(source_files, [])

            if n == default_channel
                default_channel += 1
            end

            available_source[n] = false
        end
    end

    for index in 1:length(source_files[default_channel])
        for channel in 1:3
            if rand(0:1:100) >= 2
                cp(joinpath(source_channels[channel], source_files[channel][index]), joinpath(destination_channels[channel], source_files[channel][index]), force=true)
                if rand(0:1:100) >= 98
                    sleep(rand(0:0.001:0.5))
                end

            end
        end

        sleep(1/31.25)
        if rand(0:1:100) >= 98
            sleep(rand(0:0.001:0.5))
        end
    end
else
    println("Folder error")
end



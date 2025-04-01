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
    return match(r"c(.*?).sdt", filename) !== nothing ? parse(Int64, match(r"c(.*?).sdt", filename).captures[1]) : nothing
end

source_folderpath = get_folder_path("source")
destination_folderpath = get_folder_path("destination")


source_files = readdir(source_folderpath)

if source_files[1] == ".DS_Store"
    deleteat!(source_files, 1)
end

sort!(source_files, by=get_file_number)

for source_file in source_files
    cp(joinpath(source_folderpath, source_file), joinpath(destination_folderpath, source_file), force=true)
    sleep(0.1)
end

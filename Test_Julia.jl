using FilePathsBase

function watch_directory(directory::AbstractString)
    # Convert to Path object
    dir_path = Path(directory)
    
    # Previous state of the directory
    prev_files = Set(readdir(string(dir_path)))
    
    while true
        sleep(0.1)  # Check every 1 second
        
        current_files = Set(readdir(string(dir_path)))
        
        # Detect newly created files
        new_files = setdiff(current_files, prev_files)
        if !isempty(new_files)
            for new_file in new_files
                println("New file created: ", joinpath(string(dir_path), new_file))
            end
        end
        
        # Update the previous state
        prev_files = current_files
    end
end

# Directory to be monitored
directory_to_watch = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/Fichier_test"

# Start watching the directory
watch_directory(directory_to_watch)






# path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/Fichier_test"
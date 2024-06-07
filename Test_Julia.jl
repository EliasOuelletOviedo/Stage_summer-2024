using Gtk4

# Function to open a folder selection dialog
function select_folder_dialog(title::String="Select Folder")
    dialog = Gtk4.FileChooserDialog(
        title,
        nothing,
        Gtk4.FileChooserAction.SELECT_FOLDER,
        "Cancel" => Gtk4.ResponseType.CANCEL,
        "Select" => Gtk4.ResponseType.ACCEPT
    )
    
    # Run the dialog and capture the response
    response = Gtk4.Dialog.run(dialog)
    
    if response == Gtk4.ResponseType.ACCEPT
        folder = Gtk4.FileChooser.get_filename(dialog)
        Gtk4.Widget.destroy(dialog)
        return folder
    else
        Gtk4.Widget.destroy(dialog)
        return nothing
    end
end

# Example usage
selected_folder = select_folder_dialog("Select a folder")
println("Selected folder: ", selected_folder)

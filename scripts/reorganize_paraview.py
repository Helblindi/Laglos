import os
import sys
import shutil
import xml.etree.ElementTree as ET


assert len(sys.argv) == 3, "This file needs 2 input argument: parent_dir and target_dir, but " + str(len(sys.argv)) + " were given."
# Define the parent directory containing all the subdirectories
parent_dir = str(sys.argv[1])
# Define the target directory where the files will be moved
target_dir = str(sys.argv[2])

# Ensure the target directory exists
os.makedirs(target_dir, exist_ok=True)

# Traverse through all subdirectories
for root, dirs, files in os.walk(parent_dir):
    for file in files:
        if file == "data.pvtu":
            # Get the full path of the current .pvtu file
            pvtu_file_path = os.path.join(root, file)
            
            # Get the name of the current directory
            dir_name = os.path.basename(root)
            
            # Create the new .pvtu file name
            new_pvtu_file_name = f"{dir_name}.pvtu"
            new_pvtu_file_path = os.path.join(target_dir, new_pvtu_file_name)
            
            # Move and rename the .pvtu file
            shutil.move(pvtu_file_path, new_pvtu_file_path)
            print(f"Moved and renamed: {pvtu_file_path} -> {new_pvtu_file_path}")
            
            # Parse the .pvtu file to find referenced .vtu files
            tree = ET.parse(new_pvtu_file_path)
            root_elem = tree.getroot()
            
            # Update the paths of the referenced .vtu files
            for piece in root_elem.findall(".//Piece"):
                source = piece.get("Source")
                if source:
                    # Get the full path of the referenced .vtu file
                    vtu_file_path = os.path.join(root, source)
                    
                    # Create the new .vtu file name
                    new_vtu_file_name = f"{dir_name}_{os.path.basename(source)}"
                    new_vtu_file_path = os.path.join(target_dir, new_vtu_file_name)
                    
                    # Move and rename the .vtu file
                    shutil.move(vtu_file_path, new_vtu_file_path)
                    print(f"Moved and renamed: {vtu_file_path} -> {new_vtu_file_path}")
                    
                    # Update the reference in the .pvtu file
                    piece.set("Source", new_vtu_file_name)
            
            # Save the updated .pvtu file
            tree.write(new_pvtu_file_path)
            print(f"Updated references in: {new_pvtu_file_path}")
import os
import zipfile


# zips up the template file
def zip_template():
    # Define the directory and zip file paths
    dir_path = 'molecularnodes/assets/template/Molecular Nodes'
    zip_file_path = 'molecularnodes/assets/template/Molecular Nodes.zip'

    # Create a ZipFile object in write mode
    with zipfile.ZipFile(zip_file_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        # Walk the directory tree and add files to the zip file
        for root, dirs, files in os.walk(dir_path):
            for file in files:
                # Get the path of the file
                file_path = os.path.join(root, file)
                # Add the file to the zip file
                zipf.write(file_path, arcname=os.path.relpath(
                    file_path, start='molecularnodes/assets/template/'))


if __name__ == "__main__":
    zip_template()

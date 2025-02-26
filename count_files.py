import os


def count_files_with_fin_oscar(directory):
    count = 0
    for root, dirs, files in os.walk(directory):
        for file in files:
            if "fin.oscar" in file:
                print(file)
                count += 1
    return count


# Specify the directory to search
directory_path = "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-20-30-deuteron/"  # Change this to the target directory

# Count the files
file_count = count_files_with_fin_oscar(directory_path)
print(f"Number of files containing 'fin.oscar': {file_count}")

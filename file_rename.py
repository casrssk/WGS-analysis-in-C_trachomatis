#import libraries
import sys
import os


def file_rename(path):
    """
    renames the file in the format that is acceptable for the snakemake pipeline
    :param txtfile: .txt file path as a string
    :type str
    """
    for file in os.listdir(path):
        if not file.startswith('.') and file.endswith(".fastq.gz"):
            filepath1 = os.path.join(path, (file))
            sample_name = (os.path.splitext(file)[0]).split(".")[0]
            name_list = sample_name.split("_")
            name_list = ["1" if x == "R1" else "2" if x == "R2" else x for x in name_list]
            new_name = "_".join(name_list[:-1])
            os.rename(os.path.join(path, file), os.path.join(
                path, ''.join([new_name, '.fastq.gz'])))

# file_rename("./linux/sample/results/ct_fastq/")

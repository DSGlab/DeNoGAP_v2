import sys
import os
import sqlite3
from Bio import SeqIO

seq_dir=os.path.abspath(sys.argv[1])


def get_file_name_from_dir(directory):
	"""
	This function reads content within a directory and return name of the files 
	in the directory. key=file-name, value=full_file_path
	Input: directory name
	Return: dict
	"""
	return {os.path.splitext(name)[0]:os.path.join(directory, name) for name in os.listdir(directory)
		if os.path.isfile(os.path.join(directory, name))}
		

seq_files_dict=get_file_name_from_dir(seq_dir)

seq_files=seq_files_dict.values()

seq_index=SeqIO.index_db("../data/test_seq.idx",seq_files,"fasta")

print("%i sequences indexed" % len(seq_index))

print seq_index["sp|Q87YW7|1A1D_PSESM"].id
print seq_index["sp|Q87YW7|1A1D_PSESM"].seq

		
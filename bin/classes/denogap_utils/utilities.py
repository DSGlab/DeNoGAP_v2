#!/usr/bin/env python
import sys
import os
import subprocess
from Bio import SeqIO
from Bio import SearchIO
from collections import defaultdict

class SequenceFileUtilities:
	"""
	This is class for DeNoGAP utilities
	"""
	def get_file_name_from_dir(self,directory):
		"""
		This function reads content within a directory and return name of the files 
		in the directory. key=file-name, value=full_file_path
		Input: directory name
		Return: dict
		"""
		return {os.path.splitext(name)[0]:os.path.join(directory, name) for name in os.listdir(directory)
			if os.path.isfile(os.path.join(directory, name))}
	
	def mkfastadict(self,fastapath):

		"""
		This function takes filepath for fasta file and
		prepares 2-dimensional sequence dict where first key is file name
		and second key is sequence id. The value for the dict is string of sequence.
		""" 
		seq_dict={}
	
		fasta_sequences = SeqIO.parse(open(fastapath),'fasta')
	
		for fasta in fasta_sequences:
			seq_dict[fasta.id]=fasta.seq
		
		return(seq_dict)
	
	def getseq_forid(self,seqfile,id):
	
		"""
		This function returns fasta sequence for given ID from the sequence file
		"""	
		seq_dict={}
		
		fasta_sequences = SeqIO.parse(open(seqfile),'fasta')
		
		for fasta in fasta_sequences:
			if fasta.id==id:
				seq_dict[fasta.id]=fasta.seq
				break
				
		return(seq_dict)
		
		
	def write_seqfile(self,dirpath,seqfile_name,seqdb_dict):

		"""
		This function take path for sequence database directory and 
		dict containing sequences to write in the database as an argument.
		The output is sequence database in fasta format.
		"""
	
		with open(os.path.join(dirpath,seqfile_name),"w") as seqfile:
			for name,seq_dict in seqdb_dict.iteritems():
				for seq_id in seq_dict:		
					seqfile.write(">{}|{}\n{}\n".format(name,seq_dict[seq_id].id,seq_dict[seq_id].seq))
	
		seqfile.close()			
			
		return(os.path.join(dirpath,seqfile_name))		
	
	
	def build_msa(self,seqfile,msa_dir):
		"""
		This function builds multiple sequence alignment for the input sequence
		file and returns msa file in fasta format
		"""
		
		cluster_id=os.path.splitext(os.path.basename(seqfile))[0]
		msa_file=os.path.join(msa_dir,cluster_id+".align.fasta")
		
		mafft_cmd="mafft --auto {0} > {1}".format(seqfile,msa_file)
		
		proc=subprocess.Popen([mafft_cmd], shell=True,
		                         stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		proc.wait()
		stdout,stderr=proc.communicate()
		
		if proc.returncode!=0:
			print("Mafft execution failed, Exiting with an error {}".format(stderr))
			sys.exit()
		
		return(msa_file)                                             
		

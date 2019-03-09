#!/usr/bin/env python

from __future__ import division
import sys
import os
import re
import csv
import gc
import argparse
import shutil
import logging
import datetime
import subprocess
from Bio import SeqIO
from Bio import SearchIO
from collections import defaultdict
from argparse import RawTextHelpFormatter
from pathos.multiprocessing import ProcessPool
from classes.denogap_utils.utilities import SequenceFileUtilities as su
from classes.denogap_utils.hmmer import Hmmer
from classes.denogap_utils.dirsetup import SetDir
from classes.denogap_utils.multihsp import FixMultiHSP
from classes.denogap_utils.clustering import SequenceClustering as sc

"""
Shalabh Thakur
created October, 2018
modified November, 2018
DeNoGAP version 2.0
Guttman Lab, University of Toronto
"""
"""
Description: This script run all-by-all pairwise alignment between protein sequences 
using PHMER
"""

##### MAIN FUNCTION #######
#@profile
def main():

	"""
	Main Function
	"""
	now = datetime.datetime.now()
	
	c_args = parse_args(__file__)

	#### Define global variables ####	
	genome_info_file=os.path.abspath(c_args['genome_info'])
	project_name=c_args['project_name']
	seq_dir=os.path.abspath(c_args['seq_dir'])
	output_dir=os.path.abspath(c_args['output_dir'])
	
	filterparam={"match_cutoff":c_args["percent_match"],
				"hitcov_cutoff":c_args["hitcoverage"],
				"chicov_cutoff":c_args["chimcoverage"]}
	
	iteration=0                                     #### iteration number
	resume_iteration=None                           #### resume iteration from this number
	iter_dir=None                                   #### path for current iteration directory
	project_dir=None                                #### path for project directory 
	projectdir_dict=defaultdict()                   #### dict of project dir output
	iterdir_dict=defaultdict()                      #### dict of iter dir output
	current_timestamp=now.strftime("%Y%m%d_%H%M%S") ### current timestamp 
		
		
	if c_args["resume_iter"]!=None:
		resume_iteration=int(c_args['resume_iter'])
		
	### create project_dir paths ###
	project_dir=os.path.join(output_dir,project_name)
	
	### log file ###
	log_file=os.path.join(os.path.abspath("../log"),project_name+"_"+current_timestamp+".log")
	
	print "[{}]: Analysis Start\n".format(now.strftime("%Y:%m:%d_%H:%M:%S"))

	#### create log file for the project ####
	logging.basicConfig(filename=log_file,format="[%(asctime)s] %(levelname)-8s %(message)s",datefmt='%a, %d %b %Y %H:%M:%S',level=logging.INFO)
	logging.info("This is a run-log for DeNoGAP2 pipeline executed on {}\n".format(current_timestamp))
	#logging.info("[Command executed] {}".format(command))

	#### read/check genome information file #####
	logging.info("Reading genome information file {}".format(genome_info_file))
	print "[{}]: Reading genome information file {}\n".format(now.strftime("%Y:%m:%d_%H:%M:%S"),genome_info_file)
	
	genome_dict=read_genome_info(genome_info_file)
	
	### list of reference genome names and other genome names ###
	list_ref_genome=[genome_name for genome_name in genome_dict
						if int(genome_dict[genome_name]["REFERENCE"])==1]
	
	list_other_genome=[genome_name for genome_name in genome_dict
						if int(genome_dict[genome_name]["REFERENCE"])==0]				
	
	### read names for the fasta files in the data directory #####
	logging.info("Reading fasta file names from the sequence directory {}".format(seq_dir))
	print "[{}]: Reading fasta file names from the sequence directory {}\n".format(now.strftime("%Y:%m:%d_%H:%M:%S"),seq_dir)
	
	fasta_dict=su().get_file_name_from_dir(seq_dir)
	fasta_list=fasta_dict.values()
	
	#### check if protein sequences for all genome names specified in information file are present in the data directory ####
	for genome_name in genome_dict:
		if not genome_name in fasta_dict:
			logging.error("ERROR, Protein sequence fasta file for {} not found in {}\n".format(genome_name,data_dir))
			print "[{}]: ERROR, Protein sequence fasta file for {} not found in {}\n".format(now.strftime("%Y:%m:%d_%H:%M:%S"),genome_name,data_dir)
			sys.exit()	
	
	#### If not already present set-up project directory #####
	logging.info("Setting up project directory for DeNoGAP-HMM {}".format(project_dir))
	print "[{}]: Setting up project directory for DeNoGAP-HMM {}\n".format(now.strftime("%Y:%m:%d_%H:%M:%S"),project_dir)
	
	if not os.path.exists(project_dir):
		logging.info("Creating new project dir {}\n".format(project_dir))
		print "[{}]: Creating new project dir {}\n".format(now.strftime("%Y:%m:%d_%H:%M:%S"),project_dir)
		os.makedirs(project_dir)
		projectdir_dict=SetDir().mk_project_dirs(project_dir)
	else:
		logging.info("Using existing project dir {}\n".format(project_dir))
		print "[{}]: Using existing project dir {}\n".format(now.strftime("%Y:%m:%d_%H:%M:%S"),project_dir)
		projectdir_dict=SetDir().mk_project_dirs(project_dir)
		
	##### Create and load sequences into index database #####	
	logging.info("Creating Sequence index database {}".format(os.path.join(project_dir,"Sequence_index.idx")))
	print "[{}]: Creating Sequence index database {}\n".format(now.strftime("%Y:%m:%d_%H:%M:%S"),os.path.join(project_dir,"Sequence_index.idx"))
	
	fasta_seq_idx=SeqIO.index_db(os.path.join(projectdir_dict["TMP"],"Sequences_indexdb.idx"),
								 fasta_list,"fasta")					 
		
	#### Initial homolog mapping for reference genomes
	#### Starts with iteration 0
	#### If folder for iteration 0 is already present then program will
	#### terminate further execution to avoid any overwriting
	
	if resume_iteration==None:

		### setup dir for the iteration ###
		iterdir_dict=SetDir().mk_hmm_iter_dirs(project_dir,iteration)
	
		logging.info("Starting pairwise sequence comparision step\n")
		print "[{}]: Starting pairwise sequence comparision step\n".format(now.strftime("%Y:%m:%d_%H:%M:%S"))
		
		#### make sequence database ###
		logging.info("Creating sequence database for pairwise comparision\n")
		print "[{}]: Creating sequence database for pairwise comparision\n".format(now.strftime("%Y:%m:%d_%H:%M:%S"))

		seqdb_dict=defaultdict(dict)
		
		for genome_name in list_ref_genome:
			seq_dict=SeqIO.index(fasta_dict[genome_name],'fasta')
			seqdb_dict[genome_name]=seq_dict
			seq_dict.close()
				
		seqdb_dir=projectdir_dict["SEQ_DB"]
		seqdb_file="DBSEQ.fasta"
		
		logging.info("Saved sequence database {} for pairwise comparision at {}\n".format(seqdb_file,seqdb_dir))
		print "[{}]: Saved sequence database {} for pairwise comparision {}\n".format(now.strftime("%Y:%m:%d_%H:%M:%S"),seqdb_file,seqdb_dir)		
				
		seqdb_path=su().write_seqfile(seqdb_dir,seqdb_file,seqdb_dict)
		
		#### Align query sequences against sequence database ####
		
		logging.info("Aligning Sequences\n")
		print "[{}]: Aligning Sequences\n".format(now.strftime("%Y:%m:%d_%H:%M:%S"))
		
		parsed_alignment_dict=defaultdict(dict)		
		
		for genome_name in list_ref_genome:
			
				logging.info("QUERY: {}\n".format(genome_name))
				print "[{}]: QUERY: {}\n".format(now.strftime("%Y:%m:%d_%H:%M:%S"),genome_name)
					
				outpath=os.path.join(iterdir_dict["ALL_MATCH"],genome_name+".hmmalign.txt")
				domtabpath=os.path.join(iterdir_dict["ALL_MATCH"],genome_name+".domtab.txt")
				bestmatchpath=os.path.join(iterdir_dict["BEST_MATCH"],genome_name+".bhh.txt")
				partialmatchpath=os.path.join(iterdir_dict["PARTIAL_MATCH"],genome_name+".phh.txt")
				chimeramatchpath=os.path.join(iterdir_dict["CHIMERIC_MATCH"],genome_name+".chh.txt")                                                             
					                        
						
				hmmer_proc_returncode,hmmerstderr=Hmmer("phmmer").run_hmmer(c_args,
				                                                            seqdb_path,
				                                                            fasta_dict[genome_name],
				                                                            outpath,
				                                                            domtabpath)
					
				if hmmer_proc_returncode!=0:
					print "[{}]: HMMER Execution Failed, Exiting with an error {}, "\
						"returncode: {}\n".format(now.strftime("%Y:%m:%d_%H:%M:%S"),
												  hmmerstderr,
					                              hmmer_proc_returncode)
					logging.error("HMMER Execution Failed, Exiting with an error {}, "\
						"returncode: {}\n".format(hmmerstderr,
												  hmmer_proc_returncode))
					sys.exit()
						
				domhmmer_dict=Hmmer("phmmer").parse_hmmer_domtab(domtabpath)
				simhmmer_dict=Hmmer("phmmer").parse_hmmer_similarity(outpath)
				hmmer_result=Hmmer("phmmer").add_hmmer_stats(domhmmer_dict,simhmmer_dict)
					
				domhmmer_dict=None
				simhmmer_dict=None
				gc.collect()
					
				fixed_hmmer_result=FixMultiHSP("phmmer").fix_hsps(hmmer_result,
																  c_args['avg_accuracy'])
				filtered_hmmer_result=Hmmer("phmmer").filter_hits(filterparam,
																  fixed_hmmer_result)
					
				hmmer_result=None
				fixed_hmmer_result=None
				gc.collect()
					
				write_hmmer_result(genome_name,
					                bestmatchpath,
					                filtered_hmmer_result["BEST"])
				write_hmmer_result(genome_name,
					                partialmatchpath,
					                filtered_hmmer_result["PARTIAL"])
				write_hmmer_result(genome_name,
					                chimeramatchpath,
					                filtered_hmmer_result["CHIMERA"])
		
		
		logging.info("MCL Clustering Best-Hits\n")
		print "[{}]: MCL Clustering Best-Hits\n".format(now.strftime("%Y:%m:%d_%H:%M:%S"))
		
		mcl_abc_file=os.path.join(projectdir_dict["TMP"],"seq.abc")
			                   	
		os.system("cat {}/* | cut -f 1,2,11 | sed '1d' > {}".format(
				iterdir_dict["BEST_MATCH"],mcl_abc_file))
				
		mcl_cluster=sc().mcl_clustering(mcl_abc_file,c_args["mcl_inflation"])
		
		logging.info("Adding Group-IDs\n")
		print "[{}]: Adding Group-ID\n".format(now.strftime("%Y:%m:%d_%H:%M:%S"))
		
		grouped_cluster=sc().add_cluster_ids(mcl_cluster,start_at=1000)
		
		hmm_cluster_file=os.path.join(iterdir_dict['CLUSTER'],
		                              "Hmmcluster_"+now.strftime("%Y%m%d_%H%M%S")+".txt")
		
		logging.info("Saving Cluster file at {}\n".format(hmm_cluster_file))
		print "[{}]: Saving Cluster file at {}\n".format(now.strftime("%Y:%m:%d_%H:%M:%S"),
		                                                 hmm_cluster_file)
		
		cp_cmd="cp {0} {1}".format(grouped_cluster,hmm_cluster_file)
		proc=subprocess.Popen([cp_cmd],shell=True,stdout=subprocess.PIPE,
		                      stderr=subprocess.PIPE)
		proc.wait()
		stdout,stderr=proc.communicate()
		
		if proc.returncode!=0:
			print stderr
			sys.exit()
		else:	                     
			logging.info("Pairwise sequence comparision completed successfully\n")
			print "[{}]: Pairwise sequence comparision completed successfully\n"\
				  .format(now.strftime("%Y:%m:%d_%H:%M:%S"))
				  
	elif resume_iteration!=None:
		"""
		Execute iterative block for each genome
		read cluster file from iteration to resume from
		build hmm models/ sequence database
		scan hmm models/ sequence database iteratively against each genome
		"""
		
		while len(list_other_genome)!=0:
		
			iteration=int(resume_iteration)+1
			iterdir_dict=SetDir().mk_hmm_iter_dirs(project_dir,iteration)
			
			### getting model cluster directory #####
			model_cluster_dir=os.path.join(projectdir_dict["BASE"],
										   "iter_{0}".format(resume_iteration),"CLUSTER")
			
			### getting model cluster file ####							   
			model_cluster_file=[os.path.join(model_cluster_dir,file_name) 
			                     for file_name in os.listdir(model_cluster_dir)
			                        if file_name.startswith("Hmmcluster")]
			
			### Call function to build hmm models and if singleton then add 
			### to sequence database                        
			#hdb().build_hmmdb(model_cluster_file[0],seq_dir,
			#								outdir_dict["HMM_DB"],
			#								outdir_dict["SEQ_DB"])
			model_dir=os.path.join(projectdir_dict["HMMER_DB"],"HMM_MODELS")
			os.makedirs(model_dir)
			
			with open(model_cluster_file[0],"r") as clustfile:
			
				pool=ProcessPool(nodes=c_args["cpu"])
				[pool.apply(build_hmm,args=(clustline,fasta_dict,model_dir))
					     for clustline in clustfile]	     
											                 
		
##### FUNCTIONS #######
def parse_args(desc):

	"""
	Command line option parser
	
	:param desc: Short program description.
	:type desc: str
	
	:return arg: dict of command line arguments with key=command line
        argument and val=argument value
        
	:return_type: dict
	"""

	parser = argparse.ArgumentParser(description=desc, formatter_class=RawTextHelpFormatter)
	parser._optionals.title="Argument for DeNoGAP2 pairwise sequence comparision"
	
	#### Mandatory arguments ####
	mandatoryArguments = parser.add_argument_group('Mandatory arguments for DeNoGAP-HMM')
	mandatoryArguments.add_argument("--genome_info", metavar="<FILE PATH>",
									help="specify name and path of the tab-delimited "\
	                                     "genome information file", 
	                                required=True)
	mandatoryArguments.add_argument("--project_name", metavar="<STR>", 
									help="specify name of the project", required=True)
	mandatoryArguments.add_argument("--seq_dir", metavar="<DIR PATH>", 
									help="specify name and path to the directory "\
										 "containing sequences in fasta format", 
									required=True)
	mandatoryArguments.add_argument("--output_dir", metavar="<DIR PATH>", 
									help="specify path to the output directory", 
									required=True)
	mandatoryArguments.add_argument("--resume_iter",metavar="<INT>",type=int,
									help="specify last iteration number from which "\
										"analysis should be resumed") 
	
	seq_comparision=parser.add_mutually_exclusive_group(required=True)
	seq_comparision.add_argument("--pairwise",action="store_true",
								help="perform pairwise sequence comparision")
	seq_comparision.add_argument("--iterative",action="store_true",
								help="perform iterative sequence comparision")
	
	#### required arguments for hmmer iterative alignment ####
	mandatoryIterativeArguments= parser.add_argument_group('Arguments for iterative alignment')
	mandatoryIterativeArguments.add_argument("--inclust",metavar="<FILE PATH>", 
											default="NA", 
											help="use specified file as initial "\
												"cluster file for iterative hmmer alignment")	
	mandatoryIterativeArguments.add_argument("--seqdb",metavar="<FILE PATH>", 
											default="NA", 
											help="use specified file as sequence database "\
											"during pairwise or iterative hmmer alignment")	
	mandatoryIterativeArguments.add_argument("--hmmdb",metavar="<FILE PATH>", 
											default="NA", 
											help="use specified file as hmm database "\
											"during iterative hmmer alignment")
		
	#### Optional arguments for hmmer ####
	OptionalArguments= parser.add_argument_group('Optional arguments for DeNoGAP-HMM')
	OptionalArguments.add_argument("--evalue",metavar="<FLOAT>", default=10.0,type=float, 
								help="report sequences <= this E-value threshold in output "\
									"[default: 10.0] (x>0)")
	OptionalArguments.add_argument("--cpu",metavar="<INT>", default=2,type=int, 
								help="number of parallel CPU workers to use for "\
									"multithreads [default: 2]")
	OptionalArguments.add_argument("--mcl_inflation",metavar="<FLOAT>", default=1.5,
									type=float, 
									help="inflation value for MCL clustering, "\
										"varying this parameter effects "\
										"clustering granularity [default: 1.5]")

	#### Optional arguments for parsing alignment output ###
	hmmerparseArguments = parser.add_argument_group('Parsing arguments for DeNoGAP-HMM result')
	hmmerparseArguments.add_argument("--avg_accuracy", metavar="<FLOAT>",default=0.7,
									choices=range(0, 1),type=float,
									help="parse hits >= this accuracy threshold as "\
										"significant hit [default: 0.7] (x>0)")
	hmmerparseArguments.add_argument("--percent_match", metavar="<FLOAT>",default=70.0,
									choices=range(0, 100),type=float,
									help="parse hits with percent sequence match between "\
										"query and target >= this threshold as "\
										"signficant hit [default: 70.0] (x>0)")
	hmmerparseArguments.add_argument("--hitcoverage", metavar="<FLOAT>",default=70.0,
									choices=range(0, 100),type=float,
									help="parse hits >= this query/target coverage "\
										"threshold as significant hit "\
										"[default: 70.0] (x>0)")
	hmmerparseArguments.add_argument("--chimcoverage", metavar="<FLOAT>",default=20.0,
									choices=range(0, 100),type=float,
									help="minimum coverage threshold for query/target "\
									"to be reported as chimera alignment "\
									"[default: 20.0] (x>0)")		
	args=parser.parse_args()

	return vars(args)
	
#### Function to read genome information file #####
def read_genome_info(genome_info_file):

	"""
	This function read tab-delimited genome information file.
	Parse the file and store each column into dict
	"Note: On Mac OS save genome_info_file in Windows-formatted tab-delimited text file"
	"""
	### Initialize dict to store genome information ####
	genome_dict=defaultdict(dict)
	
	genome_info_dict=csv.DictReader(open(genome_info_file,"rU"),dialect="excel-tab")
	for dict_line in genome_info_dict:
		genome_abbreviation=dict_line['GENOME_ABBREVIATION']
		genome_dict[genome_abbreviation]=dict_line
			
	return(genome_dict)
	
def write_hmmer_result(genome_name,outfile,result_hmmer):
	"""
	This function writes parsed hmmer result to the output file
	Input: output file path
           result_dict
	"""
	
	with open(outfile,"w") as result_out:
	
		result_out.write("QUERY_ID\tHMM_ID\tQUERY_LEN\tHMM_LEN\tDOMAIN_INDEX\t"\
		                 "QUERY_START\tQUERY_END\tHMM_START\tHMM_END\t"\
		                 "AVERAGE_ACCURACY\tEVALUE\tDOM_BITSCORE\tDOM_IDENTICAL_COUNT"\
		                 "DOM_SIMILAR_COUNT\tDOM_PIDENT\tDOM_PSIM\tDOM_QCOVERAGE\t"\
		                 "DOM_HMMCOVERAGE\tSEQ_BITSCORE\tSEQ_IDENTICAL_COUNT\t"\
		                 "SEQ_SIMILAR_COUNT\tSEQ_PIDENT\tSEQ_PSIM\tSEQ_QCOVERAGE\t"\
		                 "SEQ_HMMCOVERAGE\n")
		
		for feature in result_hmmer:
			
			feature["QUERY_ID"]=genome_name+"|"+feature["QUERY_ID"]
				
			result_line="{QUERY_ID}\t{HMM_ID}\t{QUERY_LEN}\t{HMM_LEN}\t"\
					    "{DOMAIN_INDEX}\t{QUERY_START}\t{QUERY_END}\t{HMM_START}\t"\
					    "{HMM_END}\t{AVERAGE_ACCURACY}\t{EVALUE}\t{BITSCORE}\t"\
					    "{IDENTICAL_COUNT}\t{SIMILAR_COUNT}\t{PERCENT_IDENTITY}\t"\
					    "{PERCENT_SIMILARITY}\t{QUERY_COVERAGE}\t{HMM_COVERAGE}\t"\
					    "{TOTAL_BITSCORE}\t{TOTAL_IDENTICAL_COUNT}\t{TOTAL_SIMILAR_COUNT}\t"\
					    "{TOTAL_PERCENT_IDENTITY}\t{TOTAL_PERCENT_SIMILARITY}\t"\
					    "{TOTAL_QUERY_COVERAGE}\t{TOTAL_HMM_COVERAGE}\n"\
					    .format(**feature)
					    
			#print result_line		    

			result_out.write(result_line)
	
	result_out.close()


def build_hmm(cluster_line,fasta_files,model_dir):
	"""
	This functions builds hmm model for each cluster
	"""
	# remove new-line character and split cluster line into columns
	cluster_line=cluster_line.rstrip("\n")
	cluster_column=cluster_line.split("\t")
	cluster_id=cluster_column.pop(0)
	cluster_id=cluster_id.rstrip("\:")
	
	print cluster_id
	print mp.current_process()
	
	os.makedirs(os.path.join(model_dir,cluster_id))
	msa_fasta_file=os.path.join(model_dir,cluster_id,cluster_id+".fasta")
	
	with open(msa_fasta_file,"w") as msa_fasta:
		for seq_id in cluster_column:
			split_seqid=seq_id.split("|")
			genome_name=split_seqid.pop(0)
			id="|".join(split_seqid)
			seq_dict=su().getseq_forid(fasta_files[genome_name],id)
			msa_fasta.write(">{}|{}\n{}\n".format(cluster_id,id,seq_dict[id]))
	msa_fasta.close()
	
	#### build msa from multi gene msa_fasta ###
	if len(cluster_column)>1:
		msa_file=su().build_msa(msa_fasta_file,os.path.join(model_dir,cluster_id))
		hmm_file=Hmmer("hmmbuild").run_hmmbuild(msa_file,model_dir)
		
	
	

	
	
if __name__ == '__main__':
	main()
#!/usr/bin/env python
from __future__ import division
import os
import sys
import re
import subprocess
from Bio import SeqIO
from Bio import SearchIO
from collections import defaultdict

class Hmmer:
	"""
	This class is for running Hmmer programs
	"""
	
	def __init__(self,hmmer_program):
		self.hmmer_program=hmmer_program
	
	def run_hmmer(self,hmmer_args,dbpath,querypath,outpath,domtabpath):
		"""
		This function runs sequence-sequence,sequence-profile or profile-sequence
		alignment using HMMER programs
		"""
		### prepare argument for command line ###
		hmm_args=None
		first=None
		second=None
		out=outpath
		domtabout=domtabpath
		cpu=hmmer_args["cpu"]
		evalue=hmmer_args["evalue"]
		
		if self.hmmer_program=="phmmer" or\
		   self.hmmer_program=="hmmsearch":
			first=querypath
			second=dbpath
		elif self.hmmer_program=="hmmscan":
			first=dbpath
			second=querypath
		
		hmm_args="--cpu {0} -E {1} --domE {1} --incE {1} --incdomE {1} -o {2} "\
		         "--domtblout {3} {4} {5}".format(cpu,evalue,out,domtabout,first,second) 		
		
		proc=subprocess.Popen([self.hmmer_program+" "+hmm_args], shell=True,
								stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		proc.wait()
		stdout,stderr=proc.communicate()
		
		return(proc.returncode,stderr)
		
	def parse_hmmer_domtab(self,domtab_file):
		"""
	This function parse result from hmmer domtab alignment output
		"""
		### define dictionary to store parsed domtab output
		domtab_all_dict=defaultdict(dict)
		
		with open(domtab_file,"r") as domtab_result:
			seq_records=SearchIO.parse(domtab_result, 
										"{}3-domtab".format(self.hmmer_program))
			
			for QueryResult in seq_records:
					
				for Hit in QueryResult:
	
					domain_dict=defaultdict(dict)
					
					for HSP in Hit:
						domain_index=HSP.domain_index
						domain_dict[domain_index]["QUERY_ID"]=QueryResult.id
						domain_dict[domain_index]["QUERY_LEN"]=QueryResult.seq_len
						domain_dict[domain_index]["HMM_ID"]=Hit.id
						domain_dict[domain_index]["HMM_LEN"]=Hit.seq_len
						domain_dict[domain_index]["EVALUE"]=Hit.evalue
						domain_dict[domain_index]["BITSCORE"]=Hit.bitscore						
						domain_dict[domain_index]["AVERAGE_ACCURACY"]=HSP.acc_avg
						domain_dict[domain_index]["DOMAIN_INDEX"]=HSP.domain_index
						domain_dict[domain_index]["HMM_END"]=HSP.hit_end
						domain_dict[domain_index]["HMM_START"]=HSP.hit_start
						domain_dict[domain_index]["QUERY_END"]=HSP.query_end
						domain_dict[domain_index]["QUERY_START"]=HSP.query_start
						
						if HSP.query_start==0:
							domain_dict[domain_index]["QUERY_START"]=1
						if HSP.hit_start==0:
							domain_dict[domain_index]["HMM_START"]=1		
					
					domtab_all_dict[QueryResult.id][Hit.id]=domain_dict	
					
		return(domtab_all_dict)				
		
	
	def parse_hmmer_similarity(self,alignment_file):
		"""
		This function parse result from hmmer alignment and 
		fetch similarity/identity values
		"""	
		
		align_dict=defaultdict(dict)
		
		with open(alignment_file,"r") as alignment_result:
			
			flag_block=0
			flag_sim=0
			query_id=""
			target_id=""
			domain_index=""
			first_regx=""
			second_regx=""
			dom_sim_line=defaultdict(dict)
			    
			for line in alignment_result:
			
				if line.startswith("Query:"):
					query_row=re.split(" +",line)
					query_id=query_row[1]
					dom_sim_line=defaultdict(dict)
					continue
				elif line.startswith(">>"):
					if len(dom_sim_line)!=0:
						align_dict[query_id][target_id]=dom_sim_line
					target_row=re.split(" +",line)
					target_id=target_row[1]
					dom_sim_line=defaultdict(dict)
					continue
				elif line.startswith("  == domain"):
					domain_row=re.split(" +",line)
					domain_index=int(domain_row[3])
					flag_block=1
					continue
				elif line.startswith("Internal pipeline"):
					align_dict[query_id][target_id]=dom_sim_line
				
				else:
					if self.hmmer_program=="phmmer" or self.hmmer_program=="hmmsearch":
						first_regx=r"^(\s+)"+ re.escape(query_id)
						second_regx=r"^(\s+)"+ re.escape(target_id)
					elif self.hmmer_program=="hmmscan":
						first_regx=r"^(\s+)"+ re.escape(target_id)
						second_regx=r"^(\s+)"+ re.escape(query_id)
        			
				if flag_block==1 and\
        		   		re.match(first_regx,line):
        		   		flag_sim=1
				elif flag_sim==1:
        				if not domain_index in dom_sim_line:
        					sim_line=line.replace(" ","").rstrip("\n")
        					identical_count=len(sim_line)-sim_line.count("+")
        					similar_count=len(sim_line)
        					dom_sim_line[domain_index]["IDENTICAL_COUNT"]=identical_count
        					dom_sim_line[domain_index]["SIMILAR_COUNT"]=similar_count
        					
        				else:
        					sim_line=line.replace(" ","").rstrip("\n")
        					identical_count=len(sim_line)-sim_line.count("+")
        					similar_count=len(sim_line)
        					dom_sim_line[domain_index]["IDENTICAL_COUNT"]+=identical_count
        					dom_sim_line[domain_index]["SIMILAR_COUNT"]+=similar_count
        				flag_sim=0

		return(align_dict)
		
		
	def add_hmmer_stats(self,domhmmer_dict,simhmmer_dict):

		"""
		This function adds stats to hmmer result after calculation of
		percent identity, similarity and sequence coverage
		"""
		for query_id in domhmmer_dict:
					
			for target_id in domhmmer_dict[query_id]:
						
				for domain_index in domhmmer_dict[query_id][target_id]:
				
					hmmdom_dict=domhmmer_dict[query_id][target_id][domain_index]
					hmmdom_dict.update(simhmmer_dict[query_id][target_id][domain_index])
				
					#### get similar/identical residue counts ##			
					identical_count=hmmdom_dict["IDENTICAL_COUNT"]
					similar_count=hmmdom_dict["SIMILAR_COUNT"]
				
					### get coordinates and length for query and target					
					qlen=hmmdom_dict["QUERY_LEN"]
					hlen=hmmdom_dict["HMM_LEN"]
					qstart=hmmdom_dict["QUERY_START"]
					qend=hmmdom_dict["QUERY_END"]
					hstart=hmmdom_dict["HMM_START"]
					hend=hmmdom_dict["HMM_END"]
					avg_acc=hmmdom_dict["AVERAGE_ACCURACY"]
					qaln_len=((qend-qstart)+1)
					haln_len=((hend-hstart)+1)
				
					### use minimum of query or target to find percent identity/similarity					
					min_len=min(qaln_len,haln_len)
					percent_identity=round(float(identical_count/min_len)*100,2)
					percent_similarity=round(float(similar_count/min_len)*100,2)
				
					### get query/target sequence coverage in the alignment ###
					qcoverage=round(float(qaln_len/qlen)*100,2)
					hcoverage=round(float(haln_len/hlen)*100,2)
											
					hmmdom_dict.update({"IDENTICAL_COUNT":identical_count,
										"SIMILAR_COUNT":similar_count,
										"PERCENT_IDENTITY":percent_identity,
										"PERCENT_SIMILARITY":percent_similarity,
										"QUERY_COVERAGE":qcoverage,
										"HMM_COVERAGE":hcoverage})
					
					domhmmer_dict[query_id][target_id][domain_index]=hmmdom_dict					
				
		return(domhmmer_dict)	  


	def filter_hits(self,filter_param,hmmer_result_dict):
	
		"""
		This function filter hmmer result based on thresholds defined by user on
		command line
		"""
		
		filter_results=defaultdict(list)
		
		for query_id in hmmer_result_dict:
			
			for target_id in hmmer_result_dict[query_id]:
			
				for domain_index in hmmer_result_dict[query_id][target_id]:
					
					hmmdom_dict=hmmer_result_dict[query_id][target_id][domain_index]
					
					#print hmmdom_dict
					
					qstart=int(hmmdom_dict["QUERY_START"])
					qend=int(hmmdom_dict["QUERY_END"])
					qlen=int(hmmdom_dict["QUERY_LEN"])
					hstart=int(hmmdom_dict["HMM_START"])
					hend=int(hmmdom_dict["HMM_END"])
					hlen=int(hmmdom_dict["HMM_LEN"])
					
					### If identity or similarity and qcoverage and hcoverage >= threshold
					### then define best hit.
					### If only one of qcoverage or hcoverage >= threshold then define
					### partial hit
					### But if both qcoverage and hcoverage are less than threshold but
					### greater than chim_threshold check for chimera
					
					feature_list1=[float(hmmdom_dict["TOTAL_PERCENT_IDENTITY"]),
								   float(hmmdom_dict["TOTAL_PERCENT_SIMILARITY"])]
								   
					feature_list2=[float(hmmdom_dict["TOTAL_QUERY_COVERAGE"]),
								   float(hmmdom_dict["TOTAL_HMM_COVERAGE"])]		   
					
					if any(param>=float(filter_param["match_cutoff"])
					   for param in feature_list1):
					
						if all(param>=float(filter_param["hitcov_cutoff"]) 
						   for param in feature_list2):
						   
							filter_results["BEST"].append(hmmdom_dict)
							
							#print hmmdom_dict
						   
						elif any(param>= float(filter_param["hitcov_cutoff"])
						     for param in feature_list2):
						     
							filter_results["PARTIAL"].append(hmmdom_dict)
							
							#print hmmdom_dict
							
						elif all(param>= float(filter_param["chicov_cutoff"]) and
						     param <= float(filter_param["hitcov_cutoff"])
						     for param in feature_list2):
						     
							if qstart<=15 and qend<qlen:
								if hstart<=15 and hend<hlen:
									filter_results["CHIMERA"].append(hmmdom_dict)
									#print hmmdom_dict
							
							elif qstart>1 and qend>=(qlen-15):
								if hstart>1 and hend>=(hlen-15):
									filter_results["CHIMERA"].append(hmmdom_dict)
									#print hmmdom_dict
		
		return(filter_results)					     

	
	def run_hmmbuild(self,msa_file,model_dir):
		"""
		This function builds hmm model from multiple sequence alignment
		using hmmbuild program
		"""
		
		cluster_id=os.path.splitext(os.path.splitext(os.path.basename(msa_file))[0])[0]
		hmm_file=os.path.join(model_dir,cluster_id,cluster_id+".hmm")
		
		hmmbuild_cmd="{0} {1} {2}".format(self.hmmer_program,hmm_file,msa_file)
		
		proc=subprocess.Popen([hmmbuild_cmd], shell=True,
		                         stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		proc.wait()
		stdout,stderr=proc.communicate()
		
		if proc.returncode!=0:
			print("HmmBuild execution failed, Exiting with an error {}".format(stderr))
			sys.exit()
		
		return(hmm_file) 
							     		   
							    
					
					
					
							    				
        						
			
			
		
		
		
		
		
			
	
	
		
	
		

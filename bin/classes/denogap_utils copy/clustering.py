#!/usr/bin/env python
from __future__ import division
import os
import re
import sys
import subprocess
from collections import defaultdict


class SequenceClustering:
	"""
	This class performs clustering operation for hmmer similarity pairs
	"""
	
	def mcl_clustering(self,mcl_abc_file,mcl_inflation):
	
		#### Output files generated during mcl clustering process ###
		mcl_dir_path=os.path.dirname(mcl_abc_file)
		mcxload_mci=os.path.join(mcl_dir_path,"seq.mci")
		mcxload_tab=os.path.join(mcl_dir_path,"seq.tab")
		mcl_out=os.path.join(mcl_dir_path,"mcl_cluster")
		
		### Run mcxload #####
		mcxload_cmd="mcxload -abc {} --stream-mirror --stream-neg-log10 "\
					"-stream-tf 'ceil(200)' -o {} -write-tab {}"\
					.format(mcl_abc_file,mcxload_mci,mcxload_tab)
					
		proc_mcxload=subprocess.Popen([mcxload_cmd],shell=True,stdout=subprocess.PIPE,
										stderr=subprocess.PIPE)
										
		proc_mcxload.wait()
		stdout,stderr=proc_mcxload.communicate()
    	
		if proc_mcxload.returncode!=0:
			print stderr
    		
		mcl_cmd="mcl {} -I {} -use-tab {} "\
				"-o {}".format(mcxload_mci,mcl_inflation,mcxload_tab,mcl_out)
    	
		proc_mcl=subprocess.Popen([mcl_cmd],shell=True,stdout=subprocess.PIPE,
									stderr=subprocess.PIPE)
		proc_mcl.wait()
		stdout,stderr=proc_mcl.communicate()
    	
		if proc_mcl.returncode!=0:
			print stderr
		
		return(mcl_out)

		
	def add_cluster_ids(self,mcl_cluster,**kwargs):
	
		"""
		This function adds group ids to mcl clusters 
		"""
		
		cluster_dir=os.path.dirname(mcl_cluster)
		
		cluster_file=os.path.join(cluster_dir,"cluster.txt")
		
		clusterid_counter=None
		
		with open(cluster_file,"w") as clust_group:
		
			with open(mcl_cluster,"r") as mc:
		
				if "start_at" in kwargs:
					clusterid_counter=kwargs["start_at"]
				#else:
				#	count_line="wc -l < {0}".format(mcl_cluster)
				#	proc=subprocess.Popen([count_line],shell=True,stdout=subprocess.PIPE,
				#						  stderr=subprocess.PIPE)	
				#	proc.wait()
				#	stdout,stderr=proc.communicate()
				#	clusterid_counter=int(stdout.strip())
				#	clusterid_counter+=1
			
				for cluster_line in mc:
				
					### check if group id assigned ###
					group_regx=re.compile(r"GROUP_\d+")
					cluster_id=None
					
					if re.search(group_regx,cluster_line):
						group_ids=re.findall(group_regx,cluster_line)
						if len(group_ids)>1:
							print "Error: multiple hmm-ids clustered togather into"\
							      "same cluster"
							sys.exit()
						else:
							cluster_id=group_ids[0]
					else:
						### if group id not assigned, assign new group id  and increment
						### counter for the group id
						cluster_id="GROUP_{0}".format(clusterid_counter)
						clusterid_counter+=1	      
 
					clust_group.write("{0}:\t{1}".format(cluster_id,cluster_line))
					
		clust_group.close()
		
		return(cluster_file)		
				
				
				
				
		
					
	    	
	                                     
	    
	
#!/usr/bin/env python
import os
import sys
from collections import defaultdict

class SetDir:
	"""
	This class is for setting-up output dirs for DeNoGAP
	"""
	
	def __init__(self):
		self.outdir_dict=defaultdict()
	
	### Function to setup project dir ####	
	def mk_denogap_hmm_dirs(self,project_dir):
		"""
		This function creates project directory if not exist.
		It also create missing sub-directories for storing DeNoGAP-HMM output
		Input: project dir path
		"""
		### If not already present create sub-directories for DENOGAP HMM ####		
		if not os.path.exists(os.path.join(project_dir,"DENOGAP_HMM")):
			os.makedirs(os.path.join(project_dir,"DENOGAP_HMM"))
		self.outdir_dict["BASE"]=os.path.join(project_dir,"DENOGAP_HMM")
		
		if not os.path.exists(os.path.join(project_dir,"TMP")):
			os.makedirs(os.path.join(project_dir,"TMP"))
		self.outdir_dict["TMP"]=os.path.join(project_dir,"TMP")
		
		if not os.path.exists(os.path.join(project_dir,"DENOGAP_HMM","HMM_MODEL")):
			os.makedirs(os.path.join(project_dir,"DENOGAP_HMM","HMM_MODEL"))
		self.outdir_dict["HMM_MODEL"]=os.path.join(project_dir,"DENOGAP_HMM","HMM_MODEL")
		
		if not os.path.exists(os.path.join(project_dir,"DENOGAP_HMM","HMMER_DB")):
			os.makedirs(os.path.join(project_dir,"DENOGAP_HMM","HMMER_DB"))
		self.outdir_dict["HMMER_DB"]=os.path.join(project_dir,"DENOGAP_HMM","HMMER_DB")
			
		if not os.path.exists(os.path.join(project_dir,"DENOGAP_HMM","HMMER_DB","SEQ_DB")):
			os.makedirs(os.path.join(project_dir,"DENOGAP_HMM","HMMER_DB","SEQ_DB"))
		self.outdir_dict["SEQ_DB"]=os.path.join(project_dir,"DENOGAP_HMM","HMMER_DB","SEQ_DB")
		
		if not os.path.exists(os.path.join(project_dir,"DENOGAP_HMM","HMMER_DB","HMM_DB")):
			os.makedirs(os.path.join(project_dir,"DENOGAP_HMM","HMMER_DB","HMM_DB"))
		self.outdir_dict["HMM_DB"]=os.path.join(project_dir,"DENOGAP_HMM","HMMER_DB","HMM_DB")
		
		if not os.path.exists(os.path.join(project_dir,"DENOGAP_HMM","ALIGNMENTS")):
			os.makedirs(os.path.join(project_dir,"DENOGAP_HMM","ALIGNMENTS"))
		self.outdir_dict["ALIGNMENTS"]=os.path.join(project_dir,"DENOGAP_HMM","ALIGNMENTS")
		
		if not os.path.exists(os.path.join(project_dir,"DENOGAP_HMM","ALIGNMENTS","ALL_MATCH")):
			os.makedirs(os.path.join(project_dir,"DENOGAP_HMM","ALIGNMENTS","ALL_MATCH"))
		self.outdir_dict["ALL_MATCH"]=os.path.join(project_dir,"DENOGAP_HMM","ALIGNMENTS","ALL_MATCH")
	
		if not os.path.exists(os.path.join(project_dir,"DENOGAP_HMM","ALIGNMENTS","BEST_MATCH")):
			os.makedirs(os.path.join(project_dir,"DENOGAP_HMM","ALIGNMENTS","BEST_MATCH"))
		self.outdir_dict["BEST_MATCH"]=os.path.join(project_dir,"DENOGAP_HMM","ALIGNMENTS","BEST_MATCH")
	
		if not os.path.exists(os.path.join(project_dir,"DENOGAP_HMM","ALIGNMENTS","PARTIAL_MATCH")):
			os.makedirs(os.path.join(project_dir,"DENOGAP_HMM","ALIGNMENTS","PARTIAL_MATCH"))
		self.outdir_dict["PARTIAL_MATCH"]=os.path.join(project_dir,"DENOGAP_HMM","ALIGNMENTS","PARTIAL_MATCH")
		
		if not os.path.exists(os.path.join(project_dir,"DENOGAP_HMM","ALIGNMENTS","CHIMERA_MATCH")):
			os.makedirs(os.path.join(project_dir,"DENOGAP_HMM","ALIGNMENTS","CHIMERA_MATCH"))
		self.outdir_dict["CHIMERIC_MATCH"]=os.path.join(project_dir,"DENOGAP_HMM","ALIGNMENTS","CHIMERA_MATCH")
		
		if not os.path.exists(os.path.join(project_dir,"DENOGAP_HMM","CLUSTERS")):
			os.makedirs(os.path.join(project_dir,"DENOGAP_HMM","CLUSTERS"))
		self.outdir_dict["CLUSTER"]=os.path.join(project_dir,"DENOGAP_HMM","CLUSTERS")
		
		
		return(self.outdir_dict)	
		
		
			
	
	
		
	
		
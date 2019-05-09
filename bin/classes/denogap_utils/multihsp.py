#!/usr/bin/env python
from __future__ import division
import os
import sys
import copy
import pandas as pd
from collections import defaultdict

class FixMultiHSP:
	"""
	This class is for finding best HSPs from the Hit
	"""
	
	def __init__(self,program):
		self.program=program
		
	def fix_hsps(self,result_df):
		"""
		This function iterates over hmmer result dict and prepares list
		of hsps for each query-target pair to be given to module for computing
		best hsps
		"""
		
		fixed_hit_df=pd.DataFrame()
		
		column_header=list(result_df)
		
		grouped_hsps=result_df.groupby(['QUERY_ID','HMM_ID'])
		
		multi_hsp_df=grouped_hsps.filter(lambda x: len(x) > 1)
		single_hsp_df=grouped_hsps.filter(lambda x: len(x)==1)
		
		#### add columns to single hsp_df ####
		
		single_hsp_df.loc[:,"TOTAL_PERCENT_IDENTITY"]=single_hsp_df["PERCENT_IDENTITY"]
		single_hsp_df.loc[:,"TOTAL_PERCENT_SIMILARITY"]=single_hsp_df["PERCENT_SIMILARITY"]
		single_hsp_df.loc[:,"TOTAL_QUERY_COVERAGE"]=single_hsp_df["QUERY_COVERAGE"]
		single_hsp_df.loc[:,"TOTAL_HMM_COVERAGE"]=single_hsp_df["QUERY_COVERAGE"]
		single_hsp_df.loc[:,"TOTAL_DOM_BITSCORE"]=single_hsp_df["DOM_BITSCORE"]
		single_hsp_df.loc[:,"TOTAL_IDENTICAL_COUNT"]=single_hsp_df["IDENTICAL_COUNT"]
		single_hsp_df.loc[:,"TOTAL_SIMILAR_COUNT"]=single_hsp_df["SIMILAR_COUNT"]
		single_hsp_df.loc[:,"TOTAL_QUERY_ALN_LEN"]=(single_hsp_df["QUERY_END"]-single_hsp_df["QUERY_START"])+1
		single_hsp_df.loc[:,"TOTAL_HMM_ALN_LEN"]=(single_hsp_df["HMM_END"]-single_hsp_df["HMM_START"])+1
		
		fixed_hit_df=fixed_hit_df.append(single_hsp_df,ignore_index=True)
		
		multi_hsp_group=multi_hsp_df.groupby(['QUERY_ID','HMM_ID'])
		
		for name,hsp_group in multi_hsp_group:
			
			hsps=hsp_group.to_dict('records')
		
			compatible_hsps=self.compute_best_hsp(hsps,hsp_list=list())
			compatible_df=pd.DataFrame(compatible_hsps)
				
			hit_df=self.convert_hsp2hit(compatible_df)	
				
			fixed_hit_df=fixed_hit_df.append(hit_df,ignore_index=True)	
		
		return(fixed_hit_df)		
		
		
		
	def compute_best_hsp(self,hsps,hsp_list=list()):
  		"""
  		Memoize the function so that, if we have already computed the solution
  		for `hsps`, just fetch and return that value immediately without
  		performing below computations.

		Based on code found here:
    	http://jeff.wintersinger.org/posts/2014/07/designing-an-algorithm-to-
    	compute-the-optimal-set-of-blast-hits/
  		"""
  		
  		if len(hsps) == 0:
  			return hsp_list
  			
  		if len(hsps) == 1:
  			# Trivial solution: one HSP, so optimal solution is just itself.
  			hsp_list.append(hsps[0])
  			return hsp_list

  		# Last HSP
  		last_hsp = hsps[-1]
  		# All HSPs except last
  		previous = hsps[:-1]
  		
  		# Find subset of HSPs in `previous` that don't overlap `last_hsp`.
  		compatible = self.find_compatible(last_hsp, previous)
  		
  		without_list = copy.deepcopy(hsp_list)
  		with_list = copy.deepcopy(hsp_list)
  		with_list.append(last_hsp)
  		
  		best_without_last = self.compute_best_hsp(previous,without_list)
  		best_with_last    = self.compute_best_hsp(compatible, with_list)
  		
  		compatible_hsp=max([best_without_last, best_with_last],key=lambda x: sum([float(y['DOM_BITSCORE']) for y in x]))
  		
  		return(compatible_hsp)
		
		
	def find_compatible(self,target, hsps):
  		'''Find hsps amongst `hsps` that don't overlap or cross `target`. They
  		may not be mutually compatible, as they are only guaranteed to be
  		compatible with `target`.'''
  		compatible = []
  		
  		for hsp in hsps:
  			if hsp == target:
  				continue
  				
  			if target['QUERY_START'] <= hsp['QUERY_START']:
  				first, second = target, hsp
  			else:
  				first, second = hsp, target
  			
  			overlap = (second['QUERY_START'] <= first['QUERY_END'] or second['HMM_START'] <= first['HMM_END'])
  			
  			if not overlap:
  				compatible.append(hsp)
  				
  		return(compatible)


	def convert_hsp2hit(self,compatible_df):
		"""
		This function converts list of hsp into hit
		"""
		
		hit_df=compatible_df
		
		query_len=compatible_df["QUERY_LEN"].max()
		hmm_len=compatible_df["HMM_LEN"].max()
		
		total_dom_bitscore=compatible_df["DOM_BITSCORE"].sum()
		total_num_ident=compatible_df["IDENTICAL_COUNT"].sum()
		total_num_sim=compatible_df["SIMILAR_COUNT"].sum()
		total_query_aln=((compatible_df["QUERY_END"] - compatible_df["QUERY_START"])+1).sum()
		total_hmm_aln=((compatible_df["HMM_END"] - compatible_df["HMM_START"])+1).sum()
		
		hit_df["TOTAL_PERCENT_IDENTITY"]=round(float(int(total_num_ident)/int(min(total_query_aln,total_hmm_aln)))*100,2)
			
		hit_df["TOTAL_PERCENT_SIMILARITY"]=round(float(int(total_num_sim)/int(min(total_query_aln,total_hmm_aln)))*100,2)
		
		hit_df["TOTAL_QUERY_COVERAGE"]=round(float(int(total_query_aln) / int(query_len))*100,2)
			
		hit_df["TOTAL_HMM_COVERAGE"]=round(float(int(total_hmm_aln) / int(hmm_len))*100,2)
		
		hit_df["TOTAL_DOM_BITSCORE"]=total_dom_bitscore
		hit_df["TOTAL_IDENTICAL_COUNT"]=total_num_ident
		hit_df["TOTAL_SIMILAR_COUNT"]=total_num_sim
		hit_df["TOTAL_QUERY_ALN_LEN"]=total_query_aln
		hit_df["TOTAL_HMM_ALN_LEN"]=total_hmm_aln
	
		return(hit_df)							 					
		
		
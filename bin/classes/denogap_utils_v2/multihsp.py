#!/usr/bin/env python
from __future__ import division
import os
import sys
import copy
from collections import defaultdict

class FixMultiHSP:
	"""
	This class is for finding best HSPs from the Hit
	"""
	
	def __init__(self,program):
		self.program=program
		
	def fix_hsps(self,result_dict,accuracy):
		"""
		This function iterates over hmmer result dict and prepares list
		of hsps for each query-target pair to be given to module for computing
		best hsps
		"""
		
		fixed_hsp_dict=defaultdict(dict)
		
		for query_id in result_dict:
			for target_id in result_dict[query_id]:
				
				hsps=[domain_line for domain_index,domain_line in result_dict[query_id][target_id].iteritems()
					  if domain_line['AVERAGE_ACCURACY']>=accuracy]
				
				if len(hsps)!=0:	
					compatible_hsps=self.compute_best_hsp(hsps,hsp_list=list())
					fixed_hsp_dict[query_id][target_id]=compatible_hsps
		
		hit_dict=self.convert_hsp2hit(fixed_hsp_dict)
		
		#print hit_dict		
		
		return(hit_dict)		
		
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

		compatible_hsp=max([best_without_last, best_with_last],key=lambda (x): sum([float(y['BITSCORE']) for y in x]))
		
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


	def convert_hsp2hit(self,fixed_hsp_dict):
		"""
		This function converts list of hsp into hit
		"""
	
		hit_dict=defaultdict(dict)
		
		for query_id in fixed_hsp_dict:
		
			for target_id in fixed_hsp_dict[query_id]:

				hsp_dict=defaultdict(dict)
				
				hsp_list=fixed_hsp_dict[query_id][target_id]
				
				total_pident=hsp_list[0]['PERCENT_IDENTITY']
				total_psim=hsp_list[0]['PERCENT_SIMILARITY']
				total_bitscore=hsp_list[0]['BITSCORE']
				total_qcov=hsp_list[0]['QUERY_COVERAGE']
				total_hcov=hsp_list[0]['HMM_COVERAGE']
				total_num_ident=hsp_list[0]['IDENTICAL_COUNT']
				total_num_sim=hsp_list[0]['SIMILAR_COUNT']
				query_len=hsp_list[0]['QUERY_LEN']
				hmm_len=hsp_list[0]['HMM_LEN']
				num_hsp=len(hsp_list)
				
				if num_hsp>1:
					total_num_ident=sum(int(k['IDENTICAL_COUNT']) for k in hsp_list)
					total_num_sim=sum(int(k['SIMILAR_COUNT']) for k in hsp_list)
					sum_query_aln=sum(int((k['QUERY_END']-k['QUERY_START'])+1) for k in hsp_list)
					sum_hmm_aln=sum(int((k['HMM_END']-k['HMM_START'])+1) for k in hsp_list)					
					total_bitscore=sum(float(k['BITSCORE']) for k in hsp_list)
					total_pident=round(float(int(total_num_ident)/int(min(sum_query_aln,sum_hmm_aln)))*100,2)
					total_psim=round(float(int(total_num_sim)/int(min(sum_query_aln,sum_hmm_aln)))*100,2)
					total_qcov=round(float(int(sum_query_aln) / int(query_len))*100,2)
					total_hcov=round(float(int(sum_hmm_aln) / int(hmm_len))*100,2)
				
				for hsp in hsp_list:
					hsp_dict[hsp['DOMAIN_INDEX']]=hsp
					hsp_dict[hsp['DOMAIN_INDEX']].update({'TOTAL_IDENTICAL_COUNT':total_num_ident,
									 'TOTAL_SIMILAR_COUNT':total_num_sim,
									 'TOTAL_BITSCORE':total_bitscore,
									 'TOTAL_PERCENT_IDENTITY':total_pident,
									 'TOTAL_PERCENT_SIMILARITY':total_psim,
									 'TOTAL_QUERY_COVERAGE':total_qcov,
									 'TOTAL_HMM_COVERAGE':total_hcov})
				
				hit_dict[query_id][target_id]=hsp_dict
				
		return(hit_dict)							 					
		
		
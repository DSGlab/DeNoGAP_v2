#!/usr/bin/env python
import os
import sys
import subprocess
from Bio import SeqIO
from Bio import SearchIO
from collections import defaultdict

class SequenceAlignment:
	"""
	This class comprises of functions for sequence alignment
	"""
	
	def pairwise_alignment(self,**kwargs):
		"""
		This function performs pairwise hmmer alignment between sequence file
		and sequence database
		"""
		
		

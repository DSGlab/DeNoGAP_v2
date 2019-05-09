#!/usr/bin/env python
import os
import sys
import re
import time
import argparse
import pandas as pd
from Bio import SeqIO
from Bio import SearchIO
from collections import defaultdict
from argparse import RawTextHelpFormatter
from classes.denogap_utils.multihsp import FixMultiHSP

#### CMD ARGS ######
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
    parser._optionals.title="Argument for HMMER Result Parser"

    parser.add_argument("--hmmer_program",metavar="<STR>",
                        choices=["phmmer","hmmscan","hmmsearch"],
                        help="Name of the hmmer program",required=True)
    parser.add_argument("--hmmer_dom",metavar="<STR>",
                        help="Dom-tab result file for HMMER",required=True)
    parser.add_argument("--hmmer_aln",metavar="STR>",
                        help="Alignment result file for HMMER",required=True)
    parser.add_argument("--hmmer_result_dir",metavar="<STR>",
                        help="Directory path for storing parsed hmmer results",required=True)                    
    parser.add_argument("--accuracy_thresh", metavar="<FLOAT>",default=0.7,
									choices=range(0, 1),type=float,
									help="parse hits >= this accuracy threshold as "\
										"significant hit [default: 0.7] (x>0)")
    parser.add_argument("--ident_thresh", metavar="<FLOAT>",default=70.0,
									choices=range(0, 100),type=float,
									help="parse hits with percent sequence match between "\
										"query and target >= this threshold as "\
										"signficant hit [default: 70.0] (x>0)")
    parser.add_argument("--hit_cov", metavar="<FLOAT>",default=70.0,
									choices=range(0, 100),type=float,
									help="parse hits >= this query/target coverage "\
										"threshold as significant hit "\
										"[default: 70.0] (x>0)")
    parser.add_argument("--chimera_cov", metavar="<FLOAT>",default=20.0,
									choices=range(0, 100),type=float,
									help="minimum coverage threshold for query/target "\
									"to be reported as chimera alignment "\
									"[default: 20.0] (x>0)")
    
    args=parser.parse_args()

    return vars(args)
    

### parsing function for dom-tab file ####
def parse_hmmer_domtab(hmmer_program,domtab_file):
    """
    This function parse result from hmmer domtab alignment output
    """
    ### define dictionary to store parsed domtab output
    domtab_list=list()
		
    with open(domtab_file,"r") as domtab_result:
        seq_records=SearchIO.parse(domtab_result,"{}3-domtab".format(hmmer_program))
			
        for QueryResult in seq_records:	
            for Hit in QueryResult:
                for HSP in Hit:
                    domain_dict=defaultdict()
                    domain_dict["QUERY_ID"]=QueryResult.id
                    domain_dict["QUERY_LEN"]=QueryResult.seq_len
                    domain_dict["HMM_ID"]=Hit.id
                    domain_dict["HMM_LEN"]=Hit.seq_len
                    domain_dict["DOMAIN_INDEX"]=HSP.domain_index
                    domain_dict["SEQ_EVALUE"]=Hit.evalue
                    domain_dict["SEQ_BITSCORE"]=Hit.bitscore
                    domain_dict["DOM_EVALUE"]=HSP.evalue
                    domain_dict["DOM_BITSCORE"]=HSP.bitscore
                    domain_dict["AVERAGE_ACCURACY"]=HSP.acc_avg
                    domain_dict["HMM_END"]=HSP.hit_end
                    domain_dict["HMM_START"]=HSP.hit_start
                    domain_dict["QUERY_END"]=HSP.query_end
                    domain_dict["QUERY_START"]=HSP.query_start
                    
                    if HSP.query_start==0:
                        domain_dict["QUERY_START"]=1
                    if HSP.hit_start==0:
                        domain_dict["HMM_START"]=1
  
                    domtab_list.append(domain_dict)
		
    dom_df=pd.DataFrame(domtab_list)			

    return(dom_df)				
		
### parsing function for hmm alignment file ####
def parse_hmmer_similarity(hmmer_program,alignment_file):
    """
    This function parse result from hmmer alignment and 
    fetch similarity/identity values
    """	
    
    align_list=list()
    
    with open(alignment_file,"r") as alignment_result:
        flag_block=0
        flag_sim=0
        query_id=""
        target_id=""
        domain_index=""
        first_regx=""
        second_regx=""
        dom_sim=defaultdict()
        
        for line in alignment_result:
            if line.startswith("Query:"):
                query_row=re.split(" +",line)
                query_id=query_row[1]
                dom_sim=defaultdict(dict)
                continue
            elif line.startswith(">>"):
                if len(dom_sim)!=0:
                    align_list.append(dom_sim)
                target_row=re.split(" +",line)
                target_id=target_row[1]
                dom_sim=defaultdict(dict)
                continue
            elif line.startswith("  == domain"):
                domain_row=re.split(" +",line)
                domain_index=int(domain_row[3])
                if len(dom_sim)!=0:
                    align_list.append(dom_sim)
                dom_sim=defaultdict(dict)
                dom_sim["QUERY_ID"]=query_id
                dom_sim["HMM_ID"]=target_id
                dom_sim["DOMAIN_INDEX"]=domain_index
                flag_block=1
                continue
            elif line.startswith("Internal pipeline"):
                align_list.append(dom_sim)
            else:
                if hmmer_program=="phmmer" or hmmer_program=="hmmsearch":
                    first_regx=r"^(\s+)"+ re.escape(query_id)
                    second_regx=r"^(\s+)"+ re.escape(target_id)
                elif hmmer_program=="hmmscan":
                    first_regx=r"^(\s+)"+ re.escape(target_id)
                    second_regx=r"^(\s+)"+ re.escape(query_id)
                    
                if flag_block==1 and re.match(first_regx,line):
                    flag_sim=1
                elif flag_sim==1:
                    if not "IDENTICAL_COUNT" in dom_sim:
                        sim_line=line.replace(" ","").rstrip("\n")
                        identical_count=len(sim_line)-sim_line.count("+")
                        similar_count=len(sim_line)
                        dom_sim["IDENTICAL_COUNT"]=identical_count
                        dom_sim["SIMILAR_COUNT"]=similar_count
                    else:
                        sim_line=line.replace(" ","").rstrip("\n")
                        identical_count=len(sim_line)-sim_line.count("+")
                        similar_count=len(sim_line)
                        dom_sim["IDENTICAL_COUNT"]+=identical_count
                        dom_sim["SIMILAR_COUNT"]+=similar_count
                    flag_sim=0
    
    sim_df=pd.DataFrame(align_list)
    
    return(sim_df)


### parsing function for add hmmer statistics ####	
def add_hmmer_stats(merge_df1):
    """
    This function adds stats to hmmer result after calculation of
    percent identity, similarity and sequence coverage
    """
    new_column_dict=defaultdict(list)
    merge_df2=merge_df1
    
    for index, rows in merge_df1.iterrows():
                
                #### get similar/identical residue counts ##
                identical_count=rows["IDENTICAL_COUNT"]
                similar_count=rows["SIMILAR_COUNT"]
                
                ### get coordinates and length for query and target
                qlen=rows["QUERY_LEN"]
                hlen=rows["HMM_LEN"]
                qstart=rows["QUERY_START"]
                qend=rows["QUERY_END"]
                hstart=rows["HMM_START"]
                hend=rows["HMM_END"]
                avg_acc=rows["AVERAGE_ACCURACY"]
                
                qaln_len=((qend-qstart)+1)
                haln_len=((hend-hstart)+1)
                
                ### use minimum of query or target to find percent identity/similarity
                min_len=min(qaln_len,haln_len)
                percent_identity=round(float(identical_count/min_len)*100,2)
                percent_similarity=round(float(similar_count/min_len)*100,2)
                
                ### get query/target sequence coverage in the alignment ###
                qcoverage=round(float(qaln_len/qlen)*100,2)
                hcoverage=round(float(haln_len/hlen)*100,2)
                
                new_column_dict["PERCENT_IDENTITY"].append(percent_identity)
                new_column_dict["PERCENT_SIMILARITY"].append(percent_similarity)
                new_column_dict["QUERY_COVERAGE"].append(qcoverage)
                new_column_dict["HMM_COVERAGE"].append(hcoverage)
       
    merge_df2["PERCENT_IDENTITY"]=new_column_dict["PERCENT_IDENTITY"]
    merge_df2["PERCENT_SIMILARITY"]=new_column_dict["PERCENT_SIMILARITY"]
    merge_df2["QUERY_COVERAGE"]=new_column_dict["QUERY_COVERAGE"]
    merge_df2["HMM_COVERAGE"]=new_column_dict["HMM_COVERAGE"]   					
				
    return(merge_df2)	  

### parsing function for filtering hits ####
def filter_hits(filter_param,hmmer_result_dict):
    """
    This function filter hmmer result based on thresholds defined by user on
    command line
    """
    
    filter_results=defaultdict(list)
    
    for query_id in hmmer_result_dict:
        for target_id in hmmer_result_dict[query_id]:
            for domain_index in hmmer_result_dict[query_id][target_id]:
                hmmdom_dict=hmmer_result_dict[query_id][target_id][domain_index]
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
                        
                    elif any(param>= float(filter_param["hitcov_cutoff"])
                        for param in feature_list2):
                        
                        filter_results["PARTIAL"].append(hmmdom_dict)
                        
                    elif all(param>= float(filter_param["chicov_cutoff"]) and 
                        param <= float(filter_param["hitcov_cutoff"])
					    for param in feature_list2):

                        if qstart<=15 and qend<qlen:
                            if hstart<=15 and hend<hlen:
                                filter_results["CHIMERA"].append(hmmdom_dict)
					            
                        elif qstart>1 and qend>=(qlen-15):
                            if hstart>1 and hend>=(hlen-15):
                                filter_results["CHIMERA"].append(hmmdom_dict)
					            
    return(filter_results)	
    

##### MAIN ####
def main():

    """
    Main Function
    """
    start=time.time()
    c_args = parse_args(__file__)
    
    hmmer_program=c_args["hmmer_program"]
    hmmer_dom_file=os.path.abspath(c_args["hmmer_dom"])
    hmmer_aln_file=os.path.abspath(c_args["hmmer_aln"])
    hmmer_result_dir=os.path.abspath(c_args["hmmer_result_dir"])
   
    ### parse domain-tab file ###
    dom_df=parse_hmmer_domtab(hmmer_program,hmmer_dom_file)
    ### parse full-alignment file ###
    sim_df=parse_hmmer_similarity(hmmer_program,hmmer_aln_file)
    
    ### merge dom and sim dataframes ####
    merge_df1=pd.merge(dom_df,sim_df,on=["QUERY_ID","HMM_ID","DOMAIN_INDEX"])
    
    ### additional stats to dataframe
    merge_df2=add_hmmer_stats(merge_df1)
    
    merge_df2_filtered=merge_df2[(merge_df2.PERCENT_IDENTITY>=c_args["ident_thresh"]) &
                                 (merge_df2.AVERAGE_ACCURACY>=c_args["accuracy_thresh"])]
    
    fixed_hmmer_result_df=FixMultiHSP(hmmer_program).fix_hsps(merge_df2_filtered)
    
    fixed_hmmer_result_df.to_csv(os.path.join(hmmer_result_dir,"test_fixed_hmmer.csv"),sep="\t")
    
    end=time.time()
    print(end-start)
	
if __name__ == '__main__':
	main()  
	

    							
	
	
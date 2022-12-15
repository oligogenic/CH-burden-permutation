import argparse

import time
import numpy as np
from scipy import stats
import random
print("------------------------------------------")
print("Starting...")
start_time = time.time()



# Get commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument('-v','--variants', help='Location of variants file', required=True)
parser.add_argument('-g','--genes', help="Location of genes file", required=True)
parser.add_argument('-n','--number_of_runs', help="# runs", required=True)
parser.add_argument('-ng','--number_of_genes', help="# genes in subset", required=True)
parser.add_argument('-p','--patients', help='Location of patients file')

args = vars(parser.parse_args())


variants_fl = args['variants']
genes_fl = args['genes']

a_genes =  set([line.strip() for line in open(genes_fl)])

patients_fl = args['patients']


exl_p = set()
if patients_fl is not None:
    f = open(patients_fl)
    for l in f:
        exl_p.add(l.strip())

number_of_runs = int(args['number_of_runs'])

number_of_genes = int(args['number_of_genes'])
#['patient', 'pathology', 'gene_symbol', 'hgvs_protein', 'chr', 'pos', 'exon_intron_rank', 'exon_intron_total', 'cdna_pos', 'cdna_length', 'protein_pos', 'protein_length', 'reference', 'alternative', 'filters', 'confidence', 'variant_confidence_by_depth', 'read_depth', 'allelic_depth_proportion_alt', 'allelic_depth_proportion_ref', 'zygosity', 'genotype_likelihood_het', 'genotype_likelihood_hom_alt', 'genotype_likelihood_hom_ref', 'detected_as_error', 'change_num_germline_filt', 'change_num_germline_filt_num_groups', 'change_num_germline_filt_details_groups', 'short_tandem_repeat', 'the1000G_AF', 'the1000G_AC', 'gonl_af', 'gonl_ac', 'exac_af', 'exac_ac', 'ESP6500_EA_AF', 'lof_tolerant_or_recessive_gene', 'biotype', 'snpeff_effect', 'snpeff_impact', 'is_scSNV_Ensembl', 'is_scSNV_RefSeq', 'splicing_ada_pred', 'splicing_rf_pred', 'aggregation_pred_lr', 'aggregation_pred_radial_svm', 'reliability_index', 'consensus_prediction', 'cadd_phred', 'vest_score', 'num_genes', 'check_segregation', 'of_interest_gene', 'private_comments_gene', 'public_comments_gene', 'of_interest_variant', 'private_comments_variant', 'public_comments_variant', 'check_insilico', 'check_validated_change', 'check_segregation_username', 'evaluation', 'evaluation_comments', 'hgvs_dna', 'transcript_uniprot_id', 'transcript_refseq_mrna\n']


control_paths= ["MEDINT","TOL","FERTILITY","TOL_FULL","THYROID","NEPHRO","CTRL","CAR_ARR"]
case_paths = ["HYDROCEPH"]




for i in range(number_of_runs):
    random.seed(random.random())
    genes = set(random.sample(list(a_genes), k=number_of_genes))

    variants_f = open(variants_fl)
    
    variants_f.readline() # header
    
    n_vars = 0
    cases = {}
    controls = {}
    
    
    for l in variants_f:
        n_vars+=1
        l = l.strip().split("\t")
        pt = l[0]
        if pt in exl_p:
            continue
        pa = l[1]
        zyg = l[20]
        zyg_v = None
        if zyg == 'Heterozygous':
            zyg_v = 1
        elif zyg == 'Homozygous':
            zyg_v = 2    
        g = l[2]
        
    
        gt = "NH"     # check if hydro gene
        if g in genes:
            gt = "H"
        
    
    
            if pa in control_paths:
                if pt in controls:
                    controls[pt] += zyg_v
                else:
                    controls[pt] = zyg_v       
            elif pa in case_paths:
                if pt in cases:
                    cases[pt] += zyg_v
                else:
                    cases[pt] = zyg_v       
            else:
             raise Exception('Not case or control pathology')    
    
    
    cases_v = list(cases.values())
    controls_v = list(controls.values())
    
    cases_m = np.mean(cases_v)
    cases_std = np.std(cases_v)
    cases_med = np.median(cases_v)
    
    controls_m = np.mean(controls_v)
    controls_std = np.std(controls_v)
    controls_med = np.median(controls_v)
    
    wr = stats.ranksums(cases_v,controls_v)

    #if len(cases)==0 or len(controls)==0:
    #    continue
    
    print("# genes: " + str(len(genes)))
    print("# vars: " + str(n_vars))
    print("\t\t cases \t | \tcontrols")
    print("# patients \t "+str(len(cases)) + "\t | \t"+str(len(controls)) )
    print("# vars \t\t "+str(sum(cases_v)) + "\t | \t"+str(sum(controls_v)) )
    print("# vars/pats \t "+str(sum(cases_v)/len(cases)) + "\t | \t"+str(sum(controls_v)/len(controls)))
    print("# mean vars \t "+str(round(cases_m, 2)) + "\t | \t"+str(round(controls_m, 2)) )
    print("# stdev  \t "+str(round(cases_std,2)) + "\t | \t"+str(round(controls_std,2)) )
    print("# median  \t "+str(cases_med) + "\t | \t"+str(controls_med) )
    print()
    print("wilcox: " + str(wr))
    
    print()
    print()
    print("Analysis completed. Ran %s seconds " % (round(time.time() - start_time, 1)))
    print("------------------------------------------")
    
    
    
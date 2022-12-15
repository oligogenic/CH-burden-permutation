import argparse

import time
import numpy as np
from scipy import stats
import random
import sys
import seaborn as sns

print("------------------------------------------")
print("Starting...", end = " ")
print(" ".join(sys.argv))
start_time = time.time()



# Get commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument('-p','--patients', help='Location of patients file')
parser.add_argument('-v','--variants', help='Location of variants file', required=True)




args = vars(parser.parse_args())


variants_fl = args['variants']
patients_fl = args['patients']



exl_p = set()
if patients_fl is not None:
    f = open(patients_fl)
    for l in f:
        exl_p.add(l.strip())




#patient pathology   gene_symbol hgvs_protein    chr pos exon_intron_rank    exon_intron_total   cdna_pos    cdna_length protein_pos protein_length  reference   alternative filters confidence  variant_confidence_by_depth read_depth  allelic_depth_proportion_alt    allelic_depth_proportion_ref    zygosity    genotype_likelihood_het genotype_likelihood_hom_alt genotype_likelihood_hom_ref detected_as_error   change_num_germline_filt    change_num_germline_filt_num_groups change_num_germline_filt_details_groups short_tandem_repeat the1000G_AF the1000G_AC gonl_af gonl_ac exac_af exac_ac ESP6500_EA_AF   lof_tolerant_or_recessive_gene  biotype snpeff_effect   snpeff_impact   is_scSNV_Ensembl    is_scSNV_RefSeq splicing_ada_pred   splicing_rf_pred    aggregation_pred_lr aggregation_pred_radial_svm reliability_index   consensus_prediction    cadd_phred  vest_score  num_genes   check_segregation   of_interest_gene    private_comments_gene   public_comments_gene    of_interest_variant private_comments_variant    public_comments_variant check_insilico  check_validated_change  check_segregation_username  evaluation  evaluation_comments hgvs_dna
#['patient', 'pathology', 'gene_symbol', 'hgvs_protein', 'chr', 'pos', 'exon_intron_rank', 'exon_intron_total', 'cdna_pos', 'cdna_length', 'protein_pos', 'protein_length', 'reference', 'alternative', 'filters', 'confidence', 'variant_confidence_by_depth', 'read_depth', 'allelic_depth_proportion_alt', 'allelic_depth_proportion_ref', 'zygosity', 'genotype_likelihood_het', 'genotype_likelihood_hom_alt', 'genotype_likelihood_hom_ref', 'detected_as_error', 'change_num_germline_filt', 'change_num_germline_filt_num_groups', 'change_num_germline_filt_details_groups', 'short_tandem_repeat', 'the1000G_AF', 'the1000G_AC', 'gonl_af', 'gonl_ac', 'exac_af', 'exac_ac', 'ESP6500_EA_AF', 'lof_tolerant_or_recessive_gene', 'biotype', 'snpeff_effect', 'snpeff_impact', 'is_scSNV_Ensembl', 'is_scSNV_RefSeq', 'splicing_ada_pred', 'splicing_rf_pred', 'aggregation_pred_lr', 'aggregation_pred_radial_svm', 'reliability_index', 'consensus_prediction', 'cadd_phred', 'vest_score', 'num_genes', 'check_segregation', 'of_interest_gene', 'private_comments_gene', 'public_comments_gene', 'of_interest_variant', 'private_comments_variant', 'public_comments_variant', 'check_insilico', 'check_validated_change', 'check_segregation_username', 'evaluation', 'evaluation_comments', 'hgvs_dna', 'transcript_uniprot_id', 'transcript_refseq_mrna\n']


control_paths= ["MEDINT","TOL","FERTILITY","TOL_FULL","THYROID","NEPHRO","CTRL","CAR_ARR"]
case_paths = ["HYDROCEPH"]



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

print()
print()
print()
print("cases counts")
for e in cases:
    print(e, end = ":")
    print(cases[e], end=", ")

print()
print()
print("controls counts")

for e in controls:
    print(e, end = ":")
    print(controls[e], end=", ")

print()
print()
print()


xs = []

for _ in cases_v:
    xs.append(1)

for _ in controls_v:
    xs.append(5)

cc = cases_v + controls_v

cases_m = np.mean(cases_v)
controls_m = np.mean(controls_v)

pi=sns.stripplot(xs,cc, jitter=0.2)
pi.axhline(cases_m, color= "blue")
pi.axhline(controls_m, color= "green")

p = pi.get_figure()
fn = variants_fl+variants_fl.replace("/","_").split(".")[0]+".png"
print(fn)
p.savefig(fn)




print("Analysis completed. Ran %s seconds " % (round(time.time() - start_time, 1)))
print("------------------------------------------")
time.sleep(1)    
    
    
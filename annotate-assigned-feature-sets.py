import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-out", type=str, required=True,
                    help="out file name")
parser.add_argument("-mrna", type=str, required=True,
                    help="file with peaks assigned to mRNA features")
parser.add_argument("-exon", type=str, required=True,
                    help="file with peaks assigned to exon features")
parser.add_argument("-cds", type=str, required=True,
                    help="file with peaks assigned to cds features")
parser.add_argument("-rel_cds_limits", type=str, required=True,
                    help="file mapping mRNA feature IDs to relative CDS start and end")

args = parser.parse_args()

assigned_mrna_file = args.mrna
assigned_exon_file = args.exon
assigned_cds_file = args.cds
relative_cds_limits_file = args.rel_cds_limits
outfile_path = args.out

##### function definitions

# given a text key and a gff description string,
# return a the corresponding value
def GetGffDescValue(gene_desc, key):
    for tag in gene_desc.strip().split(';'):
        k, v = tag.split('=', maxsplit=1)
        if k == key:
            return v

# Applied to a pandas dataframe row by row, returns a series for each row
# with the summit distance to the feature's true (stranded) start and true
# (stranded) end as a signed int. Positive dist-to-start and negative
# dist-to-end implying the summmit falls within the feature
def CalcSummitDistToFeature(s):
    # s is a series from a pandas dataframe of assigned gene (or other) features
    s['abs_peak_summit'] = s['peak_summit'] + s['peak_start']
    feature_strand = s['feature_strand']
    summit_loc = s['abs_peak_summit']
    if feature_strand == '+':
        stranded_feature_start = s['feature_start']
        stranded_feature_end = s['feature_end']
        summit_to_transcript_start = summit_loc - stranded_feature_start
        summit_to_transcript_end = summit_loc - stranded_feature_end
    elif feature_strand == '-':
        stranded_feature_start = s['feature_end']
        stranded_feature_end = s['feature_start']
        summit_to_transcript_start = stranded_feature_start - summit_loc
        summit_to_transcript_end = stranded_feature_end - summit_loc
        
    # returns the same series plus added fields
    s['summit_to_transcript_start'] = summit_to_transcript_start
    s['summit_to_transcript_end'] = summit_to_transcript_end
    return s

# Filters to remove any peaks where the summit is not overlapping the features
def FilterForOverlappingSummits(df):
    df['abs_peak_summit'] = df['peak_summit'] + df['peak_start']
    df = df[df['abs_peak_summit'] >= df['feature_start']]
    df = df[df['abs_peak_summit'] <= df['feature_end']]
    return df

# For each peak, annotate the type of region where the summit is located
def AnnotatePeakRegion(s):
    if s['summit_to_transcript_start'] >= 0 and s['summit_to_transcript_end'] <= 0:
        # this is within the mRNA feature
        if s['within_exon'] and s['within_cds']:
            region = 'cds'
        elif s['within_exon'] and not s['within_cds']:
            # this is within one of the UTRs
            if abs(s['summit_to_transcript_start']) < abs(s['summit_to_transcript_end']):
                # closer to the start, so assigned 5prime utr
                region = 'utr5prime'
            else:
                # otherwise must be 3prime UTR
                region = 'utr3prime'
        else:
            # if in gene but not in a UTR or exon, must be in intron
            region = 'intron'
    elif s['summit_to_transcript_start'] < 0:
        region = 'upstream'
    elif s['summit_to_transcript_end'] > 0:
        region = 'downstream'
    else:
        region = 'undefined'
    return region

# Assign a boolean based on whether peak that was assigned to multiple
# features is valid. Features marked False are those that were annotated as
# 'upstream' or 'downstream', but also fall within another mRNA feature
def FilterMultiAnnotated(peak_group):
    
    peak_group['pass_multi_assign_filter'] = True
    
    has_gene_overlap = False
    overlap_options = ['cds', 'utr5prime', 'utr3prime', 'intron']
    for o in overlap_options:
        if o in peak_group['annotated_peak_region'].values:
            has_gene_overlap = True
            break

    if has_gene_overlap:
        peak_group.loc[peak_group['annotated_peak_region'].isin(['upstream', 'downstream']), 'pass_multi_assign_filter'] = False
        
    peak_group = peak_group[peak_group['pass_multi_assign_filter']]
    
    return peak_group

# Takes in a row of the assigned_mrna dataframe and
# the entire overlapping_exon dataframe
# and marks whether the given peak is contained within an exon
# of the the corresponding assigned mRNA (ignoring possible other overlapping genes)
def MarkExonicFeatures(assigned_mrna_series, overlapping_exon):
    mrna_id = assigned_mrna_series['feature_id']
    peak_name = assigned_mrna_series['peak_name']
    relevant_exons = overlapping_exon[overlapping_exon['feature_parent'] == mrna_id]
    if peak_name in relevant_exons['peak_name'].values:
        within_exon = True
    else:
        within_exon = False
    assigned_mrna_series['within_exon'] = within_exon
    return assigned_mrna_series

# Takes in a row of the assigned_mrna dataframe and
# the entire overlapping_exon dataframe
# and marks whether the given peak is contained within an exon
# of the the corresponding assigned mRNA (ignoring possible other overlapping genes)
def MarkCdsFeatures(assigned_mrna_series, overlapping_cds):
    mrna_id = assigned_mrna_series['feature_id']
    peak_name = assigned_mrna_series['peak_name']
    relevant_cds = overlapping_cds[overlapping_cds['feature_parent'] == mrna_id]
    if peak_name in relevant_cds['peak_name'].values:
        within_cds = True
    else:
        within_cds = False
    assigned_mrna_series['within_cds'] = within_cds
    return assigned_mrna_series

# Main wrapper function
def AnnotateAssignedGenes(assigned_mrna, assigned_exon, assigned_cds):

    assigned_mrna = assigned_mrna.apply(CalcSummitDistToFeature, axis='columns')
    assigned_mrna['feature_id'] = assigned_mrna['feature_desc'].apply(lambda x: GetGffDescValue(x, 'ID'))
    
    overlapping_exon = FilterForOverlappingSummits(assigned_exon)
    overlapping_cds = FilterForOverlappingSummits(assigned_cds)
    
    overlapping_exon['feature_parent'] = overlapping_exon['feature_desc'].apply(lambda x: GetGffDescValue(x, 'Parent'))
    overlapping_cds['feature_parent'] = overlapping_cds['feature_desc'].apply(lambda x: GetGffDescValue(x, 'Parent'))
    
    # assigned_mrna['within_exon'] = assigned_mrna['peak_name'].isin(overlapping_exon['peak_name'])
    # assigned_mrna['within_cds'] = assigned_mrna['peak_name'].isin(overlapping_cds['peak_name'])

    assigned_mrna = assigned_mrna.apply(lambda x: MarkExonicFeatures(x, overlapping_exon), axis='columns')
    assigned_mrna = assigned_mrna.apply(lambda x: MarkCdsFeatures(x, overlapping_cds), axis='columns')
    
    assigned_mrna['annotated_peak_region'] = assigned_mrna.apply(AnnotatePeakRegion, axis='columns')
    
    assigned_mrna = assigned_mrna.groupby('peak_name').apply(FilterMultiAnnotated)
    
    return assigned_mrna

assigned_features_cols = ['peak_chr', 'peak_start', 'peak_end', 'peak_name', 'peak_display', 'unused', 'peak_foldch', 'peak_pscore',
                          'peak_qscore', 'peak_summit', 'feature_chr', 'feature_source', 'feature_type', 'feature_start', 'feature_end', 'unused2',
                          'feature_strand', 'feature_phase', 'feature_desc']

assigned_mrna = pd.read_csv(assigned_mrna_file, sep='\t', header=None, names=assigned_features_cols)
print(f'Read assigned mRNA feautures: {assigned_mrna.shape[0]}') 
assigned_exon = pd.read_csv(assigned_exon_file, sep='\t', header=None, names=assigned_features_cols)
print(f'Read assigned exon feautures: {assigned_exon.shape[0]}')
assigned_cds = pd.read_csv(assigned_cds_file, sep='\t', header=None, names=assigned_features_cols)
print(f'Read assigned CDS feautures: {assigned_cds.shape[0]}')

relative_cds_limits = pd.read_csv(relative_cds_limits_file, sep='\t').set_index('feature_id', drop=True)

print('Annotating assigned mRNA features.')
assigned_mrna_annotated = AnnotateAssignedGenes(assigned_mrna, assigned_exon, assigned_cds)
columns_to_drop = ['pass_multi_assign_filter', 'peak_display', 'unused', 'unused2', 'feature_phase', 'within_exon', 'within_cds']
assigned_mrna_annotated = assigned_mrna_annotated.drop(columns=columns_to_drop)

# map cds locations relative to mRNA start
assigned_mrna_annotated['feature_id'] = assigned_mrna_annotated['feature_desc'].apply(lambda x: GetGffDescValue(x, 'ID'))
assigned_mrna_annotated['relative_cds_start'] = assigned_mrna_annotated['feature_id'].map(relative_cds_limits['relative_cds_start'])
assigned_mrna_annotated['relative_cds_end'] = assigned_mrna_annotated['feature_id'].map(relative_cds_limits['relative_cds_end'])

# add columns for peak summit to cds start and end
assigned_mrna_annotated['summit_to_cds_start'] = assigned_mrna_annotated['summit_to_transcript_start'] - assigned_mrna_annotated['relative_cds_start']
assigned_mrna_annotated['summit_to_cds_end'] = assigned_mrna_annotated['summit_to_transcript_end'] - assigned_mrna_annotated['relative_cds_end']

print(f'Write annotated mRNA features: {assigned_mrna_annotated.shape[0]}')
assigned_mrna_annotated.to_csv(outfile_path, sep='\t', index=False)
print('Done')

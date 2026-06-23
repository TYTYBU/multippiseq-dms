import os
import pandas as pd
import pysam
import re
import argparse
import logging

import warnings
warnings.filterwarnings("ignore")

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

def parse_args():
    parser = argparse.ArgumentParser(
        description='This script calculates edit fractions for different chromatin accessibility bins.',
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Required parameters
    parser.add_argument('bamFile', 
                        type=str, 
                        default=None,
                        help=('BAM file containing the aligned reads. The filename needs to be structured as "*_[gene]_DMS_[block]_Dictionary.*.bam"\n'
                              'Default: None'))
    
    parser.add_argument('fastaFile', 
                        type=str, 
                        default=None,
                        help=('Fasta file of the genome. \n'
                              'Default: None'))

    parser.add_argument('varLib', 
                        type=str, 
                        default=None,
                        help=('DMS variant library with mutated oligos as a CSV file. \n'
                              'File must contain the columns: gene_aa_str, block_number, and variant_seq.\n'
                              'Default: None'))
        
    parser.add_argument('--out',
                        type=str,
                        default=None,
                        help=('Output CSV path; compression is inferred from the extension (e.g. .csv.gz → gzip). \n'
                              'Default: None'))

    return parser.parse_args()


# load arguments
args = parse_args()
bam_path=args.bamFile
fasta_path=args.fastaFile
lib_path=args.varLib


# load reference oligo sequences
logging.info(f'Loading reference oligo sequences: {fasta_path}')
fasta = pysam.FastaFile(fasta_path)

# load bam file
logging.info(f'Loading align reads bam file: {bam_path}')
if not os.path.exists(bam_path + '.bai'):
    pysam.index(bam_path)
bam = pysam.AlignmentFile(bam_path,'rb')


# extract gene and block info from bam filename.
# Parse the basename only: a directory in the path may contain underscores
# (e.g. .../Dual_NEBuilder/), which would otherwise truncate the gene token.
pat = re.compile('([^_]+)_DMS_([\\d]+)_Dictionary')
bam_basename = os.path.basename(bam_path)
gene = pat.search(bam_basename).group(1)
blk = pat.search(bam_basename).group(2)
ref_oligo_id = gene + '_blk_' + blk
ref_oligo_seq = fasta.fetch(ref_oligo_id)
df_r_seq = pd.DataFrame({'r_pos': range(len(ref_oligo_seq)), 'ref': list(ref_oligo_seq)})
logging.info(f'Parsed from bam filename, gene: {gene}, block: {blk}')



# load variant library
logging.info(f'Loading variant library file: {lib_path}')
df_variants = pd.read_csv(lib_path)

# Check for required columns
required_columns = {'gene_aa_str', 'block_number', 'variant_seq'}
missing_columns = required_columns - set(df_variants.columns)
if missing_columns:
    raise ValueError(f"Variant library file is missing required columns: {', '.join(missing_columns)}")

df_variants['gene'] = df_variants['gene_aa_str'].apply(lambda x: x.split('---')[0])
df_variants = df_variants.loc[(df_variants['gene'] == gene) & (df_variants['block_number'] == int(blk))]
logging.info(f'Number of variant oligos for {gene} block {blk}: {len(df_variants)}')



# identify mismatches for each variant
logging.info(f'Identify mismatches in variant oligos ...')
gene_aa_str = []
oligo_var = []
for index, row in df_variants.iterrows():
    var_oligo_seq = row['variant_seq']

    # trim the var_oligo_seq to the length of the ref_oligo_seq
    var_oligo_seq = var_oligo_seq[:len(ref_oligo_seq)]

    df_align = pd.DataFrame({'r_pos': range(len(ref_oligo_seq)), 'ref': list(ref_oligo_seq), 'alt': list(var_oligo_seq)})
    df_align = df_align.dropna().loc[df_align['ref'] != df_align['alt']]

    if len(df_align):
        oligo_var.append((df_align['ref'] + df_align['r_pos'].apply(str) + df_align['alt']).tolist())
        gene_aa_str.append(row['gene_aa_str'])
    
df_library_oligo_var = pd.DataFrame({'gene_aa_str': gene_aa_str, 'oligo_var': oligo_var})
df_library_oligo_var['n_oligo_var'] = df_library_oligo_var['oligo_var'].apply(lambda x: len(x))

# df_library_oligo_var.to_csv("df_library_oligo_var.csv")



# identify mismatches for each read
logging.info(f'Identify mismatches in reads ...')
bc_pair = []
oligo_var = []
ct = 0
for read in bam.fetch():
    aligned_pairs = read.get_aligned_pairs(with_seq=True, matches_only=True)
    
    if read.reference_name == ref_oligo_id:
        q_seq = read.query_sequence.upper()
        df_q_seq = pd.DataFrame({'q_pos': range(len(q_seq)), 'alt': list(q_seq)})
        df_align = pd.DataFrame(aligned_pairs, columns=['q_pos', 'r_pos', 'ref'])
        df_align['ref'] = df_align['ref'].str.upper()
        df_align = df_r_seq.merge(df_align, how='left', on=['r_pos', 'ref']).merge(df_q_seq, how='left', on='q_pos')
        df_align = df_align.dropna().loc[df_align['ref'] != df_align['alt']]

        if len(df_align):
            bc_pair.append(read.query_name)
            oligo_var.append((df_align['ref'] + df_align['r_pos'].apply(str) + df_align['alt']).tolist())

    ct += 1
    if ct % 10000 == 0:
        logging.info(f'Processed {ct} reads ...')

        # df_read_oligo_var = pd.DataFrame({'bc_pair': bc_pair, 'oligo_var': oligo_var})
        # df_read_oligo_var['n_oligo_var'] = df_read_oligo_var['oligo_var'].apply(lambda x: len(x))
        # df_read_oligo_var.to_csv("df_read_oligo_var.csv")
        
logging.info(f'Processed {ct} reads in total.')


df_read_oligo_var = pd.DataFrame({'bc_pair': bc_pair, 'oligo_var': oligo_var})
df_read_oligo_var['n_oligo_var'] = df_read_oligo_var['oligo_var'].apply(lambda x: len(x))





# match variants to reads
logging.info(f'Matching variants to reads ...')
gene_aa_str = []
ct = 0
for index, row in df_read_oligo_var.iterrows():
    read_oligo_var = row['oligo_var']
    gene_aa_str.append(df_library_oligo_var.loc[df_library_oligo_var['oligo_var'].apply(lambda x: set(x).issubset(set(read_oligo_var)))]['gene_aa_str'].tolist())

    ct += 1
    if ct % 10000 == 0:
        logging.info(f'Processed {ct} reads ...')
logging.info(f'Processed {ct} reads in total.')

df_read_oligo_var['gene_aa_str'] = gene_aa_str
df_read_oligo_var['n_gene_aa_str'] = df_read_oligo_var['gene_aa_str'].apply(lambda x: len(x))

df_read_oligo_var2 = df_read_oligo_var.loc[df_read_oligo_var['n_gene_aa_str'] == 1]
df_read_oligo_var2['gene_aa_str'] = df_read_oligo_var2['gene_aa_str'].apply(lambda x: x[0])

logging.info(f'Writing reads mapped with variants to {args.out}')
df_read_oligo_var2[['bc_pair', 'gene_aa_str']].to_csv(
    args.out,
    header=False,
    index=False,
    compression='infer',
)













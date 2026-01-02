#!/usr/bin/env python3 
import sys 
def parse_sample_info(<OUTPUT_DIR_OF_LNCPIPE.SBATCH>/sample_info):
   parts = sample_info.split('|')
   if len(parts) >= 7:
      gene_id = parts[0].split(':')[1] if ':' in parts[0] else parts[0]
      transcript_id = parts[1]
      exon_count = parts[2]
      length = parts[6]
      return gene_id, transcript_id, exon_count, length
      return None, None, None, None

def map_class_code(class_code): 
   mapping = {'u': 'Intergenic', 'i': 'Intronic', 'x': 'Antisense', 'o': 'Overlap', '=': 'Known', 'c': 'Contained', 'j': 'Novel_isoform', 'e': 'Extension' }
   return mapping.get(class_code, 'Unknown')

def determine_type(ref_gene, class_code):
   return 'known' if ref_gene != '-' and class_code == '=' else 'novel'

# Parse CPAT results
cpat_scores = {}
with open('<PATH_TO_CPAT_RESULTS>/cpat_results_full.txt', 'r') as f:
   for line in f:
      fields = line.strip().split('\t')
      if len(fields) >= 6:
         transcript_id = fields[0]
         cpat_score = float(fields[5])
         cpat_scores[transcript_id] = cpat_score

# Read tracking file
with open('<OUTPUT_DIR_OF_LNCPIPE.SBATCH>/gffcompare/stringtie_merged/stringtie_merged.tracking', 'r') as fin:
   with open('<OUTPUT_DIR_OF_LNCPIPE.SBATCH>/basic_charac.txt', 'w') as fout:
   # Write header
   fout.write("Gene\tTranscript\tType\tPotential\tLength\tExons\tClass\n")

   for line in fin:
      fields = line.strip().split('\t')
      if len(fields) >= 5:
      ref_gene = fields[2]
      class_code = fields[3]
      sample_info = fields[4]
      gene_id, transcript_id, exon_count, length = parse_sample_info(sample_info)
      if not gene_id:
         continue
      # Get CPAT score
      cpat_score = cpat_scores.get(transcript_id, 0.0)

      lnc_type = determine_type(ref_gene, class_code)
      bio_class = map_class_code(class_code)

      output_line = '\t'.join([ gene_id, transcript_id, lnc_type, str(cpat_score), length, exon_count, bio_class ])
      fout.write(output_line + '\n')

print("✓ Created basic_charac.txt with CPAT scores")

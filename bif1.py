pip install biopython

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

# DNA sequence
dna_sequence = Seq("ATGCGTACGTAGCTAGCTGATCGTAGCTAGTACGATCGTACGTAGCTGAC")

# 1. Calculate GC Content using gc_fraction
# [GC content = (Number of Gs + Number of Cs)/Total Number of Bases] * 100
gc_content = gc_fraction(dna_sequence)

# 2. Find specific motifs
motif1 = "ATG"
motif2 = "TAG"
motif1_positions = []
motif2_positions = []

for i in range(len(dna_sequence) - len(motif1) + 1):
    if dna_sequence[i:i+len(motif1)] == motif1:
        motif1_positions.append(i)

for i in range(len(dna_sequence) - len(motif2) + 1):
    if dna_sequence[i:i+len(motif2)] == motif2:
        motif2_positions.append(i)

# 3. Identify Coding Regions (ORFs)
start_codon = "ATG"
stop_codons = ["TAA", "TAG", "TGA"]

# checks for the occurrence of the specified motif in the sequence
# searches for ORFs, starting from the start codon ("ATG") and ending at any of the stop codons ("TAA", "TAG", "TGA")
def find_orfs(sequence):
    orfs = []
    for i in range(len(sequence) - 2):
        codon = str(sequence[i:i+3])
        if codon == start_codon:
            for j in range(i, len(sequence) - 2, 3):
                stop_codon = str(sequence[j:j+3])
                if stop_codon in stop_codons:
                    orfs.append((i, j+3))
                    break
    return orfs

coding_regions = find_orfs(dna_sequence)

# 4. Create a Report
report = "DNA Sequence Analysis Report\n\n"
report += f"Provided DNA Sequence:\n{dna_sequence}\n\n"
report += "Analysis 1: Finding Motifs\n"
report += f"Motif 1 ({motif1}) found at positions: {motif1_positions}\n"
report += f"Motif 2 ({motif2}) found at positions: {motif2_positions}\n\n"
report += "Analysis 2: Calculating GC Content\n"
report += f"GC Content: {gc_content:.2%}\n\n"
report += "Analysis 3: Identifying Coding Regions\n"
if len(coding_regions) > 0:
    report += "Coding regions found:\n"
    for start, stop in coding_regions:
        report += f"Start: {start}\nStop: {stop}\n"
else:
    report += "No coding regions found in the sequence."

# 5. Save the Report to a File
with open("dna_sequence_analysis_report.txt", "w") as report_file:
    report_file.write(report)

print("Report generated and saved as 'dna_sequence_analysis_report.txt'")

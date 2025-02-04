python

# Topics
"""
1. Sequence Handling and Manipulation (DNA to RNA)
2. Alignments (pairwise, multiple alignment)
3. Phylogenetics tree
4. BioSQL and Database Handling (Relational Database)
5. BioPandas (tabular format)
6. Query PubMed, Gene Ontology, Motifs
"""

# Import necessary libraries from Biopython
from Bio import Entrez, SeqIO, SwissProt, ExPASy, Geo
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo import PhyloXML
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqRecord
from Bio.ExPASy import ScanProsite
from urllib.request import urlopen
import pandas as pd

# Set your email for Entrez
Entrez.email = 'nsangani@iu.edu'

# --------------------------
# Section 1: Sequence Handling & Manipulation
# --------------------------

# Example: Translating a DNA sequence to a protein sequence
sequence = Seq("ATGAGTACGTA")
protein = sequence.translate()
print("Translated Protein Sequence:", protein)

# Example: Reverse complement of a DNA sequence
reverse_complement = sequence.reverse_complement()
print("Reverse Complement:", reverse_complement)

# Example: Creating and working with a SeqRecord (used for annotation)
record = SeqRecord(sequence, id="seq1", description="Sample DNA sequence")
print(record)

# --------------------------
# Section 2: Sequence Alignments
# --------------------------

# Reading a sequence alignment (e.g., from ClustalW or FASTA file)
alignment = SeqIO.parse("alignment.fasta", "fasta")
for seq in alignment:
    print(seq.id, seq.seq)

# Example: Pairwise sequence alignment (using Biopython's pairwise2)
from Bio import pairwise2
alignments = pairwise2.align.globalxx("AGCT", "AGT")
for alignment in alignments:
    print(alignment)

# --------------------------
# Section 3: Phylogenetics (Working with Phylogenetic Trees)
# --------------------------

# Reading a phylogenetic tree in Newick format
tree = Phylo.read("tree.nwk", "newick")
tree.rooted = True
print("Rooted Tree:")
Phylo.draw(tree)  # Visualize the tree

# Example: Working with tree objects (getting clade names)
print("Clade Names:", [clade.name for clade in tree.get_terminals()])

# --------------------------
# Section 4: BLAST Searches (Querying Sequence Databases)
# --------------------------

# Perform a BLAST search against the NCBI database
seq = Seq("AGCT")
result_handle = NCBIWWW.qblast("blastn", "nt", seq)
blast_records = NCBIXML.parse(result_handle)
for record in blast_records:
    print(f"Query: {record.query_id}")
    for alignment in record.alignments:
        print(f"Alignment title: {alignment.title}")
        print(f"Alignment length: {alignment.length}")

# --------------------------
# Section 5: Working with Protein Data (SwissProt and Prosite)
# --------------------------

# Retrieve and display a SwissProt record for a protein
handle = ExPASy.get_sprot_raw('P12345')  # Example protein ID
record = SwissProt.read(handle)
print("SwissProt Record Description:", record.description)

# Fetch Prosite motifs for a protein
handle = ExPASy.get_prosite_raw('P12345')
prosite_data = handle.read()
print("Prosite Data:", prosite_data)

# --------------------------
# Section 6: Gene Ontology (GO) and Sequence Annotations
# --------------------------

# Retrieve Gene Ontology (GO) annotations for a gene
# (For simplicity, here we use mock data or can pull data from specific GO databases)
handle = Entrez.efetch(db="gene", id="672", rettype="xml", retmode="text")
record = Entrez.read(handle)
go_annotations = record[0]["Entrezgene_properties"]["Gene-commentary_text"]
print("Gene Ontology Annotations:", go_annotations)

# --------------------------
# Section 7: BioSQL - Store and Retrieve Data from Databases
# --------------------------

# Store and query biological data using BioSQL (requires setup of a database)
# The following is a simplified example of storing and retrieving data using BioSQL
from BioSQL import BioSeqDatabase
from Bio import SeqIO
db = BioSeqDatabase.open_database("my_database")
seq_record = SeqIO.read("sequence.fasta", "fasta")
db.write(seq_record)  # Store a sequence
retrieved_record = db.lookup("seq1")  # Retrieve sequence by ID
print("Retrieved Sequence:", retrieved_record.seq)

# --------------------------
# Section 8: Working with BioPandas (DataFrames for Biological Data)
# --------------------------

# Example: Load gene expression data into a pandas DataFrame
data = {
    "Gene": ["Gene1", "Gene2", "Gene3"],
    "Expression_Level": [10.5, 12.3, 7.8]
}
df = pd.DataFrame(data)
print("Gene Expression Dataframe:")
print(df)

# --------------------------
# Section 9: Fetching and Analyzing Protein Sequences with ExPASy and ScanProsite
# --------------------------

# Example of scanning a sequence with ScanProsite for functional motifs
sequence = 'MALWMRLLPLLALLALWGPDPHMAEGLC'
handle = ScanProsite.scan(seq=sequence)
scan_result = ScanProsite.read(handle)
print("ScanProsite Result:", scan_result)

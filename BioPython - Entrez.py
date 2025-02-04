#!/usr/bin/env python

# Import necessary libraries from Biopython
from Bio import Entrez, SwissProt, ExPASy, Geo
from Bio.ExPASy import ScanProsite
from urllib.request import urlopen

# Set your email for Entrez (this is required by NCBI for identification purposes)
Entrez.email = 'nsangani@iu.edu'

# Step 1: EINFO - Retrieve information about available NCBI databases
# This command gets metadata about available databases
handle = Entrez.einfo()
result = handle.read()
handle.close()

# Print the raw result to see the structure
print("Result from Entrez.einfo():")
print(result)

# We can examine the type of result and access its keys
print("\nType of result:", type(result))  # It should be a string
print("Keys of the result (after reading into a dictionary):")
record = Entrez.read(handle)
print(record)  # This will be a dictionary after reading the result

# Step 2: Get detailed information about PubMed using EINFO
# Let's specifically query the PubMed database to get information about it
handle = Entrez.einfo(db='pubmed')
record = Entrez.read(handle)
print("\nRecord for PubMed database:")
print(record)

# Step 3: Search PubMed for articles related to "biopython" in the title
# We can use esearch to search for publications with a specific term
handle = Entrez.esearch(db='pubmed', term='biopython[title]')
record = Entrez.read(handle)

print("\nSearch result for 'biopython[title]':")
print(record)

# We can check the IDs of the articles retrieved from the search
print("\nPubMed IDs retrieved from search:")
print(record['IdList'])

# Step 4: Retrieve a summary of a PubMed article using its ID
# Here we fetch the article details based on a PubMed ID (replace with actual ID for real use)
handle = Entrez.esummary(db='pubmed', id='19304878')
record = Entrez.read(handle)
print("\nSummary of article with ID 19304878:")
print(record)

# Step 5: Fetch a GenBank nucleotide sequence using its ID
# Using efetch to retrieve sequence data from the Nucleotide database
handle = Entrez.efetch(db='nucleotide', id='EU490707', rettype='gb', retmode='text')
print("\nGenBank record for EU490707:")
print(handle.read())

# Step 6: Handling BioGeo data (e.g., GEO datasets)
# Use Entrez to search for GEO datasets (e.g., GSE)
handle = Entrez.esearch(db='gds', term='GSE')
record = Entrez.read(handle)

print("\nGEO dataset search result:")
print(record['Count'])  # Total number of results
print("Dataset IDs from GEO search:")
print(record['IdList'])

# Fetch each dataset's details
for rec_ID in record['IdList']:
    handle = Entrez.efetch(db='gds', id=rec_ID)
    data = Geo.parse(handle)  # Parse the GEO data
    # We could process the data here, like extracting relevant information

# Step 7: Using URL to fetch data (e.g., SwissProt data)
# Fetch a file from a URL to demonstrate downloading sequence information
url = 'https://raw.githubusercontent.com/biopython/biopython/master/Tests/SwissProt/F2CXE6.txt'
handle = urlopen(url)
sequence_data = handle.read().decode('utf-8')  # Read and decode to string
print("\nSwissProt record fetched from URL:")
print(sequence_data)

# Step 8: Retrieve a SwissProt entry using ExPASy (for protein data)
handle = ExPASy.get_sprot_raw('Q99527')  # Example accession number
record = SwissProt.read(handle)

print("\nSwissProt Record for Q99527:")
print(record.description)  # Print description of the protein

# Step 9: Get sequence length and other attributes of the protein
print("\nProtein sequence length:", record.sequence_length)

# Step 10: Fetch Prosite information for a protein (using the same SwissProt ID)
handle = ExPASy.get_prosite_raw('Q99527')
text = handle.read()
print("\nProsite information for Q99527:")
print(text)  # Raw data from Prosite

# Step 11: Fetch Prosite data using a Prosite ID
handle = ExPASy.get_prosite_raw('PDOC00001')  # Example Prosite ID
text = handle.read()
print("\nProsite information for PDOC00001:")
print(text)

# Step 12: Save Prosite data to an HTML file
# We can save raw data to a file for later use or inspection
with open('new.html', 'w') as out_handle:
    out_handle.write(text)
    print("\nProsite data saved to new.html.")

# Step 13: Scan a protein sequence using ScanProsite
# Example protein sequence for scanning with Prosite motifs
sequence = 'MEHKEVVLLLLLFLKSGQGEPLDDYVNTQGASLFSVTKKQLGAGSIEECAAKCEEDEEFTCRAFQYHSKEQQCVIMAENRKSSIIIRMRDVVLFEKKVYLSECKTGNGKNYRGTMSKTKN'
handle = ScanProsite.scan(seq=sequence, mirror='https://www.expasy.org/')
result = ScanProsite.read(handle)

# Output the result of the scan
print("\nScanProsite result:")
print(result)

result.n_seq
result_n.match
result[:]


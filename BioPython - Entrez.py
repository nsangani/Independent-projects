#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio import Entrez


# In[2]:


Entrez.email = 'nsangani@iu.edu'
handle = Entrez.einfo()
result = handle.read()
handle.close()


# In[3]:


type(result)
print(result)


# In[4]:


handle = Entrez.einfo()
record = Entrez.read(handle)


# In[5]:


print(record)
print(type(record))


# In[6]:


record.keys()
record['DbList']


# In[7]:


handle = Entrez.einfo(db='pubmed')
record = Entrez.read(handle)
print(record)
print(type(record))


# In[8]:


record['DbInfo'].keys()


# In[9]:


record['DbInfo']['Count']


# In[10]:


handle = Entrez.esearch(db='pubmed', term = 'biopython[title]')
record = Entrez.read(handle)
record.keys()
print(record)


# In[11]:


record['IdList']


# In[12]:


handle = Entrez.esummary(db='pubmed', id = '19304878')
record = Entrez.read(handle)
print(record)


# In[13]:


handle = Entrez.efetch(db='nucleotide', id = 'EU490707', rettype = 'gb', retmode='text')
print(handle.read())


# In[ ]:


handle = Entrez.esummary(db = 'pubmed', id = ID)
record = Entrez.read(handle)
for info in record:


# In[14]:


from Bio import Geo


# In[29]:


handle = Entrez.esearch(db='gds', term = 'GSE')
record = Entrez.read(handle)


# In[30]:


print(record['Count'])
record['IdList']


# In[24]:


for rec_ID in record['IdList']:
    handle = Entrez.efetch('gds', id = rec_ID)
    data = Geo.parse(handle)


# In[31]:


from urllib.request import urlopen


# In[32]:


url = 'https://raw.githubusercontent.com/biopython/biopython/master/Tests/SwissProt/F2CXE6.txt'
handle = urlopen(url)


# In[35]:


from Bio import ExPASy
handle = ExPASy.get_sprot_raw('Q99527')


# In[38]:


from Bio import SwissProt
record = SwissProt.read(handle)
print(record.description)


# In[39]:


dir(record)


# In[43]:


record.sequence_length


# In[44]:


handle = ExPASy.get_prosite_raw('Q99527')
text = handle.read()
print(text)


# In[46]:


handle = ExPASy.get_prosite_raw('PDOC00001')
text = handle.read()
print(text)


# In[47]:


with open ('new.html', 'w') as out_handle:
    out_handle.write(text)
    out_handle.close()


# In[50]:


sequence = 'MEHKEVVLLLLLFLKSGQGEPLDDYVNTQGASLFSVTKKQLGAGSIEECAAKCEEDEEFTCRAFQYHSKEQQCVIMAENRKSSIIIRMRDVVLFEKKVYLSECKTGNGKNYRGTMSKTKN'
from Bio.ExPASy import ScanProsite
handle = ScanProsite.scan(seq = sequence, mirror = 'https://www.expasy.org/')
result = ScanProsite.read(handle)


# In[ ]:


result.n_seq
result_n.match
result[:]


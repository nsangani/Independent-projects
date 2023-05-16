#!/usr/bin/env python
# coding: utf-8

# In[1]:


from bioinfokit import analys, visuz
import pandas as pd


# In[2]:


df = pd.read_csv('A_H_genes_gtf.csv')


# In[5]:


df['chr']=df['chr'].str.replace('chr', "")
df.head()


# In[7]:


visuz.marker.mhat(df=df, chr='chr', pv = 'pvalue', show=True, gwas_sign_line=True, markeridcol = 'gene', gstyle=2, dim=(10,2))


# In[ ]:





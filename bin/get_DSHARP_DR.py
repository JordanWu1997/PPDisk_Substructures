#!/usr/bin/env python
# coding: utf-8

# ## DSHARP Data Release Webpage
# - https://almascience.eso.org/almadata/lp/DSHARP/

# In[1]:


import wget
from os import chdir, mkdir
from time import sleep
from urllib.request import urlopen
from bs4 import BeautifulSoup


# In[2]:


def find_label(input_list, label):
    output_list = [inp for inp in input_list if label in inp]
    return output_list
def add_data_archive(input_list, data_archive):
    output_list = [data_archive + inp for inp in input_list]
    return output_list
def find_target(input_list, label='continuum', sep='_', target_id=0):
    output_list = [inp.split(sep)[target_id] for inp in input_list]
    return output_list


# In[3]:


data_archive = 'https://almascience.eso.org/almadata/lp/DSHARP/images/'
html  = urlopen(data_archive)
bsObj = BeautifulSoup(html, 'html.parser')
abras = bsObj.find_all('a')
hrefs = [abra.get('href') for abra in abras]


# In[4]:


href_cont = find_label(hrefs, 'continuum')
href_CO21 = find_label(hrefs, 'CO')
wget_cont = add_data_archive(href_cont, data_archive)
wget_CO21 = add_data_archive(href_CO21, data_archive)
target_ls = find_target(href_cont)


# In[5]:


print('{:3}: {:10} {:25} {:25}'.format('ID', 'Target', 'Continuum', 'COJ2-1'))
for i, target in enumerate(target_ls):
    print('{:>3d}: {:10} {:25} {:25}'.format(i, target, href_cont[i], href_CO21[i]))


# In[6]:


get_ipython().run_cell_magic('time', '', "#chdir('DSHARP_DR')\nfor i, target in enumerate(target_ls):\n    print(target)\n    mkdir(target)\n    chdir(target)\n    print('Cont: {}'.format(href_cont[i]))\n    wget.download(wget_cont[i])\n    sleep(15)\n    print('CO21: {}'.format(href_CO21[i]))\n    wget.download(wget_CO21[i])\n    sleep(15)\n    chdir('../')")


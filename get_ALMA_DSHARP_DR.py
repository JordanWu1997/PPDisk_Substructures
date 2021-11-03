#!/usr/bin/env python
# coding: utf-8

# ############################################################################
# :Author: Kuan-Hsien Wu
# :Email: jordankhwu@gapp.nthu.edu.tw
# :Date: 2021-03-28
# :Description: This code is to download ALMA DSHARP project observation data
#       DSHARP data release: https://almascience.eso.org/almadata/lp/DSHARP/
# ############################################################################

from os import chdir, mkdir
from time import sleep
from urllib.request import urlopen

from bs4 import BeautifulSoup

import wget


def find_label(input_list, label):
    """

    Args:
      input_list(1D list): input list
      label(str): label to filter

    Returns:
      1D list: label filtered input list

    """
    return [inp for inp in input_list if label in inp]


def add_data_archive(input_list, data_archive):
    """

    Args:
      input_list(1D list): input list
      data_archive(str): location of storage

    Returns:
      1D list: data archive added input list

    """
    return [data_archive + inp for inp in input_list]


def find_target(input_list, label='continuum', sep='_', target_id=0):
    """

    Args:
      input_list(1D list): input list
      sep(str): separation (Default value = '_')
      target_id: target index (Default value = 0)
      label: label of data (Default value = 'continuum')

    Returns:
      1D list: joint target list

    """
    return [inp.split(sep)[target_id] for inp in input_list]


def main():
    """Main function in get_ALMA_DSHARP_DR.py"""
    data_archive = 'https://almascience.eso.org/almadata/lp/DSHARP/images/'
    html = urlopen(data_archive)
    bsObj = BeautifulSoup(html, 'html.parser')
    abras = bsObj.find_all('a')
    hrefs = [abra.get('href') for abra in abras]

    href_cont = find_label(hrefs, 'continuum')
    href_CO21 = find_label(hrefs, 'CO')
    wget_cont = add_data_archive(href_cont, data_archive)
    wget_CO21 = add_data_archive(href_CO21, data_archive)
    target_ls = find_target(href_cont)

    print('{:3}: {:10} {:25} {:25}'.format('ID', 'Target', 'Continuum',
                                           'COJ2-1'))
    for i, target in enumerate(target_ls):
        print('{:>3d}: {:10} {:25} {:25}'.format(i, target, href_cont[i],
                                                 href_CO21[i]))


if __name__ == '__main__':
    main()

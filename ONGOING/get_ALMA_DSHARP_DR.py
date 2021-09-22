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

    :param input_list:
    :param label:

    """
    return [inp for inp in input_list if label in inp]


def add_data_archive(input_list, data_archive):
    """

    :param input_list:
    :param data_archive:

    """
    return [data_archive + inp for inp in input_list]


def find_target(input_list, label='continuum', sep='_', target_id=0):
    """

    :param input_list:
    :param label:  (Default value = 'continuum')
    :param sep:  (Default value = '_')
    :param target_id:  (Default value = 0)

    """
    return [inp.split(sep)[target_id] for inp in input_list]


def main():
    """ """
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

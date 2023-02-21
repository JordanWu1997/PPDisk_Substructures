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

import wget
from bs4 import BeautifulSoup


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


def query_ALMA_DSHARP_DR(
        data_archive='https://almascience.eso.org/almadata/lp/DSHARP/images/'):
    """

    Args:
      data_archive: root URL of data archive (here is DSHARP)

    Returns:
      1D list: target list
      nD list: fits file list of targets

    """

    # Connect to ALMA DSHARP archive
    html = urlopen(data_archive)
    bsObj = BeautifulSoup(html, 'html.parser')
    abras = bsObj.find_all('a')
    hrefs = [abra.get('href') for abra in abras]
    # Find continuum and CO 2-1
    href_cont = find_label(hrefs, 'continuum')
    href_CO21 = find_label(hrefs, 'CO')
    wget_cont = add_data_archive(href_cont, data_archive)
    wget_CO21 = add_data_archive(href_CO21, data_archive)
    target_list = find_target(href_cont)
    # Show found target fits files
    print('\nTarget found in {}:\n'.format(data_archive))
    print('{:3}: {:10} {:25} {:25}'.format('ID', 'Target', 'Continuum',
                                           'CO J2-1'))
    for i, target in enumerate(target_list):
        print('{:>3d}: {:10} {:25} {:25}'.format(i, target, href_cont[i],
                                                 href_CO21[i]))
    return target_list, [href_cont, href_CO21]


def get_ALMA_DSHARP_DR(target_list, fits_list):

    # Enter target id to download, default is to download all
    target_id_list = input(
        '\nEnter IDs of target that you want to download (default: all): '
    ).split() or [i for i in range(len(target_ls))]
    # Download files of targets
    for target_id in target_id_list:
        target_id = int(target_id)
        # Download fits files of target (e.g. continuum, CO J2-1)
        for fits in fits_list:
            print('\nDownloading: {:10} {:25}'.format(target_list[target_id],
                                                      fits[target_id]))
            wget.download(
                'https://almascience.eso.org/almadata/lp/DSHARP/images/{}'.
                format(fits[target_id]))


def main():
    """Main function in get_ALMA_DSHARP_DR.py"""

    target_ls, fits_list = query_ALMA_DSHARP_DR()
    get_ALMA_DSHARP_DR(target_ls, fits_list)


if __name__ == '__main__':
    main()

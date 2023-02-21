#!/usr/bin/env python
# coding: utf-8

# ############################################################################
# :Author: Kuan-Hsien Wu
# :Email: jordankhwu@gapp.nthu.edu.tw
# :Date: 2021-03-28
# :Description: This code is to download ALMA observation data from Pinte
# ############################################################################

import wget

# HD_97048_13CO32_briggs_selfcal_nocontsub.image.fits
wget.download(https://figshare.com/ndownloader/files/15457643)

# HD_97048_b7_continuum_centered_selfcal_manmask_briggs.image.fits
wget.download('https://figshare.com/ndownloader/files/15457640')

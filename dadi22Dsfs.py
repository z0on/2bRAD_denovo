#!/usr/bin/env python

import moments
import numpy as np
from numpy import array
from moments import Misc,Spectrum
import sys

infile=sys.argv[1]
pop_ids=[sys.argv[2],sys.argv[3]]
projections=[int(sys.argv[4]),int(sys.argv[5])]
outname=pop_ids[0]+pop_ids[1]+"_sfs"

dd = Misc.make_data_dict(infile)
data = Spectrum.from_data_dict(dd, pop_ids,projections,polarized=True)
Spectrum.to_file(data,outname,precision=4,foldmaskinfo=False)

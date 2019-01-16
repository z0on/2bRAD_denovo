# fuzzy calculator of probability that heterozygotes constitute 
# the majority of calls for a given site
# (for filtering out lumped paralogs)
# by Nathaniel "Nate" S. Pope (nspope@utexas.edu)

# uses output by ANGSD run with options "-doPost 2 -doGeno 11"
# output: contig, position, number non-missing samples, number "hard-called" heterozygotes, 
#         expected num heterozygotes, probability of heterozygote majority 
#         (last behaves slightly differently for odd/even numbers of samples)
# requires numpy, poibin, scipy
# put provided poibin.py into your PYTHONPATH location

from sys import stdin, stdout
from numpy import array, array_split, exp
from scipy.misc import logsumexp
from poibin import PoiBin
for line in stdin:
  line = line.strip().split("\t")
  chrom, pos, gl = line[0], line[1], array_split(array(line[4:], "float"), (len(line)-4)/4)
  pr_heteroz = array([x[2] for x in gl if not x[0]<0.])
  num_heteroz = sum([int(x[0]) for x in gl if x[0]==1])
  h_expected = sum(pr_heteroz)/len(pr_heteroz)
  try:
    pois_binom = PoiBin(pr_heteroz)
    utail_prob = pois_binom.pval(len(pr_heteroz)/2+1)
  except:
    utail_prob = 'NaN'
  stdout.write("\t".join([chrom, pos, str(len(pr_heteroz)), str(num_heteroz), str(h_expected), str(utail_prob)]) + "\n")

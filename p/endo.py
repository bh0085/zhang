#!/usr/bin/env python
'''
A utility to look up enzyme information from NEB.
'''
import sys, itertools as it, textwrap as tw

from pyquery import PyQuery as pq
import compbio.config as cfg
import subprocess as spc

def usage():
    print '\n'.join(list(it.chain(*[tw.wrap(l,80) for l in '''
usage: endo [enzyme_name]
mm
Finds the product number and opens the NEB website for the given enzyme.
Ignores case.
 
'''.splitlines()])))
    return 1

def enzyme_link(name):
    index = open(cfg.dataPath('zhang/neb_products.html')).read()
    d = pq(index)
    print d('a')
    named_elt = d('a').filter(lambda x: name.lower() in pq(this).text().lower())[0]
    product_link = named_elt.attrib['href']
    return product_link

    return 0

if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) != 1:
        exit(usage())

    product_link = enzyme_link(args[0])
    spc.Popen('open http://www.neb.com/nebecomm/products/{0}'.\
                  format(product_link),
              shell = True)
    exit(0)
    #link = re.search(re.compile(

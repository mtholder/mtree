#!/usr/bin/env python
import subprocess
import itertools
import sys
import os
import re
TOL = 0.0001
exe, arg = sys.argv[1:3]
inp_fn = os.path.split(arg)[-1]
expected = [float(i) for i in sys.argv[3:]]
fn = '.test_output'
with open(fn, 'w') as tout:
    p = subprocess.Popen([exe, arg],
                         stdout=tout,
                         stderr=None)
    rc = p.wait()
if rc != 0:
    sys.stdout.write(open(fn, 'r').read())
    sys.stderr.write('\nExited with err code = {}\n'.format(rc))
    sys.exit(1)
LN_PAT = re.compile(r'^\s*lnL\s*=\s*(-[0-9.]+)\s*\n')
lnL_list = []
with open(fn, 'r') as fo:
    for line in fo:
        m = LN_PAT.match(line)
        if m:
            gs = float(m.group(1))
            lnL_list.append(gs)
if len(expected) != len(lnL_list):
    sys.exit('''Expecting {} values to check the lnl values
{}
but only found:
{}
'''.format(len(lnL_list),
           " ".join([str(i) for i in lnL_list]),
           expected))
n = 0
for e, o in itertools.izip(expected, lnL_list):
    if abs(e - o) > TOL:
        sys.exit('lnL #{n} from "{f}" incorrect: {o} != {e}\n'.format(
            n=n,
            f=inp_fn,
            o=o,
            e=e))
    n += 1


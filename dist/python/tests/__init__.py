import sys

import os
import re

index_pcon = [i for i, word in enumerate(sys.path) if re.search("pcon-.+.egg$", word)][0]

sys.path.insert(0, sys.path.pop(index_pcon))

import pcon

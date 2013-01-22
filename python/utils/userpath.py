#! /usr/bin/env python
# userpath.py

"""
This module checks the user and adjusts the git directory accordingly.  Also adds git directory (and python subdirectory) to the path!

Usage (assumes userpath.py must be in the python path):
In [1]: from userpath import userpath
In [2]: file = userpath + '/dat/EOS/something.file'
"""

import os

user = os.environ['USER']

if user == 'youd':
    userpath = '/Users/youd/Documents/Work/Repositories/core_accretion'
elif user == 'ana-mariapiso':
    userpath ='/Users/ana-mariapiso/Documents/core_accretion'
elif user == 'apiso':
    userpath = '/home/apiso/repos/core_accretion'
    

def cduser(s = '/python'):
    """cd to the users python directory, or other directory in repository as specified by string"""
    os.chdir(userpath+s)

    

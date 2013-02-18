#!/usr/bin/env python
#--------------------------------------------------------------------------
#	Copyright 2013 Wolfgang Friederich
#
#	This file is part of Gemini II.
#
#	Gemini II is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 2 of the License, or
#	any later version.
#
#	Gemini II is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with Gemini II.  If not, see <http://www.gnu.org/licenses/>.
#----------------------------------------------------------------------------
#-----------------------------------------------------------
#  Python module to create dependencies from use and include
#  statements in Fortran programs
#
#  Usage: makeDepFromUseInclude [[dir] [dir]....]
#
#  Search all f90 files in . and given directories and
#  create a list of 
#  basename.o: deps.o lines
#  to be included into a Makefile
#-----------------------------------------------------------
import sys
import glob
import re
import os.path
#-----------------------------------------------------------
def search_text_file(datei, pattern):
    """ Search each line of a text file for a pattern.

    Start searching at the beginning of the line.
    Returns the group(1) match.
    """
    matchlist = []
    f = open(datei,'r')
    for line in f:
        m = re.match(pattern,line)
        if m is not None:
            matchlist.append(m.group(1))
    f.close()
    return matchlist
#-----------------------------------------------------------
def build_makefile_entry(datei, matchlist, depext):
    """Build a makefile entry.

    Form: datei_without_path_and_extension.o: dependencies.
    matchlist: list of dependencies
    depext: extension to be appended to the dependencies
    """
    target = str.split(os.path.basename(datei),'.')[0]+'.o:'
    entry = target
    for el in set(matchlist):
         entry = entry+' '+el+depext
    return entry
#-----------------------------------------------------------
dirs = ['.']+sys.argv[1:]                                 # list of directories including current one and those given on command line
#
#  search for use statements in f90 files
#
for folder in dirs:                                       # loop through folders
    for datei in glob.glob(folder+'/*.f90'):              # loop through Fortran 90 files
        matchlist = search_text_file(datei,'\tuse (\w+)')
        entry = build_makefile_entry(datei, matchlist,'.o')
        print entry
#
#  search for include statements in f-files
#
    for datei in glob.glob(folder+'/*.f'):                # loop through Fortran 77 files
        matchlist = search_text_file(datei,'\tinclude \'(\w+)')
        entry = build_makefile_entry(datei, matchlist,'.h')
        print entry
#
#  search for include statements in .f.m4-files
#
    for datei in glob.glob(folder+'/*.f.m4'):              # loop through m4-Fortran 77 files
        matchlist = search_text_file(datei,'\tinclude \'(\w+)')
        entry = build_makefile_entry(datei, matchlist,'.h')
        print entry

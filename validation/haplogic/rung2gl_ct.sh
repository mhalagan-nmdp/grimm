#!/bin/bash

DONFILE_IN=./data/recip_dpre.in.txt
DONFILE_OUT=./data/don.gl.txt
./g2gl_ct.pl -i $DONFILE_IN -o $DONFILE_OUT -t D

PATFILE_IN=./data/recip_dpre.in.txt
PATFILE_OUT=./data/pat.gl.txt
./g2gl_ct.pl -i $PATFILE_IN -o $PATFILE_OUT -t R


#!/usr/bin/env python
import subprocess
import argparse
import GFA

parser = argparse.ArgumentParser(description="Plots a gfa file from canu, generated both a png and a gml file")
parser.add_argument("gfaFile", help="Reads a GFA 1 file generated by Canu, and then plots it" )
args = parser.parse_args()
gfaFile = args.gfaFile


GFA.GFA( gfaFile  )



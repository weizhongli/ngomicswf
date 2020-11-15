#!/usr/bin/env python
import os
import re
import glob
import argparse
from argparse import RawTextHelpFormatter
import textwrap
import csv
from xlsxwriter.workbook import Workbook

banner = '''
==================================================================
Conert multiple tsv files into a single xlsx file
==================================================================
'''

if __name__ == "__main__":
  parser = argparse.ArgumentParser(formatter_class = RawTextHelpFormatter,
                                   description     = textwrap.dedent(banner))

  parser.add_argument('-i', '--input',       help="input tsv files, required", required=True)
  parser.add_argument('-o', '--output',      help="output xlsx file, required", required=True)
  parser.add_argument('-l', '--label',       help="optional labels for each sheet, optional")

  args = parser.parse_args()

  labs = re.split(',', args.input)
  for i in range(len(labs)):
    labs[i] = re.sub('.tsv','',labs[i])
    labs[i] = re.sub('.txt','',labs[i])

  if (not args.label is None):
    labs = re.split(',', args.label)

  workbook = Workbook(args.output, {'strings_to_numbers': True} )
  i = 0;
  for tsvfile in re.split(',', args.input):
    worksheet = workbook.add_worksheet(labs[i])
    with open(tsvfile, 'rt') as f:
      reader = csv.reader(f, delimiter='\t')
      for r, row in enumerate(reader):
        for c, col in enumerate(row):
          worksheet.write(r, c, col)
    i = i+1
  workbook.close()


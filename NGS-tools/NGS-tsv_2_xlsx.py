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
  args = parser.parse_args()

  workbook = Workbook(args.output, {'strings_to_numbers': True} )
  for tsvfile in re.split(',', args.input):
    worksheet = workbook.add_worksheet( re.sub('.tsv', '', tsvfile))
    with open(tsvfile, 'rt') as f:
      reader = csv.reader(f, delimiter='\t')
      for r, row in enumerate(reader):
        for c, col in enumerate(row):
          worksheet.write(r, c, col)
  workbook.close()


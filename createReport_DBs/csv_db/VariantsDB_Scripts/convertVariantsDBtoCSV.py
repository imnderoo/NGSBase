#!/usr/bin/python

import csv
import sqlite3

from glob import glob
from os.path import expanduser

conn = sqlite3.connect("../variants.db")
cursor = conn.cursor()
cursor.execute("select * from variants;")

with open("variants_export.csv", "wb") as csv_file:
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow([i[0] for i in cursor.description]) # write headers
    csv_writer.writerows(cursor)

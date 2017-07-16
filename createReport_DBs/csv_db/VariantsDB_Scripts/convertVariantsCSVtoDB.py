#!/usr/bin/python

import sqlite3, csv
import os, time
import argparse

def main():
	
	parser = argparse.ArgumentParser(description='Iterates through variant folder')
	parser.add_argument('variantsCSV', help='variants CSV')
	parser.add_argument('variantsDB', help='variants DB file with .db extension: [existing or new]')

	args = parser.parse_args()
	
	conn = sqlite3.connect (args.variantsDB)
	conn.text_factory = str
	cur = conn.cursor()
	cur.execute("""CREATE TABLE IF NOT EXISTS variants (chr VARCHAR, start_pos VARCHAR, ref VARCHAR, alt VARCHAR, gene VARCHAR, transcript VARCHAR, cdna_change VARCHAR, aa_change VARCHAR, cdna_type VARCHAR, aa_type VARCHAR, mutation_type VARCHAR, acmg VARCHAR, omit VARCHAR, omit_comment VARCHAR, date_last_updated VARCHAR, user_updated VARCHAR, PRIMARY KEY (chr, start_pos, transcript, ref, alt))""")

	reader = csv.DictReader(open(args.variantsCSV, "rb"))


	col_header = ",".join(reader.fieldnames)
	#print (col_header)

	rowCount = 0
	rowBatch = 2

	#cur.execute('BEGIN TRANSACTION')

	for row in reader:
		print ("-".join([row['chr'],row['start_pos'],row['gene'],row['transcript'],row['cdna_change'],row['mutation_type']]))

		for key in reader.fieldnames:
			row[key] = row[key].replace("\"", "")
			row[key] = row[key].replace("\'", "")

		col_value = ",".join([ "'" + row[key] + "'" for key in reader.fieldnames])
		rowCount = rowCount + 1
		
		sqlStatement = "INSERT OR REPLACE INTO variants (" + col_header + ") VALUES (" + col_value + ")"
		#print sqlStatement
		cur.execute(sqlStatement)

	conn.commit()
main()

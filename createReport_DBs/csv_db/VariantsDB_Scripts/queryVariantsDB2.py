#!/usr/bin/python

import sqlite3 as lite
import itertools

def main():
	conn = lite.connect('/home/miseq/miseq_pipe/resources/csv_db/variants.db')
	conn.row_factory = lite.Row
	variantsDict = getVariantsInfo(conn, '17', '41209079', '-', 'G', 'NM_007294')
	print variantsDict

def getVariantsInfo(conn,chr,start_pos,ref,alt,transcript):
	variantsDict = {}

	cur = conn.cursor()
	
	sql = """SELECT * FROM variants WHERE chr='{0}' AND start_pos={1} AND ref='{2}' AND alt='{3}' AND transcript='{4}'""".format(chr,start_pos,ref,alt,transcript)
	cur.execute(sql)
	
	if cur:
		for row in cur:

			variantsDict = dict(itertools.izip(row.keys(), row))
		
	return variantsDict
	
main()

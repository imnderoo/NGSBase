#!/usr/bin/python

import sqlite3 as lite

conn = lite.connect('/home/miseq/miseq_pipe/resources/csv_db/run_metrics.db')
cur = conn.cursor()

def get_posts():
	cur.execute("SELECT * FROM run_qc")
	varDB = cur.fetchall()

	fieldNames = ['ID', 'UNIFORMITY', 'PCT_DUPLICATE', 'TI_TV_RATIO', 'TOTAL_PF_READS', 'PCT_ALIGNED', 'MEDIAN_FRAGMENT_LENGTH', 'MEAN_COVERAGE_DEPTH', 'PCT_Q30', 'PCT_ON_MANIFEST', 'R1_MEAN_ERROR', 'R2_MEAN_ERROR', 'GC_DROPOUT', 'AT_DROPOUT']

	print ','.join(fieldNames)

	for entry in varDB:
		print ','.join(str(e) for e in list(entry))
get_posts()

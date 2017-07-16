#!/usr/bin/python

import sqlite3 as lite

conn = lite.connect('/home/miseq/miseq_pipe/resources/csv_db/variants.db')
cur = conn.cursor()

print ("chr,start_pos,ref,alt,gene,transcript,cdna_change,aa_change,cdna_type,aa_type,mutation_type,acmg,omit,omit_comment,date_last_updated,user_updated")

def get_posts():
	cur.execute("SELECT * FROM variants")
	varDB = cur.fetchall()

	for entry in varDB:
		print (",".join(list(entry)))
get_posts()

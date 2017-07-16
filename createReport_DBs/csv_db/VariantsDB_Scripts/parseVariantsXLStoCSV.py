#!/usr/bin/python

import xlrd, csv
import os, time
import argparse

def main():
	
	parser = argparse.ArgumentParser(description='Iterates through variant folder')
	parser.add_argument('varFolder', help='variant folder')
	parser.add_argument('outCSV', help='output csv file')
	args = parser.parse_args()
	
	variantsList = []
	variantKeys = ['chr','start_pos','ref','alt','gene','transcript','cdna_change','aa_change','cdna_type','aa_type','mutation_type','acmg','omit','omit_comment','date_last_updated','user_updated']

	# Extract Excel Files into VCF
	#for xlsFile in os.listdir(args.varFolder):
	for root, subdirs, files in os.walk(args.varFolder):
		for xlsFile in files:
			if ".xls" in xlsFile:
				print xlsFile
				xlsFile = os.path.join(root, xlsFile)
				#xlsFile = os.path.join(args.varFolder, xlsFile)
				parseExcelToVarDB(xlsFile, variantsList)
	
	with open(args.outCSV, 'wb') as f:
	#with open('/home/miseq/miseq_pipe/resources/csv_db/variants_new.csv', 'wb' ) as f:
	
		dict_writer = csv.DictWriter(f, fieldnames=variantKeys)
	
		dict_writer.writeheader()

		for variant in variantsList:
			try:
				dict_writer.writerow(variant)
			except: 
				continue

def parseExcelToVarDB(xlsFile, variantsList):
	
	variantsDict = createVariantsDict()
	
	try:
		workbook = xlrd.open_workbook(xlsFile, on_demand=True)
	except:
		return

	try:
		worksheet = workbook.sheet_by_name('Alamut.Common')
	except:
		return
	try:
		worksheetSummary = workbook.sheet_by_name('Summary')
	except:
		return
	
	print xlsFile

	num_rows = worksheet.nrows - 1
	num_cells = 3 #1st col: Heading, 2nd col: Alamut.Heading, 3rd col: value
	curr_row = -1
	
	alamutDict = createAlamutToVarDBDict()
	acmgDict = createACMGDict()
	
	while curr_row < num_rows:
		curr_row += 1
		row = worksheet.row(curr_row)
		
		alamutCat = worksheet.cell_value(curr_row, 1)
		value = worksheet.cell_value(curr_row, 2)
		valueType = worksheet.cell_type(curr_row, 2)

		# Remove the decimal from transcript
		if alamutCat == 'Alamut.transcript':
			value = value.split('.')[0]
		
		# Remove the decimal from numerical values
		if valueType == 2:
			value = int(value)
		
		# Keep the values
		if alamutCat in alamutDict:
			variantsDictCat = alamutDict[alamutCat]
			variantsDict[variantsDictCat] = value
			
	workbook.unload_sheet('Alamut.Common')
	
	num_rows = worksheetSummary.nrows - 1
	num_cells = 3 #1st col: Heading, 2nd col: Alamut.Heading, 3rd col: value
	curr_row = -1
	
	while curr_row < num_rows:
		curr_row += 1
		row = worksheetSummary.row(curr_row)
		
		cat = worksheetSummary.cell_value(curr_row, 0)

		if 'Fellow Interpretation' in cat:
			value = worksheetSummary.cell_value(curr_row, 1)
			
			if value:
				value = value.split('-')[0].strip().title().replace("Likely", "Predicted").replace("Unknown Significance", "VUS").replace("Vus", "VUS")
				variantsDictCat = alamutDict[cat]
				variantsDict[variantsDictCat] = value
	
				try:
					variantsDict['acmg'] = acmgDict[value]
				except:
					variantsDict['acmg'] = 0
					variantsDict[variantsDictCat] = value
	
		if 'Variant Interpretation' == cat:
			value = worksheetSummary.cell_value(curr_row, 1)
		
			if value:	
				value = value.split('-')[0].strip().title().replace("Likely", "Predicted").replace("Unknown Significance", "VUS").replace("Vus", "VUS")

				variantsDictCat = alamutDict[cat]
				variantsDict[variantsDictCat] = value
	
				try:
					variantsDict['acmg'] = acmgDict[value]
				except:
					variantsDict['acmg'] = 0
					variantsDict[variantsDictCat] = value
	
	variantsDict['date_last_updated'] = time.strftime('%m/%d/%Y', time.gmtime(os.path.getmtime(xlsFile)))
	#os.path.getmtime(xlsFile)
	
	variantsList.append(variantsDict)
	
def createAlamutToVarDBDict():
	# chr	start_pos	ref	alt	gene	transcript	cdna_change	aa_change	cdna_type	aa_type	mutation_type	acmg	omit	omit_comment

	dict = {
		'Alamut.chrom':'chr',
		'Alamut.gDNAstart':'start_pos',
		'Alamut.wtNuc':'ref',
		'Alamut.varNuc':'alt',
		'Alamut.gene':'gene',
		'Alamut.transcript':'transcript',
		'Alamut.cNomen':'cdna_change',
		'Alamut.alt_pNomen':'aa_change',
		'Alamut.varType':'cdna_type',
		'Alamut.codingEffect':'aa_type',
		'Fellow Interpretation & Rationale':'mutation_type',
		'Variant Interpretation':'mutation_type'
	}	
	
	return dict

def createACMGDict():
	
	dict = {
		'Not Assessed':'0',
		'Pathogenic':'1',
		'Predicted Pathogenic':'2',
		'VUS':'3',
		'Predicted Benign':'4',
		'Benign':'5'
	}
	
	return dict
	
def createVariantsDict():
	
	dict = {
		'chr':'-',
		'start_pos':'-',
		'ref':'-',
		'alt':'-',
		'gene':'-',
		'transcript':'-',
		'cdna_change':'-',
		'aa_change':'-',
		'cdna_type':'-',
		'aa_type':'-',
		'mutation_type':'Not Assessed',
		'acmg':'0',
		'omit':'N',
		'omit_comment':'-',
		'date_last_updated':'-',
		'user_updated':'Script'
	}
	
	return dict
	
main()

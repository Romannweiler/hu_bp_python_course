# -*- coding: utf-8 -*-
import pandas as pd
import sqlite3 as lite
import cPickle as pickle
import collections as col

class KnowledgeBase(object):
	"""
	Databaseobject for organizing data in a database and being able to export those tables (in .sql-files).

	@ivar availableTables: a list of the available tables of the object. Note that editing this will not change any method and thus not change the output in any way. This variable solely exists as a help for the user.
	@type availableTables: list of strings
	@authors: Steffen Miels & Kora Schmitt
	"""
	def __init__(self):
		self.availableTables=['genes', 'protein_hl','metabolites','conc','stopcodon', 'genome']
		
	def datainput(self):
		"""
		Much rather a script than a method that takes the data from the csv-files, imports them as a dataframe-object and then exports the newly formated data to sql-files in the export folder. Furthermore a database gets set up using SqLite, in which tables can be found (specifically the ones named in).

		Note that the export-folder and its content is NOT necessary for the database to work.
		"""
		####Input aus Excel-Tabellen
		#Einlesen der Gen-, Metabolit- und Mediumsspezifikation
		genes = pd.read_csv('./csvfastadata/genes.csv',header=4, index_col=None)
		metabolites =pd.read_csv('./csvfastadata/metabolite.csv',header=4, index_col=0, na_values=['NaN'], keep_default_na=False)
		media = pd.read_csv('./csvfastadata/media.csv',header=3,index_col=0, na_values=['NaN'],keep_default_na=False)
		protein_hl = pd.read_csv('./csvfastadata/mRNA.csv',header=4, index_col=0, na_values=[' '],keep_default_na=False)
		stopcodon = pd.read_csv('./csvfastadata/stopcodon.csv',header=2, index_col=None)

		#Editieren der Dataframes
		#Herauswerfen aller unwichtigen Spalten
		metabolites = metabolites['Name']
		genes = genes[['BioCyc','SwissProt','GenBank','Sequence','Coordinate']]
		media = media['Concentration mM']
		protein_hl = protein_hl[['Symbol', 'Protein Halflife']]


		#Konkatenieren der Listen an entsprechender Stelle und Herauswerfen der NaN Werte
		conc = pd.concat([metabolites, media], axis=1)
		conc = conc[pd.notnull(conc['Concentration mM'])]

		##Laden des vollstaendigen Genoms in einen Objekt
		with open('./csvfastadata/sequence.fasta','r') as f:
			genome = f.readlines()
		f.close()
		#Herauswerfen der ersten Zeile und zusammenhaengen der Sequenzen in einen String
		genome = [line[:-1] for line in genome if line[0] != '>']
		genome = ''.join(genome)
	#Ueberschreiben der csv files

	##Datenbank
		def writeData(data,name):
			f = open(name+'.sql','w')
			with f:
				f.write(data)

		#Schreiben in die Datenbank
		con = lite.connect('./database/knowledgebase.db')
		con.text_factory = str



		#Writing sql-files for the database
		#helplist
		helplist1 = ['genes', 'protein_hl','metabolites','conc','stopcodon']
		helplist2 = [genes, protein_hl, metabolites, conc, stopcodon]
		helper = dict((helplist1[t], helplist2[t]) for t in xrange(len(helplist1)))
		for j in helper.keys():
			with con:
				cur = con.cursor()
				cur.execute("DROP TABLE IF EXISTS "+j )
				wert=helper[j].to_sql(j,con)
				data = '\n'.join(con.iterdump())
				writeData(data, "./database/export/"+j)

		with con:
			cur = con.cursor()
			cur.execute("DROP TABLE IF EXISTS genome" )
			cur.execute("CREATE TABLE genome(Sequence TEXT)")
			cur.execute("INSERT INTO genome VALUES('"+genome+"')")
			data = '\n'.join(con.iterdump())
			writeData(data, "./database/export/genome")
		return None



	###############################################################################

	def getfromdb(self,namestring):
		"""
		Extracts wanted tables from the database.
		
		@param namestring: this specifies, which table will be extracted.
		@type namestring: str
		@return: if namestring is 'genome'
					- returns a str
				 else 
					- a dataframe-object
		"""
		helplist1 = ['genes', 'protein_hl','metabolites','conc', 'genome','stopcodon']
		if namestring not in helplist1:
			raise ValueError('Wanted Table not in database')
		con = lite.connect('./database/knowledgebase.db')
		with con:
			j=namestring
			cur = con.cursor()    
			cur.execute('SELECT * FROM '+j)
			data = cur.fetchall()
			indexhelp = [data[i][0] for i in xrange(len(data))]
			#Naming columns and rows properly
			if j=='genes':
				dframe = pd.DataFrame(data, index=indexhelp, columns=['dump','BioCyc', 'SwissProt','GenBank','Sequence', 'Coordinate'])
				dframe = dframe.drop('dump',axis=1)
			if j=='protein_hl':
				indexhelp = [data[i][0] for i in xrange(len(data))]
				dframe = pd.DataFrame(data,index = indexhelp, columns=['dump','Name','Halflife min'], dtype=float)
				dframe = dframe.drop('dump',axis=1)
			if j=='metabolites':
				indexhelp = [data[i][0] for i in xrange(len(data))]
				dframe = pd.DataFrame(data, index=indexhelp, columns=['dump','Name'])
				dframe = dframe.drop('dump',axis=1)
			if j=='conc':
				dframe = pd.DataFrame(data, index=indexhelp, columns=['dump','Name', 'Concentration mM'], dtype=float)
				dframe = dframe.drop('dump',axis=1)
			if j=='stopcodon':
				dframe = pd.DataFrame(data, index=indexhelp, columns=['dump','RNA','DNA'])
				dframe = dframe.drop('dump',axis=1)
			if j=='genome':
				dframe = data[0][0]
		return dframe


	###############################################################################
	def putindb(self,name, dframe):
		"""
		Stores a dataframe table in the database.

		@param name: name under which the dframe gets stored in the database
		@type name: str or unicode

		@param dframe: dataframe with the data
		@type dframe: dataframe
		"""
		con = lite.connect('./database/knowledgebase.db')
		con.text_factory = str
		cur = con.cursor()
		cur.execute("DROP TABLE IF EXISTS "+name )
		dframe.to_sql(name,con)
		data = '\n'.join(con.iterdump())
		writeData(data, "./database/export/"+name)
	


	def deletefromdb(self,name):
		"""
		Deletes a table from the db.

		@param name: name of the table that is supposed to be dropped.
		@type name: str or unicode
		"""
		con = lite.connect('./database/knowledgebase.db')
		cur = con.cursor()
		cur.execute("DROP TABLE IF EXISTS "+name )
	###############################################################################
	def get_sequence(self,name):
		"""
		Solely returns the dna-sequence of the gene specified by 'name'.

		Admissable names are the ones of SwissProt, GenBank, BioCyc

		@param name: Name of the gene in SwissProt or GenBank or BioCyc
		@type name: str

		@return: 'unicode'
		"""
		return self.get_gene(name)[0]

	def get_gene(self,name):
		"""
		Returns the dna-sequence and the coordinate of the gene specified by 'name'.

		Admissable names are the ones of SwissProt, GenBank, BioCyc. Can be found in ./database/genenames.ods

		@param name: Name of the gene in SwissProt or GenBank or BioCyc
		@type name: str

		@return: tuple of 'unicode' and 'numpy.int64'
		"""
		genes = self.getfromdb('genes')
		if name not in map(str,genes['BioCyc'])and name not in map(str,genes['SwissProt'])and name not in map(str,genes['GenBank']):
			raise AttributeError('Attribute '+str(name)+' does not exist.')
		else:
			try:
				stelle=map(str,genes['BioCyc']).index(name)
				return genes.iloc[stelle][-2],genes.iloc[stelle][-1]
			except:
				try:
					stelle=map(str,genes['SwissProt']).index(name)
					return genes.iloc[stelle][-2],genes.iloc[stelle][-1]
				except:
					stelle=map(str,genes['GenBank']).index(name)
					return genes.iloc[stelle][-2],genes.iloc[stelle][-1]
	def get_protein_hl(self,name):
		"""
		Returns the protein half-life of a specified protein.

		Admissable names range from MG_001 to MG_470.

		@param name: Name of the protein
		@type name: str

		@return: 'unicode'
		"""
		protein_hl = self.getfromdb('protein_hl')
		if name not in protein_hl.index:
			raise AttributeError('Attribute '+str(name)+' does not exist.')
		else:
			return protein_hl.loc[name]['Halflife min']

	def get_concentration(self,name):
		"""
		Returns the concentration of <name> in the cell.
		Admissable names can be found in ./database/medianames.ods

		@param name: name of the medium
		@type name: str

		@returns: unicode
		"""
		conc = self.getfromdb('conc')

		if name not in conc.index and name not in map(str,conc['Name']):
			raise AttributeError('Attribute '+str(name)+' does not exist.')
		else:
			try:
				return conc.loc[name]['Concentration mM']
			except:
				stelle=map(str,conc['Name']).index(name)
				return conc.iloc[stelle][1]
	def get_metabolites(self):
		"""
		Returns metabolites table from the database.

		@returns: dataframe
		"""
		#returns a dataframe consisting of contained metabolites
		return self.getfromdb('metabolites')
	def get_genome(self):
		"""
		Returns the complete genome of Mycoplasma genitalium.

		@return: unicode	
		"""
		##returns a string
		return self.getfromdb('genome')
	def get_stopcodon(self):
		"""
		Returns a table with the stopcodons for DNA and RNA stored in a dataframe.

		@return: dataframe
		"""
		return self.getfromdb('stopcodon')
	###############################################################################

if __name__=='__main__':
	print 'hello world'	
	#import KnowledgeBase
	kb=KnowledgeBase()
	print kb.get_sequence('MG_001')

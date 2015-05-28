# -*- coding: utf-8 -*-
import pandas as pd
import sqlite3 as lite
import cPickle as pickle
import collections as col


class KnowledgeBase(object):
	"""
	Epydoc
	@param
	@type

	@return

	"""
	def __init__(self):
		pass
	
	def datainput(self):
		"""
			doku
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
				writeData(data, "./database/"+j)

		with con:
			cur = con.cursor()
			cur.execute("DROP TABLE IF EXISTS genome" )
			cur.execute("CREATE TABLE genome(Sequence TEXT)")
			cur.execute("INSERT INTO genome VALUES('"+genome+"')")
			data = '\n'.join(con.iterdump())
			writeData(data, "./database/genome.sql")
		return None




	#Reading from db
	def getfromdb(self,namestring):
		"""
		Doku
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


	###########


	#Uebergabefunktionen
	def get_sequence(self,name):
		"""
		doku
		"""
		return self.get_genes(name)[0]

	def get_genes(self,name):
		"""
		doku
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
		doku
		"""
		protein_hl = self.getfromdb('protein_hl')
		if name not in protein_hl.index:
			raise AttributeError('Attribute '+str(name)+' does not exist.')
		else:
			return protein_hl.loc[name]['Halflife min']

	def get_concentration(self,name):
		"""
		doku
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
			doku
		"""
		#returns a dataframe consisting of contained metabolites
		return self.getfromdb('metabolites')
	def get_genome(self):
		"""
			doku
		"""
		##returns a string
		return self.getfromdb('genome')
	def get_stopcodon(self):
		"""
			doku
		"""
		return self.getfromdb('stopcodon')


if __name__=='__main__':
	print 'hello world'	

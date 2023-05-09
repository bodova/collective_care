import subprocess
import numpy as np
import utils.auxSolomonAnalysis as aSA

colorabbrv = {'o':'orange',
			  'b':"blue",
			  'g':'green',
			  'y':'yellow',
			  'r':'red',
			  'p':'purple'}

class pelletFileReader:

	def __init__(self, filename,separator=','):
		self.dropping_dict=dict()
		with open(filename) as fin:
			num_read_lines = 0;
			for line in fin:
				num_read_lines += 1
				if num_read_lines < 2:
					continue
				lines = line.strip().split(separator)
				treat=lines[0]
				rep = int(lines[1])
				color = colorabbrv[lines[3]]
				key = (treat, rep, color)
				if key not in self.dropping_dict.keys():
					self.dropping_dict[key] = []
				frame=int(lines[5])
				self.dropping_dict[key].append(frame)

	def get_droppings(self,treatment,replicate,color):
		key = (treatment, replicate, color)
		if key not in self.dropping_dict.keys():
			return []
		return self.dropping_dict[key]

class SolomonReader:

	paramTranslator = { 'CBCount':'numButtons',
					'SubjCount':'numAnimals',
					'TimeResolution':'timeResolution'}

	def wccount(self,filename):
	   out = subprocess.Popen(['wc', '-l', filename],
	                                                stdout=subprocess.PIPE,
	                                                stderr=subprocess.STDOUT
	                                                ).communicate()[0]
	   return int(out.partition(b' ')[0])



	#Initializes an object by reading a solomon annotation file from fileName.
	# It only considers the frames between initFrame and endFrame.
	# If initFrame is None, the starting frame is the first one, and if endFrame is none the ending frame is the last.
	#
	#This object contains:
	# 0. numAnts : the number of ants
	# 1. theMatr :  a matrix of number of frames times number of ants which contains the recorded time-series::
	#      The first numAnts entries are interpreted as follows:
	#         if the i'th entry is j, then it means ant i+1 is grooming ant j
	#      The next numAnimals entries are interepreted as follows:
	#         if the (i-numAnts)'th entry is y it means ant i+1 is doing action y
	#	   The next numAnimals entries are interpreted as follows:
	#  		 if the (i+2*numAnts)'th entry is j, it means ant i+1 and ant j are doing trophallaxis
	# 2. buttonNames : which contains the names of the solomon buttons,
	#   the first numAnts entries are the identifiers of the ants (colors),
	#   and the next are the name of the actions they can make,
	#   such that the j'th entry is the name of action j+1, as encoded in theMatr.
	#   That is, if theMatr[f,numAnts+i] = j, it means that in frame initFrame+f, ant i+1 is 
	#   performing the action whose name is given by buttonNames[j-1]
	#
	# Furthermore, this object has the functions:
	#   getGroomingMatrix   and  getTrophallaxisMatrix
	#        both of which return numAnts x numAnts matrices
	#       whose i,j entries count the number of frames in which  ant i+1 is doing grooming or trophallaxis, respectively, to ant j+1

	def __init__(self,fileName,initFrame=None,endFrame=None,enforceSymetryInTrophallaxis=True):
		totNumLines = self.wccount(fileName)
		fopen = open(fileName,'r');
		self.params = dict();
		self.buttonNames = [];
		lineNum = 0;
		print(fileName+" "+str(initFrame)+" "+str(endFrame))
		#First we read the headers to get the parameters
		for line in fopen:
			lineNum += 1;
			lstrp  = line.strip();
			lsplit = lstrp.split("\t")
			if lsplit[0] == "Time":
				break
			self.readParams(lsplit);

		self.numAnts = self.params['numAnimals']
		#We now know the size of the  time series matrix
		totNumFrames = totNumLines - lineNum;
		headerLines = lineNum;

		otherCode = [self.buttonNames.index('other')+1,
		             self.buttonNames.index('rest')+1]  
		self.theMatr  = np.zeros((totNumFrames,3*self.numAnts),dtype=np.uint8)
		self.enforceSymetryInTrophallaxis=enforceSymetryInTrophallaxis;

		#We can now read the rest of the data to fill up the matrix
		t = 0;
		for line in fopen:
			lineNum += 1;
			if ',' in line:
				continue
			lstrp  = line.strip("\n\r");
			lsplit = lstrp.split("\t")
			if len(lsplit)>2:
				lsplit = lsplit[-2:]
			
			vect = self.readDataLine(lsplit,otherCode = otherCode)
			self.theMatr[t,:] = vect;
			t+=1;

		fopen.close()

		self.theMatr = self.theMatr[initFrame:endFrame,:];

		#we clean up the all-zero rows
		numMarks = self.theMatr.sum(axis=1);
		withMark = np.nonzero(numMarks)[0];
		tInit = min(withMark);
		tEnd  = max(withMark)+1;
		self.theMatr = self.theMatr[tInit:tEnd,:];


		
	#Returns a vector of size:
	#     3 * numAnimals
	#     The first numAnimals entries are interpreted as follows:
	#        if the i'th entry is j, then it means ant i+1 is grooming ant j
	#     The next numAnimals entries are interepreted as follows:
	#        if the (i-numAnimals)'th entry is y it means ant i+1 is doing action y
	#	  The next numAnimals entries are interpreted as follows:
	# 		 if the (i-2*numAnimals)'th entry is j, it means ant i+1 and ant j are doing trophallaxis
	def readDataLine(self,lsplt,otherCode=[8,9]): 

		# Note: trophallaxis was encoded in a similar way to grooming:
		# with a key for a performer, a key for a receiver and a key to indicate 
		# the behaviour is trophallaxis: 'E' or by 'O'. 
		# These keys are used to code 'rest' and 'other' behaviors when used by themselves. 

		if len(lsplt)==0:
			return None

		vector = np.zeros(3*self.numAnts,dtype=np.uint8);


		groomingPart = lsplt[0].replace(',','.');
		groomingDestination = groomingPart.split(";");
		if groomingPart != "":
			for groomingSource in range(len(groomingDestination)):
				posDes = groomingDestination[groomingSource];
				if posDes != "":
					destination = int(posDes);
					vector[groomingSource]=destination;
					if destination>self.numAnts:
						print("\t ** "+str(groomingPart)+": ")

		if len(lsplt)==1:
			return vector

		selfPart = lsplt[1].replace(',','.');
		selfAction = selfPart.split(";");
		ammendment = np.ones_like(vector);
		if selfPart!="":
			for selfActor in range(len(selfAction)): #between 0 and numAnimals-1
				posAction = selfAction[selfActor];
				if posAction != "":
					action = int(posAction);
					if not(action in otherCode and vector[selfActor]>0):
						vector[self.numAnts+selfActor] = action
					else:
						vector[2*self.numAnts+selfActor] = vector[selfActor];
						ammendment[selfActor] = 0;

		return vector*ammendment


	def readParams(self,lsplt):
		if len(lsplt) < 2:
			return

		configKey = lsplt[0];
		if configKey in SolomonReader.paramTranslator.keys():
			toSave = lsplt[1];
			if configKey == 'TimeResolution':
				toSave = float(toSave.replace(',','.'))
			else:
				toSave = int(toSave);
			self.params[SolomonReader.paramTranslator[configKey]] = toSave;

		if configKey == 'CBName':
			self.buttonNames.append( lsplt[1] );


	def getTrophallaxisMatrix(self,timeSeries=None):
		if timeSeries == None:
			timeSeries = self.theMatr;
		nAnimals = self.params['numAnimals'];
		tMatrix = np.zeros((nAnimals,nAnimals));
		for a1 in range(nAnimals):
			for a2 in range(1,nAnimals+1):
				a1TROPHSa2 = np.nonzero(timeSeries[:,2*nAnimals+a1]==a2)[0];
				if self.enforceSymetryInTrophallaxis:
					a2TROPHSa1 = np.nonzero(timeSeries[:,2*nAnimals+a2-1]==a1+1)[0];
					tMatrix[a1,a2-1] = len(set(a1TROPHSa2) & set(a2TROPHSa1));
				else:
					tMatrix[a1,a2-1] = len(a1TROPHSa2);

		if self.enforceSymetryInTrophallaxis:
			return (tMatrix+tMatrix.T)/2;
		return tMatrix


	#returns a list of tuples of the form:
	#  (start,end,actor,groomedAnt)
	#  for each grooming event
	# The order of this events is by actor number, then time
	def getGroomingEvents(self,timeSeries=None):
		if timeSeries == None:
			timeSeries = self.theMatr;

		nAnimals = self.params['numAnimals'];
		result = [];
		for actor in range(nAnimals):
			thisActorsTS = timeSeries[:,actor];
			nzi = aSA.nonzero_intervals(thisActorsTS);
			for evStart,evEnd in nzi:
				groomed = thisActorsTS[evStart]-1;
				result.append([evStart,evEnd,actor,groomed])

		return result



	def getGroomingMatrix(self,timeSeries=None):
		if timeSeries == None:
			timeSeries = self.theMatr;
		nAnimals = self.params['numAnimals'];
		gMatrix = np.zeros((nAnimals,nAnimals));
		for a1 in range(nAnimals):
			for a2 in range(1,nAnimals+1):
				a1GROOMSa2 = np.nonzero(timeSeries[:,a1]==a2)[0];
				gMatrix[a1,a2-1] = len(a1GROOMSa2);

		return gMatrix


	#returns a list of tuples of the form
	# (start,end,actor,type)
	def getSelfGroomingEvents(self,typesOfBehaviour = ['self','head grooming','gastertip']):
		numAnts = self.params['numAnimals'];
		allEvs = []
		for bt in typesOfBehaviour:
			if bt=='self':
				evs = self.getGroomingEvents();
				evs = [(s,e,a,bt) for (s,e,a,ga) in evs if ga==a]
			else:
				beVal = self.buttonNames.index(bt)+1
				evs = [];
				for an in range(numAnts):
					nzi = aSA.nonzero_intervals(self.theMatr[:,numAnts+an] == beVal)
					for evStart,evEnd in nzi:
						evs.append((evStart,evEnd,an,bt))
			allEvs = allEvs + evs
		return allEvs



	def getEntireTimeSeries(self):
		return self.theMatr

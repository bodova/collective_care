import numpy as np
import networkx as nx
import pydot
import glob

import matplotlib.pyplot as plt
import random as ran
import matplotlib.patches as patches
import math

# The following dishes were swapped under the cameras. We correct the mistake here to ensure that pre-treatment and post-treatment 
# obsevations correspond to the same group of ants. 
swapcorrections = {('HL',13):('HH',14),
				   ('HH',14):('HL',13),
				   ('LL',15):('LTx',15),
				   ('LTx',15):('LL',15) }

#For an ant with a given treatment, level and counts, in a video with given levels of R and G
#Normalizes the counts in the following way:
#If the ant is a nestmate, it normalzes its headG count by dividing it by the mean count of G in the heads of the exposure controls of level videoLevelG
#                         it normalzes its headR count by dividing it by the mean count of R in the heads of the exposure controls of level videoLevelR
#                         it normalzes its bodyR count by dividing it by the mean count of R in the bodies of the exposure controls of level videoLevelR
#                         it normalzes its bodyG count by dividing it by the mean count of G in the bodies of the exposure controls of level videoLevelG
#If the ant was exposed, it is first normalized the same way as nestmates, resulting in a normalized level x, then the final level is -(1-x)
#
#In the case of nestmates, this normalization indicates what percentage of the mean level of each spore in each bodypart they gained
#In the case of exposed, this normalization indicates what percentage of the mean leve of each spore in each bodypart they lost

def deltaSporeCount(treatment,level,headRFP,headGFP,bodyRFP,bodyGFP,videoLevelR,videoLevelG,useNegatives=True):
	maximalSporeCounts = { 'R': {'headL':13902.0,'headH':27155.0,'bodyL':92224.0,'bodyH':153539.0},
						   'G': {'headL':16490.0,'headH':23316.0,'bodyL':93996.0,'bodyH':152606.0} }

	factorHeadR = 1.0;
	factorBodyR = 1.0;
	factorHeadG = 1.0;
	factorBodyG = 1.0;

	if videoLevelR in ['L','H']:
		factorBodyR = maximalSporeCounts['R']['body'+videoLevelR]
		factorHeadR = maximalSporeCounts['R']['head'+videoLevelR]

	if videoLevelG in ['L','H']:
		factorBodyG = maximalSporeCounts['G']['body'+videoLevelG]
		factorHeadG = maximalSporeCounts['G']['head'+videoLevelG]

	if useNegatives:
		if treatment == 'R':
			headRFP = -(maximalSporeCounts['R']['head'+level]-headRFP);
			bodyRFP = -(maximalSporeCounts['R']['body'+level]-bodyRFP);
		if treatment == 'G':
			headGFP = -(maximalSporeCounts['G']['head'+level]-headGFP);
			bodyGFP = -(maximalSporeCounts['G']['body'+level]-bodyGFP);

	return headRFP/factorHeadR, headGFP/factorHeadG, bodyRFP/factorBodyR, bodyGFP/factorBodyG


def myColorMap(ax,wi=800,he=70,exp=1,useNegatives=True):
	dimToGrow = 0 if wi>he else 1;
	delta = float(2)/float(wi);
	ima  = np.zeros([wi,he,3])

	x = -1;
	pos = 0;
	yticpos = []
	yticlab = []
	while pos < wi:
		color1 =  floatToColor(x,'R',exp=exp)
		color2 =  floatToColor(x,'G',exp=exp)
		if dimToGrow == 0:
			ima[wi-pos-1,:he/2,:] = color1;
			ima[wi-pos-1,he/2:,:] = color2;
		if pos%100 == 0:
			yticpos.append(pos)
			if np.abs(x) < 0.000001:
				x=0
			yticlab.append(str(int(x*100))+"%")
		pos+=1
		x+=delta

	if not useNegatives:
		ima = np.flipud(ima)



	ax.imshow(ima)
	ax.set_xticks([])
	yticpos = yticpos+[wi-1]
	yticlab = yticlab+['100%']
	yticlab.reverse()
	if not useNegatives:
		ytp = [wi/8 + yy for yy in yticpos[int(len(yticpos)/2):]];

		ax.set_yticks(ytp)
		ytl = yticlab[:int(len(yticlab)/2)];
		ytl.reverse()
		ax.set_yticklabels(ytl)
		ax.set_ylim((wi/2,wi))
	else:
		ax.set_yticks(yticpos)
		ax.set_yticklabels(yticlab)
	#plt.yticklabels(yticlab)
	#print(plt.yticks())
	#plt.show()





def floatToColor(fl,col,exp=1):
	if fl > 1:
		fl = 1
	if fl < -1:
		fl = -1

	if fl < 0:
		fl = -math.pow(np.abs(fl),exp)
	else:
		fl = math.pow(np.abs(fl),exp)
	if col=='G':
		if fl > 0 and fl < 0.5:
			return (1-2*fl,1,1-2*fl  )#white(1,1,1) to green (0,1,0)
			#return (0,2*fl,0) #black to green
		if fl > 0 and fl >= 0.5:
			#return (2*fl-1),1,2*fl-1)   #green to white
			return  (0, 2*(1-fl),0)   #green(0,1,0) to black (0,0,0)
		else:
			return  (1,1+fl,1) #white to magenta
			#return (-fl,0,-fl)    #black (0,0,0) to magenta(1,0,1)
	else:
		if fl > 0 and fl < 0.5:
			#return (fl*2,0,0)    #black to red
			return (1,1-2*fl,1-2*fl)#white to red
		if fl > 0 and fl >= 0.5:
			#return (1,2*fl-1,2*fl-1)   #red to white
			return (2*(1-fl),0,0)#red to black
		else:
			#return (0,-fl,-fl)  #black to cyan
			return (1+fl,1,1)  #white to cyan

#counts is of the form:  headRed,headGreen,bodyRed,bodyGreen
def plotGraphWithPies(fileName,counts,treatments,gm,widthFactor = 10,withOutsideRings=True):
	numNodes = 6;
	#angles for self arrows. counterclockwise starting from East
	angles={0:[-(1/6.0)*np.pi,(1/6.0)*np.pi],
			1:[-(1/6.0)*np.pi,(1/6.0)*np.pi],
			5:[-(1/6.0)*np.pi,(1/6.0)*np.pi],
			2:[np.pi-(1/6.0)*np.pi,np.pi+(1/6.0)*np.pi],
			3:[np.pi-(1/6.0)*np.pi,np.pi+(1/6.0)*np.pi],
			4:[np.pi-(1/6.0)*np.pi,np.pi+(1/6.0)*np.pi]}

	G = nx.DiGraph();
	for g in range(numNodes):
		G.add_node(g);
		#for g2 in range(numNodes):
		#	G.add_edge(g,g2)
	posO=nx.circular_layout(G)
	nums = range(numNodes)
	nums = nums+nums
	nums.reverse()
	off = +3
	pos = dict()
	for a in posO.keys():
		pos[nums[a+off%numNodes]]=1.3*posO[a]-np.array([0.15,0.1])

	fig=plt.figure(figsize=(6,6),frameon = False)
	ax=plt.axes([0,0,1,1])
	ax.set_axis_off()
	ax.set_aspect('equal')

	plt.xlim(-0.5,1.5)
	plt.ylim(-0.5,1.5)

	trans=ax.transData.transform
	trans2=fig.transFigure.inverted().transform


	#The background

	piesize=0.25
	p2=piesize/2.0
	for n in G:
		xx,yy=trans(pos[n]) # figure coordinates
		xa,ya=trans2((xx,yy)) # axes coordinates
		#print("node: "+str(xa)+","+str(ya))

		a = plt.axes([xa-p2,ya-p2, piesize, piesize])
		a.set_aspect('equal')
		fracs = np.array([1])
		if withOutsideRings:
			thisColor = 'w'
			if 'R' in treatments[n]:
				thisColor = 'r'
			if 'G' in treatments[n]:
				thisColor = 'g'
			if 'T' in treatments[n] :
				thisColor = 'b'
		else:
			thisColor = '';
			if 'R' in treatments[n] or 'G' in treatments[n] or 'T' in treatments[n]:
				thisColor = treatments[n]

		if withOutsideRings:
			a.pie(100*fracs/fracs.sum(),colors=[thisColor])
		else:
			a.set_axis_off()
			a.text(0.5,1,thisColor)



	#The node's piecharts
	piesize=0.2 if withOutsideRings else 0.25
	p2=piesize/2.0
	for n in G:
		xx,yy=trans(pos[n]) # figure coordinates
		xa,ya=trans2((xx,yy)) # axes coordinates
		#print("node: "+str(xa)+","+str(ya))
		a = plt.axes([xa-p2,ya-p2, piesize, piesize])
		a.set_aspect('equal')
		#counts is of the form:  headRed,headGreen,bodyRed,bodyGreen
		#colors is headG headR bodyR bodyG
		colors=[(0,counts[n][1],0),(counts[n][0],0,0),(counts[n][2],0,0),(0,counts[n][3],0)]
		colors=[floatToColor(counts[n][1],'G') , floatToColor(counts[n][0],'R') , floatToColor(counts[n][2],'R'), floatToColor(counts[n][3],'G')]
		fracs = np.array([0.25 for i in range(4)])
		a.pie(100*fracs/fracs.sum(),colors=colors)

	piesize=0.25
	p2=piesize/2.0
	s = 1*piesize;
	alpha = np.pi/18;

	for ex in range(numNodes):
		for ey in range(numNodes):
			xx1,yy1=pos[ex] # figure coordinates
			xx2,yy2=pos[ey] # figure coordinates
			ew = gm[ex,ey]
			if ew==0:
				continue
			if (ex==ey):
				bx = nums[ex+off+1]
				alpha1 = angles[nums[bx]][0]
				alpha2 = angles[nums[bx]][1]
				ox = xx1+0.8*s*math.cos(alpha1);
				oy = yy1+0.8*s*math.sin(alpha1);
				fx = xx1+0.8*s*math.cos(alpha2);
				fy = yy1+0.8*s*math.sin(alpha2);

				pa = patches.FancyArrowPatch((ox,oy),(fx,fy),
				connectionstyle='arc3, rad=1.2',linewidth=widthFactor*ew,edgecolor='k',
				facecolor='k',arrowstyle='simple,head_width=8,head_length=8')

			else:
				dx = xx2-xx1;
				dy = yy2-yy1;
				nd = np.linalg.norm((dx,dy));
				x1 = xx1 + s* (dx)/nd;
				y1 = yy1 + s* (dy)/nd;
				x2 = xx1 + (nd-0.8*s)* (dx)/nd;
				y2 = yy1 + (nd-0.8*s)* (dy)/nd;
				tx = -(yy2-y2);
				ty = (xx2 - x2);
				nt = np.linalg.norm((tx,ty));
				b  = (s/2)*math.sin(alpha)
				fx = x2+b*tx/nt;
				fy = y2+b*ty/nt;
				ox = xx1+b*tx/nt
				oy = yy1+b*ty/nt

				pa = patches.FancyArrow(ox,oy,
						fx-ox,fy-oy,length_includes_head=True,linewidth=widthFactor*ew,edgecolor='k',facecolor='k')
			ax.add_patch(pa)

	if fileName == None:
		print("showing!")
		fig.show()
	else:
		fig.savefig(fileName, dpi=200	 )
		plt.close()



'''
Find islands of non-zeros in the vector vec
'''
def nonzero_intervals(vec):
	# Code released under creative commons license (http://stackoverflow.com/users/851699/peter) 
	# found in http://stackoverflow.com/a/27642744/3074835

    if len(vec)==0:
        return []
    elif not isinstance(vec, np.ndarray):
        vec = np.array(vec)

    edges, = np.nonzero(np.diff((vec==0)*1))
    edge_vec = [edges+1]
    if vec[0] != 0:
        edge_vec.insert(0, [0])
    if vec[-1] != 0:
        edge_vec.append([len(vec)])
    edges = np.concatenate(edge_vec)
    return zip(edges[::2], edges[1::2])
#
# For pre, use:
#   fileNamePRE = aSA.findSolomonFile(treatment,replicate,solomonDir,preCorrection=True);
def findSolomonFile(treatment,replicate,dir, preCorrection=False):
	if preCorrection and (treatment,replicate) in swapcorrections.keys():
		treatment,replicate = swapcorrections[(treatment,replicate)]
	if dir[-1] != '/':
		dir = dir+"/"
	query = dir+treatment+"_rep"+str(replicate)+"_*.arch"
	allFiles = glob.glob(query)
	if len(allFiles)==0:
		print(" !E: "+str(treatment)+" "+str(replicate)+" Not found in "+str(dir))
	return allFiles[0];

def findTimes(treatment,replicate,fileLocation):
	startContColumm		 = 6
	endContColumm   	 = 7
	startTreatmentColumn = 8;
	endTreatmentColumn	 = 9;
	fi = open(fileLocation,'r');
	numLines = 0;
	for line in fi:
		numLines+=1;
		#Skip the header
		if numLines==1:
			continue
		lstrp = line.strip();
		lsplt = lstrp.split(",")
		if lsplt[0]==treatment and int(lsplt[1])==replicate:
			initControl =  int(lsplt[startContColumm])
			endControl = int(lsplt[endContColumm])
			initTreatment = int(lsplt[startTreatmentColumn])
			endTreatment =int(lsplt[endTreatmentColumn])
			fixed_duration=15*60*15;
			#initControl=endControl-fixed_duration; #last15
			#initTreatment=endTreatment-fixed_duration; #last15
			#endControl=initControl+fixed_duration; #first15
			#endTreatment=initTreatment+fixed_duration; #first15
			return initControl, endControl, initTreatment, endTreatment

	print(" !E: "+str(treatment)+" "+str(replicate)+" Not found in "+str(fileLocation))


def makeGraph(adjMatrixO,duration,useLog=True):

	adjMatrix = adjMatrixO / duration;
	normAdjMatrix = (adjMatrix-adjMatrix.min())/(adjMatrix.max()-adjMatrix.min())
	G = nx.DiGraph();
	numAnimals = adjMatrix.shape[0];
	for a1 in range(numAnimals):
		G.add_node(a1)
		for a2 in range(numAnimals):
			Eweight = np.log10(adjMatrix[a1,a2]) if useLog else adjMatrix[a1,a2];
			if Eweight == -np.inf:
				Eweight = 0;
			if Eweight>0:
				G.add_edge(a1,a2,weight=Eweight);


	return G

def get_circular_layout(numNodes,scale):
	t=np.arange(0,2.0*np.pi,2.0*np.pi/numNodes,dtype=np.float32)
	pos=np.transpose(np.array([scale*np.cos(t),scale*np.sin(t)]))

	return dict(zip(range(numNodes),pos))

def makeGraphForPyDot(adjMatrixO,
					duration,
					pos,
					animalNames,
					nodeLabels,
					useLog=True,
					weightRatio = 5,
					posOffset=2,
					weightNorm = lambda x:x):

	adjMatrix = adjMatrixO / duration;
	numNodes = adjMatrix.shape[0]
	print("AA\n"+str(animalNames)+"\n")
	#print(pos)

	graph = pydot.Dot(graph_type='digraph')
	graph.set_ranksep(1)
	#add nodes
	nodeList = [];

	for a in range(numNodes):
		#Coloring the nodes:
		if nodeLabels[animalNames[a]][0:2]=='NM':
			colorThisNode='white';
		if nodeLabels[animalNames[a]][0]=='R':
			colorThisNode='red3'
			if 'L' in nodeLabels[animalNames[a]]:
				colorThisNode='tomato'
		if nodeLabels[animalNames[a]][0]=='G':
			colorThisNode='green4'
			if 'L' in nodeLabels[animalNames[a]]:
				colorThisNode='greenyellow'
		if nodeLabels[animalNames[a]][0]=='T':
				colorThisNode='skyblue'

		node = pydot.Node(animalNames[a],style='filled',fillcolor=colorThisNode);
		#this trick is just to make the first two nodes be ontop
		Apos = range(numNodes)[a-posOffset]
		print(animalNames[a])
		print(nodeLabels[animalNames[a]])
		print(str(Apos)+" "+str(int(pos[Apos].tolist()[0]))+","+str(int(pos[Apos].tolist()[1])))
		print(".")
		node.set_pos(int(pos[Apos].tolist()[0]),int(pos[Apos].tolist()[1]));
		node.set_label(nodeLabels[animalNames[a]]);
		node.set_width(0.1)
		#node.set_height()
		#node.set_fixedsize(1.2)
		graph.add_node(node);
		nodeList.append(node);

	for a1 in range(numNodes):
		for a2 in range(numNodes):
			Eweight = np.log10(adjMatrix[a1,a2]) if useLog else adjMatrix[a1,a2];
			if Eweight == -np.inf:
				Eweight = 0;
			if Eweight>=0:
				edge = pydot.Edge(nodeList[a1],nodeList[a2])
				edge.set_penwidth(weightRatio*Eweight)
				graph.add_edge(edge)

	return graph

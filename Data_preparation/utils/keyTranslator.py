
def toTwoLetter(di):
	v1 = list(di.values())[0];
	v2 = list(di.values())[1];

	code = [];
	for treat in ['H','L','Tx']:
		if treat in v1:
			code.append(treat[0]);
		if treat in v2:
			code.append(treat[0]);
	code.sort();
	return code[0]+code[1];

def toRoman(i):
	if i == 1:
		return "I"
	if i == 2:
		return "II"
	if i == 3:
		return "III"
	if i == 4:
		return "IV"
	if i == 5:
		return "V"
	if i == 6:
		return "VI"

def onlySporeOrTx(s):
	if 'g' in s:
		return 'G'
	if 'r' in s:
		return 'R'
	if 'T' in s:
		return 'T'

	return s;

def onlyLevel(s):
	if 'H' in s:
		return 'H'
	if 'L' in s:
		return 'L'
	return 'X';

def orderByTreatmentLevel(a):
	sortIndex=['RH','GH','RL','GL','TX','NX']
	return sorted(range(len(a)),key=lambda x:sortIndex.index(a[x]))

#Replicate is a POSITVE integer
#Treatent is one of "HL" "LTx" "HTx" "HH" "LL" "TxTx"
#Returns a dictionary of the form 'color':'treatment' where treatment can be either one of the above or N for none
def kTrans(replicate,treatmentO):
	d = [[{"purple":"Hr", "green":"gH"}, {"yellow":"Hr", "red":"Tx"}, {"orange":"Lr", "purple":"gH"}, {"red":"gL", "yellow":"Tx"}, {"yellow":"Tx", "purple":"Tx"}, {"blue":"Lr", "orange":"gL"}], [{"blue":"gH", "orange":"Tx"}, {"red":"Lr", "blue":"gH"}, {"yellow":"Tx", "purple":"Tx"}, {"green":"Lr", "blue":"Tx"}, {"yellow":"gL", "blue":"Lr"}, {"orange":"gH", "purple":"Hr"}], [{"orange":"gL", "yellow":"Tx"}, {"purple":"Lr", "green":"gL"}, {"yellow":"gH", "blue":"Hr"}, {"orange":"gL", "red":"Hr"}, {"orange":"Tx", "blue":"Tx"}, {"red":"Hr", "green":"Tx"}], [{"orange":"gH", "yellow":"Tx"}, {"yellow":"gH", "orange":"Hr"}, {"green":"gL", "blue":"Lr"}, {"yellow":"gL", "blue":"Tx"}, {"green":"Tx", "purple":"Tx"}, {"blue":"gL", "purple":"Hr"}], [{"yellow":"Lr", "purple":"gH"}, {"blue":"gL", "red":"Lr"}, {"blue":"Hr", "green":"gH"}, {"red":"Hr", "green":"Tx"}, {"green":"Tx", "blue":"Tx"}, {"blue":"Lr", "purple":"Tx"}], [{"green":"Tx", "red":"Tx"}, {"purple":"gL", "green":"Lr"}, {"yellow":"gH", "blue":"Hr"}, {"red":"Lr", "green":"Tx"}, {"orange":"Hr", "green":"Tx"}, {"blue":"gL", "yellow":"Hr"}], [{"orange":"Lr", "blue":"Tx"}, {"yellow":"Hr", "purple":"Tx"}, {"yellow":"gH", "green":"Hr"}, {"purple":"Lr", "blue":"gL"}, {"purple":"Tx", "green":"Tx"}, {"green":"Lr", "blue":"gH"}], [{"yellow":"Tx", "purple":"Tx"}, {"yellow":"Lr", "red":"gH"}, {"purple":"Lr", "red":"gL"}, {"orange":"gH", "red":"Tx"}, {"red":"Hr", "orange":"gH"}, {"green":"Lr", "orange":"Tx"}], [{"purple":"gL", "yellow":"Tx"}, {"red":"Lr", "orange":"gH"}, {"green":"Lr", "red":"gL"}, {"purple":"Hr", "yellow":"gH"}, {"purple":"gH", "orange":"Tx"}, {"blue":"Tx", "green":"Tx"}], [{"red":"gL", "blue":"Tx"}, {"purple":"Tx", "orange":"Tx"}, {"orange":"gL", "yellow":"Lr"}, {"yellow":"gH", "green":"Tx"}, {"blue":"gL", "purple":"Hr"}, {"blue":"Hr", "green":"gH"}], [{"orange":"Lr", "yellow":"Tx"}, {"blue":"Tx", "yellow":"Tx"}, {"red":"gL", "orange":"Hr"}, {"yellow":"Hr", "blue":"Tx"}, {"yellow":"Lr", "orange":"gL"}, {"blue":"Hr", "orange":"gH"}], [{"red":"Hr", "purple":"Tx"}, {"green":"Tx", "purple":"Tx"}, {"green":"gH", "blue":"Hr"}, {"blue":"Lr", "purple":"gL"}, {"yellow":"gL", "blue":"Hr"}, {"blue":"Lr", "purple":"Tx"}], [{"orange":"gL", "blue":"Tx"}, {"red":"Lr", "green":"gL"}, {"yellow":"gH", "red":"Tx"}, {"purple":"Tx", "orange":"Tx"}, {"orange":"Hr", "red":"gH"}, {"blue":"gL", "red":"Hr"}], [{"yellow":"Hr", "blue":"gH"}, {"yellow":"gH", "purple":"Tx"}, {"red":"Tx", "yellow":"Tx"}, {"red":"gL", "purple":"Tx"}, {"green":"gL", "orange":"Lr"}, {"blue":"gL", "green":"Hr"}], [{"red":"gL", "yellow":"Tx"}, {"purple":"Lr", "green":"gL"}, {"blue":"Tx", "orange":"Tx"}, {"yellow":"gH", "red":"Tx"}, {"yellow":"Hr", "purple":"gH"}, {"purple":"Lr", "yellow":"gH"}], [{"red":"Tx", "yellow":"Tx"}, {"yellow":"gL", "blue":"Tx"}, {"yellow":"Lr", "blue":"gL"}, {"orange":"Hr", "red":"gH"}, {"green":"Lr", "orange":"gH"}, {"green":"Hr", "blue":"Tx"}], [{"purple":"Hr", "green":"Tx"}, {"green":"gL", "blue":"Hr"}, {"purple":"gL", "red":"Tx"}, {"orange":"Lr", "purple":"gL"}, {"orange":"Hg", "purple":"Hr"}, {"purple":"Tx", "red":"Tx"} ], [{"purple":"gH", "red":"Tx"}, {"orange":"Hr", "purple":"gH"}, {"orange":"gL", "purple":"Lr"}, {"purple":"Tx", "blue":"Tx"}, {"yellow":"Lr", "orange":"gH"}, {"yellow":"gL", "orange":"Tx"}]];

	#We sort treatment alfabetically
	treatmentSorted = [xx  for xx in treatmentO]
	treatmentSorted.sort();
	treatment = ''.join(treatmentSorted)[:2];

	thisRep = d[replicate-1];
	twoLetThisRep = [toTwoLetter(xx) for xx in thisRep]
	indTreatment = twoLetThisRep.index(treatment);
	modif = thisRep[indTreatment]; #is a dictionary that contains two pairs of color:treatmentLevel, corresponding to two ants treated in each petri

	resultTreatment = {'yellow':'N',
 				'orange':'N',
 				'red':'N',
 				'purple':'N',
 				'blue':'N',
 				'green':'N'}

	resultTreatment[list(modif.keys())[0]] = onlySporeOrTx(list(modif.values())[0])
	resultTreatment[list(modif.keys())[1]] = onlySporeOrTx(list(modif.values())[1])

	resultLevel = {'yellow':'X',
 				'orange':'X',
 				'red':'X',
 				'purple':'X',
 				'blue':'X',
 				'green':'X'}

	resultLevel[list(modif.keys())[0]] = onlyLevel(list(modif.values())[0])
	resultLevel[list(modif.keys())[1]] = onlyLevel(list(modif.values())[1])

	dishNumber = toRoman(indTreatment+1);

	return resultTreatment,resultLevel,dishNumber

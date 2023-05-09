import numpy as np
from matplotlib import pylab as P
import pickle
import itertools
import scipy.ndimage as ndimage
import os

import utils.keyTranslator as kT
import utils.auxSolomonAnalysis as aSA
import utils.reader as r

behaviourTypes = ["groomIn", "groomOut", "groomTot"] #"groomedByAtLeastOne", "acidopore"]


keysForOutput = ['color', 'treatment', 'level', 'headRFP', 'headGFP', 'bodyRFP', 'bodyGFP', 'gaster', 'other',
                 'antenna', 'head', 'selfG', 'groomIn', 'groomOut', 'totGroomIn', 'totGroomOut', 'groomToR', 'groomToG',
                 'groomToT', 'rankGroomInGeneral', 'rankGroomInParticular', 'rankGroomOutGeneral',
                 'rankGroomOutParticular', 'groomedByAtLeastOne'];
headOfWindowedOutput = ['color', 'treatment', 'level', 'headRFP', 'headGFP', 'bodyRFP', 'bodyGFP']


class TableBuilder():

    def __init__(self,
                 solomonDir,
                 outFolder,
                 ddpcrFilePath,
                 metadataFilePath,
                 fps=15,
                 symmetricTrophallaxis = False):
        self.solomonDir = solomonDir
        self.outFolder = outFolder
        self.ddpcrFilePath = ddpcrFilePath
        self.metadataFilePath = metadataFilePath
        self.fps=fps
        self.symmetricTrophallaxis = symmetricTrophallaxis

        if self.outFolder[-1] != '/':
            self.outFolder += '/';

    # behaviourType can be groomIn, groomOut, groomTot
    def buildWindowedTables(self, treatmentReplicatePairs ,
                            behaviourType='groomIn',
                            windowSizeInFrames=None,
                            separator="," ,
                            eventStartCount=False,
                            eventEnds=False):
        if windowSizeInFrames is None:
            windowSizeInFrames=30 * self.fps
        countp = ""
        if eventStartCount:
            countp = "_EventStartCounts"
        if eventEnds:
            countp = "_EventEndCounts"
        fileSufix = "_windowed_" + str(windowSizeInFrames) + "_" + behaviourType + countp +".csv"
        controlOutFile = open(self.outFolder + "PRE" + fileSufix, "wt")
        treatmentOutFile = open(self.outFolder + "POST" + fileSufix, "wt")
        cDs = []
        tDs = []
        for treatment, replicate in treatmentReplicatePairs:
            print(treatment, replicate, " : ")
            controlDict, treatmentDict = self.buildWindowedTimeSeries(treatment,
                                                                replicate,
                                                                behaviourType=behaviourType,
                                                                windowSizeInFrames=windowSizeInFrames,
                                                                eventCount=eventStartCount,
                                                                eventEnds=eventEnds)
            antColors = controlDict['colors']
            controlDict['treatment_replicate_period'] = (treatment,replicate,'PRE')
            treatmentDict['treatment_replicate_period'] = (treatment, replicate,'POST')
            cDs.append(controlDict)
            tDs.append(tDs)
            for color in antColors:
                for outFile, inDict in [(controlOutFile, controlDict), (treatmentOutFile, treatmentDict)]:
                    antDescription = inDict[color]['antDescription']
                    antSpores = inDict[color]['antSpores']
                    antTimeSeries = inDict[color]['antTimeSeries']
                    antRow = [treatment, replicate] \
                             + antDescription \
                             + antSpores.tolist() \
                             + antTimeSeries.tolist()
                    outFile.write(separator.join([str(x) for x in antRow]))
                    outFile.write("\n")
            controlOutFile.flush()
            treatmentOutFile.flush()
        controlOutFile.close()
        treatmentOutFile.close()

        return cDs,tDs


    def buildAllTables(self,treatmentReplicatePairs, writePickle=False):
        bothDicts = []
        for useGorT in ['G']: #, 'T']:
            wholeDict = dict();
            fileNames = [self.outFolder + "header.csv"]
            hf = open(self.outFolder + "header.csv", 'w')
            hf.write("treatmentgroup" + "," + "replicate" + "," + "period" + ",")
            for aa in keysForOutput:
                hf.write(str(aa) + ",")
            hf.write("\n")
            hf.close()

            for treatment, replicate in treatmentReplicatePairs:
                print(treatment, replicate, " : ")
                controlDict, treatmentDict, outFiles = self.buildTable(treatment, replicate,  useGorT=useGorT)
                wholeDict[(treatment, replicate, 'PRE')] = controlDict;
                wholeDict[(treatment, replicate, 'POST')] = treatmentDict;
                if useGorT == 'T':
                    print(str(treatment) + " " + str(replicate) + " Done")
                fileNames = fileNames + outFiles

            if self.outFolder != None:
                if writePickle:
                    wholeDictFileName = self.outFolder + 'wholeDict_' + useGorT + '.pkl'
                    pickle.dump(wholeDict, open(wholeDictFileName, 'wb'));
                print(fileNames)
                catCommand = 'cat '
                rmCommand = 'rm '
                for fn in fileNames:
                    catCommand += fn + " "
                    rmCommand += fn + " "
                catCommand += " >> " + self.outFolder + useGorT + "_ALL.csv"
                os.system(catCommand)
                os.system(rmCommand)

            bothDicts.append(wholeDict)

        return bothDicts


    # rnum, gnum and tnum store the 0-based indices of the ants that have read, green or tritonx treatment, respectively, or [] if no ant has such treatment.
    def findNumbersOfRedAndGreenAnts(self,colorTreatmentDictionary, antColors):
        antTreatments = [colorTreatmentDictionary[xx] for xx in antColors[:len(colorTreatmentDictionary)]];

        rnum = [i for i, at in enumerate(antTreatments) if at == 'R']
        gnum = [i for i, at in enumerate(antTreatments) if at == 'G']
        tnum = [i for i, at in enumerate(antTreatments) if at == 'T']

        return rnum, gnum, tnum


    # behaviour can be:  groomIn, groomOut, groomTot, groomedByAtLeastOne, or acidopore
    def buildWindowedVectorFromSolomon(self,A, duration, color, colorTreatmentDictionary, outName, durationLimitsDict, useGorT='G',
                        behaviour='groomTot', windowSizeInFrames=None, eventCount=False, eventEnds=False, TotIncludesHead=False):
        if windowSizeInFrames is None:
            windowSizeInFrames=30 * self.fps
        if eventCount or eventEnds:
            return self.buildWindowedEventVectorFromSolomon(A,duration,color,colorTreatmentDictionary,outName,durationLimitsDict,
                                                       useGorT=useGorT,behaviour=behaviour,windowSizeInFrames=windowSizeInFrames,
                                                       eventEnds=eventEnds)
        antColors = A.buttonNames;
        rnum, gnum, tnum = self.findNumbersOfRedAndGreenAnts(colorTreatmentDictionary, antColors);
        numAnts = A.params['numAnimals'];

        colGrooming = antColors.index(color);
        indThisAnt = colGrooming + 1;

        timeLimit = durationLimitsDict[outName] if outName in durationLimitsDict.keys() else 30 * 60 * self.fps
        theMatr = A.theMatr[:min(timeLimit, A.theMatr.shape[0]), :]
        numWindows = int(timeLimit / windowSizeInFrames)
        bins = [i * windowSizeInFrames for i in range(numWindows + 1)]

        # groomOut = allogrooming performed towards other ants (i.e. excluding body selfgrooming)
        # groomIn  = allogrooming received from other ants (i.e excluding body selfgrooming)
        # groomTot = allogrooming received + body selfgrooming

        if behaviour == 'groomOut':
            groomOutCol = theMatr[:, colGrooming]
            y = np.zeros_like(groomOutCol)
            y[groomOutCol > 0] += 1;
            y[groomOutCol != indThisAnt] += 1
            framesToBin = np.nonzero(y == 2)[0]
        if behaviour in ['groomIn','groomTot']:
            allIndices = list(range(numAnts))
            if behaviour == 'groomIn':
                allIndices.remove(colGrooming)
            groomInMatrix = theMatr[:, allIndices]            
            framesToBin = np.nonzero(groomInMatrix == indThisAnt)[0]            
            if TotIncludesHead:
                #print("A.buttonames=",A.buttonNames)
                action_value = A.buttonNames.index('head grooming') + 1
                col_head = colGrooming + numAnts
                framesToBin2 = np.nonzero(theMatr[:,col_head]==action_value)[0]
                #print(f"{color}\t1: {framesToBin.shape}  2:{framesToBin2.shape}")
                framesToBin = np.hstack((framesToBin,framesToBin2))
        if behaviour == 'groomedByAtLeastOne':
            framesToBin = np.nonzero(self.find_groomIn_collapsed(A, color))[0]
        if behaviour == 'acidopore':
            action_value = A.buttonNames.index('gastertip') + 1
            col_acidopore = colGrooming + numAnts
            framesToBin = np.nonzero(theMatr[:,col_acidopore]==action_value)

        theVector, theBins = np.histogram(framesToBin, bins=bins)
        return theVector



    # behaviour can be:  groomIn, groomOut, groomTot, groomedByAtLeastOne, or acidopore
    def buildWindowedEventVectorFromSolomon(self, A, duration, color, colorTreatmentDictionary, outName, durationLimitsDict,
                                       useGorT='G',
                                       behaviour='groomTot', windowSizeInFrames=None, eventEnds=False):
        if windowSizeInFrames is None:
            windowSizeInFrames=30 * self.fps
        antColors = A.buttonNames;
        rnum, gnum, tnum = self.findNumbersOfRedAndGreenAnts(colorTreatmentDictionary, antColors);
        numAnts = A.params['numAnimals'];

        colGrooming = antColors.index(color);
        indThisAnt = colGrooming + 1;

        timeLimit = durationLimitsDict[outName] if outName in durationLimitsDict.keys() else 30 * 60 * self.fps
        theMatr = A.theMatr[:min(timeLimit, A.theMatr.shape[0]), :]
        numWindows = int(timeLimit / windowSizeInFrames)
        bins = [i * windowSizeInFrames for i in range(numWindows + 1)]

        if behaviour == 'groomOut':
            groomOutCol = theMatr[:, colGrooming]
            y = np.zeros_like(groomOutCol)
            y[groomOutCol > 0] += 1;
            y[groomOutCol != indThisAnt] += 1
            framesToBin = np.nonzero(y == 2)[0]
        if behaviour == 'groomedByAtLeastOne':
            framesToBin = np.nonzero(self.find_groomIn_collapsed(A, color))[0]
        if behaviour == 'acidopore':
            action_value = A.buttonNames.index('gastertip') + 1
            col_acidopore = colGrooming + numAnts
            framesToBin = np.nonzero(theMatr[:,col_acidopore]==action_value)[0]

        if behaviour in ['groomIn','groomTot']:
            allIndices = list(range(numAnts))
            if behaviour == 'groomIn':
                allIndices.remove(colGrooming)
            groomInMatrix = theMatr[:, allIndices]
            evmatrix = np.zeros_like(groomInMatrix)
            evmatrix[groomInMatrix==indThisAnt] = 1
            tmpmatrix = np.zeros_like(evmatrix)
            if eventEnds:
                framesToBin = np.nonzero((evmatrix[:-1,:]-evmatrix[1:,:])==1)[0]
            else:
                framesToBin = np.nonzero((evmatrix[1:,:]-evmatrix[:-1,:])==1)[0]

        if len(framesToBin)>0:
            #print(behaviour,framesToBin.shape)
            tmpvect = np.zeros(framesToBin.max()+1)
            tmpvect[framesToBin] = 1
            if eventEnds:
                framesToBin = np.nonzero((tmpvect[:-1]-tmpvect[1:])==1)[0]
            else:
                framesToBin = np.nonzero((tmpvect[1:]-tmpvect[:-1])==1)[0]

        theVector, theBins = np.histogram(framesToBin, bins=bins)
        return theVector

    def find_groomIn_collapsed(self, A, color_reciever):
        """

        :param A: a solomonReader object
        :param color_reciever:  a color of the ant of interest e.g. 'orange'
        :return:  A binary vector, with as many entries a frames in the video. The j'th entry is 1 iff the ant whose color
        is color_reciever is recieveing grooming by any other ant at frame j.
        """
        numAnts = A.params['numAnimals']
        antColors = A.buttonNames
        col_self_grooming = antColors.index(color_reciever)
        other_ants = list(range(numAnts))[:col_self_grooming] + list(range(numAnts))[col_self_grooming + 1:]
        groom_events = A.theMatr[:, other_ants] == col_self_grooming + 1
        sum_groom = groom_events.sum(axis=1)
        return sum_groom >= 1


    def buildVectorFromSolomon(self, A, duration, color, colorTreatmentDictionary, useGorT='G', outName=None):
        antColors = A.buttonNames;
        rnum, gnum, tnum = self.findNumbersOfRedAndGreenAnts(colorTreatmentDictionary, antColors);
        numAnts = A.params['numAnimals'];

        colGrooming = antColors.index(color);
        colSelfs = numAnts + colGrooming;

        # print("antColors:",antColors,
        # 		"color:",color,
        # 		"colGrooming",colGrooming,
        # 		"treatment:",colorTreatmentDictionary[color])

        if useGorT == 'G':
            gm = A.getGroomingMatrix() / duration;
        else:
            gm = A.getTrophallaxisMatrix() / duration;
        # print(antColors)
        # print("\n\n")

        gasterCode = antColors.index('gastertip') + 1;
        otherCode = antColors.index('other') + 1;
        antennaCode = antColors.index('antenna stroke') + 1;
        headCode = antColors.index('head grooming') + 1;

        gaster = len(np.nonzero(A.theMatr[:, colSelfs] == gasterCode)[0]) / duration;
        other = len(np.nonzero(A.theMatr[:, colSelfs] == otherCode)[0]) / duration;
        antenna = len(np.nonzero(A.theMatr[:, colSelfs] == antennaCode)[0]) / duration;
        head = len(np.nonzero(A.theMatr[:, colSelfs] == headCode)[0]) / duration;
        selfG = gm[colGrooming, colGrooming];

        invdic = {v: k for k, v in colorTreatmentDictionary.items()}

        if 'R' in colorTreatmentDictionary.values():
            redAntColor = invdic['R']
            redAntIndex = antColors.index(redAntColor);
            # print("redAntIndex "+str(redAntIndex))
            groomToR = gm[colGrooming, redAntIndex];
        else:
            groomToR = 0;

        if 'G' in colorTreatmentDictionary.values():
            greenAntColor = invdic['G']
            greenAntIndex = antColors.index(greenAntColor)
            # print("greenAntIndex "+str(greenAntIndex))
            groomToG = gm[colGrooming, greenAntIndex];
        else:
            groomToG = 0;

        if 'T' in colorTreatmentDictionary.values():
            tritonAntColor = invdic['T']
            tritonAntIndex = antColors.index(tritonAntColor);
            # print("tritonAntIndex "+str(tritronAntIndex))
            groomToT = gm[colGrooming, tritonAntIndex];
        else:
            groomToT = 0;

        totGroomIn = gm[:, colGrooming].sum();
        totGroomOut = gm[colGrooming, :].sum();
        groomIn = totGroomIn - selfG
        groomOut = totGroomOut - selfG

        allIndices = set(range(len(colorTreatmentDictionary)))
        treatedIndices = set(rnum + gnum + tnum)
        untreatedIndices = allIndices - treatedIndices
        if colorTreatmentDictionary[color] in ['R', 'T', 'G']:
            particularIndices = list(treatedIndices)
        else:
            particularIndices = list(untreatedIndices)

        # The rank of an ant is only with respect to grooming recieved and performed to other ants, no self grooming must be counted.
        gm_noSelf = gm.copy()
        for i in range(gm_noSelf.shape[0]):
            gm_noSelf[i, i] = 0
        sumPerformed = gm_noSelf.sum(axis=1)
        sumRecieved = gm_noSelf.sum(axis=0)
        ranksIn = list(np.argsort(sumRecieved))
        rankGroomInGeneral = len(ranksIn) - ranksIn.index(colGrooming)  # make it 1-based
        groomedByAtLeastOne = self.find_groomIn_collapsed(A, color).sum() / float(duration)

        ranksOut = list(np.argsort(sumPerformed))
        rankGroomOutGeneral = len(ranksOut) - ranksOut.index(colGrooming)
        # print(colGrooming, "treated:", treatedIndices,
        # "untreated:",untreatedIndices,
        # "particular:",particularIndices)
        if len(particularIndices) > 0:
            thisInParticular = particularIndices.index(colGrooming)

            ranksIn = list(np.argsort(sumRecieved[particularIndices]))
            rankGroomInParticular = len(ranksIn) - ranksIn.index(thisInParticular)

            ranksOut = list(np.argsort(sumPerformed[particularIndices]))
            rankGroomOutParticular = len(ranksOut) - ranksOut.index(thisInParticular)
        else:
            rankGroomInParticular = -1
            rankGroomOutParticular = -1

        theVector = [gaster, other, antenna, head, selfG, groomIn, groomOut, totGroomIn, totGroomOut, groomToR, groomToG,
                     groomToT, rankGroomInGeneral, rankGroomInParticular, rankGroomOutGeneral, rankGroomOutParticular,
                     groomedByAtLeastOne]

        return np.array(theVector);


    def getHeadAndBodyCounts(self, replicate, dish, color, useLog=False):
        replicateColumn = 0;
        dishColumn = 1;
        colorColumn = 2;
        partColumn = 3;
        greenCountColumn = 5;  
        redCountColumn = 4;  

        separator = ","

        colsOfInterest = [replicateColumn, dishColumn, colorColumn, partColumn, greenCountColumn, redCountColumn];

        headRFP = None;
        headGFP = None;
        bodyRFP = None;
        bodyGFP = None;

        lineNum = 0;
        fin = open(self.ddpcrFilePath, "r");
        # print("\tOPENING FILE:  "+self.ddpcrFilePath+ddpcrFileName+" looking for "+dish+" "+str(replicate)+" "+color+"\n")

        for line in fin:
            lineNum += 1;
            if lineNum < 2:
                continue
            lstrp = line.strip();
            lsplt = lstrp.split(separator);
            if len(lsplt) < max(colsOfInterest) + 1:
                print("breaking read of ddpcr data in line " + str(lineNum) + " because it is malformed")
                break
            if all([lsplt[x] != "" for x in colsOfInterest]):
                replicateVal = int(lsplt[replicateColumn]);
                dishVal = lsplt[dishColumn];
                colorVal = lsplt[colorColumn];

                if not (dishVal == dish and replicateVal == replicate and colorVal == color):
                    continue

                if lsplt[partColumn] == "head":
                    headRFP = np.nan if lsplt[redCountColumn] == '-' else int(lsplt[redCountColumn]);
                    headGFP = np.nan if lsplt[greenCountColumn] == '-' else int(lsplt[greenCountColumn]);
                if lsplt[partColumn] == "body":
                    bodyRFP = np.nan if lsplt[redCountColumn] == '-' else int(lsplt[redCountColumn]);
                    bodyGFP = np.nan if lsplt[greenCountColumn] == '-' else int(lsplt[greenCountColumn]);

        fin.close()
        if useLog:
            return np.log10(headRFP + 1), np.log10(headGFP + 1), np.log10(bodyRFP + 1), np.log10(bodyGFP + 1)
        else:
            return headRFP, headGFP, bodyRFP, bodyGFP


    # Replicate is a POSITVE integer
    # Treatent is one of "HL" "LTx" "HTx" "HH" "LL" "TxTx"

    '''
    Output columns:
    	treatmentgroup
    	replicate
    	period
    	color
    	treatment   	{R,G,T,N}
    	level			{H,L,X}
    	headRFP			fromddpcr
    	headGFP			fromddpcr
    	bodyRFP			fromddpcr
    	bodyGFP			fromddpcr
    	gaster			fromSolomon
    	other			fromSolomon
    	antenna			fromSolomon
    	head			fromSolomon
    	selfG			fromSolomon
    	totGroomIn		fromSolomon
    	totGroomOut		fromSolomon
    	groomToR		fromSolomon
    	groomToG		fromSolomon
    	groomToT		fromSolomon
    	rankGroomInGeneral
    	rankGroomInParticular
    	rankGroomOutGeneral
    	rankGroomOutParticular
    '''


    # This filter will be used for dilations
    def myFilter(X, y):
        if X.sum() < y:
            return 0
        a = X != 0;
        if X[a].sum() == y * a.sum():
            return y
        return 0


    # Writes a csv file with 3 x numAnts columns. Each row is a frame in the video
    #
    #     The first numAnts entries are interpreted as follows:
    #        if the i'th entry is j, then it means ant i+1 is grooming ant j
    #     The next numAnts entries are interepreted as follows:
    #        if the (i-numAnts)'th entry is y it means ant i+1 is doing action y
    #	  The next numAnts entries are interpreted as follows:
    # 		 if the (i-2*numAnts)'th entry is j, it means ant i+1 and ant j are doing trophallaxis

    #

    def buildWholeTimeSeries(self, dicto, treatment, replicate, prefix, dilateRadius=-1):
        # header = "yellow(1)Groom orange(2)Groom red(3)Groom purple(4)Groom blue(5)Groom green(6)Groom gastertipGroom rest other antennaStroke headGrooming yellowTro orangeTro redTro purpleTro blueTro greenTro"

        header = "yG oG rG pG bG gG yA oA rA pA bA gA yT oT rT pT bT gT"
        # where middle columns can have entries {7-11}, {R,E,O,A,H}

        A = dicto['A']
        colorTreatmentDictionary, colorLevelDictionary, dishNumber = kT.kTrans(replicate, treatment);
        outName = self.outFolder + "TIMESERIES_" + prefix + "_" + treatment + "_" + str(replicate) + ".csv"
        np.savetxt(outName, A.theMatr, fmt='%d', header=header, comments="% ")


    # A is a reader object
    # Write in self.outFolder a file with all the events that are recorder in the reader object. An event is a line of the format:
    #  start,duration,actor,treatmentActor,levelActor,typeOfBehaviour,reciever,treatmentReciever,levelReciever
    # the actor and receiver are one of the ant colors ('yellow', 'orange', ...)
    # type is one '
    # 'G':  grooming
    # 'T':  trophallaxis
    # 'S':  selfgrooming
    # 'R':  gastergrooming
    # 'A':  antenna stroke
    # 'H':  headgrooming
    # 'O':  other
    # 'E':  rest
    def buildEventRecords(self, dicto,  treatment, replicate, prefix, dilateRadius=-1):
        A = dicto['A']
        colorTreatmentDictionary, colorLevelDictionary, dishNumber = kT.kTrans(replicate, treatment);
        outName = self.outFolder + "EVENTS_" + prefix + "_" + treatment + "_" + str(replicate) + ".csv"
        oF = open(outName, 'wt')
        antColors = A.buttonNames[:A.numAnts];
        for col in range(A.theMatr.shape[1]):
            column = A.theMatr[:, col]
            newCol = np.zeros_like(column);
            if dilateRadius > 0:
                st = np.ones(dilateRadius)
                # For every value in the column, that is <7 (i.e. only grooming or trophallaxis), we dilate it, so that "holes" in the column of size smaller than dilateRadius are "filled".
                for v in [vu for vu in np.unique(column).tolist() if vu != 0 and vu < 7]:
                    X = ndimage.generic_filter(column,
                                               myFilter,
                                               dilateRadius,
                                               mode='constant',
                                               extra_arguments=tuple([v]));
                    mask = np.zeros_like(X);
                    mask[X != 0] = 1;
                    nm = ndimage.binary_erosion(mask, structure=st);
                    X[nm == 0] = 0;
                    newCol += X;
                newCol[np.nonzero(column)] = column[np.nonzero(column)]

                if newCol.max() != column.max():
                    print("!!!!!!!!1")
                column = newCol;
            streaks = self.findConsecutiveStreaks(column.tolist())
            for (val, start, dur) in streaks:
                actor, reciever, typeOfBehaviour = self.singleEventRecord(col, val, antColors);
                treatmentActor = colorTreatmentDictionary[actor];
                levelActor = colorLevelDictionary[actor];
                treatmentReciever = colorTreatmentDictionary[reciever];
                levelReciever = colorLevelDictionary[reciever]
                toWrite = ",".join([str(x) for x in
                                    [treatment, replicate, start, dur, actor, treatmentActor, levelActor, typeOfBehaviour,
                                     reciever, treatmentReciever, levelReciever]]) + "\n"
                oF.write(toWrite);

        oF.close()


    # Given a vector, returns a list triplets of the form (value,start,duration).
    def findConsecutiveStreaks(self,a):
        V = itertools.groupby(a)
        vals = [x for x, y in itertools.groupby(a)]
        lens = [len(list(y)) for x, y in itertools.groupby(a)]
        starts = np.hstack(([0], np.cumsum(lens)))

        return [(vals[i], starts[i], starts[i + 1] - starts[i]) for i in range(len(vals)) if vals[i] != 0]


    def singleEventRecord(self,colNumber, value, antColors):
        actionTranslator = {7: 'R', 8: 'E', 9: 'O', 10: 'A', 11: 'H'}
        numAnts = len(antColors)
        # First  numAnts Columns (Grooming)
        if colNumber < numAnts:
            actor = antColors[colNumber]
            reciever = antColors[value - 1]
            typeOfBehaviour = 'G'
            if colNumber + 1 == value:
                typeOfBehaviour = 'S'

        # Last numAnts Columns (trophallaxis)
        if colNumber >= 2 * numAnts:
            typeOfBehaviour = 'T'
            actor = antColors[colNumber - 2 * numAnts]
            reciever = antColors[value - 1]

        # Middle numAnts columns (actions)
        if colNumber >= numAnts and colNumber < 2 * numAnts:
            actor = antColors[colNumber - numAnts];
            reciever = actor;
            typeOfBehaviour = actionTranslator[value]

        return actor, reciever, typeOfBehaviour


    def buildTable(self, treatment, replicate, useGorT='G', buildTimeSeries=False):
        controlDict, treatmentDict = self.buildDictionaries(treatment, replicate, useGorT=useGorT);
        antColors = controlDict['colors'];

        # -------- Write the event records ---
        if useGorT == 'G':
            self.buildEventRecords(controlDict,  treatment, replicate, 'PRE')
            self.buildEventRecords(treatmentDict,  treatment, replicate, 'POST')
            if buildTimeSeries:
                self.buildWholeTimeSeries(controlDict,  treatment, replicate, 'PRE')
                self.buildWholeTimeSeries(treatmentDict,  treatment, replicate, 'POST')

        outFiles = []

        if self.outFolder == None:
            controlOutFile = '/dev/stdout';
            treatmentOutFile = '/dev/stdout';
        else:
            if self.outFolder[-1] != '/':
                self.outFolder += '/';
            controlOutFile = self.outFolder + treatment + "_" + str(replicate) + "_PRE" + useGorT + "_wRanks.csv"
            treatmentOutFile = self.outFolder + treatment + "_" + str(replicate) + "_POST" + useGorT + "_wRanks.csv"
            outFiles = [controlOutFile, treatmentOutFile]

        cf = open(controlOutFile, 'w');
        tf = open(treatmentOutFile, 'w');

        # Control (pre)
        for color in antColors:
            cf.write(str(treatment) + "," + str(replicate) + ",PRE,")
            for key in keysForOutput:
                cf.write(str(controlDict[color][key]) + ",")
            cf.write("\n")
        cf.close()

        # Treatment

        for color in antColors:
            tf.write(str(treatment) + "," + str(replicate) + ",POST,")
            for key in keysForOutput:
                tf.write(str(treatmentDict[color][key]) + ",")
            tf.write("\n")
        tf.close()

        return controlDict, treatmentDict, outFiles


    def buildWindowedTimeSeries(self,treatment,
                                replicate,
                                behaviourType='totGroomIn',
                                windowSizeInFrames=None,
                                eventCount=False,
                                eventEnds=False):
        if windowSizeInFrames is None:
            windowSizeInFrames=30 * self.fps
        limitsDict = {'controlOutput': 30 * 60 * self.fps, 'treatmentOutput': 90 * 60 * self.fps}
        vectorFunction = lambda A, du, co, ctd, useGorT, outName: self.buildWindowedVectorFromSolomon(A, du, co, ctd, outName,
                                                                                                 limitsDict, useGorT,
                                                                                                 behaviourType,
                                                                                                 windowSizeInFrames,
                                                                                                 eventCount=eventCount,
                                                                                                 eventEnds=eventEnds)

        controlDict, treatmentDict = self.buildDictionaries(treatment=treatment,
                                                       replicate=replicate,useGorT='G',
                                                       vectorFunction=vectorFunction)
        return controlDict, treatmentDict


    def buildDictionaries(self,treatment,
                          replicate,
                          useGorT='G',
                          vectorFunction=None):
        if vectorFunction is None:
            vectorFunction = self.buildVectorFromSolomon
        # --------- The different translations -------------
        colorTreatmentDictionary, colorLevelDictionary, dishNumber = kT.kTrans(replicate, treatment);
        initControl, endControl, initTreatment, endTreatment = aSA.findTimes(treatment, replicate, self.metadataFilePath);
        fileName = aSA.findSolomonFile(treatment, replicate, self.solomonDir);
        fileNamePRE = aSA.findSolomonFile(treatment, replicate, self.solomonDir, preCorrection=True);
        # fixed_duration=15*60*15; #for 15 minutes only
        # endControl=initControl+fixed_duration;
        # endTreatment=initTreatment+fixed_duration;

        # --------- Read the Solomon files
        Acontrol = r.SolomonReader(fileNamePRE,
                                   initFrame=initControl,
                                   endFrame=endControl,
                                   enforceSymetryInTrophallaxis=self.symmetricTrophallaxis)
        Atreatment = r.SolomonReader(fileName,
                                     initFrame=initTreatment,
                                     endFrame=endTreatment,
                                     enforceSymetryInTrophallaxis=self.symmetricTrophallaxis)
        if useGorT == 'G':
            gmControl = (Acontrol.getGroomingMatrix()) / float(endControl - initControl);
            gmTreatment = (Atreatment.getGroomingMatrix()) / float(endTreatment - initTreatment);
        else:
            gmControl = (Acontrol.getTrophallaxisMatrix()) / float(endControl - initControl);
            gmTreatment = (Atreatment.getTrophallaxisMatrix()) / float(endTreatment - initTreatment);

        numAnts = Acontrol.params['numAnimals'];
        antColors = Acontrol.buttonNames[:numAnts];

        # printOutControl (before)
        outControlDict = {'timeSeries': Acontrol,
                          'matrix': gmControl,
                          'duration': float(endControl - initControl),
                          'colors': antColors,
                          'GorT': useGorT,
                          'A': Acontrol}
        outTreatmentDict = {'timeSeries': Atreatment,
                            'matrix': gmTreatment,
                            'duration': float(endTreatment - initTreatment),
                            'colors': antColors,
                            'GorT': useGorT,
                            'A': Atreatment}
        allVars = [('controlOutput', Acontrol, float(endControl - initControl), outControlDict)
            , ('treatmentOutput', Atreatment, float(endTreatment - initTreatment), outTreatmentDict)]
        for (outName, A, duration, outDict) in allVars:
            # print("\n******\n"+outName+"\n*******\n")
            for color in antColors:
                thisAntsTreatment = colorTreatmentDictionary[color];  # <-
                thisAntsLevel = colorLevelDictionary[color];  # <-
                headRFP, headGFP, bodyRFP, bodyGFP = self.getHeadAndBodyCounts(
                                                                          replicate,
                                                                          dishNumber,
                                                                          color);
                thisAntsddpcrVector = np.array([headRFP, headGFP, bodyRFP, bodyGFP]);  # <-

                # By default this points to buildVectorFromSolomon
                thisAntsSolomonVector = vectorFunction(A,
                                                       duration,
                                                       color,
                                                       colorTreatmentDictionary, useGorT=useGorT, outName=outName)  # <-
                outputV = [color, thisAntsTreatment, thisAntsLevel];

                dictThisAnt = dict()
                # If this is a table output, the ant's dict contains as keys keysForOutput
                if vectorFunction == self.buildVectorFromSolomon:
                    for av in thisAntsddpcrVector:
                        outputV.append(av)
                    for av in thisAntsSolomonVector:
                        outputV.append(av)
                    for i in range(len(keysForOutput)):
                        dictThisAnt[keysForOutput[i]] = outputV[i];
                else:
                    # Otherwise, it contains three keys: antDescription antSpores and antTimeSeries
                    dictThisAnt['antDescription'] = outputV
                    dictThisAnt['antSpores'] = thisAntsddpcrVector
                    dictThisAnt['antTimeSeries'] = thisAntsSolomonVector

                outDict[color] = dictThisAnt

        return outControlDict, outTreatmentDict

    # behaviourType can be groomIn, groomOut, groomTot
    def buildWindowedTablesPairwise(self,treatmentReplicatePairs,
                            windowSizeInFrames=None,
                            separator="," ,eventStartCount=False,
                            secondOnlyTreated=True):
        if windowSizeInFrames is None:
            windowSizeInFrames=30 * self.fps

        countp = "_EventStartCounts" if eventStartCount else ""
        fileSufix = "_windowed_" + str(windowSizeInFrames) + "_outPairwise"  + countp +".csv"
        controlOutFile = open(self.outFolder + "PRE" + fileSufix, "wt")
        treatmentOutFile = open(self.outFolder + "POST" + fileSufix, "wt")
        cDs = []
        tDs = []
        for treatment, replicate in treatmentReplicatePairs:
            print(treatment, replicate, " : ")
            controlDict, treatmentDict = self.buildWindowedTimeSeries_pairwise(treatment, replicate,
                                                                 windowSizeInFrames=windowSizeInFrames,
                                                                 eventCount=eventStartCount,
                                                                 secondOnlyTreated=secondOnlyTreated)
            antColors = controlDict['colors']
            controlDict['treatment_replicate_period'] = (treatment,replicate,'PRE')
            treatmentDict['treatment_replicate_period'] = (treatment, replicate,'POST')
            cDs.append(controlDict)
            tDs.append(tDs)

            for outFile, inDict in [(controlOutFile, controlDict), (treatmentOutFile, treatmentDict)]:
                for color in inDict.keys():
                    if type(color) is not tuple or color[0] not in antColors or color[1] not in antColors:
                        continue
                    indicolor = inDict[color]
                    #print("\t\t"+str(color)+" ty"+str(type(indicolor)))
                    ant1Description = indicolor['ant1Description']
                    ant1Spores = indicolor['ant1Spores']
                    ant2Description = indicolor['ant2Description']
                    ant2Spores = indicolor['ant2Spores']
                    antPairTimeSeries = indicolor['antPairTimeSeries']
                    antRow = [treatment, replicate] \
                             + ant1Description \
                             + ant1Spores.tolist() \
                             + ant2Description \
                             + ant2Spores.tolist() \
                             + [" "+str(x) for x in antPairTimeSeries.tolist()]
                    outFile.write(separator.join([str(x) for x in antRow]))
                    outFile.write("\n")
            controlOutFile.flush()
            treatmentOutFile.flush()
        controlOutFile.close()
        treatmentOutFile.close()

        return cDs,tDs



    def buildVectorForPairPairwise(self,A, duration, color, color2, colorTreatmentDictionary, outName,
                        durationLimitsDict, useGorT='G',
                        windowSizeInFrames=None,
                        eventCount=False):
        if windowSizeInFrames is None:
            windowSizeInFrames=30 * self.fps
        antColors = A.buttonNames;
        rnum, gnum, tnum = self.findNumbersOfRedAndGreenAnts(colorTreatmentDictionary, antColors);
        numAnts = A.params['numAnimals'];

        colGrooming = antColors.index(color);
        indThisAnt = colGrooming + 1;

        indAnt2 = antColors.index(color2) + 1

        timeLimit = durationLimitsDict[outName] if outName in durationLimitsDict.keys() else 30 * 60 * self.fps
        theMatr = A.theMatr[:min(timeLimit, A.theMatr.shape[0]), :]
        numWindows = int(timeLimit / windowSizeInFrames)
        bins = [i * windowSizeInFrames for i in range(numWindows + 1)]

        # the column that contains to whom is ant1 applying grooming
        groomOutCol = theMatr[:, colGrooming]
        y = np.zeros_like(groomOutCol)
        # Only the frames in which the receptor is ant2 we mark as 1
        y[groomOutCol == indAnt2] += 1;
        if eventCount:
            # we take only the frames in which the grooming event started
            framesToBin = np.nonzero(y[1:]-y[:-1]==1)[0]
            lf = framesToBin.shape[0]
            if len(framesToBin)>0:
                tmpvect = np.zeros(max(framesToBin)+1)
                tmpvect[framesToBin] = 1
                framesToBin = np.nonzero((tmpvect[1:]-tmpvect[:-1])==1)[0]
        else:
            framesToBin = np.nonzero(y)[0]

        theVector, theBins = np.histogram(framesToBin, bins=bins)
        #print("\t\t"+color+","+color2+" : "+str(y.sum())+" == "+str(len(framesToBin))+" == "+str(theVector.sum()))
        return theVector

    def buildWindowedTimeSeries_pairwise(self, treatment,
                                replicate,
                                windowSizeInFrames=None,
                                eventCount=False,
                                secondOnlyTreated=True):
        if windowSizeInFrames is None:
            windowSizeInFrames=30 * self.fps
        limitsDict = {'controlOutput': 30 * 60 * self.fps, 'treatmentOutput': 90 * 60 * self.fps}
        vectorFunction = lambda A, du, co, co2, ctd, useGorT, outName: self.buildVectorForPairPairwise(A, du, co, co2, ctd, outName,
                                                                                                 limitsDict, useGorT,
                                                                                                 windowSizeInFrames,
                                                                                                 eventCount=eventCount)

        controlDict, treatmentDict = self.buildDictionariesPairwise(treatment=treatment,
                                                       replicate=replicate,
                                                       vectorFunction=vectorFunction,
                                                       useGorT='G',
                                                       secondOnlyTreated=secondOnlyTreated)
        return controlDict, treatmentDict


    def buildDictionariesPairwise(self, treatment,
                          replicate,
                          vectorFunction,
                          useGorT='G',
    					  secondOnlyTreated=True):
        # --------- The different translations -------------
        colorTreatmentDictionary, colorLevelDictionary, dishNumber = kT.kTrans(replicate, treatment);
        initControl, endControl, initTreatment, endTreatment = aSA.findTimes(treatment, replicate, self.metadataFilePath);
        fileName = aSA.findSolomonFile(treatment, replicate, self.solomonDir);
        fileNamePRE = aSA.findSolomonFile(treatment, replicate, self.solomonDir, preCorrection=True);
        # fixed_duration=15*60*15; #for 15 minutes only
        # endControl=initControl+fixed_duration;
        # endTreatment=initTreatment+fixed_duration;

        # --------- Read the Solomon files
        Acontrol = r.SolomonReader(fileNamePRE,
                                   initFrame=initControl,
                                   endFrame=endControl,
                                   enforceSymetryInTrophallaxis=self.symmetricTrophallaxis)
        Atreatment = r.SolomonReader(fileName,
                                     initFrame=initTreatment,
                                     endFrame=endTreatment,
                                     enforceSymetryInTrophallaxis=self.symmetricTrophallaxis)
        if useGorT == 'G':
            gmControl = (Acontrol.getGroomingMatrix()) / float(endControl - initControl);
            gmTreatment = (Atreatment.getGroomingMatrix()) / float(endTreatment - initTreatment);
        else:
            gmControl = (Acontrol.getTrophallaxisMatrix()) / float(endControl - initControl);
            gmTreatment = (Atreatment.getTrophallaxisMatrix()) / float(endTreatment - initTreatment);

        numAnts = Acontrol.params['numAnimals'];
        antColors = Acontrol.buttonNames[:numAnts];

        # printOutControl (before)
        outControlDict = {'timeSeries': Acontrol,
                          'matrix': gmControl,
                          'duration': float(endControl - initControl),
                          'colors': antColors,
                          'GorT': useGorT,
                          'A': Acontrol}
        outTreatmentDict = {'timeSeries': Atreatment,
                            'matrix': gmTreatment,
                            'duration': float(endTreatment - initTreatment),
                            'colors': antColors,
                            'GorT': useGorT,
                            'A': Atreatment}
        allVars = [('controlOutput', Acontrol, float(endControl - initControl), outControlDict)
            , ('treatmentOutput', Atreatment, float(endTreatment - initTreatment), outTreatmentDict)]
        for (outName, A, duration, outDict) in allVars:
            # print("\n******\n"+outName+"\n*******\n")
            for color in antColors:
                thisAntsTreatment = colorTreatmentDictionary[color];  # <-
                thisAntsLevel = colorLevelDictionary[color];  # <-
                headRFP, headGFP, bodyRFP, bodyGFP = self.getHeadAndBodyCounts(replicate,
                                                                          dishNumber,
                                                                          color);
                thisAntsddpcrVector = np.array([headRFP, headGFP, bodyRFP, bodyGFP]);  # <-
                outputV1 = [color, thisAntsTreatment, thisAntsLevel];

                for color2 in antColors:
                    ant2Treatment = colorTreatmentDictionary[color2]
                    ant2Level = colorLevelDictionary[color2]
                    if ant2Treatment == "N" and secondOnlyTreated:
                        continue
                    #print("\tDoing: "+color+","+color2)
                    outputV2 = [color2, ant2Treatment, ant2Level];
                    ant2headRFP, ant2headGFP, ant2bodyRFP, ant2bodyGFP = self.getHeadAndBodyCounts(replicate,
                                                                              dishNumber,
                                                                              color2);
                    ant2ddpcrVector = np.array([ant2headRFP, ant2headGFP, ant2bodyRFP, ant2bodyGFP]);
                    # By default this points to buildVectorFromSolomon
                    thisAntPairSolomonVector = vectorFunction(A,
                                                           duration,
                                                           color,
                                                           color2,
                                                           colorTreatmentDictionary, useGorT=useGorT, outName=outName)  # <-



                    # Otherwise, it contains five keys: ant1Description ant1Spores ant2Description ant2Spores and antTimeSeries
                    dictThisAnt = dict()
                    dictThisAnt['ant1Description'] = outputV1
                    dictThisAnt['ant2Description'] = outputV2
                    dictThisAnt['ant1Spores'] = thisAntsddpcrVector
                    dictThisAnt['ant2Spores'] = ant2ddpcrVector
                    dictThisAnt['antPairTimeSeries'] = thisAntPairSolomonVector

                    outDict[(color,color2)] = dictThisAnt
                    #print("buildingDict: "+outName+" "+color+","+color2+" ty:"+str(type(dictThisAnt)))

        return outControlDict, outTreatmentDict

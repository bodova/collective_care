import utils.auxSolomonAnalysis as aSA
import utils.keyTranslator as kT
import utils.tableBuilder as tB
import utils.reader as r

import numpy as np
from matplotlib import pylab as P
import networkx as nx
import pydot
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib


def writeSelfFile(solomonDir,
                  metadataFilePath, treatmentReplicatePairs, resultsDir, sep=",",
                  typesOfBehaviour=['self', 'head grooming', 'gastertip'],
                  ):
    catCommand = "cat "

    for treatment, replicate in treatmentReplicatePairs:
        colorTreatmentDictionary, colorLevelDictionary, dishNumber = kT.kTrans(replicate, treatment);
        initControl, endControl, initTreatment, endTreatment = aSA.findTimes(treatment, replicate, metadataFilePath);
        fileName = aSA.findSolomonFile(treatment, replicate, solomonDir);
        fileNamePRE = aSA.findSolomonFile(treatment, replicate, solomonDir, preCorrection=True);
        APRE = r.SolomonReader(fileNamePRE, initFrame=initControl, endFrame=endControl)
        APOST = r.SolomonReader(fileName, initFrame=initTreatment, endFrame=endTreatment)
        numAnts = APRE.params['numAnimals'];
        antColors = APRE.buttonNames[:numAnts];

        commonToAllRows = treatment + sep + str(replicate);

        outFileName = resultsDir + "Self_" + str(treatment) + str(replicate) + ".csv"
        fout = open(outFileName, 'w')
        catCommand += outFileName + " "

        # returns a list of tuples of the form
        # (start,end,actor,type)
        evPRE = APRE.getSelfGroomingEvents(typesOfBehaviour=typesOfBehaviour)
        evPOST = APOST.getSelfGroomingEvents(typesOfBehaviour=typesOfBehaviour)
        toProcess = [(evPRE, 'PRE'), (evPOST, 'POST')]
        for evl, inst in toProcess:

            for (start, end, actor, behTypeL) in evl:
                behType = behTypeL.split(" ")[0]
                antColor = antColors[actor]
                antTreatment = colorTreatmentDictionary[antColor];
                antLevel = colorLevelDictionary[antColor];
                if antTreatment == 'N':
                    antLevel = 'M'
                beg = str(start);
                dur = str(end - start);
                toWrite = commonToAllRows + sep + inst + sep + antColor + sep + antTreatment + sep + antLevel + sep + behType + sep + beg + sep + dur;
                fout.write(toWrite + "\n")

        fout.close()

    catCommand += " > " + resultsDir + "Self_ALL.csv"
    import os
    os.system(catCommand)


def beforeAndAfterForSelfBehaviours(treatmentReplicatePairs, metadataFilePath, solomonDir, useNumberOfEvents=False):
    behaviourTypes = ['gastertip', 'head grooming']
    pointsPRE = dict()
    pointsPOST = dict()
    for bt in behaviourTypes:
        pointsPOST[bt] = {'X': [], 'H': [], 'L': []};
        pointsPRE[bt] = {'X': [], 'H': [], 'L': []};

    for treatment, replicate in treatmentReplicatePairs:
        colorTreatmentDictionary, colorLevelDictionary, dishNumber = kT.kTrans(replicate, treatment);
        initControl, endControl, initTreatment, endTreatment = aSA.findTimes(treatment, replicate, metadataFilePath);
        fileName = aSA.findSolomonFile(treatment, replicate, solomonDir);
        fileNamePRE = aSA.findSolomonFile(treatment, replicate, solomonDir, preCorrection=True);
        APRE = r.SolomonReader(fileNamePRE, initFrame=initControl, endFrame=endControl)
        APOST = r.SolomonReader(fileName, initFrame=initTreatment, endFrame=endTreatment)
        numAnts = APRE.params['numAnimals'];
        antColors = APRE.buttonNames[:numAnts];
        durationPRE = float(endControl - initControl);
        durationPOST = float(endTreatment - initTreatment);

        for an in range(numAnts):
            thisColor = antColors[an];
            thisLevel = colorLevelDictionary[thisColor];
            for bt in behaviourTypes:
                beVal = APRE.buttonNames.index(bt) + 1
                if useNumberOfEvents:
                    evPRE = aSA.nonzero_intervals(APRE.theMatr[:, numAnts + an] == beVal)
                    evPOST = aSA.nonzero_intervals(APOST.theMatr[:, numAnts + an] == beVal)
                    countPRE = len(evPRE);
                    countPOST = len(evPOST);
                else:
                    countPRE = len(np.nonzero(APRE.theMatr[:, numAnts + an] == beVal)[0])
                    countPOST = len(np.nonzero(APOST.theMatr[:, numAnts + an] == beVal)[0])
                pointsPRE[bt][thisLevel].append(countPRE / durationPRE)
                pointsPOST[bt][thisLevel].append(countPOST / durationPOST)

    styDict = {'X': 'o', 'H': '^', 'L': 'v'}
    xpos = 1;
    xticks = []
    xlabels = [];
    for bt in behaviourTypes:
        for le in ['X', 'H', 'L']:
            for dd in [pointsPRE, pointsPOST]:
                Y = np.array(dd[bt][le]);
                X = xpos * np.ones_like(Y);
                xlabel = bt + " " + le + " " + ('PRE' if xpos % 2 == 1 else 'POST');
                P.plot(X, Y, styDict[le], markerfacecolor='none');
                xticks.append(xpos);
                xlabels.append(xlabel);
                xpos += 1;

    P.xticks(xticks, xlabels, rotation=90)
    P.xlim(0, max(xticks) + 1)
    P.show()


def selfGroomingVsIngrooming(treatmentReplicatePairs, solomonDir,
                             ddpcrDir,
                             metadataFilePath):
    pairsHighPRE = [];
    pairsLowPRE = [];
    pairsTritonPRE = [];
    pairsNestmatePRE = [];
    pairsHighPOST = [];
    pairsLowPOST = [];
    pairsTritonPOST = [];
    pairsNestmatePOST = [];

    totalMinI = np.Inf;
    totalMaxI = -np.Inf;
    totalMinS = np.Inf;
    totalMaxS = -np.Inf;

    # For every file
    for treatment, replicate in treatmentReplicatePairs:
        colorTreatmentDictionary, colorLevelDictionary, dishNumber = kT.kTrans(replicate, treatment);
        initControl, endControl, initTreatment, endTreatment = aSA.findTimes(treatment, replicate, metadataFilePath);
        fileName = aSA.findSolomonFile(treatment, replicate, solomonDir);
        fileNamePRE = aSA.findSolomonFile(treatment, replicate, solomonDir, preCorrection=True);
        APRE = r.SolomonReader(fileNamePRE, initFrame=initControl, endFrame=endControl)
        APOST = r.SolomonReader(fileName, initFrame=initTreatment, endFrame=endTreatment)
        numAnts = APRE.params['numAnimals'];
        antColors = APRE.buttonNames[:numAnts];

        # We read the pre and post grooming matrices
        gmPRE = APRE.getGroomingMatrix() / float(endControl - initControl)
        gmPOST = APOST.getGroomingMatrix() / float(endTreatment - initTreatment)

        # And extract the self grooming and the sum of the ingrooming
        selfPRE = np.diag(gmPRE);
        inPRE = (gmPRE - np.diag(selfPRE)).sum(axis=0);
        selfPOST = np.diag(gmPOST);
        inPOST = (gmPOST - np.diag(selfPOST)).sum(axis=0);

        totalMinI = min([totalMinI] + inPRE.tolist() + inPOST.tolist());
        totalMaxI = max([totalMaxI] + inPRE.tolist() + inPOST.tolist());
        totalMinS = min([totalMinS] + selfPRE.tolist() + selfPOST.tolist())
        totalMaxS = max([totalMaxS] + selfPRE.tolist() + selfPOST.tolist())

        # store them into the commulative lists
        for an in range(numAnts):
            level = colorLevelDictionary[antColors[an]];
            treat = colorTreatmentDictionary[antColors[an]];
            thisPairPRE = (selfPRE[an], inPRE[an])
            thisPairPOST = (selfPOST[an], inPOST[an])

            if level == 'H':
                pairsHighPRE.append(thisPairPRE)
                pairsHighPOST.append(thisPairPOST)
            elif level == 'L':
                pairsLowPRE.append(thisPairPRE)
                pairsLowPOST.append(thisPairPOST)
            elif treat == 'T':
                pairsTritonPRE.append(thisPairPRE)
                pairsTritonPOST.append(thisPairPOST)
            elif treat == 'N':
                pairsNestmatePRE.append(thisPairPRE)
                pairsNestmatePOST.append(thisPairPOST)

    P.subplot(1, 2, 1)
    P.plot([x[0] for x in pairsHighPRE], [x[1] for x in pairsHighPRE], '^', label='High', markersize=15)
    P.plot([x[0] for x in pairsLowPRE], [x[1] for x in pairsLowPRE], 'v', label='Low', markersize=15)
    P.plot([x[0] for x in pairsTritonPRE], [x[1] for x in pairsTritonPRE], 'x', label='Triton X', markersize=15)
    P.plot([x[0] for x in pairsNestmatePRE], [x[1] for x in pairsNestmatePRE], 'ko', label='Nest Mate', markersize=15)
    P.xlabel("Self grooming\n (percentage of total time)")
    P.ylabel("Alorooming received\n (percentage of total time)")
    P.legend(numpoints=1)
    P.title("Pre-exposure")
    P.xlim(totalMinS - 0.01, totalMaxS)
    P.ylim(totalMinI - 0.01, totalMaxI)

    P.subplot(1, 2, 2)
    P.plot([x[0] for x in pairsHighPOST], [x[1] for x in pairsHighPOST], '^', label='High', markersize=15)
    P.plot([x[0] for x in pairsLowPOST], [x[1] for x in pairsLowPOST], 'v', label='Low', markersize=15)
    P.plot([x[0] for x in pairsTritonPOST], [x[1] for x in pairsTritonPOST], 'x', label='Triton X', markersize=15)
    P.plot([x[0] for x in pairsNestmatePOST], [x[1] for x in pairsNestmatePOST], 'ko', label='Nest Mate', markersize=15)
    P.xlabel("Self grooming\n (percentage of total time)")
    P.ylabel("Allogrooming received\n (percentage of total time)")
    P.legend(numpoints=1)
    P.title("Post-exposure")
    P.xlim(totalMinS - 0.01, totalMaxS)
    P.ylim(totalMinI - 0.01, totalMaxI)

    P.show()


def treatmentLevel2Plotstyle(treatment, level):
    style = ""
    if treatment == 'T':
        return 'bx'
    if treatment == 'N':
        return 'ko'
    if treatment == 'R':
        style += 'r'
    if treatment == 'G':
        style += 'g'
    if level == 'H':
        style += '^'
    if level == 'L':
        style += 'v'

    return style


def plotstyle2Treatmentlevel(plotstyle):
    color2treatment = {'b': 'T', 'k': 'N', 'r': 'R', 'g': 'G'}
    marker2level = {'x': 'X', 'o': 'M', '^': 'H', 'v': 'L'}
    return color2treatment[plotstyle[0]] + marker2level[plotstyle[1]];


# get duration of grooming events to infected
# plot duration vs time
#     duration vs grooming time from all ants
#     duration vs grooming time from each ant

def plotEventDurations(treatmentReplicatePairs,
                       solomonDir,
                       metadataFilePath):
    allTreatments = set([tr[0] for tr in treatmentReplicatePairs]);

    points_VsTimePRE = dict();
    points_VsGroomingAllPRE = dict();
    points_VsGroomingEachPRE = dict();
    points_VsTimePOST = dict();
    points_VsGroomingAllPOST = dict();
    points_VsGroomingEachPOST = dict();
    styles_PRE = dict();
    styles_POST = dict();
    for treatment in allTreatments:
        points_VsTimePRE[treatment] = [];
        points_VsGroomingAllPRE[treatment] = [];
        points_VsGroomingEachPRE[treatment] = [];
        points_VsTimePOST[treatment] = [];
        points_VsGroomingAllPOST[treatment] = [];
        points_VsGroomingEachPOST[treatment] = [];
        styles_PRE[treatment] = [];
        styles_POST[treatment] = [];

    for treatment, replicate in treatmentReplicatePairs:
        colorTreatmentDictionary, colorLevelDictionary, dishNumber = kT.kTrans(replicate, treatment);
        initControl, endControl, initTreatment, endTreatment = aSA.findTimes(treatment, replicate, metadataFilePath);
        fileName = aSA.findSolomonFile(treatment, replicate, solomonDir);
        fileNamePRE = aSA.findSolomonFile(treatment, replicate, solomonDir, preCorrection=True);
        APRE = r.SolomonReader(fileNamePRE, initFrame=initControl, endFrame=endControl)
        APOST = r.SolomonReader(fileName, initFrame=initTreatment, endFrame=endTreatment)
        numAnts = APRE.params['numAnimals'];
        antColors = APRE.buttonNames[:numAnts];

        # These are list of tuples of the form (start,end,actor,groomedAnt)
        eventsPRE = APRE.getGroomingEvents();
        eventsPOST = APOST.getGroomingEvents();

        allAntNumbers = set(range(numAnts))

        for (start, end, actor, groomedAnt) in eventsPRE:
            level = colorLevelDictionary[antColors[groomedAnt]];
            if not (level in ['L', 'H']):  # We consider only events to groomed ants
                continue
            treat = colorTreatmentDictionary[antColors[groomedAnt]];
            style = treatmentLevel2Plotstyle(treat, level)
            dur = end - start;

            groomingThisActor = APRE.theMatr[0:start, actor]
            timeEach = len(np.nonzero(groomingThisActor == groomedAnt + 1)[0])

            antsExceptGroomed = list(allAntNumbers - set([groomedAnt]));
            groomingAllExceptGroomed = APRE.theMatr[0:start, antsExceptGroomed]
            timeAll = len(np.nonzero(groomingAllExceptGroomed == groomedAnt + 1)[0])

            points_VsGroomingAllPRE[treatment].append((dur, timeAll, style))
            points_VsGroomingEachPRE[treatment].append((dur, timeEach, style))
            points_VsTimePRE[treatment].append((dur, start, style))

        for (start, end, actor, groomedAnt) in eventsPOST:
            level = colorLevelDictionary[antColors[groomedAnt]];
            if not (level in ['L', 'H']):  # We consider only events to groomed ants
                continue
            treat = colorTreatmentDictionary[antColors[groomedAnt]];
            style = treatmentLevel2Plotstyle(treat, level)

            dur = end - start;

            groomingThisActor = APOST.theMatr[0:start, actor]
            timeEach = len(np.nonzero(groomingThisActor == groomedAnt + 1)[0])

            antsExceptGroomed = list(allAntNumbers - set([groomedAnt]));
            groomingAllExceptGroomed = APOST.theMatr[0:start, antsExceptGroomed]
            timeAll = len(np.nonzero(groomingAllExceptGroomed == groomedAnt + 1)[0])

            points_VsGroomingAllPOST[treatment].append((dur, timeAll, style))
            points_VsGroomingEachPOST[treatment].append((dur, timeEach, style))
            points_VsTimePOST[treatment].append((dur, start, style))

    f1 = plotScatter1(points_VsTimePRE, points_VsTimePOST, "Grooming event start time", allTreatments)

    f2 = plotScatter1(points_VsGroomingAllPRE, points_VsGroomingAllPOST, "Time ant was groomed by all ants before",
                 allTreatments)
    f3 = plotScatter1(points_VsGroomingEachPRE, points_VsGroomingEachPOST, "Time ant was groomed by this ant before",
                 allTreatments)

    return [f1,f2,f3]


def plotScatter1(pointsPRE, pointsPOST, ylab, allTreatments):
    fig = P.figure()
    subplotNum = 0;
    maxY = max([max([x[0] for x in evdata]) for evdata in pointsPRE.values() if len(evdata) > 0]
               + [max([x[0] for x in evdata]) for evdata in pointsPOST.values() if len(evdata) > 0])
    for ti, treatment in enumerate(allTreatments):
        subplotNum += 1;
        P.subplot(2, len(allTreatments), subplotNum)
        P.ylim(0, maxY)
        P.title(treatment)
        P.xticks([])
        if ti > 0:
            P.yticks([])
        allevents = [ss[0] for ss in pointsPRE[treatment]] + [ss[0] for ss in pointsPOST[treatment]]
        if len(allevents) == 0:
            continue
        maxX = max(allevents)

        b = pointsPRE[treatment];
        allStyles = set([ss[2] for ss in b])
        for style in allStyles:
            thisStyleX = [ss[0] for ss in b if ss[2] == style]
            thisStyleY = [ss[1] for ss in b if ss[2] == style]
            mfc = 'none' if style[0] == 'k' else None
            P.plot(thisStyleX, thisStyleY, style, label=plotstyle2Treatmentlevel(style), markerfacecolor=mfc)

        P.legend(numpoints=1, title="PRE")
        P.xlim(0, maxX)
        P.ylim(0, maxY)

        P.subplot(2, len(allTreatments), subplotNum + len(allTreatments))
        b = pointsPOST[treatment];
        allStyles = set([ss[2] for ss in b])
        for style in allStyles:
            thisStyleX = [ss[0] for ss in b if ss[2] == style]
            thisStyleY = [ss[1] for ss in b if ss[2] == style]
            mfc = 'none' if style[0] == 'k' else None
            P.plot(thisStyleX, thisStyleY, style, label=plotstyle2Treatmentlevel(style), markerfacecolor=mfc)

        P.legend(numpoints=1, title="POST")
        P.xlim(0, maxX)
        P.ylim(0, maxY)
        if ti > 0:
            P.yticks([])

    fig.text(0.5, 0.02, "Grooming event duration", ha='center')
    fig.text(0.04, 0.5, ylab, va='center', rotation='vertical')
    fig.subplots_adjust(bottom=0.11, right=0.97, top=0.92, wspace=0.11, hspace=0.03)
    return fig


def plotNetworksAndMatricesSeveral(treatmentReplicatePairs,
                                   solomonDir,
                                   ddpcrDir,
                                   resultsDir,
                                   metadataFilePath,
                                   useGorT='G', colorMax=None, inTime=False,
                                   withOutsideRings=True,
                                   useNegatives=True,
                                   imageDir=None, doShow=True):
    for treatment, replicate in treatmentReplicatePairs:
        maxes = plotNetworksAndMatricesSingle(treatment,
                                              replicate,
                                              solomonDir=solomonDir,
                                              ddpcrDir=ddpcrDir,
                                              metadataFilePath=metadataFilePath,
                                              useGorT=useGorT,
                                              colorMax=colorMax,
                                              inTime=inTime,
                                              resultsDir=resultsDir,
                                              withOutsideRings=withOutsideRings,
                                              useNegatives=useNegatives,
                                              imageDir=imageDir)

        print(str(treatment) + str(replicate) + "\t" + str(maxes))

    if doShow:
        P.show()
    else:
        for treatment, replicate in treatmentReplicatePairs:
            P.close()
            P.close()


def getPerVideoColorLevels(colorTreatmentDictionary, colorLevelDictionary):
    levelR = 'X'
    levelG = 'X'
    for acolor in colorTreatmentDictionary.keys():
        if colorTreatmentDictionary[acolor] == 'R':
            levelR = colorLevelDictionary[acolor]
        if colorTreatmentDictionary[acolor] == 'G':
            levelG = colorLevelDictionary[acolor]

    return levelR, levelG


def plotNetworksAndMatricesSingle(treatment,
                                  replicate,
                                  solomonDir,
                                  ddpcrDir,
                                  resultsDir,
                                  metadataFilePath,
                                  useGorT='G',
                                  colorMax=None,
                                  inTime=False,
                                  withOutsideRings=True,
                                  useNegatives=True,
                                  imageDir=None):
    # --------- The different translations -------------
    colorTreatmentDictionary, colorLevelDictionary, dishNumber = kT.kTrans(replicate, treatment);
    initControl, endControl, initTreatment, endTreatment = aSA.findTimes(treatment, replicate, metadataFilePath);
    fileName = aSA.findSolomonFile(treatment, replicate, solomonDir);
    fileNamePRE = aSA.findSolomonFile(treatment, replicate, solomonDir, preCorrection=True);
    videoLevelR, videoLevelG = getPerVideoColorLevels(colorTreatmentDictionary, colorLevelDictionary);
    print(str(replicate) + " " + str(dishNumber))
    # print(str(treatment)+"\t vG:  "+videoLevelG+"\t vR:  "+videoLevelR)
    # --------- Read the Solomon files
    Acontrol = r.SolomonReader(fileNamePRE, initFrame=initControl, endFrame=endControl)
    Atreatment = r.SolomonReader(fileName, initFrame=initTreatment, endFrame=endTreatment)

    numAnts = Acontrol.params['numAnimals'];
    antColors = Acontrol.buttonNames[:numAnts];
    listAntTypes = [colorTreatmentDictionary[x] + colorLevelDictionary[x] for x in antColors]
    # Permute antColors by orderAnts
    orderAnts = kT.orderByTreatmentLevel(listAntTypes)
    antColors = [antColors[x] for x in orderAnts]
    countsDictionary = dict()
    countsDictionaryNumbers = dict()
    maxR = -np.inf;
    maxG = -np.inf;
    minR = np.inf;
    minG = np.inf;
    for color in antColors:
        thisAntsTreatment = colorTreatmentDictionary[color];
        thisAntsLevel = colorLevelDictionary[color];
        headRFP, headGFP, bodyRFP, bodyGFP = tB.getHeadAndBodyCounts(ddpcrDir,
                                                                     replicate,
                                                                     dishNumber,
                                                                     color);
        headRFP, headGFP, bodyRFP, bodyGFP = aSA.deltaSporeCount(thisAntsTreatment, thisAntsLevel, headRFP, headGFP,
                                                                 bodyRFP, bodyGFP, videoLevelR, videoLevelG,
                                                                 useNegatives=useNegatives);

        strheadRFP = str(round(headRFP, 3))
        strheadGFP = str(round(headGFP, 3))
        strbodyRFP = str(round(bodyRFP, 3))
        strbodyGFP = str(round(bodyGFP, 3))

        countsDictionary[color] = "h:" + strheadRFP + "/" + strheadGFP + "\nb:" + strbodyRFP + "/" + strbodyGFP;
        thisAntNum = antColors.index(color);
        # counts is of the form:  headRed,headGreen,bodyRed,bodyGreen
        countsDictionaryNumbers[thisAntNum] = [headRFP, headGFP, bodyRFP, bodyGFP]
        maxR = max([maxR, headRFP, bodyRFP])
        maxG = max([maxG, headGFP, bodyGFP])
        minR = min([minR, headRFP, bodyRFP])
        minG = min([minG, headGFP, bodyGFP])
    # print(str([headRFP,headGFP,bodyRFP,bodyGFP]))
    # print("maxes:  "+str([minR,maxR,minG,maxG]))

    for thisAntNum in range(len(antColors)):
        cv = countsDictionaryNumbers[thisAntNum];
        hR = cv[0]
        hG = cv[1]
        bR = cv[2]
        bG = cv[3]
        countsDictionaryNumbers[thisAntNum] = [hR, hG, bR, bG]

    if inTime:
        doPlotsByTimesteps(fileName,
                           fileNamePRE,
                           initControl,
                           endControl,
                           initTreatment,
                           endTreatment,
                           countsDictionary,
                           aControl=Acontrol,
                           aTreatment=Atreatment,
                           antTreatments=colorTreatmentDictionary,
                           antLevels=colorLevelDictionary,
                           colorMax=colorMax,
                           useGorT=useGorT,
                           orderAnts=orderAnts);

    else:
        return doPlots1(fileName, fileNamePRE,
                        initControl,
                        endControl,
                        initTreatment,
                        endTreatment,
                        countsDictionary,
                        antColors=antColors,
                        orderAnts=orderAnts,
                        aControl=Acontrol,
                        aTreatment=Atreatment,
                        makePlots=True,
                        makeGraphs=True,
                        antTreatments=colorTreatmentDictionary,
                        antLevels=colorLevelDictionary,
                        colorMax=colorMax,
                        useGorT=useGorT,
                        resultsDir=resultsDir,
                        countsDictionaryNumerical=countsDictionaryNumbers,
                        withRings=withOutsideRings, useNegatives=useNegatives,
                        imageDir=imageDir);


def doPlotsByTimesteps(fileName,
                       fileNamePRE,
                       initControl,
                       endControl,
                       initTreatment,
                       endTreatment,
                       countsDictionary,
                       stepSize=None,
                       aControl=None,
                       aTreatment=None,
                       antTreatments=[],
                       antLevels=[],
                       colorMax=None,
                       useGorT='G',
                       orderAnts=None):
    if stepSize == None:
        stepSize = 13500;
    maxes = [];
    mins = [];
    gms = [];
    timeStepLabels = [];
    numTimeSteps = 0;

    offset = (endControl - initControl) % stepSize - 1;

    A = r.SolomonReader(fileName, initFrame=initControl, endFrame=endTreatment)

    initFrame = offset;
    endFrame = initFrame + stepSize;
    switchedToTreatmentSegment = False;

    print("initPRE" + str(initControl) + " endPRE" + str(endControl))
    print("initPOST" + str(initTreatment) + " endPOST" + str(endTreatment))

    while (endFrame <= endTreatment - initControl):

        timeSeries = A.theMatr[initFrame:endFrame, :]
        print(str(initFrame + initControl) + " : " + str(endFrame + initControl) + " " + str(
            switchedToTreatmentSegment) + "  " + str(timeSeries.shape))
        print(str((endFrame - initFrame) / (60 * 15))),
        if useGorT == 'G':
            gm = (A.getGroomingMatrix(timeSeries)) / float(stepSize);
        else:
            gm = (A.getTrophallaxisMatrix(timeSeries)) / float(stepSize);

        gm = gm[orderAnts, :]
        gm = gm[:, orderAnts]

        # This is only required once, but it's ok to do it in every time step
        numAnts = A.params['numAnimals'];
        maxes.append(gm.max());
        mins.append(gm.min());
        gms.append(gm);
        numTimeSteps += 1;

        if switchedToTreatmentSegment:
            timeStepLabels.append('POST')
            print("POST")
        else:
            timeStepLabels.append('PRE')
            print("PRE")

        initFrame = endFrame + 1;
        endFrame = initFrame + stepSize;
        if endFrame > endControl - initControl and not (switchedToTreatmentSegment):
            switchedToTreatmentSegment = True;
            initFrame = initTreatment - initControl;
            endFrame = initFrame + stepSize;

    antColors = A.buttonNames[:numAnts];
    antColors = [antColors[x] for x in orderAnts]
    nodeLabels = dict();
    for aname in antColors:
        nodeLabels[aname] = antTreatments[aname] + antLevels[aname];
        if nodeLabels[aname] == 'NX':
            nodeLabels[aname] = "NM";
        nodeLabels[aname] += " (" + aname[0:2] + ")\n" + countsDictionary[aname];

    totalMax = max(maxes) if colorMax == None else colorMax;
    totalMin = min(mins);
    labs = [nodeLabels[xx] for xx in antColors]
    labsX = [nodeLabels[xx].replace("\n", " ") for xx in antColors]

    btype = 'Trophallaxis' if useGorT == 'T' else 'Grooming'

    fig, axes = plt.subplots(nrows=3, ncols=3)

    # for timeStep in range(numTimeSteps):
    subplotNum = -1;
    for ax in axes.flat:
        subplotNum += 1;
        # P.subplot(1,numTimeSteps,timeStep+1)
        if subplotNum == 2:
            ax.set_yticks([0, 1])
            ax.set_xticks([0, 1])
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            ax.text(0.09, 0.8, fileName.split("/")[-1].split('.')[0] + " " + useGorT)
            continue
        timeStep = subplotNum if subplotNum < 2 else subplotNum - 1;
        print(str(subplotNum) + " ts:" + str(timeStep))

        # P.subplot(3,3,subplotNum)
        im = ax.imshow(gms[timeStep], interpolation='none', vmin=totalMin, vmax=totalMax)

        if (subplotNum % 3) == 0:
            ax.set_yticks(range(numAnts))
            ax.set_yticklabels(labsX)
        else:
            ax.set_yticks([]);
        ax.set_xticks([]);
        ax.set_title(timeStepLabels[timeStep])
    # if timeStep == numTimeSteps-1:
    # P.colorbar();

    fig.colorbar(im, ax=axes.flat[2])
    fig.subplots_adjust(left=0.16, right=1.0, top=0.96, bottom=0.04, wspace=0.43)


def doPlots1(fileName, fileNamePRE,
             initControl,
             endControl,
             initTreatment,
             endTreatment,
             countsDictionary,
             resultsDir,
             antColors,
             orderAnts,
             countsDictionaryNumerical,
             aControl=None,
             aTreatment=None,
             makePlots=True,
             makeGraphs=True,
             antTreatments=[],
             antLevels=[],
             colorMax=None,
             useGorT='G',
             withRings=True, useNegatives=True, imageDir=None):
    vname = fileName.split("/")[-1].split('.')[0]
    if imageDir != None:
        fig1Name = imageDir + vname + '_Matrix.png'
        fig2Name = imageDir + vname + '_Graph.png'

    maxes = [];
    mins = [];
    gms = [];

    # ----- First we get the matrices ------
    # Control ----
    if aControl == None:
        A = r.SolomonReader(fileNamePRE, initFrame=initControl, endFrame=endControl)
    else:
        A = aControl;
    # get grooming matrix
    if useGorT == 'G':
        gm = (A.getGroomingMatrix()) / float(endControl - initControl);
    else:
        gm = (A.getTrophallaxisMatrix()) / float(endControl - initControl);

    # Sort specified by orderAnts vector
    gm = gm[orderAnts, :]
    gm = gm[:, orderAnts]
    numAnts = A.params['numAnimals'];
    if len(antLevels) == 0:
        antLevels = ['' for i in range(numAnts)]
    if len(antTreatments) == 0:
        antTreatments = ['' for i in range(numAnts)]
    # Save PRE interaction matrix
    fileName_noExt = fileName.split('/')[-1].split('.')[0]
    np.savetxt(resultsDir + 'PRE' + '_' + useGorT + 'Matrix' + fileName_noExt + '.csv', gm)
    maxes.append(gm.max());
    mins.append(gm.min());
    gms.append(gm);
    # Treatment ----
    # Load the data
    if aTreatment == None:
        A = r.SolomonReader(fileName, initFrame=initTreatment, endFrame=endTreatment)
    # get grooming matrix
    else:
        A = aTreatment;
    if useGorT == 'G':
        gm = (A.getGroomingMatrix()) / float(endTreatment - initTreatment);
    else:
        gm = (A.getTrophallaxisMatrix()) / float(endTreatment - initTreatment);
    # Sort specified by orderAnts vector
    gm = gm[orderAnts, :]
    gm = gm[:, orderAnts]
    # Save POST interaction matrix
    np.savetxt(resultsDir + 'POST' + '_' + useGorT + 'Matrix' + fileName_noExt + '.csv', gm)
    maxes.append(gm.max());
    mins.append(gm.min());
    gms.append(gm);

    # antColors = A.buttonNames[:numAnts];
    nodeLabels = dict();
    treatmentNumerical = dict()
    for aname in antColors:
        nodeLabels[aname] = antTreatments[aname] + antLevels[aname];
        if nodeLabels[aname] == 'NX':
            nodeLabels[aname] = "NM";
        nodeLabels[aname] += " (" + aname[0:2] + ")\n" + countsDictionary[aname];
        thisAntNum = antColors.index(aname)
        treatmentNumerical[thisAntNum] = antTreatments[aname] + antLevels[aname]

    totalMax = max(maxes) if colorMax == None else colorMax;
    totalMin = min(mins);
    # labs = [nodeLabels[xx] for xx in  A.buttonNames[:numAnts ]]
    labs = [nodeLabels[xx] for xx in antColors]

    if makePlots:  #
        # P.figure()
        # ----- Now we plot --------
        fig = doPlots2(labs, totalMin, totalMax, gms, numAnts, useGorT, vname)
        fig.subplots_adjust(left=0.12, right=0.91, top=0.97, bottom=0.02, wspace=0.11)
        if imageDir != None:
            plt.savefig(fig1Name)

    if makeGraphs:
        maxEdgeWeight = max(gms[0].max(), gms[1].max());
        minEdgeWeight = min(gms[0].min(), gms[1].min());
        edgeWeightNormalizations = lambda ew: (ew - minEdgeWeight) / (maxEdgeWeight - minEdgeWeight)
        # print(str(minEdgeWeight)+" "+str(maxEdgeWeight))
        normGMs = [(gm - minEdgeWeight) / (maxEdgeWeight - minEdgeWeight) for gm in gms]

        # print(str(countsDictionaryNumerical))
        # print(str([gm.shape for gm in gms]))
        aSA.plotGraphWithPies('control.png', countsDictionaryNumerical, treatmentNumerical, gms[0],
                              withOutsideRings=withRings)
        aSA.plotGraphWithPies('treatment.png', countsDictionaryNumerical, treatmentNumerical, gms[1],
                              withOutsideRings=withRings)
        # print("PRE "+str(gms[0].min())+" "+str(gms[0].max()))
        # print("POST "+str(gms[1].min())+" "+str(gms[1].max()))
        # print(str(minEdgeWeight)+" ->"+str(maxEdgeWeight))

        # Now the display
        fig = plt.figure(figsize=(21, 11), dpi=100)
        gs = gridspec.GridSpec(1, 3, width_ratios=[5, 5, 1])
        ax0 = plt.subplot(gs[0])
        imageC = P.imread('control.png')
        ax0.imshow(imageC)
        ax0.set_xticks([])
        ax0.set_yticks([])
        btype = 'Trophallaxis' if useGorT == 'T' else 'Grooming'
        ax0.set_title(btype + '\n PRE')

        imageT = P.imread('treatment.png')
        ax0 = plt.subplot(gs[1])
        ax0.imshow(imageT)
        ax0.set_xticks([])
        ax0.set_yticks([])
        ax0.set_title(vname + "\n" + 'POST')

        ax = plt.subplot(gs[2])
        aSA.myColorMap(ax, useNegatives=useNegatives)

        fig.subplots_adjust(left=0.02, right=0.99, top=0.97, bottom=0.02)
        if imageDir != None:
            plt.savefig(fig2Name)

    return maxes


def doPlots2(labs,
             totalMin,
             totalMax,
             gms,
             numAnts,
             useGorT,
             fileName,
             subplotRows=1,
             subPlotCols=2,
             subplotNums=[1, 2]):
    # Control ----
    fig = plt.figure(figsize=(21, 11), dpi=100)
    gs = gridspec.GridSpec(1, 3, width_ratios=[8, 8, 1])

    labs2 = [la.split("\n")[0] for la in labs]

    ax = plt.subplot(gs[0])
    # P.subplot(1,2,1)
    im = ax.imshow(gms[0], interpolation='none', vmin=totalMin, vmax=totalMax)
    ax.set_xticks(range(numAnts))
    ax.set_xticklabels(labs2)
    ax.set_yticks(range(numAnts))
    ax.set_yticklabels(labs)
    btype = 'Trophallaxis' if useGorT == 'T' else 'Grooming'
    ax.set_title(btype + '\n PRE')
    # P.colorbar();

    # Treatment ----
    # P.subplot(1,2,2)
    ax = plt.subplot(gs[1])
    im = ax.imshow(gms[1], interpolation='none', vmin=totalMin, vmax=totalMax)
    ax.set_xticks(range(numAnts))
    ax.set_xticklabels(labs2)
    ax.set_yticks([])
    ax.set_yticklabels([])
    btype = 'Trophallaxis' if useGorT == 'T' else 'Grooming'
    ax.set_title(fileName + '\n POST')

    plt.colorbar(im, cax=plt.subplot(gs[2]), orientation='vertical')

    return fig

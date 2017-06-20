# Author:      Robby Boney, robby.boney@wsu.edu
# License:     ...
# Last Update: May 23, 2017
# Start Date:  Feb 20, 2016
# Version:     alpha 1.1

import os
import shutil
import time
from time import gmtime, strftime

from Bio.Align.Applications import MuscleCommandline
from io import StringIO as strIO
from Bio import AlignIO

# ////////////////////////////////////////////////////////////////
# Utility Functions //////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////
    
def TextToString(path):
    """ Takes a file path and outputs the described file as string """
    return open(path, 'r').read()

def SortProtein(alignmentfile, outputfile, orderCmdList):
    """
    OVERVIEW:  
        This function will take a protein sequence and
        then reorders it based on the order of the list (orderCmdList)
    USE of FUNCTION:
        alignment file ==> sortProtein ==> outputfile
    """

    AliF = open(alignmentfile, "r")
    AliF = AliF.read()
    OutF = open(outputfile, "w")

    temp = AliF.split(">")
    tempList = []
    output = ""

    for j in range(len(orderCmdList)):
        for i in range(1, len(temp)):
            if orderCmdList[j] in temp[i].splitlines()[0]:
                if temp[i] not in tempList:
                    tempList.append("\n" + temp[i])
                break
    output = output.join(tempList)
    return output

def RemoveSameGene(file):
    """
    OVERVIEW: 
        This function takes any output text file from MainPC_analysis or DomainHits_analysis
        and returns the same output excluding any duplicate genes leaving the gene with the highest 
        score on the far left.
    """

    file = open(file, "r")
    file = file.read()
    file = file.splitlines()

    genecounter = ""
    tophit = 0
    bestrows = "-"
    indnum = 0
    L = file[0].split()

    for i in range(len(L)):
        if ("DARG" in L[i]):
            indnum = i
            break

    for i in file:        
        L = i.split()
        gene = L[indnum].split("_")[0]

        if (gene != genecounter):
            genecounter = gene
            if (genecounter not in bestrows):    
                bestrows = genecounter + bestrows
                tophit = 0
                print(i)

        elif (float(L[0]) > tophit):
            tophit = float(L[0])
    return 0

def stats(f):
    """
    OVERVIEW:
        stats, is a decorator function for timing other functions
    NOTES:
        this decorator function has one package dependancy
        - time
    """
    def caller(*args, **kwds):
        t0 = time.time()
        ret = f(*args, **kwds)
        t1 = time.time()
        out = str(
        '''function title:    ''' + str(f.__name__) +
        '''\n run time(sec):    ''' + str(round(t1-t0, 5)))
        print (out)
        return ret
    return caller


# ////////////////////////////////////////////////////////////////
# Core Functions    //////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////
def MainPC_defineHit(RegenL, NonRegenL):
    """
     OVERVIEW: this function checks two groups of amino acids RL and NRL. and then determines if
        the amino acid site classifies as a hit or not. the function returns 1(yes) / 0(no)
     USE of FUNCTION:
        helper function for MainPC
    """
    return 0

def MainProteinCompare(alignmentF, groupsfile, outputtxt):
    """
    DATE_Started: Jan 6, 2017
    OVERVIEW: 
        this function is a refactor of the function "ProteinCompare".
        This will allow for simpler use of the function rather then using it
        as a helper function.
        - the goal of the Protein compare function is to compare amino acids from aligned seq's and
        find hits as defined by the user, in our case, hits where regenerator species are all similar
        and non-regenerators are different from regenerators.
        NOTES:
        -alignmentF = folder path of aligned protein seq's
        -groupsfile = .txt file path with every line containing (SpeciesID R/NR)
        -outputtxt = .txt file path where output will be written
        -stats = bool for printing stats on output (recommended: keep false)

    USE of FUNCTION:
        folder path -----> alignmentF --\
                                         |
        txt with type ---> groupsfile ---|----> outputtxt
                                         |
        [1/0] ---> print stats [T/F] ---/
    """

    hitDefinition = ""
    outputtxt = open(outputtxt, "w+")

    speID = open(groupsfile, "r")
    speID = speID.read()
    speID = speID.splitlines()

    for file in os.listdir(alignmentF):
        if file.endswith(".txt"):
            uafloc = alignmentF + file
            fileName = "" + file
            tempF = open(uafloc, "r")
            tempFText = tempF.read()
            tempFText_split = tempFText.split(">")
            t = "".join(tempFText_split[1].splitlines()[1:])  # seperates seq from id line
            aliLen = len(t)
            hits = 0

            for site in range(aliLen):
                regens = []
                nonregens = []
                for seq in range(1,len(tempFText_split)):
                    seqName = tempFText_split[seq].splitlines()[0]
                    seqSite = "".join(tempFText_split[seq].splitlines()[1:])[site]
                    for file in speID:
                        if file.split()[0] in seqName:
                            if file.split()[1] == "R":
                                regens.append(seqSite)
                                break
                            elif file.split()[1] == "NR":
                                nonregens.append(seqSite)
                                break
                #  ======== ======== ======== analysis Section ======== ======== ========
                # add species count in read out
                # definition of hit inside functional domain

                regenLjoin = "".join(regens)
                hitTest = 1

                if len(regenLjoin) > 0 and len(nonregens) > 0:

                    # checks if R sites are all the same and "-" chars
                    if regenLjoin == len(regenLjoin) * regenLjoin[0] and regenLjoin[0] is "-":
                        hitTest = 0

                    non_re = "".join(nonregens)
                    if non_re == len(non_re) * non_re[0] and non_re[0] is "-":
                        hitTest = 0
                    
                    RLsp = "".join( regenLjoin.split("-") )
                    if len(RLsp) < 2:
                        hitTest = 0
                    if len(RLsp) > 0 and RLsp != len(RLsp) * RLsp[0]:
                        hitTest = 0

                    for i in nonregens:                           # i = amino acid at site
                        if i in regenLjoin and i is not "-":      # definition of NOT a hit
                            hitTest = 0
                            break
                        
                    if hitTest == 1:
                        hits += 1
                        output1 = "\n|  site = " + str(site) + "  R: "+ str(regenLjoin) + " NR: "+ "".join(nonregens)
                        print(output1)
                        outputtxt.write(output1)

            output2 = "\nL__ GENE:"+ str(fileName)+\
                    "  ---total Length:" +\
                    str(aliLen) +\
                    " ---#hits ="+\
                    str(hits) +\
                    " R: " + str(len(regens)) + \
                    " NR: " + str(len(nonregens))
            outputtxt.write(output2)
            tempF.close()

def MainPC_analysis(hitstxt, outtxt):
    """
    DATE_Started: Jan 6 2017
    OVERVIEW:
        this function analyzes the output(txt file) from mainproteincompare
        as specified by user.
    USE of FUNCTION:
    MainProteinCompare ==> MainPC_analysis ==> textfile
    MainProteinCompare ==> MainPC_analysis ==> markdown
    """

    data = open(hitstxt, "r")
    outputFile = open(outtxt, "w")
    outputFile.write("GENERATED from MainPC_analysis on " + strftime("%Y-%m-%d %H:%M:%S", gmtime()) )
    data = data.read()
    datasplit = data.splitlines()
    tempL = []
    outL = []
    for line in datasplit:
        if "L__" in line:
            tempL.append(line)
    for line in tempL:
        item = line.split()  
        perc = round(int(item[5][1:]) /int(item[3][7:]), 3)

        # this is what gets printed to the output file
        outL.append( str("\n" +str(perc) + " = " + item[5][1:] + "/" + item[3][7:]+ " " + item[1][5:] +\
                         " " + "".join(item[6:8]) + " " + "".join(item[8:])   ))
    outL.sort()
    outL.reverse()
    for i in outL:
        if ( int(i.split()[5][3:]) > 0 and int(i.split()[4][2:]) > 0):         #filter for output
            outputFile.write(i)

def SeqFuncDomain_Fast(alignPath="none", combinationPath="none", funcDomPath="none", outtxt=""):
    """
    OVERVIEW: this function analyses a folder of alignments
    USE of FUNCTION:
    funcDomPath ==== > folder of folders containing functional domains for each species

    alignPath ======== > folder of alignment txt files
    combinationPath == > txt file of combinations (from CombinOfProID output)
    """
    outtxt = open(outtxt, "w")
    if alignPath != "none":
        alignP = TextToString(alignPath)
    if combinationPath != "none":
        combP = TextToString(combinationPath)
        combSplit = combP.splitlines()
    curGene = ""
    takeout = "[]',_"
    domDataBase = []
    domCountDataBase = []
    # compiles all domain ID's for each sp into 1 list of lists
    for folder in os.listdir(funcDomPath):
        newSp = []
        if not folder.endswith(".DS_Store"):
            folderLoc = funcDomPath + folder + "/"
            newSp.append(folder)
            for file in os.listdir(folderLoc):
                if file.endswith(".txt"):
                    fileloc = folderLoc + file
                    fileName = "" + file[:-4]  # gets file name without .txt
                    tempF = open(fileloc, "r")
                    tempFText = tempF.read()
                    tempFText = tempFText.splitlines()
                    for line in range(1, len(tempFText)):
                        newSp.append(tempFText[line])
        domDataBase.append(newSp)
    # goes through the compiled domains and sorts and counts the # of occurences
    for sp in domDataBase:
        if len(sp) > 0:
            newSp = []
            curdomain = ""
            curID = ""
            counter = 1
            spSort = sorted(sp)
            prntTrue = 0
            for ID in spSort:
                if len(ID.split()) > 1:
                    if ID != curID:
                        curID = ID
                        prntTrue = 1
                    if ID.split()[1] != curdomain:
                        curdomain = ID.split()[1]
                        prntTrue = 1
                    if prntTrue:
                        inputString = ID.split()[0] +" " + ID.split()[1] + " cnt=" + str(counter)
                        newSp.append(inputString)
                        prntTrue = 0
                        counter = 1
                    elif ID.split()[1] == curdomain and ID == curID:
                        counter += 1

        domCountDataBase.append(newSp)

    if alignPath != "none":
        alignP = TextToString(alignPath)
    if combinationPath != "none":
        combP = TextToString(combinationPath)
        combP = combP.splitlines()
    combinationDataBase = []
    for cb in combP:
        spComb = []
        if cb[0:3] == "ENS":
            curGene = cb.split()[0]
        else:
            # this takes out  artifacts in string
            for i in range(len(cb)):        
                if cb[i] in takeout:
                    cb = cb[:i] + " " + cb[i + 1:]
            cbSplit = cb.split()
            spComb.append(curGene + " " + cbSplit[0])
            for scp in cbSplit[1:]:
                for spdom in domCountDataBase:
                    #print(spdom[0][3:7], " --- ", scp[3:7])     # = GALP  ---  SSCP
                    #print(spdom[0].split()[0], " --- ", scp)    # = ENSGALP00000000209  ---  ENSSSCP00000022814
                    if spdom[0][3:7] == scp[3:7]:
                        for line in spdom:
                            if line.split()[0] == scp:
                                spComb.append(line.split()[1] +" " + line.split()[0] +" " + line.split()[2])
                            pass
            combinationDataBase.append(spComb)
    return combinationDataBase

def SeqFuncDomain_FastAnalysis(L, SPID):
    """
    OVERVIEW: this is a helper function for SeqFuncDomain_Fast. input L should be the list var
        returned from that function.
    """

    SPID = open(SPID, "r")
    SPID = SPID.read()
    SPIDsplit = SPID.splitlines()

    for i in L:
        isort = sorted(i[1:])
        curdom = ""
        curSP = ""
        NRL = []
        RL = []
        NRLid = []
        RLid = []
        for j in isort:
            curSP = j.split()[1]
            if curdom != j.split()[0]:
                if len(NRL) > 0 and len(RL) > 0:
                    if max(NRL) < min(RL):
                        print(i[0], " ", curdom, " NRL = ", NRL, "  RL = ", RL)
                curdom = j.split()[0]
                NRL = []
                RL = []
                NRLid = []
                RLid = []
            for line in SPIDsplit:
                if line.split()[0] in curSP:
                    if line.split()[1] == "R":
                        RL.append(int(j.split()[2][4:]))
                        RLid.append(j.split()[1][3:7])
                    elif line.split()[1] == "NR":
                        NRL.append(int(j.split()[2][4:]))
                        NRLid.append(j.split()[1][3:7])

def SeqFuncDomain(alignPath="none", combinationPath="none", funcDomPath="none", outtxt=""):
    """
    OVERVIEW: this function analyses a folder of alignments
    USE of FUNCTION:
    funcDomPath ==== > folder of folders containing functional domains for each species

    alignPath ======== > folder of alignment txt files
    combinationPath == > txt file of combinations (from CombinOfProID output)
    """
    outtxt = open(outtxt, "w")
    if alignPath != "none":
        alignP = TextToString(alignPath)
    if combinationPath != "none":
        combP = TextToString(combinationPath)
        combP = combP.splitlines()

    # === cycles through combinations  === ########
    curGene = ""
    takeout = "[]',_"
    outputLine = ""
    domains = []        #a list of the domain names as in the folder

    for cb in combP:
        if cb[0:3] == "ENS":
            curGene = cb.split()[0]
        else:
            for i in range(len(cb)):
                if cb[i] in takeout:
                    cb = cb[:i] + " " + cb[i+1:]
            curCombN = cb.split()[0]
            outputLine += curGene + " " + curCombN + "  "
            cb = cb.split()[1:]
            for sp in cb:                   # for each species in the combination
                outputLine += ">"+sp[3:7] + ": "
                spDomList = []
                if len(sp) < 10:
                    outputLine += "-- "
                else:
                    # === cycles through functional domains folders === ########
                    for folder in os.listdir(funcDomPath):
                        if not folder.endswith(".DS_Store") and sp[3:7] in folder:
                            folderLoc = funcDomPath + folder + "/"
                            for file in os.listdir(folderLoc):
                                spIdCount = 0
                                if file.endswith(".txt"):
                                    fileloc = folderLoc + file
                                    fileName = "" + file[:-4]         #gets file name without .txt
                                    if fileName not in "".join(domains):
                                        domains.append(fileName + " ")
                                    tempF = open(fileloc, "r")
                                    tempFText = tempF.read()
                                    tempFText = tempFText.splitlines()
                                    for line in tempFText:
                                        if sp in line and len(line.split()) > 1:
                                            spIdCount += 1
                                    outputLine += str(spIdCount) + " "
                            break
            cLine = outputLine.split(">")
            print (cLine[0])
            for sp in cLine:
                outtxt.write(sp + "\n")

        outputLine = ""
                # end === cycles through functional domains folders === ########
    # end === cycles through combinations === ########

def SeqFuncDomain_analysis(SFD, SPID, outtxt=""):
    """
    OVERVIEW: this is a helper function for SFD (seqFuncDomain) it takes the output from SFD and
        looks for hits by showing R to NR comparisons. it refines hit definition by filtering out
        instances where R's have more function domains
    USE of FUNCTION:
    seqFuncDomain output ==== > SFD ===> outtxt
    """

    SFD = open(SFD, "r")
    SFD = SFD.read()
    SFD = SFD.split("ENS")

    SPID = open(SPID, "r")
    SPID = SPID.read()
    SPIDsplit = SPID.splitlines()

    NRL = []
    RL = []
    domainsL = SFD[0].split()
    outtxt = open(outtxt, "w")
    # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= gets content of output
    curGene = ""
    c = 0
    for cmb in SFD:                # line = new combination
        
        if cmb == SFD[0]:
            print(domainsL)
            outtxt.write(str(domainsL) + "\n")
            for i in SPIDsplit:
                if i.split()[1] == "R":
                    RL.append(i.split()[0] + ":")
                elif i.split()[1] == "NR":
                    NRL.append(i.split()[0] + ":")
            final = str("R=="+ "".join(RL)+ "  NR=="+ "".join(NRL))
            print(final)
            outtxt.write(final + "\n")
        else:
            cmbL = cmb.splitlines()[1:]

            for dom in range(len(domainsL)):
                RL = []
                NRL = []
                for cmbSp in cmbL:
                    for SPline in SPIDsplit:
                        if cmbSp.split()[0][:4] in SPline.split()[0]:
                            if SPline.split()[1] == "R":

                                if len(list(cmbSp.split()[1:])) == len(domainsL):
                                    RL.append(list(cmbSp.split()[1:])[dom] + " ")
                                elif len(list(cmbSp.split()[1:])) == 1:
                                    RL.append(list(cmbSp.split()[1:])[0] + " ")

                            elif SPline.split()[1] == "NR":

                                if len(list(cmbSp.split()[1:])) == len(domainsL):
                                    NRL.append(list(cmbSp.split()[1:])[dom] + " ")
                                elif len(list(cmbSp.split()[1:])) == 1:
                                    NRL.append(list(cmbSp.split()[1:])[0] + " ")
                #########       define HITS       #########
                if "".join(NRL).replace(" ", "") != "".join(NRL).replace(" ", "")[0] * len("".join(NRL).replace(" ", "")):

                    thrshd = 1
                    thrshdtest = 0
                    NRLmax = 0
                    RLlow = 0
                    newRL = "".join(RL).replace("--", "0").split()
                    for i in range(len(newRL)):
                        newRL[i] = int(newRL[i])
                    RLsort = sorted(newRL)
                    for i in NRL:
                        if i != "-- " and int(i) > NRLmax:
                            NRLmax = int(i)

                    for i in range(len(RLsort)):
                        if RLsort[i] != "-- ":
                            if int(RLsort[i]) > int(RLlow) and thrshdtest < thrshd:
                                RLlow = RLsort[i]
                                thrshdtest += 1
                            if thrshdtest == thrshd:
                                break
                    if int(NRLmax) < int(RLlow):
                        comparison = str("NRLmax = "+ str(NRLmax) + " RLlow = "+ str(RLlow))
                        final = str(str(cmb.splitlines()[0])+ str(domainsL[dom])\
                                    +"   R===" + "".join(RL) + " NR===" + "".join(NRL)+\
                                    " "+ comparison)
                        print(final)
                        outtxt.write(final + "\n")
    return 0

def FunctionDomainLen(ensemblHits):
    """ FORMAT DEMO          Database     site of interest          reg.Start       reg.End
    ENSDARP00000062559     2.60.40.10       174                        155             217
    ENSDARP00000062559     2.60.40.10       310                        218             312
    ENSDARP00000062559     2.60.40.10       383                        313             408
    :param ensemblHits: String of ALL databases and ALL functional domains
    :return: The length of the functional domain
    """

    j = ensemblHits.splitlines()
    currentId = ""
    length = 0
    funcPoints = []
    output = []
    tempOut = ""
    for i in range(len(j)):
        temp = j[i].split()
        if currentId != temp[0]:
            currentId = temp[0]
            length = 0
            funcPoints.clear()
        if currentId == temp[0]:
            for i in range(int(temp[3]), int(temp[4])+1):
                if i not in funcPoints:
                    funcPoints.append(i)

        output.append(temp[0] + '  ' + temp[3] + "  " + temp[4] + "  len: " + str(len(funcPoints)))

    out = []
    for i in range(len(output)):
        tempAhead = output[i].split()

        if i > 0:
            tempAt = output[i - 1].split()

        elif i == 0:
            tempAt = output[0].split()

        if tempAt[0] != tempAhead[0]:
            out.append(output[i-1])

        if i == len(output)-1:
            out.append(output[len(output)-1])

    for i in range(len(out)):
        print(out[i])

    return length

def EnsemblHits(strToAnalyze, hitList):
    """
    This Function takes a string to analyze (ensembl database output) and a hit list from the protein compare function
    output. It then compares the two and outputs a list of results
    """

    lineiterator = strToAnalyze.splitlines()

    STAtemp = lineiterator
    HLtemp = hitList
    temp = []
    STAcomTem = ""
    j = 1

    for i in range(1, len(STAtemp)):

        STAcomTem = STAtemp[i].split()
        temp = HLtemp[j].split()

        if temp[0] not in STAcomTem[0]:
            j += 1
            temp = HLtemp[j].split()

        for ii in range(1, len(temp)):
            if len(STAtemp[i].split()) > 1:
                if int(temp[ii]) >= int(STAtemp[i].split()[2]) and int(temp[ii]) <= int(STAtemp[i].split()[3]):
                    print(STAtemp[i].split()[0], "   ", STAtemp[i].split()[1],"  ",  temp[ii], "  ", "   ", STAtemp[i].split()[2], "  ", STAtemp[i].split()[3])
                    break

def DataList(dir, idList):
    """
    OVERVIEW: his function will take a folder of text files and find all files with the same prefix,
        then compile those files into 1 string.
    Loc = /Volumes/HDD/-Apps/projects_Python/ProteinConvergence2016/Resources/Protein_Domains-2016-07-26/Protein/Domains

    """

    #first compile all .txt's with the same prefix
    print(idList)
    fileList = []
    outlist = []

    for file in os.listdir(dir):
        if file.endswith(".txt") and file:
            place = (dir + file)
            fo = open(place, 'r')
            fileList.append(fo.read())

    for j in range(len(fileList)):
        for i in range(len(fileList[j].splitlines())):
            if (len(fileList[j].splitlines()[i].split()) == 3):
                tempTitle =  fileList[j].splitlines()[i].split()[1].upper()
                tempString = tempTitle + "  " + fileList[j].splitlines()[i].split()[2] + "  " + fileList[j].splitlines()[i].split()[0]
                outlist.append(tempString)

    ##### SORT THE OUTPUT ####
    outlist = sorted(outlist)

    # analyze output 
    # this loops through each gene and compares regenerators then compares non regenerators
    # then cross compares the hits from both comparisons together

    geneListSearched = []
    proteinListSearched = []
    
    #locators will be defined when a new instance is encountered
    geneLocator = ""
    proteinLocator = ""
    numberCountList = [0] * len(idList)
    for gene in range(len(outlist)):
        if (str(geneLocator) != str(outlist[gene].split()[0])):
            geneLocator = outlist[gene].split()[0]

            for pro in range(gene, len(outlist)):
                if (str(geneLocator) != str(outlist[pro].split()[0])):
                    break

        #########################   Check for a hit, then +1 the row of that species ####################
                for xp in range(len(numberCountList)):
                    if idList[xp] in str(outlist[pro].split()[2]):
                        numberCountList[xp] += 1
        #################### Make list of regenerators, then count to check that at least 4 have the same count ####################
                regenCount = 0
                for num in range(2, len(numberCountList)):
                    if numberCountList[0] < numberCountList[num] and numberCountList[1] < numberCountList[num]:
                        regenCount += 1
                if (pro < len(outlist)-1 and str(outlist[pro].split()[1]) != str(outlist[pro+1].split()[1])):     
                    # print  out selective comparisons
                    if regenCount > 3:
                        print(outlist[pro].split()[0], "  ",outlist[pro].split()[1], "  ", numberCountList)
                        True
                    numberCountList = [0] * len(idList)
    # if the gene is not in geneListSearched
    return True

def comb_TrPr(id_dir, spgr):
    """
    OVERVIEW
        combines a folder of transcript -- protein associations
        into a long list.
    INPUTS
        alignments_dir = path to folder with alignments
        spgr           = string for species id to look for (i.e. "DARP")
    """
    tr_pr_L = []
    for file in os.listdir(id_dir):
        
        if file.endswith(".txt"):
            fileloc = id_dir + file
            fileName = "" + file
            fileString = open(fileloc, "r").read()
            fileString_split = fileString.splitlines()
            
            for line in fileString_split:
                ln = line.split()
                if len(ln) == 2:
                    tr_pr_L.append(str(ln[0] + "   " + ln[1]))
    return tr_pr_L

def comb_rm_dups(tr_pr_L, names):
    """
    OVERVIEW
        re-sorts the list (tr_pr_L) and removes duplicate entries
    INPUTS
        tr_pr_L = list of ( transcipts - protein id's) 
                  from comb_TrPr()
        names   = list of ordered speices
    """
    spacer = "                  "
    spacer_cnt = 14
    dataCol = []
    curgene = ""
    masterLForG = []
    testString = ""

    for i in tr_pr_L:
        if curgene != i.split()[0]:
            curgene = i.split()[0]
            for msti in masterLForG:
                dataCol.append(list(msti))
            masterLForG.clear()
            masterLForG.append(list())
        for row in masterLForG:
            if len(row) == 0 and masterLForG.index(row) == 0:
                row.append(curgene)
            elif len(row) == 0 and masterLForG.index(row) != 0:
                row.append(spacer)
            if i.split()[1][:7] in testString.join(row) and masterLForG.index(row) == len(masterLForG)-1:
                masterLForG.append(list())
            elif i.split()[1][:7] not in testString.join(row):
                row.append(i.split()[1])
                break

    for i in dataCol:
        print(i)
        for nme in names:
            if ''.join(i).find(nme) == -1:
                if "DART" in ''.join(i):
                    i.append(nme + " "*spacer_cnt)
                else:
                    i.append(nme +     "-"*spacer_cnt)

        i[1:] = sorted(i[1:])
    return dataCol

def comb_gen_combs(mstrList, out_fl):
    """
    OVERVIEW
        This Function prints out the combinations of the output from SortDuplicates.
        using data not outputed from sortedduplicates will most likely not return valid
        output.
    INPUTS
        mstrList = List object from comb_rm_dups() output
    """

    out_fl = open(out_fl, "w")
    reg = mstrList[0][0]
    mstL = []
    tmpL = []
    for i in mstrList:
        if "DART" in i[0] and i[0] != reg :
            reg = i[0]
            mstL.append(list(tmpL))
            tmpL.clear()
        tmpL.append(i[1:])
        if "DART" in i[0] and mstrList.index(i) == len(mstrList)-1:
            reg = i[0]
            mstL.append(list(tmpL))
            tmpL.clear()
            
    startN = 0
    n = 0
    for num in range(0,len(mstL)):

        new_egg = zip(*mstL[num])
        new_egg = ([x for x in tup if "-" not in x] for tup in new_egg)
        combb = 1
        for i in new_egg:
            if len(i) != 0:
                combb *= len(i)
        print(mstrList[n][0], "  combs: " , combb)

        new_egg = zip(*mstL[num])
        new_egg = ([x for x in tup if "-" not in x] for tup in new_egg)
        combsCount = []
        each = [[]]
        print(str(mstrList[n][0]  + "  combs: " + str(combb) + "\n"))
        #out_fl.write(str(mstrList[n][0]  + "  combs: " + str(combb) + "\n"))
        n += 1
        while n < len(mstrList) and "DART" not in mstrList[n][0]:
            n += 1

        for l in new_egg:
            neach = [x + [t] for t in l for x in each]
            each = neach
        for e in each:
            print(str(str(each.index(e)) + "   " + str(e) + "\n"))
            #out_fl.write(str(str(each.index(e)) + "   " + str(e) + "\n"))
    out_fl.close()
    return True

def DomainHits_GeneProtein(alignments_dir, SPID):
    """
    OVERVIEW: Aquires the proteinID that associates with its GeneID

    USE of FUNC:
        alignments_dir ===> path to folder with alignments
        SPID =============> string for species id to look for (i.e. "DARP")

    RETURNS:    list of "geneID   ProID"
    """

    g_pL = []    #gene_protein list
    curGene = " "
    for file in os.listdir(alignments_dir):
        if file.endswith(".txt") and curGene not in str(file):
            fileloc = alignments_dir + file
            fileName = "" + file
            fileString = open(fileloc, "r")
            fileString = fileString.read()
            fileString_split = fileString.splitlines()
            
            curGene = fileName.split(sep="_")[0]
            curPro = ""
            for line in fileString_split:
                if ">" in line and SPID in line:
                    g_pL.append(str(curGene + "   " + line[1:]))
    return g_pL

def DomainHits(domain_dir, alignments_dir, mainpc_file, output_file):
    """
    OVERVIEW: 
        this function takes a folder of text files with domain content (gene, database, startsite, endsite)
    and an output txt from MainPCProteinCompare. an procedure finds seq hits which fall inside functional
    domains and outputs the results to a text file.

    USE of FUNC:
        domain_dir =====> folder path containing ensembl domain output
        mainpc_file ====> txt path from MainPC_proteinCompare output
    """

    output_file = open(output_file, "w")

    # STEP 1 ---> create the domainList
    domainList = []
    for file in os.listdir(domain_dir):
        if file.endswith(".txt"):
            fileloc = domain_dir + file
            fileName = "" + file
            fileString = open(fileloc, "r")
            fileString = fileString.read()
            fileString_split = fileString.splitlines()
            
            for line in fileString_split:
                if len(line.split()) > 1 and "ENS" in line:
                    domainList.append(line)
    domainList = sorted(domainList)

    # STEP 2 ---> collect and rearrange gene and hits

    mainpc = open(mainpc_file, "r")
    mainpc = mainpc.read()
    mainpc_split = mainpc.split("GENE")
    db_GP = DomainHits_GeneProtein(alignments_dir, SPID="DARP")
    
    curGene = ""
    curPro = ""
    gp_cntr = 0
    gpTF = 0
    tempC = 0
    tempC_bkmk = tempC
    
    for i in mainpc_split:
        psl = i.splitlines()
        if (len(psl) > 1):
            curGene = psl[0].split()[0][1:-4]
            gpTF = 0
            for gp in range(tempC, len(db_GP)):
                gene = db_GP[gp].split()[0]      
                if (gene in curGene):
                    print(gene, " ", curGene)
                    curGene = gene
                    curPro = db_GP[gp].split()[1]
                    gpTF = 1
                    tempC_bkmk = tempC
                    break
                tempC += 1

            if (gpTF == 0):
                curPro = "ERROR: no corresponding protein"
                tempC = tempC_bkmk

    # STEP 3 ---> check if sites are inside functional domains
            
            site = 0
            for line in psl:  
                if (line == psl[0]):
                    # check for at least 2R and 2NR
                    if (int(line.split()[6]) < 2 or int(line.split()[8]) < 2):     
                        break
                    output_file.write(str(line+ " "+ curPro)+"\n")
                    print (line, " ", curPro)
                    pass
                if (not gpTF):
                    output_file.write(str(line+ " "+ curPro)+"\n")
                    break    
                if "|" == line[0]:
                    site = int(line.split()[3])
                    for num in domainList:
                        if (num.split()[0] == curPro):
                            lower = int(num.split()[2])
                            upper = int(num.split()[3])
                            if (lower < site and site < upper):
                                output_file.write(str(curPro + " "+ str(site) +"  " + str(lower) + "-" + str(upper))+"\n")
                                break        
    return 0

def DomainHits_analysis(dm_output_file):
    """
    OVERVIEW: This function takes the file outputed by the core function and organizes the data
        to be more readable.
    """
    dm = open(dm_output_file, "r")
    dm = dm.read()
    dm = dm.splitlines()
    curGene = ""
    prevGene = ""
    hitsCnt = 0
    collection = []
    rangesL = []
    regionLen = 0
    for line in dm:
        
        if(line):
            if(line[0] == ":" or line[0] == "["):
                testg = line.split()[0][1:-4]
                if (testg != curGene):
                    
                    # ================ get length of dom region ================ start
                    if (rangesL):                   # rangesL[0] = [1,2]
                        rng_low = rangesL[0][0]     # where 1 = lower, and 2 = upper
                        rng_hi = rangesL[0][1]
                        regionLen = 0
                        for i in rangesL:
                            if (i[0] > rng_hi):         # RNGlow is > than highest ranger return data and start range
                                regionLen += rng_hi - rng_low
                                rng_hi = i[1]
                                rng_low = i[0]
                            if (i[0] < rng_hi < i[1]):
                                rng_hi = i[1]
                        regionLen += rng_hi - rng_low

                    # ================ get length of dom region ================ end
                    if (regionLen > 0):
                        collection.append([round(hitsCnt/regionLen, 5), str("DMhits=" + str(hitsCnt) ), str("DomLen=" + str(regionLen) ), "Gene/Pro: ",curGene, prevGene])
                    else:
                        collection.append([regionLen, str("DMhits=" + str(hitsCnt) ), str("DomLen=" + str(regionLen) ), "Gene/Pro: ",curGene, prevGene])
                    hitsCnt = 0
                    curGene = testg
                    rangesL.clear()

                prevGene = " ".join(line.split()[2:])
            elif (line and line[0] == "E"):
                hitsCnt += 1
                rng = line.split()[2].split("-")
                rangesL.append( [int(rng[0]), int(rng[1])] )       
    genecounter = ""
    tophit = 0
    bestrows = "-"
    for i in sorted(collection, reverse=1):        
        if (i[4].split("_")[0] != genecounter):
            genecounter = i[4].split("_")[0]
            if (genecounter not in bestrows):    
                bestrows = genecounter + bestrows
                tophit = 0
                print("   ".join( [str(k) for k in i] ) )
        elif (i[0] > tophit):
            tophit = i[0]
    return 0

def SortDuplicates(listGP, names):
    """
    OVERVIEW
        re-sorts the list input arranging the non-zebrafish proteins together,
        accounting for when zebrafish proteins and genes are duplicates.

    INPUTS
        listGP = list of ( DAR gene - Species pro ID) from DomainHits_GeneProtein()
        names  = list of ordered species id's

    NOTES
        DARG DARP species proteins......
    """

    spacer = "                  "
    dataCol = []
    curgene = ""
    masterLForG = []
    testString = ""

    for i in listGP:
        if curgene != i.split()[0]:
            curgene = i.split()[0]
            for msti in masterLForG:
                dataCol.append(list(msti))
            masterLForG.clear()
            masterLForG.append(list())
        for row in masterLForG:
            if len(row) == 0 and masterLForG.index(row) == 0:
                row.append(curgene)
            elif len(row) == 0 and masterLForG.index(row) != 0:
                row.append(spacer)
            if i.split()[1][:7] in testString.join(row) and masterLForG.index(row) == len(masterLForG)-1:
                masterLForG.append(list())
            elif i.split()[1][:7] not in testString.join(row):
                row.append(i.split()[1])
                break

    for i in dataCol:
        for nme in names:
            if ''.join(i).find(nme) == -1:
                if "DARG" in ''.join(i):
                    i.append(nme + "           ")
                else:
                    i.append(nme +     "-----------")

        i[1:] = sorted(i[1:])
    return dataCol

def CombinOfProID(mstrList):
    """
    This Function prints out the combinations of the output from SortDuplicates.
     using data not outputed from sortedduplicates will most likely not return valid
     output.

     /Volumes/HDD/-Apps/projects_Python/ProteinConvergence2016/Resources/CombOutPut.txt
    """

    combOut = open("/Volumes/HDD/-Apps/projects_Python/ProteinConvergence2016/Resources/CombOutPut.txt", "w")
    reg = mstrList[0][0]
    mstL = []
    tmpL = []
    for i in mstrList:

        if "DARG" in i[0] and i[0] != reg :
            reg = i[0]
            mstL.append(list(tmpL))
            tmpL.clear()
        tmpL.append(i[1:])
        if "DARG" in i[0] and mstrList.index(i) == len(mstrList)-1:
            reg = i[0]
            mstL.append(list(tmpL))
            tmpL.clear()
    startN = 0
    n = 0
    for num in range(0,len(mstL)):

        new_egg = zip(*mstL[num])
        new_egg = ([x for x in tup if "-" not in x] for tup in new_egg)
        combb = 1
        for i in new_egg:
            if len(i) != 0:
                combb *= len(i)
        print(mstrList[n][0], "  combs: " , combb)

        """
        This following code writes the output to a text file for further use in the analysis

        new_egg = zip(*egg)
        new_egg = ([x for x in tup if "-" not in x] for tup in new_egg)
        each = [[]]
        for l in new_egg:
            neach = [x + [t] for t in l for x in each]
            each = neach
        for e in each:
            print(e)

        """

        new_egg = zip(*mstL[num])
        new_egg = ([x for x in tup if "-" not in x] for tup in new_egg)
        combsCount = []
        each = [[]]
        combOut.write(str(mstrList[n][0]  + "  combs: " + str(combb) + "\n"))
        n += 1
        while n < len(mstrList) and "DARG" not in mstrList[n][0]:
            n += 1

        for l in new_egg:
            neach = [x + [t] for t in l for x in each]
            each = neach
        for e in each:
            combOut.write(str(str(each.index(e)) + "   " + str(e) + "\n"))
    combOut.close()
    return True
 
def combGeneCompile(combfile, unalignedfiles, outputFolder):
    """
    OVERVIEW: this function takes a .txt file for combinations from CombinOfProID and makes a new .txt
     containing the associated unaligned protein sequences of that combination
    USE of FUNCTIONS:
        SortDuplicates => CombinOfProID => ( combGeneCompile ) => bioMuscleAlign => proteinCompare

    NOTES:
        # iterate through all combinations
        # make new .txt for combination to output to. file title should be <geneID><#combination>.txt
        # ∞∞ iterate through IDs of combination
            # if ID is included, go to folder of UnalignedSeq's
                # iterate through seq's till correct ID is found.
                # then iterate through that ID's file till correct gene is found
                # copy species ID, gene ID, and seq to current output file
            # ∞∞ move on to next spec. ID in combination
        # after all combination spec. info have been copied to outputgene file, go to next combination
    """

    combfile = TextToString(combfile)           # here we define vars
    combfile = combfile.splitlines()
    currentGene = ""
    errorLogFileLoc = outputFolder + "1_errorLog" + ".txt"
    errorLog = open(errorLogFileLoc, "w+")

    num = 0

    for gene in combfile:                       # here we iterate through combination file
        gene = gene.replace(",","")
        gene = gene.replace("[", "")
        gene = gene.replace(",", "")
        gene = gene.replace("]", "")
        gene = gene.replace("'", "")
        gene = gene.replace("_", "")
        temp = gene.split()
        if "ENS" in temp[0]:
            currentGene = temp[0]
            print(currentGene)
        else:
            outloc = outputFolder + currentGene + "_" + temp[0] + ".txt"
            outfile = open(outloc, "w+")
            for sp in temp:
                if len(sp) > 10:
                    for file in os.listdir(unalignedfiles):
                        if file.endswith(".txt"):
                            uafloc = unalignedfiles + file
                            tempF = open(uafloc, "r")
                            tempFText = tempF.read()
                            curFileName = os.path.basename(tempF.name).split(".txt")[0]
                            if curFileName in sp:
                                tempFsplit = tempFText.split(">")
                                found = 0
                                for item in tempFsplit:
                                    if sp in item:
                                        found = 1
                                        item = ">" + item
                                        outfile.write(item)
                                if found == 0:
                                    errorLog.write("\n"+currentGene + "_" + temp[0] + " :: " + sp + " unable to locate in seq list")
            num += 1

def bioMuscleAlign(inputF, musclePath, outputF=""):
    """
    OVERVIEW: returns aligned .fa files of the unaligned .fa files from "inputF" folder in the same
        folder using the MUSCLE algorithm for peptide alignment.
    USE of FUNCTION:
        1. copy unaligned fasta (.fa) files to new folder
        2. use that folder path as inputF
        3. function will overwrite all .fa files with aligned seq's
    """

    if (outputF == ""):
        outputF = inputF
    else:
    # Copy files from source folder to output folder
        src_files = os.listdir(inputF)
        for file_name in src_files:
            full_file_name = os.path.join(inputF, file_name)
            if (os.path.isfile(full_file_name)):
                shutil.copy(full_file_name, outputF)

    # Run Alignment        
    muscExe = open(musclePath)
    for file in os.listdir(inputF):
        path = outputF + file
        input_file = open(path)
        cline = MuscleCommandline(input=path, out=path)
        command = str(cline)
        command = command.replace("muscle", musclePath)
        os.system(command)

def checkForAlignment(folder):
    """
    OVERVIEW: this function simply looks at alignment files in "folder" and checks for hyphens to see if any
     files were/are not aligned.
    USE of FUNCTION: bioMuscleAlign ==> ( checkForAlignment )
    """
    
    for file in os.listdir(folder):
        if file.endswith(".txt"):
            uafloc = folder + file
            fileName = "" + file
            tempF = open(uafloc, "r")
            tempFText = tempF.read()
            if "-" not in tempFText:
                i = tempFText.count(">")
                if i == 1:
                    print(fileName, " may not be aligned.!!! -only 1 id line found")
                else:
                    print(fileName, " is not aligned.!!!")



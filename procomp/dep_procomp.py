# Author:      Robby Boney, robby.boney@wsu.edu
# Last Update: May 23, 2017
# Start Date:  Feb 20, 2016

__name__ = "Procomp Analysis Module"

import os
import shutil
import time
from time import gmtime, strftime
import threading
import time

from Bio.Align.Applications import MuscleCommandline
from io import StringIO as strIO
from Bio import AlignIO

def MainProteinCompare(alignmentF, groupsfile, g1, g2, outputtxt):
    """
    OVERVIEW: 
        this function is a refactor of the function "ProteinCompare".
        This will allow for simpler use of the function rather then using it
        as a helper function.
        - the goal of the Protein compare function is to compare amino acids from aligned seq's and
        find hits as defined by the user, in our case, hits where regenerator species are all similar
        and non-regenerators are different from regenerators.
    INPUTS:
        alignmentF  = (string) folder path of aligned protein seq's
        groupsfile  = (string) .txt file path with every line containing (SpeciesID R/NR)
        g1          = (string) identifier for species in groupsfile for group 1
        g2          = (string) identifier for species in groupsfile for group 2
        outputtxt   = (string) .txt file path where output will be written
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

            # Iterate over the sites in the amino acid seq for all species
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
    OVERVIEW:
        this function analyzes the output(txt file) from mainproteincompare
        as specified by user.
    USE:
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

def comb_IdPr(id_dir, name_filter="", print_computed=0):
    """
    OVERVIEW:
        combines a folder of Id -- protein associations
        into a long list. Id are typically (XXXG, XXXT, XXXP)
    INPUTS:
        alignments_dir = path to folder with alignments
        spgr           = string for species id to look for (i.e. "DARP")
        name_filter    = you can choose to use only files with the given
                         string filter, default set to accept all files.
        print_computed = will print out the files used.
    """
    tr_pr_L = []
    for file in os.listdir(id_dir):
        if file.endswith(".txt") and name_filter in file:
            if print_computed:
                print(file)
            fileloc = id_dir + file
            fileName = "" + file
            fileString = open(fileloc, "r")
            fileString_split = fileString.read().splitlines()
            
            for line in fileString_split:
                ln = line.split()
                if len(ln) == 2:
                    if (len(ln[1]) != 18):
                        n = 18-len(ln[1]) 
                        t = str(ln[1]) + "   "        
                    else:
                        t = ln[1]
                    tr_pr_L.append(str(ln[0] + "   " + t))
            fileString.close()
    print("comb_TrPr() has finished")
    return tr_pr_L

def comb_rm_dups(tr_pr_L, names, ident):
    """
    OVERVIEW:
        re-sorts the list (tr_pr_L) and removes duplicate entries
    INPUTS:
        tr_pr_L = list of ( transcipts - protein id's) 
                  from comb_TrPr()
        names   = list of ordered speices
        ident = the identifier to be used for the combinations, typically it will be 
            of the format (DARG, DART, or DARP)
    """
    spacer = "                  "
    spacer_cnt = 11
    dataCol = []
    curgene = ""
    masterLForG = []

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

            if i.split()[1][:7] in "".join(row[1:]) and masterLForG.index(row) == len(masterLForG)-1:
                masterLForG.append(list())
            elif i.split()[1][:7] not in "".join(row[1:]):
                row.append(i.split()[1])
                break

    # add in spacers if certain speices dont show up in 
    # previous processing.
    #spid_tb = spid_tb_gen(spid_path)
    for i in dataCol:
        for nme in names:
            if ''.join(i[1:]).find(nme) == -1:
                if ident in i[0]:
                    i.append(nme + " "*spacer_cnt)
                else:
                    i.append(nme + "-"*spacer_cnt)

        i[1:] = sorted(i[1:])

    print("comb_rm_dups() has finished")
    return dataCol

def comb_gen_combs(mstrList, spid, out_fl, thr_tr, ident, w=0, refac=1):
    """
    OVERVIEW:
        This Function prints out the combinations of the output from SortDuplicates.
        using data not outputed from sortedduplicates will most likely not return valid
        output.
    INPUTS:
        mstrList = List object from comb_rm_dups() output
        spid = path to the species group file
        out_fl = file where output will be written to
        thr_tr = threshold for how many combinations each transcript can have.
        thr_sp = threshold for how many combinations each species can have.
        ident = the identifier to be used for the combinations, typically it will be 
            of the format (DARG, DART, or DARP)
    """
    def _longest(L):
        """ L is list of lists """
        longest = len(L[0])
        for i in L:
            if len(i) > longest:
                longest = len(i)
        return longest

    def _get_gr_longest(L, gr, tb):
        """ returns the species with the most orthologs for a certain group """
        longest = len(L[0])
        for i in L:
            if spid_tb_get_group(i, tb) == gr:
                if len(i) > longest:
                    longest = len(i)
        return spe_id

    def _check_ens_spe(ens_id, spid_path):
        """ returns the group the ensembl id is appart of """
        for item in spid_path:
            if ens_id in item[0]:
                return item[1]
        print ("{} was not found in the group list".format(ens_id))
        return None

    def _check_R_NR_miss(comb_L):
        """ returns the number of species for each group in the combination """
        R_cnt = 0
        NR_cnt = 0
        blank_cnt = 0 
        for i_spe in comb_L:            
            check = spid_tb_get_group(i_spe[0], spid_tb)
            if " " not in i_spe[0] and "-" not in i_spe[0]:
                if check == "R":
                    R_cnt += 1
                elif check == "NR":
                    NR_cnt += 1
            else:
                blank_cnt += 1
        return R_cnt, NR_cnt, blank_cnt
    
    # Set up table for species ids
    #
    spid_tb = spid_tb_gen(spid)

    reg = mstrList[0][0]
    mstL = []
    tmpL = []
    for i in mstrList:
        if ident in i[0] and i[0] != reg :
            reg = i[0]
            mstL.append(list(tmpL))
            tmpL.clear()
        tmpL.append(i[1:])
        if ident in i[0] and mstrList.index(i) == len(mstrList)-1:
            reg = i[0]
            mstL.append(list(tmpL))
            tmpL.clear()
            
    startN = 0
    n = 0
    if w:
        out_fl = open(out_fl, "w")
    ret_val = ["Transcript, Combinations, R (#), NR (#), Missing (#), Times Refactored"]
    for num in range(0,len(mstL)):
        new_egg = zip(*mstL[num])
        new_egg = [[x for x in tup if "-" not in x] for tup in new_egg]
        combb = 1
        for i in new_egg:
            if len(i) != 0:
                combb *= len(i)
            
        # -----------------------------
        # if too many combinations exist for this transcript,
        # exclude species with large numbers of orthologs while
        # assuring that R and NR groups are not put off balance
        # 
        R_c = 0
        NR_c = 0
        higher_stack = ""
        miss_c = 0
        refac_c = 0
        spe_rm_c = 0
        while (refac and combb > thr_tr):
            comb_refac = 1
            R_c, NR_c, miss_c = _check_R_NR_miss(new_egg)
            if R_c > NR_c:
                higher_stack = "R"
            else:
                higher_stack = "NR"

            max_c = _longest(new_egg)
            for k in range(len(new_egg)):
                spe_id = new_egg[k][0]
                spe_gr = spid_tb_get_group(spe_id, spid_tb)

                # Remove particular species from the combination
                if len(new_egg[k]) == max_c:
                    new_egg[k] = [str(spe_id[:7] + " "*11)]
                    max_c = _longest(new_egg)
                    spe_rm_c += 1

                # Calculate combination count
                if len(new_egg[k]) != 0:
                    comb_refac *= len(new_egg[k])
            if comb_refac < combb:
                combb = int(comb_refac)
            refac_c += 1
            toWrite = "     {} --- Refactored {} times, removed {} species, and {} combs".format(mstrList[n][0], refac_c, spe_rm_c, combb)
            ret_val.append(toWrite)
        
        # ----------------------------
        # check to see how many species fall in R and NR
        #
        R_c, NR_c, miss_c = _check_R_NR_miss(new_egg)
        
        combsCount = []
        each = [[]]
        toWrite = "{},{},{},{},{},{}".format(mstrList[n][0], combb, \
         R_c, NR_c, miss_c, refac_c)
        ret_val.append(toWrite)
        if w:
            #toWrite = str( mstrList[n][0] + "  combs: {}".format(combb) )
            out_fl.write(toWrite + "\n")
        n += 1
        while n < len(mstrList) and ident not in mstrList[n][0]:
            n += 1

        if w:
            for l in new_egg:
                neach = [x + [t] for t in l for x in each]
                each = neach
            
            for e in each:
                toWrite = str(str(each.index(e)) + " ," + str(",".join(e)) + "\n")
                out_fl.write(toWrite)
    if w: 
        out_fl.close()
    return ret_val

def comb_gen_pro_files(combfile, seq_dir, out_dir):
    """
    OVERVIEW: 
        this function takes a .txt file for combinations from CombinOfProID and makes a new .txt
        containing the associated unaligned protein sequences of that combination
    USE:
        combfile = path to output .txt from comb_gen_combs
        seq_dir = path to dir with unaligned protein sequences
        out_dir = path to dir where combination files will be written to
        
        comb_IdPr => comb_rm_dups => comb_gen_combs => bioMuscleAlign => MainProteinCompare
    """

    # Set up hash table of species ids and sequences as well as open
    # combination files and set inital variables
    #
    seq_tb = gen_seq_hash_tb(seq_dir)
    combfile = open(combfile, 'r').read().splitlines()  
    errorLog = []
    cur_comb = ""
    num = 0

    # here we iterate through the combination file
    for comb in combfile:   

        # This is a threshold of how many combinations to iterate over   
        if num >= 100:
            errorLog.append("END PROCESS")
            return errorLog

        # convert current combination into list from CSV style
        temp = comb.split(",")

        # if current line is an id line, set current transcript id
        # to this id
        if "ENS" in temp[0]:
            cur_comb = temp[0].split()[0]
            print(cur_comb)

        # if combination line generate combination file by accessing 
        # the dictionary
        else:
            outloc = out_dir + cur_comb + "_" + temp[0] + ".txt"
            outfile = open(outloc, "w+")
            
            for sp in temp:
                if len(sp) >= 7:
                    seq = get_seq(sp, seq_tb)
                    if seq != -1:
                        outfile.write(">" + sp + "\n")
                        outfile.write(seq + "\n")
            num += 1
            outfile.close()
    return errorLog

def FunctionDomainLen(ensemblHits):
    """ 
    OVERVIEW:
        takes a file of ensembl domains and analyzes those domains for their length.
    FORMAT DEMO          Database     site of interest          reg.Start       reg.End
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
    OVERVIEW:
        This Function takes a string to analyze (ensembl database output)
        and a hit list from the protein compare function output.
        It then compares the two and outputs a list of results.
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
    OVERVIEW:
        takes a folder of text files and find all files with the same prefix,
        then compile those files into 1 string.
    NOTES:
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
            
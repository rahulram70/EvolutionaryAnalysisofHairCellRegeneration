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
import pandas
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from io import StringIO as strIO
    
def SortProtein(alignmentfile, outputfile, orderCmdList):
    """
    OVERVIEW:  
        This function will take a protein sequence and
        then reorders it based on the order of the list (orderCmdList)
    USE:
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

def RemoveSameGene(file, ident):
    """
    OVERVIEW: 
        This function takes any output text file from MainPC_analysis or DomainHits_analysis
        and returns the same output excluding any duplicate genes leaving the gene with the highest 
        score on the far left.
    USE:
        file  = path to and output file from MainPC_analysis or DomainHits_analysis
        ident = unique id to use for checking lines (i.e. "DARG" or "DART")
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
        if (ident in L[i]):
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

def spid_tb_gen(spid_path):
    """ Generates a Hash Table for all speices and their group association """
    table = {}
    spid = open(spid_path, "r")
    spid_L = spid.read().splitlines()
    spid.close()
    for i in spid_L:
        species = i.split()[0]
        group = i.split()[1]
        table[species] = group
    return table

def spid_tb_get_group(spe_id, tb):
    """ Returns the group a particular species is associated with """
    key = spe_id[:7]
    if key[-1] == "T":
        key = key[:6] + "P"
    if key == "ENST000":
        key = "ENSP000"
    try:    
        return tb[key]
    except:
        print("ERROR: cannot find ({}) in given table".format(key))
        return None

def gen_seq_hash_tb(l_seq_dir):
    """ generates a 2d dictionary for the protein sequences """
    seq_table = {}
    for file in os.listdir(l_seq_dir):
        if file.endswith(".txt"):
            uafloc = l_seq_dir + file
            tempF = open(uafloc, "r")
            tempFText = tempF.read()
            cur_sp_id =  tempFText.splitlines()[0].split()[0][1:7]
            seq_table[cur_sp_id] = {}
            for protein in tempFText.split(">"):
                L = protein.split()
                if len(L) == 2:
                    seq_table[cur_sp_id][protein.split()[0]] = protein.split()[1]
            
            tempF.close()
    return seq_table

def get_seq(seq_id, tb):
    """ returns the protein sequence for a given id """
    try:
        return tb[seq_id[:6]][seq_id]
    except:
        return -1

def algn_tb_gen(algn_path):
    """ This function takes an alignment file and generates
    a list of species and their sequences in tuples with:
    (spe_id, spe_seq) """
    with open(algn_path , "r") as fl:
        ret_val = []
        for spe in fl.read().split(">"):
            spe_lines = spe.splitlines()
            if ( len(spe_lines) > 0):
                spe_id = spe_lines[0]
                spe_seq = ""
                for line in spe_lines[1:]:
                    spe_seq += line
                ret_val.append( (spe_id, spe_seq) )
        return ret_val

def protein_tb_set(pro_path):
    """ makes a dict for a file with protein seq's from ncbi """
    """ INPUTS:
            pro_path = csv file with ncbi protein data for one species     
    """

    file = open(pro_path, "r")
    file_rd = file.read().splitlines()
    file.close()
    ret_val = {}
    cols = {}

    # Setup headers from csv file
    #
    titles_L = file_rd[0].split(",")
    for val in range( len(titles_L) ):
        cols[val] = titles_L[val]
    
    # Assemble proteins in 2d dictionary
    #
    for pr in file_rd[1:]:
        line = pr.split(",")
        entry = {}
        for i in range( len(line[1:]) ):
            entry[ titles_L[i] ] = line[i]
        ret_val[ line[0] ] = entry  
    for i in ret_val:
        print(i)

def comp_cons_seq(in_fl):
    """Returns the conserved sequence for a post-alignment file where
    conserved is defined as the most common amino acid at each site"""
    
    # setup list of tuples of species
    algn_L = algn_tb_gen(in_fl)
    
    # get length of alignment (note: all alignments have the same length)
    rng = len( algn_L[0][1] )
    
    # compare each species to get conserved sequence
    cons_seq = ""
    for pt in range(rng):
        site_L = [ i[1][pt] for i in algn_L ]
        print(site_L)
        
        com = ["",0]
        for i in site_L:
            if site_L.count(i) > com[1] and i != com[0]:
                com[0] = i
                com[1] = site_L.count(i)
        cons_seq += com[0]
        del com
    return cons_seq


def comp_for_similarity(in_fl, comparator):
    """ takes an alignment file and compares each species for percent
    similarity.
    INPUTS:
             in_fl = path to alignment file
        comparator = string id of species to be compared with (i.e ENSDAR)
        
    """
    
    # setup list of tuples of species
    algn_L = algn_tb_gen(in_fl)
    
    # get length of alignment (note: all alignments have the same length)
    rng = len( algn_L[0][1] )

    # get index for comparator species and generate log list for each species
    cm_ind = 0
    log = []
    for i in range(len(algn_L)):
        if "ENSDARP0" in algn_L[i][0]:
            cm_ind = i
        log.append( [algn_L[i][0], 1] )        
    
    # compare each species to comparator seq and log
    for pt in range(rng):
        for spe in range(len(algn_L)):         
            comp_spe = algn_L[cm_ind]
            cur_spe = algn_L[spe]
            
            if (cur_spe[1][pt] == comp_spe[1][pt]):  
                log[spe][1] += 1
    
    # divide all similarities by length of alignment
    for spe in log:
        spe[1] = round(spe[1]/rng, 3)
    return log

def comp_for_length(in_fl, comparator, thr, out_fl=""):
    """ Takes an pre-alignment file, threshold, and comparator
        and compares all species with the comparator species
        and removes the species with sequences below the threshold length."""
    pre_algn_tb = algn_tb_gen(in_fl)
    comparator_index = -1
    comparator_length = 0
    rem_proteins = []
    if(out_fl != ""):
        out_file = open(out_fl, "w+")
    else:
        out_file = open(in_fl , "w")
    for index in range(len(pre_algn_tb)):
        if(comparator in pre_algn_tb[index][0]):
            comparator_index = index
            break
    if(comparator_index != -1):
        comparator_length = len(pre_algn_tb[index][1])
        for element in pre_algn_tb:
            if(len(element[1]) >= thr * comparator_length):
                out_file.write(">" + element[0] + "\n")
                out_file.write(element[1] + "\n")
            else:
                rem_proteins.append(element[0])
        return rem_proteins
    else:
        return "No DARP found"  

def SeqFuncDomain_Fast(alignPath="none", combinationPath="none", funcDomPath="none", outtxt=""):
    """
    OVERVIEW: this function analyses a folder of alignments
    USE:
        funcDomPath ==== > folder of folders containing functional domains for each species

        alignPath ======== > folder of alignment txt files
        combinationPath == > txt file of combinations (from CombinOfProID output)
    """
    outtxt = open(outtxt, "w")
    if alignPath != "none":
        alignP = open(alignPath, 'r').read() 
    if combinationPath != "none":
        combP = open(combinationPath, 'r').read() 
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
        alignP = open(alignPath, 'r') 
    if combinationPath != "none":
        combP = open(combinationPath, 'r').read() 
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
    USE:
        funcDomPath ==== > folder of folders containing functional domains for each species

        alignPath ======== > folder of alignment txt files
        combinationPath == > txt file of combinations (from CombinOfProID output)
    """
    outtxt = open(outtxt, "w")
    if alignPath != "none":
        alignP = open(alignPath, 'r').read() 
    if combinationPath != "none":
        combP = open(combinationPath, 'r').read() 
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
    USE:
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




def DomainHits_GeneProtein(alignments_dir, SPID):
    """
    OVERVIEW: Aquires the proteinID that associates with its GeneID

    USE:
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

    USE:
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
    OVERVIEW: 
        This function takes the file outputed by the core function and organizes the data
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

def list_to_fasta(L, seq_tb, out_dir):
    """ OVERVIEW: takes a list of protein ids and generates an unaligned sequence
        file ready for alignment in the out_dir
    INPUTS:
        fl_name = name the output file will be saved as
              L = list of protein ids
         seq_tb = sequence table from gen_seq_hash_tb()
        out_dir = path to location where output will be written """
    
    with open(out_dir, "w+") as out:
        c = 0    
        # iterate over all protein ids
        for pr_id in L:
            if (" " not in pr_id and "-" not in pr_id):
                seq = get_seq(pr_id, seq_tb)
                out.write(">{}\n".format(pr_id))
                out.write(str(seq) + "\n")
                c += 1
    return "File ready for alignment with {} species".format(c)

def bioMuscleAlign(inputF, musclePath, outputF="", quiet=0):
    """
    OVERVIEW: 
        returns aligned .fa files of the unaligned .fa files from "inputF" folder in the same
        folder using the MUSCLE algorithm for peptide alignment.
    USE:
        1. copy unaligned fasta (.fa) files to new folder
        2. use that folder path as inputF
        3. function will overwrite all .fa files with aligned seq's
    """

    if (outputF == "" or outputF == inputF):
        outputF = inputF
    else:
    # Copy files from source folder to output folder
        src_files = os.listdir(inputF)
        for file_name in src_files:
            full_file_name = os.path.join(inputF, file_name)
            if (os.path.isfile(full_file_name)):
                shutil.copy(full_file_name, outputF)
    
    # Run Alignment        
    #
    muscExe = open(musclePath)
    for file in os.listdir(inputF):
        path = outputF + file
        if quiet:
            command = "{} -in {} -out {} -quiet".format(musclePath, path, path)
        else:
            command = "{} -in {} -out {}".format(musclePath, path, path)
        exit_status = os.system(command)
        print(exit_status)
        if exit_status == -1:
            print( "ERROR: could not run command\n---> {}".format(command) )
            return -1
        else:
            print("finished aligning {}".format(file))

class myThread(threading.Thread):
    def __init__(self, inputF, musclePath, outputF=""):
        threading.Thread.__init__(self)
        self.inputF = inputF
        self.musclePath = musclePath
        self.outputF = outputF
    
    def run(self):
        biomuscleAlign(self.inputF, self.musclePath, self.outputF)


def checkForAlignment(folder):
    """
    OVERVIEW: 
        this function looks at alignment files in "folder" and checks for hyphens to see if any
        files were/are not aligned.
    USE: 
        bioMuscleAlign ==> ( checkForAlignment )
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

def gen_dataframe(filePath):
    excel_file = pandas.read_excel(filePath)  
    #df = pandas.DataFrame(excel_file)
    #row_len = excel_file.iloc[0][1]
    species_list = []
    sec_keys = []
    hash_tb = {}
    outer_hash_tb = {}
    for i in range(len(excel_file.columns)):
        sec_keys.append(excel_file.iloc[0][i])
    
    row = 0
    #print(sec_keys)
    for i in range(len(excel_file.index)):
        spec = str(excel_file.iloc[i][2])
        #print(spec)
        if(i != 0 and spec not in species_list):
            species_list.append(spec)
            col = 0
            #print(spec)
            hash_tb[spec] = {}
            for element in sec_keys:

                if(element != "Species"):
                    #print(element)
                    #print(excel_file.loc[row][col])
                    hash_tb[spec][element] = excel_file.loc[row][col]
                    #outer_hash_tb[spec] = inner_hash_tb 
                col += 1
        row += 1 
    return hash_tb

def parse_file(filePath):
    ''' Creates a hashtable based on data from the blast files.

    '''

    # for i in range(row_len):
    #     spec = excel_file.iloc[i:2]
    #     if(spec not in species_list):
    #         species_list.append(spec)
    
    blast_file = open(filePath, "r")       
    hash_tb = {}
    spec_id = ""
    res_spec_id = ""
    res_spec_score = ""
    res_spec_expect = ""
    res_spec_method = ""
    res_spec_ident = ""
    res_spec_pos = ""
    res_spec_gaps = ""
    res_pro_seq = ""
    

    #hash_tb[spec_id] = {}
    found = 0
    for line in blast_file:
        if("Query=" in line):
            #print(line)
            spec_id = line.split()[1]
            
            found += 1
        
        elif(">" in line):
            res_spec_id = line.split()[0]
            res_pro_seq = ""
            found += 1
        
        
        elif("Score =" in line):
            res_spec_score = line.split(",")[0].split("=")[1]
            res_spec_expect = line.split(",")[1].split("=")[1]
            res_spec_method = line.split(",")[2].split(":")[1]
            found += 1
            
        elif("Identities =" in line):
            res_spec_ident = line.split(",")[0].split("=")[1]
            res_spec_pos = line.split(",")[1].split("=")[1]
            res_spec_gaps = line.split(",")[2].split("=")[1]
            found += 1
        
        elif("Sbjct" in line):
            res_pro_seq += line.split()[2].replace("-", "")
        

        else:
            continue
            
        hash_tb[spec_id] = {}
        hash_tb[spec_id]["Result"] = res_spec_id
        hash_tb[spec_id]["Score"] = res_spec_score
        hash_tb[spec_id]["Expect"] = res_spec_expect
        hash_tb[spec_id]["Method"] = res_spec_method[:-1]
        hash_tb[spec_id]["Identities"] = res_spec_ident
        hash_tb[spec_id]["Positives"] = res_spec_pos
        hash_tb[spec_id]["Gaps"] = res_spec_gaps[:-1]
        hash_tb[spec_id]["Sequence"] = res_pro_seq
            
            # spec_id = ""
            # res_spec_id = ""
            # res_spec_score = ""
            # res_spec_expect = ""
            # res_spec_method = ""
            # res_spec_ident = ""
            # res_spec_pos = ""
            # res_spec_gaps = ""
            

    return hash_tb

def gen_pro_list(filePath):

    hash_tb = parse_file(filePath)
    pro_list = []

    for key, value in hash_tb.items():
        #spec_id = key
        #print(key)
        pro_id_and_seq = []
        for inner_key, inner_value in value.items():
            if("Result" in inner_key):
                
                pro_id_and_seq.append(inner_value)
            
            elif('Sequence' in inner_key):
                pro_id_and_seq.append(inner_value)
            
            pro_list.append(pro_id_and_seq)
    
    return pro_list

def gen_query_list(filePath):
    hash_tb = parse_file(filePath)
    query_list = []

    for key, value in hash_tb.items():

        if(key not in query_list):    
            query_list.append(key)
    
    return query_list
        

def gen_pro_files(filePath, out_dir):
    
    pro_list = gen_pro_list(filePath)
    #pro_seq_file = open(filePath, "r")
    output_file = open(out_dir, "w+")
    
    for pair in pro_list:
        for element in pair:
            output_file.write(element +"\n")

def blast_check(filePath_ZF_X, filePath_X_ZF):
    pro_id_and_seq = gen_pro_list(filePath_ZF_X)
    query_list = gen_query_list(filePath_X_ZF)
    non_hit_pro = []

    for pair in pro_id_and_seq:
        for element in pair:
            if(element[0] not in query_list):
                non_hit_pro.append(element[0])

    return non_hit_pro

    

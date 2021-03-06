# Author:      Robby Boney, robby.boney@wsu.edu
# Last Update: May 23, 2017
# Start Date:  Feb 20, 2016

__name__ = "Procomp Analysis Module"

import os
import codecs
import shutil
import time
from time import gmtime, strftime
import threading
import time
import fnmatch
from Bio.Align.Applications import MuscleCommandline
from io import StringIO as strIO
from Bio import AlignIO

# ////////////////////////////////////////////////////////////////
# Utility Functions //////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////
    
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
#def spid_tb_is_

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

# ////////////////////////////////////////////////////////////////
# Core Functions    //////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////

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


def comp_site_by_site(alignPath, spec_group_path, func_domain_list = ""):    
    
    def _most_common(lst):
        return max(set(lst), key=lst.count)
    
    spec_L = algn_tb_gen(alignPath)
    
    group_D = spid_tb_gen(spec_group_path)
    align_file = open(alignPath, "r")
    count = 0
    hit = 0
    hit_list = []
    regen_list = []
    non_regen_list = []
    return_list = []
    amino_acid_differences = []
    for element in spec_L:
        if(spid_tb_get_group(element[0], group_D) == "R"):
            regen_list.append(element[1])
        elif(spid_tb_get_group(element[0], group_D) == "NR"):
            non_regen_list.append(element[1])
        else:
            count += 1
    # print(len(regen_list))
    # print(len(non_regen_list))
    # print(len(spec_L[0][1]))
  
    
    max_difference  = 0
    min_difference = 10
    differences = []
    max_difference_index = []
    min_difference_index = []
    index = 0
    for i in range(len(spec_L[0][1])): 
        non_regen_amino_acids = []
        regen_amino_acids = []
        
        
        
        for spec in non_regen_list:
            if(spec[i : i + 1] not in non_regen_amino_acids):
                non_regen_amino_acids.append(spec[i : i + 1])
        
        for regen_spec in regen_list:
            if(regen_spec[i : i + 1] not in regen_amino_acids):
                regen_amino_acids.append(regen_spec[i : i + 1])
            
        
        difference = len(non_regen_amino_acids) - len(regen_amino_acids)
        if(abs(difference) > max_difference):
            max_difference = abs(difference)
            
                    
        elif(abs(difference) < min_difference):
            min_difference = abs(difference)
        
        amino_acid_differences.append(difference)

        
    for element in amino_acid_differences:
        if(abs(element) == max_difference):
            max_difference_index.append(index)
        elif(abs(element) == min_difference):
            min_difference_index.append(index)

        index += 1

        
        # else:
        #     for spec in non_regen_list:
        #         if(i in func_domain_list):
        #             non_regen_amino_acids.append(spec[i : i + 1])
        
        #     for regen_spec in regen_list:
        #         if(i in func_domain_list):
        #             regen_amino_acids.append(regen_spec[i : i + 1])

    #      if(non_regen_amino_acids.count(most_common(non_regen_amino_acids)) == len(non_regen_amino_acids)):
    #         hit += 1
    #         hit_list.append(i)
    #     print(non_regen_amino_acids)

    # return_list.append(hit)
    # return_list.append(count)
    # return_list.append(hit_list)
        

    #print(amino_acid_differences)
    return amino_acid_differences

def comp_region(alignPath, spec_group_path, func_domain_list=""):
    amino_acid_differences_L = comp_site_by_site(alignPath, spec_group_path, func_domain_list)
    #print(amino_acid_differences_L)
    region_L = []
    pos_region_L = []

    for element in amino_acid_differences_L:
       
        if(element > 0):
           pos_region_L.append(element)
           
        else:
            if(pos_region_L != [] and len(pos_region_L) > 1):
                region_L.append(pos_region_L)
                pos_region_L = []
            
    
    return region_L

def comp_sort_region(alignPath, spec_group_path, func_domain_list=""):
    amino_acid_regions = comp_region(alignPath, spec_group_path, func_domain_list)
    transcript_name = alignPath.split("/")[len(alignPath.split("/")) - 1].split(".")[0]
    
    
    len_longest_amino_acid_region = 0
    largest_regional_sum = 0
    longest_amino_acid_region_L = []
    largest_regional_sum_L = []
    

    for reg in amino_acid_regions:
        reg_sum = 0
        if(len(reg) > len_longest_amino_acid_region):
            len_longest_amino_acid_region = len(reg)

        for element in reg:
            reg_sum += element
        
        if(reg_sum > largest_regional_sum):
            largest_regional_sum = reg_sum
    
    # print(len_longest_amino_acid_region)
    # print(largest_regional_sum)

    for reg in amino_acid_regions:
        reg_sum = 0
        if(len(reg) == len_longest_amino_acid_region):
            longest_amino_acid_region_L.append(reg)
        
        for element in reg:
            reg_sum += element
        
        if(reg_sum == largest_regional_sum):
            largest_regional_sum_L.append(reg)

    return transcript_name ,longest_amino_acid_region_L, largest_regional_sum_L
    
            
            

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

def DOM_gen_tb(path):

    spec_dict = {}
    for file in os.listdir(path):
        file_path = path + file
        spec_file = open(file_path, "r")
        pro_id_dict = {}
        #dom_string = ""
        prev_pro_id = ""
        dom_list = []
        count = 0
        for line in spec_file:
            if(line[0:15] == prev_pro_id):
                dom_list.append(line[17:])

            elif(count == 0):
                dom_list.append(line[17:])
                prev_pro_id = line[0:15]
                count = 1
            
            else:
                pro_id_dict[prev_pro_id] = dom_list
                dom_list = []
                dom_list.append(line[17:])
                prev_pro_id = line[15:]
        
        spec_dict[prev_pro_id[0:6]] = pro_id_dict
    
    return spec_dict 
        
def DOM_get_ranges(table, pro_id):
    spec_id = pro_id[0:6]
    pro_id_list = table[spec_id][pro_id]
    pro_id_value_list = []
    
    for element in pro_id_list:
        pro_id_string_list = element.split(",")
        pro_id_value_list.append(pro_id_string_list[0])
    
    return pro_id_value_list
       
def DOM_get_dom_map(path, pro_id):
    ''' Returns the functional domain map '''
    table = DOM_gen_tb(path)
    num_list = DOM_get_ranges(table, pro_id)
    num_tuple = tuple(num_list)
    total = []
    for i in num_tuple:
        lower = int(i.split("-")[0])
        upper = int(i.split("0")[1])
        temp_L = [k for k in range(lower, upper + 1)]
        for num in temp_L:
            if num not in total:
                total.append(num)

    return total
        

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
        test = open(path, "r")
        test_rd = test.read()
        if "-" not in test_rd: 
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
        test.close()

def bioMuscleAlignList(inputL, musclePath, outputF, quiet=0):
    """
    OVERVIEW: 
        returns aligned .fa files of the unaligned .fa files from "inputL" List in the same
        folder using the MUSCLE algorithm for peptide alignment.
    USE:
        1. copy unaligned fasta (.fa) files to new folder
        2. use that folder path as inputF
        3. function will overwrite all .fa files with aligned seq's
    """

    # Copy files from source folder to output folder
    #src_files = os.listdir(inputF)
    for index in range(len(inputL)):
        #full_file_name = os.path.join(inputL, file_name)
        if (os.path.isfile(inputL[index])):
            shutil.copy(inputL[index], outputF)
    
    # Run Alignment        
    #
    muscExe = open(musclePath)
    for element in inputL:
        file_list = element.split("/")
        #print(file_list)
        path = outputF + file_list[2]
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
            print("finished aligning {}".format(file_list[2]))


class myThread(threading.Thread):
    def __init__(self, inputL, musclePath, outputF):
        threading.Thread.__init__(self)
        self.inputL = inputL
        self.musclePath = musclePath
        self.outputF = outputF
    @stats
    def run(self):
        bioMuscleAlignList(self.inputL, self.musclePath, self.outputF)
        

def thread_creator(inputF, n, outputF, musl_path):
    # dir_path is defined as a folder of pre-alignment files and n is defined as the number 
    # of threads the user wants to create
    file_num = len(fnmatch.filter(os.listdir(inputF), '*.fasta'))
    if(file_num == 0):
        return "No Files!"

    if(file_num < n):
        n = file_num
    # Add max number of threads later
    if(n > 20):
        n = 20
    threadL = []
    file_names = []
    files_per_thread =  file_num // n
    files_remainder = file_num % n
    obj_num = 0
    index = 0
    
    for file in os.listdir(inputF):
        file_names.append(file)
    
    for thread in range(n):
        inputL = []
        for count in range(files_per_thread):
            if(index != file_num):
                inputL.append(inputF + file_names[index])
                index += 1
                
                
                
        if(files_remainder != 0):
            inputL.append(file_names[index])
            index += 1
            files_remainder -= 1
        
        threadL.append(myThread(inputL, musl_path, outputF))
        
    return threadL

    


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
    
    blast_file = codecs.open(filePath, "r",encoding='utf-8', errors='ignore')
    #blast_file = open(filePath, "r")     
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
    for line in blast_file: #.splitlines()
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
            
    #print(hash_tb)
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
            
            if(pro_id_and_seq not in pro_list):
                pro_list.append(pro_id_and_seq)
    #print(pro_list)
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
    output_file = open(out_dir, "w")
    
    for pair in pro_list:
        for element in pair:
            if len(element) > 1:
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


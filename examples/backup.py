if __name__ == '__main__':

    #================= GLOBAL VARIABLES =======================================================
    #combfile = "/Volumes/HDD/-Apps/projects_Python/ProteinConvergence2016/Resources/SeqCombinations/1-2-2017_CombOutPut.txt"
    #unalignedfiles = "/Volumes/HDD/-Apps/projects_Python/ProteinConvergence2016/Resources/unalignedPeptides/1-3-17 Cleaned up sequences/"
    #outputFolder = "/Volumes/HDD/-Apps/projects_Python/ProteinConvergence2016/Resources/*Output*/Jan-4-2017_combinationUnaligned/"
    
    #fastaFolder = "/Users/themusicman/PycharmProjects/bioinformaticsJAN2017/FastaRemain/"
    fastaComplete = "/Users/themusicman/Projects/Python/bioinformaticsJAN2017/fastaRemainCompleted/"
    #muscleExePath = "/Users/themusicman/PycharmProjects/bioinformaticsJAN2017/muscle3.8.31_i86darwin64/muscle3.8.31_i86darwin64.31_i86darwin64"
    SpeID = "/Users/themusicman/Projects/Python/bioinformaticsJAN2017/SpeciesIdentifier_Jan6-2017.txt"
    output = "/Users/themusicman/Projects/Python/bioinformaticsJAN2017/output.txt"
    #conOut = "/Users/themusicman/PycharmProjects/bioinformaticsJAN2017/condensedOutput.txt"
    compHits = "/Users/themusicman/Projects/Python/bioinformaticsJAN2017/Results/Mar16-2017_mainPC.txt"
    #combOutput = "/Users/themusicman/PycharmProjects/bioinformaticsJAN2017/CombOutPut.txt"
    #condHits = "/Users/themusicman/PycharmProjects/bioinformaticsJAN2017/Jan6-2017_condensedHits.txt"
    #domainfolder = "/Users/themusicman/PycharmProjects/bioinformaticsJAN2017/Resources/1-6-17 Domains/"
    #domainOut = "/Users/themusicman/PycharmProjects/bioinformaticsJAN2017/Resources/DomainOutput.txt"
    #example = "/Users/themusicman/PycharmProjects/bioinformaticsJAN2017/EXAMPLE/"
    domains = "/Users/themusicman/Projects/PycharmProjects/bioinformaticsJAN2017/Resources/1-23-17_ZfishDomainsStartEnd/"
    proteinHits = "/Users/themusicman/Projects/PycharmProjects/bioinformaticsJAN2017/Results/Jan17-2017_completeHitsSorted.txt"
    domOut = "/Users/themusicman/Projects/PycharmProjects/bioinformaticsJAN2017/Results/Feb4-2017_domainHits_raw.txt"
    domOut_con = "/Users/themusicman/Projects/PycharmProjects/bioinformaticsJAN2017/Results/Feb4-2017_domainHits_con_2.txt"
    file = "/Users/themusicman/Projects/PycharmProjects/bioinformaticsJAN2017/Results/Jan17-2017_condensedHits.txt"

    #=================================== FUNCTIONS ===========================
    t1 = time.clock()

    #combGeneCompile(combfile, unalignedfiles, outputFolder)
    #bioMuscleAlign(fastaFolder, muscleExePath)
    #MainProteinCompare(fastaComplete, SpeID, outputtxt=output)
    MainPC_analysis(compHits, output)
    #a = SeqFuncDomain_Fast(combinationPath=combOutput, funcDomPath=domainfolder, outtxt=output)
    #SeqFuncDomain_FastAnalysis(a, SPID=SpeID)
    #SeqFuncDomain_analysis(domainOut,SpeID, outtxt=output)

    #DomainHits(domains, fastaComplete, proteinHits, domOut_con)
    #DomainHits_analysis(domOut_con)
    #DomainHits_GeneProtein(alignments_dir=fastaComplete, SPID="DARP")
    #RemoveSameGene(file)

    #===============END PROCESS====================
    t2 = time.clock()
    print("[Finished in", round(t2-t1, 4), "sec]")
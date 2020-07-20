def addGPR2Models(model,cyc):
    '''
    This function uses the function extractGeneAndProteinAssociation to add GPR
    associations to reactions with missing GPR in the model
    Input: 1) cobra model 2) pythonCyc PGDB instance
    Output: cobra model
    Author: Sanu Shameer (sanushameer@gmail.com)
    '''
    reactions = cyc.reactions
    rxnPresentList = list()
    rxnIDed = dict()
    for CycRxn in reactions.instances:
        CycRxn_id = CycRxn.frameid
        CycRxn_id_adapted = convertCycID2sbmlID(CycRxn_id)
        tempList = list()
        for rxn in model.reactions:
            if CycRxn_id_adapted == rxn.id[0:rxn.id.rindex("_")]:
                tempList.append(rxn)
            elif CycRxn_id_adapted == rxn.id[0:rxn.id.rindex("_")].replace("_NADP","").replace("_NAD",""):
                tempList.append(rxn)
        rxnIDed[CycRxn_id]=tempList
    SoyIgnoreList = ["RXN_9650_p","2_KETO_ADIPATE_DEHYDROG_RXN_m","Phytol_biosynthesis_p" \
                ,"CYSTEINE_AMINOTRANSFERASE_RXN_m","GLYCINE_TRNA_LIGASE_RXN_c" \
                ,"RXN66_1_c","RXN_9648_p","RXN-9651","Plastidial_ATP_Synthase_p" \
                ,"GGPP_biosynthesis_p","RXN_9653_p","lycopene_biosynthesis_p" \
                ,"RXN_2141_p","SUCCINYL_COA_HYDROLASE_RXN_m","PROTON_ATPase_c" \
                ,"MDA_Fd_Ascorbate_p","MercaptoPyruvateSulfurtransferase_m" \
                ,"Phytol_degradation_p","RXN_9652_p","A_B_oxidation_x","unlProtHYPO_c" \
                ,"Mitochondrial_ATP_Synthase_m","IPP_biosynthesis_c","Mehler_Reaction_p" \
                ,"Beta_Oxidation_x","HMBPP_synthesis_p","OROTATE_REDUCTASE_NADH_RXN_p" \
                ,"Ferredoxin_Plastoquinone_Reductase_p","RXN_9651_p","NADPH_Dehydrogenase_p" \
                ,"Plastoquinol_Oxidase_p","SUCCINATE_COA_LIGASE_GDP_FORMING_RXN_m","RXN_1781_v" \
                ,"PREPHENATE_DEHYDROGENASE_NADP_RXN_p","PREPHENATEDEHYDROG_RXN_p" \
                ,"MALEYLACETOACETATE_ISOMERASE_RXN_c","RXN_9654_p","LCYSDESULF_RXN_c","RXN_9958_NAD_m" \
                ,"HEXOKINASE_RXN_MANNOSE_c","PYRUVDEH_RXN_p","PYRUVDEH_RXN_m"] #last 3 lines present in latest version of SoyCyc
    print("--------------\nThis list of metabolic reactions are ignored")
    print(SoyIgnoreList)
    print("--------------")
    IDedlist = set()
    for rxnlist in rxnIDed.values():
        IDedlist = IDedlist.union(set(rxnlist))
    for rxn in set(model.reactions) - IDedlist:
        if not("_tx" in rxn.id or "_pc" in rxn.id or \
               "_mc" in rxn.id or "_xc" in rxn.id or \
               "_im" in rxn.id or "_vc" in rxn.id or \
               "_ec" in rxn.id or "_ep" in rxn.id or \
               "_pr" in rxn.id) \
               and (not "Biomass" in rxn.id) and \
               (not "biomass" in rxn.id) and \
               (not "Protein" in rxn.id) and \
               (not "TRNA_LIGASE" in rxn.id):
               if rxn.id not in SoyIgnoreList:
                   print(rxn.id)
    for k in rxnIDed.keys():
        for v in rxnIDed.get(k):
            rxn = v
            if rxn.gene_reaction_rule == "":
                #print k
                GPR = extractGeneAndProteinAssociation(cyc,k)
                if GPR != "()":
                    GPR = GPR.replace("() or ","")
                    rxn.gene_reaction_rule = GPR
    return model

def extractGeneAndProteinAssociation(cyc,frame_id):
    '''
    This functions adds Gene Associations to cobra model from Pathway Tools via
    PythonCyc
    Input: 1) pythonCyc PGDB instance 2) Frame id of reaction from Pathway Tools
    Output: Gene-Protein-Reaction associations from a PGDB
    Author: Sanu Shameer (sanushameer@gmail.com)
    '''
    rxn = getFrame(cyc,frame_id)
    if frame_id in cyc.reactions.instances:
        print("Error check if "+frame_id+" is reaction")
        return ""
    else:
        if "enzymatic_reaction" not in dir(rxn):
            return ""
        else:
            enzrxns = cyc.get_frame_objects(rxn.enzymatic_reaction)
            GPR = "(GPR)"
            temp1 = ""
            for enzrxn in enzrxns:
                enz = getFrame(cyc,enzrxn.enzyme)
                if "names" not in dir(enz):
                    continue
                if "gene" not in dir(enz):
                    continue
                if temp1 == "":
                    temp1 = str(enz.frameid)
                else:
                    temp1 = temp1 +" or "+str(enz.frameid)
                temp2 = ""
                for gene in enz.gene:
                    gene = getFrame(cyc,gene)
                    if "accession_1" not in dir(gene):
                        temp1.replace(enz.frameid,"")
                        continue
                    if temp2 == "":
                        temp2 = gene.accession_1
                    else:
                        temp2 = temp2 +" or "+gene.accession_1
                #print temp1
                #print temp2
                temp1 = temp1.replace(enz.frameid,"("+temp2+")")
            GPR = GPR.replace("GPR",temp1)
            return GPR

def getFrame(cyc,frame_id):
    '''
    This function retrieves pythoncyc frame from a PGDB instance
    Input: 1) pythonCyc PGDB instance 2) Frame id from Pathway Tools
    Output: Python instance of a frame
    Author: Sanu Shameer (sanushameer@gmail.com)
    '''
    frame = cyc.get_frame_objects([frame_id])[0]
    return frame

def convertCycID2sbmlID(id):
    '''
    This function converts Pathway Tools IDs to one that is SBML compliant
    Input: BioCyc IDs
    Output: SBML compliant IDs
    Author: Sanu Shameer (sanushameer@gmail.com)
    '''
    new_id = id.replace(".","_PERIOD_")
    new_id = new_id.replace("%2b","_")
    new_id = new_id.replace("&#039;","_")
    new_id = new_id.replace("&amp;beta;","B")
    new_id = new_id.replace("&beta;","B")
    new_id = new_id.replace("|","")
    new_id = new_id.replace("+-","_")
    new_id = new_id.replace("--","_")
    new_id = new_id.replace("-","_")
    new_id = new_id.replace("+","_")
    new_id = new_id.replace("'","_")
    new_id = new_id.replace("(","_")
    new_id = new_id.replace(")","_")
    new_id = new_id.replace("/","_")
    new_id = new_id.replace("__","_")
    return new_id

def find_average(temp_list):
    '''
    This function calculates the average from a list of numbers
    Input: list
    Output: float
    Author:Sanu Shameer (sanushameer@gmail.com)
    '''
    return sum(temp_list)/len(temp_list)

def adjustObjectiveBasedOnDaylength(diel_leaf,daylength,obj_rxn="diel_biomass",Phloem_day = "X_Phloem_contribution_t1",Phloem_night="X_Phloem_contribution_t2"):
    '''
    This function adjusts the day/night phloem contribution
    in a diel leaf objective.
    Input:  1) diel model, 2) length of day in hours, 3)ID of objective,
            4) Metabolite ID of day-time phloem, 5) Metabolite ID of night-time
            phloem
    Output: 1)diel model
    '''
    ratio = float(3*daylength)/(24-daylength)
    rxn = diel_leaf.reactions.get_by_id(obj_rxn)
    met1 = diel_leaf.metabolites.get_by_id(Phloem_day)
    coeff = abs(rxn.metabolites.get(met1))
    new_coeff = ratio - coeff
    rxn.add_metabolites({met1:-1*new_coeff})
    return diel_leaf

def generateMetaboliteFormula(rxn):
    count = 0
    for met in rxn.metabolites:
        if met.formula=="" or met.formula=="NA" or met.formula == None:
            if met.formula == "NA" or met.formula == None:
                met.formula = ""
            count = count + 1
    if count == 1:
        unb = rxn.check_mass_balance()
        #print(unb.keys())
        for met in rxn.metabolites:
            stoich = rxn.metabolites[met]
            if met.formula == "":
                tempForm = ""
                for a in ["C","H","O"]:
                    if a in unb.keys():
                        if round(unb[a]/stoich,6)==0:
                            continue
                        tempForm = tempForm+a+str(abs(unb[a])/stoich)
                        #print(a)
                        #print(unb[a])
                        #print(stoich)
                        #print(abs(unb[a])/stoich)
                for a in unb.keys():
                    if a in ["C","H","O"]:
                        continue
                    if a=="charge" or round(unb[a]/stoich,6)==0:
                        continue
                    tempForm = tempForm+a+str(abs(unb[a])/stoich)
                met.formula = tempForm
                print(met.id)
                print(tempForm)
    else:
        print("Unable to generate missing metabolite formula")

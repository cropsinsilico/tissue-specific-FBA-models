def extractGeneAndProteinAssociation(cyc,frame_id):
    '''
    This functions adds Gene Associations to cobra model from Pathway Tools via
    PythonCyc
    Input: 1) pythonCyc PGDB instance 2) Frame id of reaction from Pathway Tools
    Output: Gene-Protein-Reaction associations from a PGDB
    Author: Sanu Shameer (sanushameer@gmail.com)
    '''
    rxn = getFrame(cyc,frame_id)
    if rxn.reaction_direction == None:
        print("Error check if "+frame_id+" is reaction")
        return
    else:
        enzrxns = cyc.get_frame_objects(rxn.enzymatic_reaction)
        GPR = "(GPR)"
        temp1 = ""
        for enzrxn in enzrxns:
            enz = getFrame(cyc,enzrxn.enzyme)
            if temp1 == "":
                temp1 = str(enz.frameid)
            else:
                temp1 = temp1 +" or "+str(enz.frameid)
                temp2 = ""
                for gene in enz.gene:
                    gene = getFrame(cyc,gene)
                    if temp2 == "":
                        temp2 = gene.accession_1
                    else:
                        temp2 = temp2 +" or "+gene.accession_1
                temp1 = temp1.replace(enz.frameid,temp2)
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
    import regex
    new_id = id.replace(".","_PERIOD_")
    new_id = new_id.replace("+-","_")
    new_id = new_id.replace("--","_")
    new_id = new_id.replace("-","_")
    new_id = new_id.replace("+","_")
    return new_id

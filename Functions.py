def extractGeneAssociations(model,cyc_id):
    '''
    This functions adds Gene Associations to cobra model from Pathway Tools via
    PythonCyc
    '''
    return

def getFrame(cyc,frame_id):
    '''
    This function retrieves pythoncyc frame from a PGDB instance
    '''
    frame = cyc.get_frame_objects([frame_id])[0]
    return frame

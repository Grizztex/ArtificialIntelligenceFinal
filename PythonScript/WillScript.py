from pymol import cmd
import get_raw_distances
import pandas as pd
import numpy as np


def PolarInteraction(filePath):
    """
    This method extract the interacting residues from docking result and stored in a table with
    Only works with hetero..mer, cannot work with homo..mer"
    return: table
    """
    cmd.load(filePath)
    receptor, ligand = list(cmd.get_names())
    
    myspace = { "listChainRec" :[],"listRec" : [],"listRecResidue": [],"listChainLig":[],"listLigResidue" : [], "listLig": [], "listDist" : []}
    
    

    cmd.distance("polar_contacts",ligand,receptor,3.2,mode=2)
    for pair in get_raw_distances.get_raw_distances("polar_contacts"):
        (file1,res1,chain1),(file2,res2,chain2),dist = pair
        if file1 == receptor and file2 == ligand:
            
            cmd.iterate(f"index {res1} and chain {chain1} and {file1}","listRecResidue.append(resn);listRec.append(resv);listChainRec.append(chain)",space = myspace)

            cmd.iterate(f"index {res2} and chain {chain2} and {file2}","listLigResidue.append(resn);listLig.append(resv);listChainLig.append(chain)", space = myspace)
        elif file1 == receptor and file2 == ligand:

            cmd.iterate(f"index {res2} and chain {chain2} and {file1}","listRecResidue.append(resn);listRec.append(resv);listChainRec.append(chain)",space = myspace)

            cmd.iterate(f"index {res1} and chain {chain1} and {file2}","listLigResidue.append(resn);listLig.append(resv);listChainLig.append(chain)",space = myspace)
        else:
            continue
        myspace["listDist"].append(dist)
    table = pd.DataFrame(myspace)
    cmd.delete("all")   
    return table

def extractFromPoses(dirName):
    """
    This method extract polar interactions from docking poses
    Input: Directory name of containings poses of a sequence result
    Output: Table
    """
    
    seqID = re.findall(pattern = ("\d+"), string = dirName.split("\\")[-1])[0]
    table =  pd.DataFrame({"Sequence ID":[],"Model":[],"listChainRec" :[],"listRec" : [],"listRecResidue": [],"listChainLig":[],"listLigResidue" : [], "listLig": [], "listDist" : []})

    ##Create table with seqID, model, and other listed features
    for pose in os.listdir(dirName):
        poseNumb = pose.split(".")[2]
        poseDir = dirName + "\\" + pose
        table_extract = script.PolarInteraction(poseDir)
        table_extract.insert(0,"Sequence ID",seqID)
        table_extract.insert(1,"Model",poseNumb)
        table = pd.concat([table_extract,table])
    
    
    
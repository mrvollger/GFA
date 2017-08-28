#!/usr/bin/env python
import subprocess
import sys
import re
import itertools 
import os
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

class GFA(object):
    
    def __init__(self, gfa_file):
        self.file  = gfa_file
        self.Graph = nx.DiGraph()
        self.readIn()
        self.drawGFA()

    def getTigID(self, tigName ):
        # remove tig and weird zeros 
        rtn = re.sub("tig0+", "", tigName)
        # remove the strand orrientation
        rtn = re.sub("_end|_start", "", rtn)
        return(rtn)
    
    
    # can add features to a node or edge despite the name, as long as it is part of the segment feature of the graph 
    def addFeatures(self, node):
        feats = {} 
        # add that all parts of this segment will have
        feats["weight"] = 5
        feats["color"] = "black" #node["tigID"]

        featFile = re.sub(".gfa$", "", self.file) + ".layout.tigInfo"
        # break if there is not feature file 
        if(not os.path.isfile(featFile) ):
            print("No file called: " + featFile)
            return(feats) 
        
        f = open(featFile)
        header = []
        values = []
        # read feature file and get features
        for line in f:
            if( line[0] == "#"):
                header = line[1:].split("\t")
            else:
                sline = line.split("\t")
                if(sline[0] == node["tigID"] ):
                    values = sline
                    break
        
        # add the features to the node
        feats.update(  dict(zip(header, values))  )
        for key, value in feats.iteritems():
            node[key] = value
        return(feats)

    # creates two connected nodes to represent a segment 
    def addSeg(self, sline):
        # create node 
        G = self.Graph 
        segment = sline[1]
        segStart = segment + "_start"
        segEnd = segment + "_end"
        G.add_node(segStart)
        G.add_node(segEnd)
        
        node1 = G.node[segStart]
        node2 = G.node[segEnd]
        # get the tig ID and make it a feature 
        tigID = self.getTigID(segStart)
        node1["tigID"] = tigID 
        node2["tigID"] = tigID 
       
        # make an enhanced tig id that includes the orrientation
        node1["eTigID"] = tigID + "+" 
        node2["eTigID"] = tigID + "-"

        # add additional features contigent on a certain file existing
        self.addFeatures(node1)
        feats = self.addFeatures(node2)
        
        # add the edge, with the features the nodes have 
        self.Graph.add_edge(segStart, segEnd, attr_dict=feats)
        
        return 
        
    def addEdge(self, sline):
        orginNode = sline[1]
        destNode = sline[3]
        
        feats = {}
        feats["orginStrand"] = sline[2]
        feats["destStrand"] = sline[4]
        feats["cigar"] = sline[5]
        feats["orginID"] = self.getTigID(orginNode) 
        feats["destID"] = self.getTigID(destNode) 
        feats["weight"] = 2
        feats["color"] = "green"
        # goes from the end of the minus strand so, which is the beginning fo the forward strand   
        if( feats["orginStrand"] == "-"):
            orginNode += "_start"
        # goes from the ned of the forward strand, so this starts at the end
        if( feats["orginStrand"] == "+" ):
            orginNode += "_end"
        
        # arrives at the beggining of the minus strand, which is the end of the forward strand   
        if( feats["destStrand"] == "-"):
            destNode += "_end"
        # arrives at the beggining of the forward strand which is the start of the segmetn
        if( feats["destStrand"] == "+" ):
            destNode += "_start"

        self.Graph.add_edge(orginNode, destNode,  attr_dict=feats)
        return

    def readIn(self):
        # possible line starts: H, S, F, E, G, O, or U
        # header, segment, fragment, edge, gap, O/U => group/subgrpah i
        # but I have GFA 1 files so L and C line have been consolidated to E lines
        # P line has been replaced with U and O lines
        # GFA 1 looks like this #   Comment; H   Header; S   Segment;  L   Link; C   Containment; P   Path
        gfafile = open(self.file)
        
        for line in gfafile:
            sline = line.split("\t")
            ID = sline[0]
            
            if(ID == "H" ):
                continue
            elif(ID == "S"):
                self.addSeg(sline)
            elif(ID in ["E","L"]):
                self.addEdge(sline)
        G = self.Graph
        print("Nodes:{} Edges:{}".format( G.number_of_nodes(), G.number_of_edges() ) )
        return
    
    def drawGFA(self):
        G = self.Graph
        self.plot = plt.figure()
        name = self.file
        
        # generate positions of the nodes
        pos = nx.spring_layout(G,  iterations=1500)
        # add the enhanced tig IDs as labels 
        labels=dict((n,d['eTigID']) for n,d in G.nodes(data=True)) 
        # add edge features as weights 
        edges = G.edges()
        weights = [G[u][v]["weight"] for u,v in edges]
        edge_color = [G[u][v]["color"] for u,v in edges]
        
        nx.draw(G, pos,labels=labels, edges=edges, width=weights, edge_color = edge_color)
        
        # add edge labels for the segments 
        edge_labels = nx.get_edge_attributes(G,"coverage")
        edge_labels2 = nx.get_edge_attributes(G,"tigLen")
        for key in edge_labels:
            edge_labels[key] = str(edge_labels[key])+"-"+str(edge_labels2[key])
        nx.draw_networkx_edge_labels(G, pos, edge_labels = edge_labels)

        plt.savefig(name + ".png")
        nx.write_gml(G, name+".gml")






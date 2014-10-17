#!/usr/bin/env python
from nose.tools import *
import networkx
"""from  betweenness  import * """

def usage():
    print " --nodes nodefile --edges edgefile --centrality_min centrality [ --Num_central_nodes size ]"
    sys.exit()

def florentine_families_graph():
    G=networkx.Graph()
    G.add_edge('Acciaiuoli','Medici')
    G.add_edge('Castellani','Peruzzi')
    G.add_edge('Castellani','Strozzi')
    G.add_edge('Castellani','Barbadori')
    G.add_edge('Medici','Barbadori')
    G.add_edge('Medici','Ridolfi')
    G.add_edge('Medici','Tornabuoni')
    G.add_edge('Medici','Albizzi')
    G.add_edge('Medici','Salviati')
    G.add_edge('Salviati','Pazzi')
    G.add_edge('Peruzzi','Strozzi')
    G.add_edge('Peruzzi','Bischeri')
    G.add_edge('Strozzi','Ridolfi')
    G.add_edge('Strozzi','Bischeri')
    G.add_edge('Ridolfi','Tornabuoni')
    G.add_edge('Tornabuoni','Guadagni')
    G.add_edge('Albizzi','Ginori')
    G.add_edge('Albizzi','Guadagni')
    G.add_edge('Bischeri','Guadagni')
    G.add_edge('Guadagni','Lamberteschi')
    return G

def weighted_G():
    G=networkx.Graph();
    G.add_edge(0,1,weight=3)
    G.add_edge(0,2,weight=2)
    G.add_edge(0,3,weight=6)
    G.add_edge(0,4,weight=4)
    G.add_edge(1,3,weight=5)
    G.add_edge(1,5,weight=5)
    G.add_edge(2,4,weight=1)
    G.add_edge(3,4,weight=2)
    G.add_edge(3,5,weight=1)
    G.add_edge(4,5,weight=4)

    return G


class TestBetweennessCentrality():
        
    def __init__(self):
        """Betweenness centrality: K5"""
        print """Betweenness centrality: K5"""
        G=networkx.complete_graph(5)
        """       b=networkx.betweenness_centrality(G, weight=False, normalized=False)  """
        b=networkx.betweenness_centrality(G, False, False)
        b_answer={0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}
        for n in sorted(G):
            print n
            assert_almost_equal(b[n],b_answer[n])

    
    def test_K5(self):
        """Betweenness centrality: K5"""
        print """Betweenness centrality: K5"""
        G=networkx.complete_graph(5)
        b=networkx.betweenness_centrality(G, False, False)
        b_answer={0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}
        for n in sorted(G):
            assert_almost_equal(b[n],b_answer[n])

    def test_K5_endpoints(self):
        """Betweenness centrality: K5 endpoints"""
        G=networkx.complete_graph(5)
        b=networkx.betweenness_centrality(G, False, False, endpoints=True)
        b_answer={0: 4.0, 1: 4.0, 2: 4.0, 3: 4.0, 4: 4.0}
        for n in sorted(G):
            assert_almost_equal(b[n],b_answer[n])


    def test_P3_normalized(self):
        """Betweenness centrality: P3 normalized"""
        G=networkx.path_graph(3)
        b=networkx.betweenness_centrality(G, False,True)
        b_answer={0: 0.0, 1: 1.0, 2: 0.0}
        for n in sorted(G):
            assert_almost_equal(b[n],b_answer[n])


    def test_P3(self):
        """Betweenness centrality: P3"""
        G=networkx.path_graph(3)
        b_answer={0: 0.0, 1: 1.0, 2: 0.0}
        b=networkx.betweenness_centrality(G, weight=False, normalized=False)
        for n in sorted(G):
            assert_almost_equal(b[n],b_answer[n])

    def test_P3_endpoints(self):
        """Betweenness centrality: P3 endpoints"""
        G=networkx.path_graph(3)
        b_answer={0: 2.0, 1: 3.0, 2: 2.0}
        b=networkx.betweenness_centrality(G, False, False, True)
        for n in sorted(G):
            print n , b[n]
            assert_almost_equal(b[n],b_answer[n])


    def test_krackhardt_kite_graph(self):
        """Betweenness centrality: Krackhardt kite graph"""
        G=networkx.krackhardt_kite_graph()
        b_answer={0: 1.667,1: 1.667,2: 0.000,3: 7.333,4: 0.000,
                  5: 16.667,6: 16.667,7: 28.000,8: 16.000,9: 0.000}
        for b in b_answer:
            b_answer[b]/=2.0
        b=networkx.betweenness_centrality(G, weight=False, normalized=False)

        for n in sorted(G):
            assert_almost_equal(b[n],b_answer[n],places=3)


    def test_krackhardt_kite_graph_normalized(self):
        """Betweenness centrality: Krackhardt kite graph normalized"""
        G=networkx.krackhardt_kite_graph()
        b_answer={0:0.023,1:0.023,2:0.000,3:0.102,4:0.000, 5:0.231,6:0.231,7:0.389,8:0.222,9:0.000}
        b=networkx.betweenness_centrality(G, weight=False, normalized=True)

        for n in sorted(G):
            assert_almost_equal(b[n],b_answer[n],places=3)


    def test_florentine_families_graph(self):
        """Betweenness centrality: Florentine families graph"""
        G=florentine_families_graph()
        b_answer=\
             {'Acciaiuoli':    0.000,
              'Albizzi':       0.212,
              'Barbadori':     0.093,
              'Bischeri':      0.104,
              'Castellani':    0.055,
              'Ginori':        0.000,
              'Guadagni':      0.255,
              'Lamberteschi':  0.000,
              'Medici':        0.522,
              'Pazzi':         0.000,
              'Peruzzi':       0.022,
              'Ridolfi':       0.114,
              'Salviati':      0.143,
              'Strozzi':       0.103,
              'Tornabuoni':    0.092}

        b=networkx.betweenness_centrality(G, weight=False, normalized=True)
        for n in sorted(G):
            assert_almost_equal(b[n],b_answer[n],places=3)


    def test_ladder_graph(self):
        """Betweenness centrality: Ladder graph"""
        G = networkx.Graph() # ladder_graph(3)
        G.add_edges_from([(0,1), (0,2), (1,3), (2,3), 
                          (2,4), (4,5), (3,5)])
        b_answer={0:1.667,1: 1.667,2: 6.667,
                  3: 6.667,4: 1.667,5: 1.667}
        for b in b_answer:
            b_answer[b]/=2.0
        b=networkx.betweenness_centrality(G, weight=False, normalized=False)
        for n in sorted(G):
            assert_almost_equal(b[n],b_answer[n],places=3)

    def test_disconnected_path(self):
        """Betweenness centrality: disconnected path"""
        G=networkx.Graph()
        G.add_path([0,1,2])
        G.add_path([3,4,5,6])
        b_answer={0:0,1:1,2:0,3:0,4:2,5:2,6:0}
        b=networkx.betweenness_centrality(G, weight=False, normalized=False)
        for n in sorted(G):
            assert_almost_equal(b[n],b_answer[n])

    def test_disconnected_path_endpoints(self):
        """Betweenness centrality: disconnected path endpoints"""
        G=networkx.Graph()
        G.add_path([0,1,2])
        G.add_path([3,4,5,6])
        b_answer={0:2,1:3,2:2,3:3,4:5,5:5,6:3}
        b=networkx.betweenness_centrality(G, weight=False, normalized=False, endpoints=True)
        for n in sorted(G):
            assert_almost_equal(b[n],b_answer[n])


    def test_directed_path(self):
        """Betweenness centrality: directed path"""
        G=networkx.DiGraph()
        G.add_path([0,1,2])
        b=networkx.betweenness_centrality(G, weight=False, normalized=False)
        b_answer={0: 0.0, 1: 1.0, 2: 0.0}
        for n in sorted(G):
            assert_almost_equal(b[n],b_answer[n])

    def test_directed_path_normalized(self):
        """Betweenness centrality: directed path normalized"""
        G=networkx.DiGraph()
        G.add_path([0,1,2])
        b=networkx.betweenness_centrality(G, weight=False, normalized=True)
        b_answer={0: 0.0, 1: 0.5, 2: 0.0}
        for n in sorted(G):
            assert_almost_equal(b[n],b_answer[n])



class TestWeightedBetweennessCentrality():
        
    def test_K5(self):
        """Weighted betweenness centrality: K5"""
        G=networkx.complete_graph(5)
        b=networkx.betweenness_centrality(G, weight=True, normalized=False)
        b_answer={0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}
        for n in sorted(G):
            assert_almost_equal(b[n],b_answer[n])

    def test_P3_normalized(self):
        """Weighted betweenness centrality: P3 normalized"""
        G=networkx.path_graph(3)
        b=networkx.betweenness_centrality(G, weight=True, normalized=True)
        b_answer={0: 0.0, 1: 1.0, 2: 0.0}
        for n in sorted(G):
            assert_almost_equal(b[n],b_answer[n])


    def test_P3(self):
        """Weighted betweenness centrality: P3"""
        G=networkx.path_graph(3)
        b_answer={0: 0.0, 1: 1.0, 2: 0.0}
        b=networkx.betweenness_centrality(G, weight=True, normalized=False)
        for n in sorted(G):
            assert_almost_equal(b[n],b_answer[n])

    def test_krackhardt_kite_graph(self):
        """Weighted betweenness centrality: Krackhardt kite graph"""
        G=networkx.krackhardt_kite_graph()
        b_answer={0: 1.667,1: 1.667,2: 0.000,3: 7.333,4: 0.000,
                  5: 16.667,6: 16.667,7: 28.000,8: 16.000,9: 0.000}
        for b in b_answer:
            b_answer[b]/=2.0

        b=networkx.betweenness_centrality(G, weight=True, normalized=False)

        for n in sorted(G):
            print n
            assert_almost_equal(b[n],b_answer[n],places=3)


    def test_krackhardt_kite_graph_normalized(self):
        """Weighted betweenness centrality: 
        Krackhardt kite graph normalized
        """
        G=networkx.krackhardt_kite_graph()
        b_answer={0:0.023,1:0.023,2:0.000,3:0.102,4:0.000,
                  5:0.231,6:0.231,7:0.389,8:0.222,9:0.000}
        b=networkx.betweenness_centrality(G, weight=True, normalized=True)

        for n in sorted(G):
            assert_almost_equal(b[n],b_answer[n],places=3)


    def test_florentine_families_graph(self):
        """Weighted betweenness centrality: 
        Florentine families graph"""
        G=florentine_families_graph()
        b_answer=\
             {'Acciaiuoli':    0.000,
              'Albizzi':       0.212,
              'Barbadori':     0.093,
              'Bischeri':      0.104,
              'Castellani':    0.055,
              'Ginori':        0.000,
              'Guadagni':      0.255,
              'Lamberteschi':  0.000,
              'Medici':        0.522,
              'Pazzi':         0.000,
              'Peruzzi':       0.022,
              'Ridolfi':       0.114,
              'Salviati':      0.143,
              'Strozzi':       0.103,
              'Tornabuoni':    0.092}

        b=networkx.betweenness_centrality(G, weight=True, normalized=True)
        for n in sorted(G):
            assert_almost_equal(b[n],b_answer[n],places=3)


    def test_ladder_graph(self):
        """Weighted betweenness centrality: Ladder graph"""
        G = networkx.Graph() # ladder_graph(3)
        G.add_edges_from([(0,1), (0,2), (1,3), (2,3), 
                          (2,4), (4,5), (3,5)])
        b_answer={0:1.667,1: 1.667,2: 6.667,
                  3: 6.667,4: 1.667,5: 1.667}
        for b in b_answer:
            b_answer[b]/=2.0
        b=networkx.betweenness_centrality(G, weight=True, normalized=False)
        for n in sorted(G):
            assert_almost_equal(b[n],b_answer[n],places=3)

    def test_G(self):
        """Weighted betweenness centrality: G"""                   
        print """Weighted betweenness centrality: G"""                   
        G = weighted_G()               
        b_answer={0: 2.0, 1: 0.1, 2: 4.3, 3: 3.2, 4: 4.0, 5: 0.0}
        for a in  G.edges():
            G[a[0]][a[1]]['weight'] = 0.3
            print a[0], a[1], G[a[0]][a[1]]['weight'] 
        b=networkx.betweenness_centrality(G, weight=True, normalized=False)
        for n in sorted(G):
            print n, b[n]
        """    assert_almost_equal(b[n],b_answer[n])  """
        print G.nodes()

       
    def test_G2(self):
        """Weighted betweenness centrality: G2"""                   
        G=networkx.DiGraph()
        G.add_weighted_edges_from([('s','u',10) ,('s','x',5) ,
                                   ('u','v',1) ,('u','x',2) ,
                                   ('v','y',1) ,('x','u',3) ,
                                   ('x','v',5) ,('x','y',2) ,
                                   ('y','s',7) ,('y','v',6)])

        b_answer={'y':5.0,'x':5.0,'s':4.0,'u':2.0,'v':2.0}

        b=networkx.betweenness_centrality(G, weight=True, normalized=False)
        for n in sorted(G):
            assert_almost_equal(b[n],b_answer[n])


class TestEdgeBetweennessCentrality():
        
    def test_K5(self):
        """Edge betweenness centrality: K5"""
        G=networkx.complete_graph(5)
        b=networkx.edge_betweenness_centrality(G, weight=False, normalized=False)
        b_answer=dict.fromkeys(G.edges(),1)
        for n in sorted(G.edges()):
            assert_almost_equal(b[n],b_answer[n])

    def test_C4(self):
        """Edge betweenness centrality: C4"""
        G=networkx.cycle_graph(4)
        b=networkx.edge_betweenness_centrality(G, weight=False, normalized=False)
        b_answer={(0, 1):2,(0, 3):2, (1, 2):2, (2, 3): 2}
        for n in sorted(G.edges()):
            assert_almost_equal(b[n],b_answer[n])

    def test_P4(self):
        """Edge betweenness centrality: P4"""
        G=networkx.path_graph(4)
        b=networkx.edge_betweenness_centrality(G, weight=False, normalized=False)
        b_answer={(0, 1):3,(1, 2):4, (2, 3):3}
        for n in sorted(G.edges()):
            assert_almost_equal(b[n],b_answer[n])

    def test_balanced_tree(self):
        """Edge betweenness centrality: balanced tree"""
        G=networkx.balanced_tree(r=2,h=2)
        b=networkx.edge_betweenness_centrality(G, weight=False, normalized=False)
        b_answer={(0, 1):12,(0, 2):12,
                  (1, 3):6,(1, 4):6,(2, 5):6,(2,6):6}
        for n in sorted(G.edges()):
            assert_almost_equal(b[n],b_answer[n])

class TestWeightedEdgeBetweennessCentrality():
        
    def test_K5(self):
        """Edge betweenness centrality: K5"""
        G=networkx.complete_graph(5)
        b=networkx.edge_betweenness_centrality(G, weight=True, normalized=False)
        b_answer=dict.fromkeys(G.edges(),1)
        for n in sorted(G.edges()):
            assert_almost_equal(b[n],b_answer[n])

    def test_C4(self):
        """Edge betweenness centrality: C4"""
        G=networkx.cycle_graph(4)
        b=networkx.edge_betweenness_centrality(G, weight=True, normalized=False)
        b_answer={(0, 1):2,(0, 3):2, (1, 2):2, (2, 3): 2}
        for n in sorted(G.edges()):
            assert_almost_equal(b[n],b_answer[n])

    def test_P4(self):
        """Edge betweenness centrality: P4"""
        G=networkx.path_graph(4)
        b=networkx.edge_betweenness_centrality(G, weight=True, normalized=False)
        b_answer={(0, 1):3,(1, 2):4, (2, 3):3}
        for n in sorted(G.edges()):
            assert_almost_equal(b[n],b_answer[n])


    def test_balanced_tree(self):
        """Edge betweenness centrality: balanced tree"""
        G=networkx.balanced_tree(r=2,h=2)
        b=networkx.edge_betweenness_centrality(G, weight=True, normalized=False)
        b_answer={(0, 1):12,(0, 2):12,
                  (1, 3):6,(1, 4):6,(2, 5):6,(2,6):6}
        for n in sorted(G.edges()):
            assert_almost_equal(b[n],b_answer[n])

    def test_weighted_graph(self):
        eList = [(0, 1, 5), (0, 2, 4), (0, 3, 3), 
                 (0, 4, 2), (1, 2, 4), (1, 3, 1), 
                 (1, 4, 3), (2, 4, 5), (3, 4, 4)]
        G = networkx.Graph()
        G.add_weighted_edges_from(eList)
        b = networkx.edge_betweenness_centrality(G, normalized=False, weight=True)
        b_answer={(0, 1):0.0,
                  (0, 2):1.0,
                  (0, 3):2.0,
                  (0, 4):1.0,
                  (1, 2):2.0,
                  (1, 3):3.5,
                  (1, 4):1.5,
                  (2, 4):1.0,
                  (3, 4):0.5}

        for n in sorted(G.edges()):
            print n
            assert_almost_equal(b[n],b_answer[n])
def  main(argvs):
     try:
         opts, args = getopt.getopt(sys.argv[1:], "n:e:c:N:",['nodes=', 'edges=', 'centrality_min=', 'Num_central_nodes='])
     except  getopt.GetoptError, err:
         usage()
         sys.exit(2)

     if len(opts) < 3:
        usage()
        sys.exit(1)


     node_file= None
     edge_file = None
     min_attachement = 0.33
     central_num = 1000
     for opt, val in opts:
         if opt == "--nodes" or opt == "-n":
            node_file = str(val)
         if opt == "--edges" or opt == "-e":
            edge_file = str(val)
         if opt == "-c" or opt =="--centrality_min":
            min_centrality = float(val)
         if opt == "-N" or opt == "--Num_central_nodes":
            central_num = float(val)
            
            
     print node_file, edge_file

     nodefile = open(node_file, 'r')
    # nodefile.readlines()
     node_id_name_map={}
     for node in nodefile.readlines():
        node=node.strip()
       
        fields = node.split("\t")
        node_id_name_map[ str(fields[0]) ] = fields[1]
        
        
       # for field in fields
       #    print field
        
   #  for key in node_id_name_map.keys():
   #      print key, node_id_name_map[key]
     nodefile.close()
     
     edge_id_map={}
     edgefile = open(edge_file, 'r')
     #allEdges= []
     for edge in edgefile.readlines():
        edge=edge.strip()
        if re.search("CORREL", edge):
            continue 

        fields = edge.split("\t")
        fields[0].strip()
        fields[1].strip()
        fields[2].strip()
     #   allEdges.append( (fields[0], fields[1], fields[2], fields[4]))
      
     #print allEdges

        edge_id_map[(str(fields[0]), str(fields[1]))]=float(fields[2])




     edgefile.close()

     G=networkx.Graph();
     for key, value in edge_id_map.items():
         #print key[0], key[1], value 
         G.add_edge( str(key[0]),str(key[1]),weight=value)

     B=networkx.betweenness_centrality(G, weight=True, normalized=True)
     sorted_nodes_by_centrality = [ x for x in B.iteritems() ]
     sorted_nodes_by_centrality.sort(key=lambda x : x[1])
     sorted_nodes_by_centrality.reverse()

#     for n in sorted(G.nodes()):
     selected_nodes = []
     count = 0
     for key, val in sorted_nodes_by_centrality:
         print key, val 
         if val > min_centrality and  count < central_num :
       #     print key, val 
            selected_nodes.append(key)
            count =  count +1
         else:
            break
         
         #print n + "\t"+  node_id_name_map[n]+"\t"+ str(B[n])

     print "Num of central nodes selected ",len(selected_nodes)
     central_edgefile = open( 'central_%1.3f_' % min_centrality + edge_file,'w')
     central_edgefile.write("FROM\tTO\tCORREL\n")
     used_nodes={}
     num_edges=0
     print "Writing the edges "
     for e in sorted(G.edges()):
         if e[0] in selected_nodes or e[1] in selected_nodes:
            used_nodes[e[0]] = 1
            used_nodes[e[1]] = 1
            num_edges = num_edges + 1
            if (e[1],e[0])  in edge_id_map.keys(): 
                central_edgefile.write(e[0]+"\t"+e[1]+"\t"+str(edge_id_map[(e[1], e[0])])+"\n") 
            else:
                central_edgefile.write(e[0]+"\t"+e[1]+"\t"+str(edge_id_map[(e[0], e[1])])+"\n") 

     central_edgefile.close()
     s= '%.2f' % min_centrality
         
     print "Writing the nodes "
     nodefile = open(node_file, 'r')
     central_nodefile = open( 'central_%1.3f_' %min_centrality + node_file,'w')
     for node in nodefile.readlines():
        node = node.rstrip('\r\n')
        fields = node.split("\t")
        fields[0].strip()
       
        if fields[0]=='ID':
             central_nodefile.write(node+"\t"+"CENTRAL"+"\n")
        if fields[0] in used_nodes.keys():
             if fields[0] in selected_nodes:
                central_nodefile.write(node+"\t"+"Y"+"\n")
             else:
                central_nodefile.write(node+"\t"+"N"+"\n")

     nodefile.close()
     central_nodefile.close()

     print "Num central nodes : %d" % len(selected_nodes)
     print "Num nodes : %d" % len(used_nodes)
     print "Num edges : %d" % num_edges





import sys
import getopt
import re


if __name__ == '__main__':
 
     main(sys.argv)
#     A=TestBetweennessCentrality()
#     print """Testing....."""
#     A.test_K5()
#     A.test_P3_endpoints()
#     B=TestWeightedBetweennessCentrality()
#     B.test_G()

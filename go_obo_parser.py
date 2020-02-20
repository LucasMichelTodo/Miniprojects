import obonet
import networkx

obo_file = '/media/lucas/Disc4T/Projects/Eli/FE_Ontologizer/go-basic.obo'
graph = obonet.read_obo(obo_file)

# N.Nodes
len(graph)

# N.Edges
graph.number_of_edges()

# Is a DAG?
networkx.is_directed_acyclic_graph(graph)

# Creating ids dict
id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}

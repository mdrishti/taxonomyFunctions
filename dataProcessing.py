import zenodo_get
import pandas as pd
import re
import gzip
import argparse
from igraph import Graph, plot
import json
import math
import sys
from collections import Counter
from rapidfuzz import fuzz


#extract taxonomic ids according to categories file
def extr(catg, fx2):
    res2 = None  # Initialize the result variable
    for i in range(len(catg)):  # Loop over rows in `catg`
        #print(catg.iloc[i, :])  # Print the current category
        # Perform regex search on the first column of `fx2`
        pattern = rf"{catg.iloc[i, 0]}[A-Z0-9]+"
        matches = [re.findall(pattern, str(val)) for val in fx2.iloc[:, 0]]
        # Extract unique matches
        #cols = pd.unique([item for sublist in matches for item in sublist]) #gives future warning
        cols = list(set(item for sublist in matches for item in sublist))
        res1 = pd.DataFrame({catg.iloc[i, 0]: [match[0] if match else None for match in matches]})
        # Remove the category prefix from the matched strings
        res1[catg.iloc[i, 0]] = res1[catg.iloc[i, 0]].str.replace(catg.iloc[i, 0], "", regex=False)
        # Combine the result with previous results
        res2 = pd.concat([res2, res1], axis=1) if res2 is not None else res1
    return res2



def extr_igraph(catg, fx2):
    g = ig.Graph(directed=True)  # Create a directed graph
    node_map = {}  # Dictionary to track added nodes
    
    for idx in range(len(fx2)):  
        main_node = str(fx2.iloc[idx, 1])  # The main node for this row

        # Ensure the main node is in the graph
        if main_node not in node_map:
            node_map[main_node] = g.add_vertex(name=main_node)

        # Process `fx2.iloc[idx, 2]`, extracting taxon IDs by prefix
        taxon_data = str(fx2.iloc[idx, 2])
        
        for i in range(len(catg)):  
            prefix = catg.iloc[i, 0]  # Example: "gbif:", "ncbi:"
            pattern = rf"{prefix}[A-Z0-9]+"

            matches = re.findall(pattern, taxon_data)  # Extract IDs with this prefix
            if matches:
                # Ensure prefix node exists
                if prefix not in node_map:
                    node_map[prefix] = g.add_vertex(name=prefix)

                # Link main node → prefix node
                g.add_edge(node_map[main_node], node_map[prefix])
            
            for match in matches:
                taxon_id = match.replace(prefix, "")  # Remove prefix
                
                # Ensure taxon ID node exists
                if taxon_id not in node_map:
                    node_map[taxon_id] = g.add_vertex(name=taxon_id)
                
                # Link prefix node → taxon ID
                g.add_edge(node_map[prefix], node_map[taxon_id])
    
    return g


def dfs_traversalOld(graph, start_node, max_id):
    dfs_result = graph.dfs(start_node)[0]
    
    # Filter nodes to stop at the specified max_id
    filtered_nodes = [node for node in dfs_result if node <= max_id]
    
    # Convert lineage to tab-separated string
    lineage_string = "\t".join(map(str, filtered_nodes))
    
    return lineage_string  # Return the full lineage as a tab-separated string


def dfs_traversal(graph, start_node, max_id):
    lineage_dict = {}  # Dictionary to store {node: parent_node}
    visited = set()  # Keep track of visited nodes
    def dfs_helper(node, parent):
        if node in visited or node > max_id:
            return  # Stop if node is already visited or exceeds max_id
        visited.add(node)
        lineage_dict[node] = parent  # Store the parent
        for neighbor in graph.neighbors(node):
            dfs_helper(neighbor, node)  # Recursively visit neighbors
    # Start DFS
    dfs_helper(start_node, None)
    return lineage_dict  # Return dictionary of {node: parent_node}



def create_graph_from_df1(df):
    edges = list(df.iloc[:, :2].itertuples(index=False, name=None))
    unique_nodes = set(df.iloc[:, 0]).union(set(df.iloc[:, 1]))
    graph = Graph(directed=True)
    graph.add_vertices(list(unique_nodes))
    graph.add_edges(edges)
    return graph



def create_graph_from_df(df):
    edges = list(df.iloc[:, :2].itertuples(index=False, name=None))  # Extract edges
    unique_nodes = set(df.iloc[:, 0]).union(set(df.iloc[:, 1]))  # Find unique nodes
    graph = Graph(directed=True)
    graph.add_vertices(list(unique_nodes))
    graph.add_edges(edges)
    node_labels = dict(zip(df.iloc[:, 0], df["name"]))  # Map node to label
    graph.vs["label"] = [node_labels.get(v["name"], "Missing") for v in graph.vs]
    return graph


def create_graph_from_df2(df):
    edges = list(df.iloc[:, :2].itertuples(index=False, name=None))  # Extract edges
    unique_nodes = set(df.iloc[:, 0]).union(set(df.iloc[:, 1]))  # Find unique nodes
    graph = Graph(directed=True)
    graph.add_vertices(list(unique_nodes))
    graph.add_edges(edges)
    node_labels = dict(zip(df.iloc[:, 0], df["name"]))  # Map node to label
    graph.vs["label"] = [node_labels.get(v["name"], "Missing") for v in graph.vs]
    node_ranks = dict(zip(df.iloc[:, 0], df["rank"]))  # Map node to label
    graph.vs["rank"] = [node_ranks.get(v["name"], "Missing") for v in graph.vs]
    return graph


def apply_dfs_to_df(df, graph, node_col, max_id):
    node_mapping = {name: idx for idx, name in enumerate(graph.vs["name"])}
    df["lineage"]= df["uid"].apply(lambda node: dfs_traversal(graph, node_mapping.get(node, -1), node_mapping.get(max_id,-1)) if node in node_mapping else "")
    return df



def export_to_sigma(lineage_dict, output_file="graph.json"):
    """Converts DFS parent mapping into Sigma.js JSON format with circular layout."""

    nodes_list = list(lineage_dict.keys())  # Get list of nodes
    num_nodes = len(nodes_list)

    # Circular Layout Calculation
    radius = 200
    center_x, center_y = 300, 300  # Center of the graph
    angle_step = (2 * math.pi) / max(1, num_nodes)  # Avoid division by zero

    nodes = []
    for i, node in enumerate(nodes_list):
        angle = i * angle_step
        x = center_x + radius * math.cos(angle)
        y = center_y + radius * math.sin(angle)
        nodes.append({"id": str(node), "label": str(node), "x": x, "y": y, "size": 1})

    # Create edges from lineage_dict (ignoring None parent nodes)
    edges = [
        {"id": f"{parent}-{node}", "source": str(parent), "target": str(node)}
        for node, parent in lineage_dict.items() if parent is not None
    ]

    sigma_graph = {"nodes": nodes, "edges": edges}

    with open(output_file, "w") as f:
        json.dump(sigma_graph, f, indent=4)

    print(f"Graph exported to {output_file} for Sigma.js.")



import math

def export_igraph_to_sigma1(graph, output_file="graph.json"):
    num_nodes = graph.vcount()
    # Circular Layout Calculation
    radius = 200
    center_x, center_y = 300, 300  # Center of the graph
    angle_step = (2 * math.pi) / max(1, num_nodes)  # Avoid division by zero
    # Nodes
    nodes = []
    for i, vertex in enumerate(graph.vs):
        angle = i * angle_step
        x = center_x + radius * math.cos(angle)
        y = center_y + radius * math.sin(angle)
        nodes.append({
            "id": str(vertex.index),
            "label": vertex["label"] if "label" in vertex.attributes() else str(vertex.index),
            "x": x,
            "y": y,
            "size": 1
        })
    edges = [
        {"id": f"{edge.source}-{edge.target}", "source": str(edge.source), "target": str(edge.target)}
        for edge in graph.es
    ]
    sigma_graph = {"nodes": nodes, "edges": edges}
    with open(output_file, "w") as f:
        json.dump(sigma_graph, f, indent=4)
    print(f"Graph exported to {output_file} for Sigma.js.")


def get_subgraph(graph, nodes, node_mapping):
    Cfull=set()
    for i in nodes:
        Cin = graph.subcomponent(node_mapping.get(i), mode="in")
        Cout = graph.subcomponent(node_mapping.get(i), mode="out")
        Call = set(Cin).union(set(Cout))
        Cfull = set(Cfull).union(set(Call))
    return Cfull


def get_recursive_branches(graph, root, branch_counts):
    children = graph.neighbors(root, mode="in")  
    branch_counts[graph.vs[root].attributes()["label"]] = 1   # Count only first-level branches
    for child in children:
        get_recursive_branches(graph, child, branch_counts)
        branch_counts[graph.vs[root].attributes()["label"]] += 1   # Count only first-level branches
    return branch_counts


def get_branch_count_select(ottDfNotInWD):
    branchesX = Counter()
    for nodeI in ottDfNotInWD["uid"]:
        subC = node_mapping.get(nodeI)
        for i in ottDfGraph.dfsiter(subC, mode="out"):  # iterates without storing all nodes as against dfs()
            branchesX[i] += 1  # optimized counting using C's counter
    return branchesX



def get_branch_count_all(ottDfNotInWD):
    branches = {}
    for nodeI in ottDfNotInWD["uid"]:
        subC = node_mapping.get(nodeI)
        while subC not in branches.keys():
            subC = ottDfGraph.neighbors(subC, mode="out")[0]
            val = ottDfGraph.degree(ottDfGraph.subcomponent(subC, mode="in")) #does not matter the mode, whatever is 1 is good
            y = sum(x for x in val if x == 1)
            branches[subC] = y


def write_counter(file, dt):
    with open(file, "w") as f:
        for key, value in dt.items():
            f.write(f"{key}\t{value}\n")  # Tab-separated format


def compare_strings(str1, str2):   #Compare two strings using fuzzy matching
    return fuzz.ratio(str1, str2)

def is_valid_for_matching(text):    #Check if a string has at least two words
    return len(str(text).split()) >= 2

def ratioNIndex(str1, str2):
    fuzzRatio = compare_strings(str1, str2)
    if fuzzRatio >= THRESHOLD:
        return pd.Series([fuzzRatio, str2])
    else:
        return pd.Series([fuzzRatio, 0])

# Apply last 3 functions like so:
#ottNotInWdDf["Similarity Score", "index"] = ottNotInWdDf.apply(lambda row: ratioNIndex(str(row['name']), str(wdDf.loc[row.name, 'WdName'])), axis=1)

#ottNotInWdDf["Similarity Score"] = ottNotInWdDf.apply(
#    lambda row: compare_strings(str(row['name']), str(wdDf.loc[row.name, 'WdName']))
#    if is_valid_for_matching(row['name']) and is_valid_for_matching(wdDf.loc[row.name, 'WdName']) else 0,
#    axis=1
#)


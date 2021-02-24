#!/usr/bin/env python3

"""
============================================================================
Author: Kouiho, Sehou Romaric
Last update: Feb 2021
email: romarickouiho@gmail.com
============================================================================

Script to retrieve the depth of specific node in a directed acyclic
                                graph given it root node.

USAGE: Node_depth_DAG.py <arguments>

============================================================================
WARNING!!! Please run "python acyclic_graph.py --help" and carefully
            read the information about the usage of this script.
============================================================================
"""

# import statements
import argparse
import logging
import numpy as np

logging.basicConfig(format="%(levelname)s : %(message)s", level=logging.INFO)


# function definitions
def get_arguments_from_cmd():
    """ Collect all the given arguments from command line

    :return: a JSON object of the mapping of arguments with their values.
    """

    parser = argparse.ArgumentParser(description="\
    Script for getting the depth of a specific node in a directed \
    acyclic graph given it root node. \n")

    parser.add_argument("-i", type=str, required=True,
                        help="Path to the CSV or tab delimited file \
                        of graph Adjacency matrix of graph")
    parser.add_argument("-r", type=int, required=True,
                        help="Root node of the graph")
    parser.add_argument("-n", type=int, required=True,
                        help="Node of interest, "
                             "for which we want to know its depth")

    return parser.parse_args()


def parse_adjacency_matrix_csv(file_path, delimiter=","):
    """Parse CSV or txt file containing adjacency matrix of graph

    :param file_path: --str, path of the file
    :param delimiter: --str, delimiter used to separate each
            element on the same line in the file. Default is  ','.

    :return: list of list, representing the graph adjacency matrix
    """

    adjacency_matrix = []

    with open(file_path) as matrix:
        for line in matrix:
            if not line.strip():
                continue

            line = line.strip().split(delimiter)
            try:
                line = list(map(int, line))
            except ValueError:
                raise ValueError("Something is wrong with the input"
                                 " adjacency matrix file. "
                                 "Make sure there is no special "
                                 "character or comma at the end"
                                 " of each line.")
            adjacency_matrix.append(line)

    # check if all the row of the matrix have the same length
    list_length = [len(vector) for vector in adjacency_matrix]
    assert len(set(list_length)) == 1, f"Different length found for matrix " \
                                       f"row. Please correct."
    return adjacency_matrix


def unique_root_node_checking(adjacency_matrix):
    """ Check if there is only one root node and find it in the graph

    :param adjacency_matrix: --list of list, representing
            the graph adjacency matrix.

    :return: a tuple of root node and boolean True indicating that
            there is only one root node or a tuple of None and False
            otherwise.
    """

    # Sum across columns to get the in-degree of each vertex
    in_degrees = np.sum(adjacency_matrix, axis=0)
    # Sum across rows to get the out-degree of each vertex
    out_degrees = np.sum(adjacency_matrix, axis=1)

    # Retrieve node with in-degree equal to 0 and out-degree not 0
    # if there is only one node with in-degree = 0 and out-degree != 0
    # then it is the root node
    in_degrees_0 = np.where(in_degrees == 0)[0]
    # number of vertex (node) with in-degree == 0
    number_vertex = len(in_degrees_0)

    if number_vertex == 1 and out_degrees[in_degrees_0[0]] != 0:
        root_node = in_degrees_0[0]
        number_root_node = 1
    else:
        number_root_node = number_vertex
        root_node = "NA"
    return root_node, number_root_node


def adjacency_matrix_to_dict_graph(adjacency_matrix):
    """ Convert adjacency matrix to dictionary graph representation

    :param adjacency_matrix: --list of list, representing
            the graph adjacency matrix.

    :return: dictionary, with keys as nodes and list of integers
            representing the list of the nodes to which it is connected.
    """

    # Check if the matrix has the right shape
    number_row_edges = len(adjacency_matrix)
    number_col_edges = len(adjacency_matrix[0])
    assert len(adjacency_matrix) == len(adjacency_matrix[0]), \
        f"Expected number of row = number of columns. {number_row_edges}" \
        f" rows and {number_col_edges} columns found."

    return {i: (np.nonzero(row)[0]).tolist() for
            i, row in enumerate(adjacency_matrix)}


def node_is_visited(graph, node, path_taken, visited_nodes):
    """ Check  if a node in graph is visited and keep track of visited nodes
    
    :param graph: --dictionary, with keys as nodes and list of integers 
            representing the list of the nodes to which it is connected.
    :param node: --int, node of interest
    :param path_taken: --list, list of nodes on the path taken
    :param visited_nodes: --list, list of nodes already visited

    :return: a tuple of a boolean indication if the node is visited,
            list of nodes in path taken and list of visited nodes.
    """

    if node in visited_nodes:
        return False, path_taken, visited_nodes
    else:
        visited_nodes.append(node)
        path_taken.append(node)

    for adjacent_node in graph.get(node, ()):
        if adjacent_node in path_taken or \
                node_is_visited(graph, adjacent_node, path_taken,
                                visited_nodes)[0]:
            return True, path_taken, visited_nodes

    path_taken.remove(node)
    return False, path_taken, visited_nodes


def is_cyclic_graph(graph):
    """ Check if there is a cycle in a directed graph

    :param graph: --dictionary, with keys as nodes and list of integers
            representing the list of the nodes to which it is connected.

    :return: a boolean indicating if there is a cycle or not in graph.
    """

    path_taken, visited_nodes, result = [], [], False
    nodes_in_graph = list(graph.keys())
    # check if the first node in list in visited
    result, path_taken, visited_nodes = node_is_visited(graph,
                                                        nodes_in_graph[0],
                                                        path_taken,
                                                        visited_nodes)
    list_results = [result]
    graph_size, results_size = len(nodes_in_graph), len(list_results)

    # iterate over each node till all nodes are visited
    while not result and graph_size == results_size:
        for node in nodes_in_graph[1:]:
            result, path_taken, visited_nodes = \
                node_is_visited(graph, node, path_taken, visited_nodes)
            list_results.append(result)
            results_size = len(list_results)
    return result


def find_incoming_nodes(graph, node):
    """ Retrieve list of incoming nodes to a node of interest

    :param graph: --dictionary, with keys as nodes and list of integers
            representing the list of the nodes to which it is connected.
    :param node: --int, node of interest

    :return: list of incoming nodes.
    """

    incoming_nodes = []
    for vertex in graph.keys():
        if node in graph.get(vertex, ()):
            incoming_nodes.append(vertex)

    return incoming_nodes


def get_node_depth(dag_graph, root_node, node_interest, path_visited=[]):
    """ Get the depth of a specific node

    :param dag_graph: --dictionary, with keys as nodes and list of integers
            representing the list of the nodes to which it is connected.
    :param root_node: --int, root node of the graph.
    :param node_interest: --int, node of interest.
    :param path_visited: --list, list of nodes on the path visited

    :return: integer indicating the depth of node.
    """

    shortest_route = None
    path_visited, node_depth = path_visited + [node_interest], 0

    # if the node of interest is == to the root node then no depth
    if node_interest == root_node:
        return path_visited, node_depth

    # for every incoming node to node of interest, backtrack till root/end
    for node in find_incoming_nodes(dag_graph, node_interest):
        if node not in path_visited:
            # get a new route
            new_route_visited = get_node_depth(dag_graph, root_node,
                                               node, path_visited)[0]
            if new_route_visited:
                if not shortest_route or \
                        len(new_route_visited) < len(shortest_route):
                    shortest_route = new_route_visited

    if shortest_route:
        node_depth = len(shortest_route) - 1
    else:
        node_depth = node_depth
    return shortest_route, node_depth


def main(args):
    """Wrapper main function for executing the script

    :param args: JSON object of the mapping of arguments with their values.

    :return: returns None. Print output to command line.S
    """

    # Take the adjacency matrix of the graph
    adjacency_mat_graph = parse_adjacency_matrix_csv(args.i)

    # Check if the number of row == number of columns of the matrix
    logging.info("Checking the matrix shape...")
    number_row_edges = len(adjacency_mat_graph)
    number_col_edges = len(adjacency_mat_graph[0])
    if number_row_edges != number_col_edges:
        raise ValueError(f"The number of nodes detected in row"
                         f"{number_row_edges} not equal to the "
                         f"number of nodes found in column "
                         f"{number_col_edges}.")

    # check and verify if graph has only one root
    logging.info("Checking if graph has only one root and verify...")
    root, number_root = unique_root_node_checking(adjacency_mat_graph)

    if root == "NA" and number_root == 0:
        raise ValueError(f"No root node found in graph")
    elif root == "NA" and number_root > 0:
        raise ValueError("More than one root node found in graph")

    # check if identify root is the same as the one provided
    if root != "NA" and root != args.r:
        raise ValueError(f"Expected {args.r} as root node but {root}"
                         f" found in graph.")

    # Convert adjacency matrix to dict representation for further steps
    logging.info("Converting matrix to dictionary graph...")
    dict_graph = adjacency_matrix_to_dict_graph(adjacency_mat_graph)

    # Check if provided node is in the graph
    logging.info("Checking if node of interest n in graph...")
    if args.n not in dict_graph:
        raise ValueError(f"Provided node of interest {args.n} "
                         f"is not found in the provided graph.")

    # Check if the provided graph is acyclic
    logging.info("Checking if graph is acyclic...")
    is_cyclic = is_cyclic_graph(dict_graph)
    if is_cyclic:
        raise ValueError("The provided graph is not acyclic. "
                         "At least one cycle was detected in graph.")
    logging.info("THE PROVIDED GRAPH IS DAG !!!")

    # Get depth
    logging.info(f"Getting the depth of node {args.n}...")
    path, depth = get_node_depth(dict_graph, args.r, args.n, [])
    assert len(path) -1 == depth, f"Length path {path} != depth {depth}"
    logging.info(f"The depth of node {args.n} is equal to {depth}")


if __name__ == '__main__':
    # get arguments from command line
    arguments = get_arguments_from_cmd()

    # run wrapper function
    main(arguments)

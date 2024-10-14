#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
from pathlib import Path
from networkx import (
    DiGraph,
    all_simple_paths,
    lowest_common_ancestor,
    has_path,
    random_layout,
    draw,
    spring_layout,
)

import networkx as nx

import matplotlib
from operator import itemgetter
import random

random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
from typing import Iterator, Dict, List

from pathlib import Path
from typing import Iterator

from itertools import combinations 



matplotlib.use("Agg")

__author__ = "Rose TCHALA SARE"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Rose TCHALA SARE"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Rose TCHALA SARE"
__email__ = "r.tchalasare@etu.u-paris.fr"
__status__ = "Developpement"


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments():  # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(
        description=__doc__, usage="{0} -h".format(sys.argv[0])
    )
    parser.add_argument(
        "-i", dest="fastq_file", type=isfile, required=True, help="Fastq file"
    )
    parser.add_argument(
        "-k", dest="kmer_size", type=int, default=22, help="k-mer size (default 22)"
    )
    parser.add_argument(
        "-o",
        dest="output_file",
        type=Path,
        default=Path(os.curdir + os.sep + "contigs.fasta"),
        help="Output contigs in fasta file (default contigs.fasta)",
    )
    parser.add_argument(
        "-f", dest="graphimg_file", type=Path, help="Save graph as an image (png)"
    )
    return parser.parse_args()

#Fonction qui retourne une liste de fichiers fastq
def read_fastq(fastq_file: Path) -> Iterator[str]:
    with open(fastq_file, "r") as file:
        sequence_count = 0
        for line in file:
            if line.startswith(("A", "T", "G", "C", "U")):
                sequence_count += 1
                yield line.strip("\n")  # Use yield to return each sequence as an iterator
        #print(f"There were {sequence_count} sequences in the fastq file.")


#Fonction qui retourne les combinaisons de 5 enchainement de nucléotides au sein de la séquence
def cut_kmer(read: str, kmer_size: int) -> Iterator[str]:
    for i in range(len(read) - kmer_size + 1):
        yield read[i:i + kmer_size]




def build_kmer_dict(fastq_file: Path, kmer_size: int) -> Dict[str, int]:
        sequences = read_fastq(fastq_file)  # Read sequences from the FASTQ file
        kmer_dict = {}

        for sequence in sequences:
            # Generate k-mers for the current sequence
            for kmer in cut_kmer(sequence, kmer_size):
                if kmer in kmer_dict:
                    kmer_dict[kmer] += 1
                else:
                    kmer_dict[kmer] = 1
        
        return kmer_dict





def build_graph(kmer_dict: Dict[str, int]) -> DiGraph:
    G = nx.DiGraph()
    for kmer, weight in kmer_dict.items():
        prefix = kmer[:-1]
        suffix = kmer[1:]
        if G.has_edge(prefix, suffix):
            G[prefix][suffix]['weight'] += weight
        else:
            G.add_edge(prefix, suffix, weight=weight)
    return G

def get_starting_nodes(graph):
    starting_nodes = []
    for node in graph.nodes():
        predecessors = list(graph.predecessors(node))
        if not predecessors: 
            starting_nodes.append(node)
    return starting_nodes

def get_sink_nodes(graph: DiGraph) -> List[str]:
    sink_nodes = []
    for node in graph.nodes():
        # Get the predecessors
        successors = list(graph.successors(node))
        if not successors:  
            sink_nodes.append(node)
    return sink_nodes


def get_contigs(
    graph: DiGraph, starting_nodes: List[str], ending_nodes: List[str]
) -> List:
    contigs = []
    for starting_node in starting_nodes :
        for ending_node in  ending_nodes :
            if nx.has_path(graph, starting_node, ending_node):
               paths =nx.all_simple_paths(graph, starting_node, ending_node)
               for path in paths:
                    contig = path[0]  
                    
                    for node in path[1:]:
                        contig += node[-1] 
                        contig_length = len(contig) 

                    if (contig, contig_length) not in contigs:
                        contigs.append((contig, contig_length))
    return contigs
  

def save_contigs(contigs_list: List[str], output_file: Path) -> None:
    with open(output_file,"w") as file:
        for i in range(0, len(contigs_list)):
            path = contigs_list[i]
            sequence = path[0]
            length = path[1]
            ligne = f">contig_{i} len={length}\n"
            ligne_sequence = textwrap.fill(sequence, width=80)
            file.write(ligne)
            file.write(ligne_sequence)




def remove_paths(
    graph: DiGraph,
    path_list: List[List[str]],
    delete_entry_node: bool,
    delete_sink_node: bool,
) -> DiGraph:
    for path in path_list:
        if delete_sink_node and delete_entry_node:
            graph.remove_nodes_from(path)
        elif delete_entry_node:
            graph.remove_nodes_from(path[:-1])  
        elif delete_sink_node:
            graph.remove_nodes_from(path[1:]) 
        else:
            graph.remove_nodes_from(path[1:-1])

    return graph

def select_best_path(
    graph: DiGraph,
    path_list: List[List[str]],
    path_length: List[int],
    weight_avg_list: List[float],
    delete_entry_node: bool = False,
    delete_sink_node: bool = False,
) -> DiGraph:
    if len(weight_avg_list) >= 2 and statistics.stdev(weight_avg_list) != 0:
        max_weight_index = weight_avg_list.index(max(weight_avg_list))

        del path_list[max_weight_index]
    elif len(path_length) >= 2 and statistics.stdev(path_length) != 0:
        max_length_index = path_length.index(max(path_length))
        del path_list[max_length_index]
    else:
        random_index = random.randint(0, len(path_list) - 1)
        del path_list[random_index]

    graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    return graph



def path_average_weight(graph: DiGraph, path: List[str]) -> float:
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean(
        [d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)]
    )


def solve_bubble(graph: DiGraph, ancestor_node: str, descendant_node: str) -> DiGraph:
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """

    #On va récuperer la liste de tous les chemins possibles entre un noeud ancètre et descendant
    paths = list(all_simple_paths(graph, ancestor_node, descendant_node))

    #On récuperer la liste de tous les poids moyens (average_weights) des chemins ainsi que leur longueur
    average_weight = []
    longueur = []

    for path in paths :
        average_weight.append(path_average_weight(graph, path))
        longueur.append(len(path)-1)

    #On retourne le choix final
    return select_best_path(graph, paths, longueur, average_weight)
    


def simplify_bubbles(graph: DiGraph) -> DiGraph:
    bubble_found = False
    ancestor_node = None

    for node in graph.nodes():
        predecessors = list(graph.predecessors(node))
        
        if len(predecessors) > 1:
            for pred_a, pred_b in combinations(predecessors, 2):
                ancestor_node = lowest_common_ancestor(graph, pred_a, pred_b)
                
                if ancestor_node is not None:
                    bubble_found = True
                    break
                    
        if bubble_found:
            break

    if bubble_found:
        graph = simplify_bubbles(solve_bubble(graph, ancestor_node, node))

    return graph



def solve_entry_tips(graph: DiGraph, starting_nodes: List[str]) -> DiGraph:
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :param starting_nodes: (list) A list of starting nodes
    :return: (nx.DiGraph) A directed graph object
    """

    import networkx as nx

def solve_entry_tips(graph: nx.DiGraph, entry_nodes: list) -> nx.DiGraph:
    """Simplify entry tips in the graph by removing undesirable entry paths.

    :param graph: (nx.DiGraph) A directed graph object
    :param entry_nodes: (list) A list of entry nodes in the graph
    :return: (nx.DiGraph) A directed graph object with simplified entry paths
    """

    nodes_to_simplify = []
    
    for node in graph.nodes():
        predecessors = list(graph.predecessors(node))
        
        connected_entry_nodes = [pred for pred in predecessors if pred in entry_nodes]
        
        if len(connected_entry_nodes) > 1:
            nodes_to_simplify.append(node)

    for node in nodes_to_simplify:
        edges_to_remove = []
        
        for entry in entry_nodes:
            if has_path(graph, entry, node):
                if entry != node: 
                    edges_to_remove.append((entry, node))
        
        for edge in edges_to_remove:
            if edge in graph.edges:
                graph.remove_edge(*edge)

    return graph



def solve_out_tips(graph: DiGraph, ending_nodes: List[str]) -> DiGraph:
    for current_node in graph.nodes:
        connected_ending_nodes = [
            successor for successor in graph.successors(current_node)
            if any(has_path(graph, successor, exit) for exit in ending_nodes) or successor in ending_nodes
        ]

        if len(connected_ending_nodes) > 1:
            potential_paths = []

            for successor in connected_ending_nodes:
                if successor in ending_nodes:
                    new_paths = [[current_node, successor]] 
                else:
                    new_paths = [
                        path for exit in ending_nodes
                        for path in all_simple_paths(graph, current_node, exit)
                        if path[1] == successor  
                    ]

                for path in new_paths:
                    if len(path) > 1:  
                        potential_paths.append(path)

            if len(potential_paths) > 1:
                average_weights = [path_average_weight(graph, path) for path in potential_paths]
                path_lengths = [len(path) for path in potential_paths]

                graph = select_best_path(
                    graph, potential_paths, path_lengths, average_weights,
                    delete_entry_node=False, delete_sink_node=True
                )

                ending_nodes = get_sink_nodes(graph)

                return solve_out_tips(graph, ending_nodes)

    return graph


def draw_graph(graph: DiGraph, graphimg_file: Path) -> None:  # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (Path) Path to the output file
    """
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d["weight"] > 3]
    # print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d["weight"] <= 3]
    # print(elarge)
    # Draw the graph with networkx
    # pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(
        graph, pos, edgelist=esmall, width=6, alpha=0.5, edge_color="b", style="dashed"
    )
    # nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file.resolve())


# ==============================================================
# Main program
# ==============================================================
def main() -> None:  # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()


    #lecture du fichier fastq
    a = read_fastq(args.fastq_file)
    #print(a)

    #Construction du dictionnaire de k-mers
    d = build_kmer_dict(a)

    #Construction du graph
    g = build_graph(d)
    #print(graph.nodes)

    e = get_starting_nodes(g) ; g = solve_entry_tips(g, e)
    f = get_sink_nodes(g) ; g = solve_out_tips(f, e)

    # Ecriture du/des contigs 
    starting_nodes = get_starting_nodes(g)
    ending_nodes = get_sink_nodes(g)
    contigs = get_contigs(g, starting_nodes, ending_nodes)
    


    save_contigs(contigs, args.output_file)

          


    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)


if __name__ == "__main__":  # pragma: no cover
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 31 12:11:18 2018

"""

from Bio import SeqIO
from collections import defaultdict
from graphviz import Digraph 
import argparse

def parse_inputs():
     parser = argparse.ArgumentParser(description='Graph de Brujin')
     parser.add_argument('-i', '--infile' , help='Input file' , metavar='Str',
                    type=str, required=True)
     parser.add_argument('-o', '--outfile' , help='Input file' , metavar='Str',
                    type=str, required=True)
     parser.add_argument('-k', '--kmersize', help = 'Kmer size', metavar = 'Int', type = int, default = 15)
     parser.add_argument('-f', '--full', help='Output = full graph', action='store_true')
     parser.add_argument('-s', '--short', help='Output = short graph', action='store_true')
     parser.add_argument('-c', '--compress', help='Output = compressed graph', action='store_true')
     
     args = parser.parse_args()
     return args.infile, args.outfile, args.kmersize, args.full, args.short, args.compress
    

class Vertex:
    
    def __init__(self, seq):
        self.seq = seq
        self.coverage = 1
        self.in_edges = {}
        self.out_edges = {}
        
    def increase_coverage(self):
        self.coverage += 1

class Edge:
    
    def __init__(self,k1,k2):
        self.seq = k1 + k2[-1]
        self.n = 2
        self.coverage = 0
    
    def calc_coverage(self, c1,c2):
        self.coverage = (c1+c2)/2


class Graph:

    def __init__(self,k):
        self.vertices = {}
        self.k = k
        self.bad_vert = []
    def add_read(self, read):
        read_lng = len(read)
        if read_lng < self.k:
            return
            
        kmer = read[:k]
        if kmer in self.vertices:
            self.vertices[kmer].increase_coverage()
        else:
            self.vertices[kmer] = Vertex(kmer)
        
        for next_kmer_index in range(1,read_lng-k+1,1):
            next_kmer = read[next_kmer_index:(next_kmer_index+k)]
            if next_kmer in self.vertices:
                self.vertices[next_kmer].increase_coverage()
            else:
                self.vertices[next_kmer] = Vertex(next_kmer)
            
            new_edge = Edge(kmer,next_kmer)
            
            self.vertices[next_kmer].in_edges[kmer]  = [new_edge]
            
            self.vertices[kmer].out_edges[next_kmer] = [new_edge]

            kmer = next_kmer
    
    def calc_init_edge_coverage(self):
        
        for current_vertex in self.vertices.keys():
            for next_vertex in self.vertices[current_vertex].out_edges.keys():
                self.vertices[current_vertex].out_edges[next_vertex][0].calc_coverage(self.vertices[current_vertex].coverage,self.vertices[next_vertex].coverage)
    
    
  
#   в узлах = конкретный камер. в ребрах - камер+1, то есть ребро перекрывается с первым и вторым 
        
    def vizualize(self):
#       создаем объект - связный граф на языке DOT
        dot = Digraph()
        
        if compress_flag:
            a= size+1
        else:
            a = size
#        если указан флаг вывода графа в полном формате    
        if full_flag:
#            все наши вершины будут вершинами графа
            for v in self.vertices:
                if len(v) == a:
                    dot.node(v)
    #                для каждой вершины мы создаем ребро со следующей за ней вершиной (которая есть в словаре исходящих вершин) . 
                    for n in self.vertices[v].out_edges:
                            dot.edge(v,n, label = self.vertices[v].out_edges[n][0].seq)
                
#       если выбран флаг вывода графа в сокращенной форме, выводим в качестве лейблов узлов покрытие каждого камера,
#       а в качестве лейблов ребер - покрытие и длину ребра (то есть длину камера + 1)
        if short_flag:
            for v in self.vertices:
                if len(v) == a:
                    dot.node(v, label = str(self.vertices[v].coverage))
                    for n in self.vertices[v].out_edges:
                            dot.edge(v,n, label = (str(self.vertices[v].out_edges[n][0].coverage) +' '+ str(self.k+1)))
                    
#       сохранение и выводод на экран     
        dot.render(out_file, view=True)
        
        
    def compact(self):
  
# проходим по графу, выписываем список вершин, которые надо сжать (один входной и один выходной узел)
           
        for v in self.vertices: 
            if (len(self.vertices[v].in_edges) == 1) and (len(self.vertices[v].out_edges) == 1):
                self.bad_vert.append(v)
           
# проходим по списку вершин, которые нужно схлопнуть
        for v in self.bad_vert:

# обозначаем новую вершину - это старая вершина + часть выходящей из нее         
            next_vert = list(self.vertices[v].out_edges)[0]
            prev_vert = list(self.vertices[v].in_edges)[0]
            new_vert = v + str(next_vert)[-1]

            self.vertices[new_vert] = Vertex(new_vert)
            
# обозначаем выходящее из новой вершины ребро - это ребро между новой вершиной и 
# выходящим из следующей за новой вершиной ребром           
            if len(list(self.vertices[next_vert].out_edges)) >= 1:
                next_out_vert = list(self.vertices[next_vert].out_edges)[0]
                new_edge1 = Edge(new_vert, next_out_vert)
                self.vertices[next_out_vert].in_edges[new_vert] = [new_edge1] 
                self.vertices[new_vert].out_edges[next_out_vert] = [new_edge1] 
    
                
# обозначаем входящее из новой вершины ребро - это ребро между новой вершиной и 
# входящим в предыдущую вершину ребро               
            if len(list(self.vertices[prev_vert].in_edges)) >= 1:
                prev_in_vert = list(self.vertices[prev_vert].in_edges)[0]
                new_edge2 = Edge(prev_in_vert, new_vert)
                self.vertices[new_vert].in_edges[prev_in_vert] = [new_edge2]
                self.vertices[prev_in_vert].out_edges[new_vert] = [new_edge2]
# добавляем информацию о покрытии  
            self.vertices[new_vert].coverage = self.vertices[v].coverage 
                
if __name__ == '__main__':
    in_file, out_file, size, full_flag, short_flag, compress_flag = parse_inputs()
    
    dataset = in_file

    k = size
    
    my_graph = Graph(k)
    

    with open(dataset, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
          
            read = str(record.seq)
            my_graph.add_read(read)
            my_graph.add_read( str(record.reverse_complement().seq) )
 
    
            
            
#   вызываем функции для подсчета покрытия и визуализации       
    my_graph.calc_init_edge_coverage()

# вызываем функцию сжатия
    my_graph.compact()
           
    for v in my_graph.vertices:
        if len(v) == k+1:
            print('Vertex: {}, coverage: {}'.format(v,my_graph.vertices[v].coverage))
            for e in my_graph.vertices[v].out_edges:
                print('-> Out edge: {}'.format(e))
            for e in my_graph.vertices[v].in_edges:
                print('-> In edge: {}'.format(e)) 
                
    my_graph.vizualize()        
            

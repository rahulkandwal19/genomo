//============================================================================================
//    GenomoKit/pathogen_network : perform operations on pathogen genetic graph network
//============================================================================================
#include<iostream>
#include<vector>
#include<tuple>
#include<queue>
#include<string>
#include<sstream>
#include<fstream>
#include"sequence_matcher.hpp"
using namespace std;

namespace GenomoKit{
    
    class Info{
        public :
            string pathogen_id;
            string pathogen_sequence;
            Info(string pathogen_id="NAN0000", string pathogen_sequence=""){
                this->pathogen_id = pathogen_id;
                this->pathogen_sequence = pathogen_sequence;
            }
    };

    class Cds{
        public :
            string cds_id;
            string cds_sequence;
            Cds(string cds_id, string cds_sequence){
                this->cds_id = cds_id;
                this->cds_sequence = cds_sequence;
            }
    };

    class Edge{
        public:
            int destination;
            float weight;
            Edge(int destination, int weight){
                this->destination = destination;
                this->weight = weight;
            }
    };

    class Node{
        public:
            int node_idx;
            Info info;
            vector<Edge*> edge_list;
            vector<Cds*> cds_list;

            Node(int node_idx, string pathogen_id, string pathogen_sequence, vector<Cds*>& cds_list){
                this->node_idx = node_idx;
                this->info = Info(pathogen_id,pathogen_sequence);
                this->cds_list = cds_list;
            }
            Node(int node_idx,string pathogen_id, string pathogen_sequence){
                this->node_idx = node_idx;
                this->info = Info(pathogen_id,pathogen_sequence);
                this->cds_list = vector<Cds*>();
            }
    };

    class Graph{
        private:
            int counter = 0;
        public:
            vector<Node*> node_list;
            void insert_node(string pathogen_id, string pathogen_sequence,vector<Cds*>&cds_list){
                Node* new_node = new Node(Graph::counter,pathogen_id,pathogen_sequence,cds_list);
                node_list.push_back(new_node);
                Graph::counter++;
            }
    
            void add_edge(int source_Idx, int destination_Idx, float similarity){
                Edge* new_edge = new Edge(destination_Idx,similarity);
                node_list[source_Idx]->edge_list.push_back(new_edge);
            }

            vector<int> get_kmer_based_matches(Node* pathogen_node, float kmer_match_threshold=0.05,int kmer=5){
                GenomoKit::SequenceMatcher sequence_matcher = GenomoKit::SequenceMatcher();
                queue<int> traversal_queue;
                vector<bool> visited (node_list.size(),false);
                
                vector<int> matched_nodes;

                for(int i=0; i<node_list.size();i++){
                    if(!visited[i]){
                        traversal_queue.push(i);
                        visited[i] = true;
                    }

                    while(!traversal_queue.empty()){
                        int current_node_index = traversal_queue.front();
                        traversal_queue.pop();

                        float match_score = sequence_matcher.find_kmer_similarity(node_list[current_node_index]->info.pathogen_sequence,
                                                                                  pathogen_node->info.pathogen_sequence,kmer);
                        
                        if(match_score >= kmer_match_threshold) matched_nodes.push_back(current_node_index);
                        
                        vector<GenomoKit::Edge*> edges = node_list[current_node_index]->edge_list;
                        for(GenomoKit::Edge* edge:edges){
                            if(!visited[edge->destination]){
                                traversal_queue.push(edge->destination);
                                visited[edge->destination] = true;
                            }   
                        }
                    }
                }
                return matched_nodes;
            }
    };



    class PathogenNetwork{        
        public:
            Graph graph = Graph();
            PathogenNetwork(){}
            
            void add_pathogen(string pathogen_id, string pathogen_sequence,vector<Cds*>cds){
                graph.insert_node(pathogen_id,pathogen_sequence,cds); 
            }
            
            void add_connection(int source_idx,int destination_idx,int similarityScore){
                graph.add_edge(source_idx,destination_idx,similarityScore);
                graph.add_edge(destination_idx,source_idx,similarityScore);
            }  

            vector<int> get_kmer_based_matches(Node* pathogen_node, float kmer_match_threshold=0.05,int kmer=5){
                return graph.get_kmer_based_matches(pathogen_node,kmer_match_threshold,kmer);
            }

            void save(string file_name){
                fstream file;
                file.open(file_name + ".gpn", ios::out);

                if (!file.is_open()) {
                    cerr << "Failed to create file for saving" << endl;
                    return;
                }

                vector<Node*> nodes = graph.node_list;
                for (Node* node : nodes) {
                    file << node->node_idx << "|";
                    file<< node->info.pathogen_id << "|";
                    file<< "["<<"data1"<<","<<"data2"<< "]|";
                    file<<node->info.pathogen_sequence<<":\n";
        
                    for (Cds* cds: node->cds_list) {
                        file <<cds->cds_id<<"|"<<cds->cds_sequence<<",";
                    }
                    file << ":\n";

                    for(Edge* edge : node->edge_list){
                        file <<"["<<edge->destination<<"-"<<edge->weight<<"]"<<",";
                    }
                    file<<";\n\n";
                }
        
                file.close();
            }  

            PathogenNetwork(string path_name) {
              // Load graph From Memory
            }
    };
}
//============================================================================================
//                    End of pathogen_network module from GenomoKit
//============================================================================================
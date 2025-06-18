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
            int threatLevel;
            Info(string pathogen_id="NAN0000", string pathogen_sequence="", int threatLevel = 3){
                this->pathogen_id = pathogen_id;
                this->pathogen_sequence = pathogen_sequence;
                this->threatLevel = threatLevel;
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
            int weight;
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

            Node(int node_idx, string pathogen_id, string pathogen_sequence, vector<Cds*>& cds_list,int threat = 4){
                this->node_idx = node_idx;
                this->info = Info(pathogen_id,pathogen_sequence,threat);
                this->cds_list = cds_list;
            }
            Node(int node_idx,string pathogen_id, string pathogen_sequence){
                this->node_idx = node_idx;
                this->info = Info(pathogen_id,pathogen_sequence);
                this->cds_list = vector<Cds*>();
            }
    };

    class Graph{
        public:
            vector<Node*> node_list;
            void insert_node(string pathogen_id, string pathogen_sequence,vector<Cds*>&cds_list, int threat = 4){
                Node* new_node = new Node(node_list.size(),pathogen_id,pathogen_sequence,cds_list, threat);
                node_list.push_back(new_node);
            }

            int insert_node(Node* pathogen_node){
                Node* new_node = new Node(node_list.size(),pathogen_node->info.pathogen_id,pathogen_node->info.pathogen_sequence,
                                         pathogen_node->cds_list,pathogen_node->info.threatLevel);
                node_list.push_back(new_node);
                return new_node->node_idx;
            }
            
    
            void add_edge(int source_Idx, int destination_Idx, int similarity){
                Edge* new_edge = new Edge(destination_Idx,similarity);
                node_list[source_Idx]->edge_list.push_back(new_edge);
            }

            vector<int> get_kmer_based_matches(Node* pathogen_node, float kmer_match_threshold,int kmer){
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
                        
                       
                        if(match_score >= kmer_match_threshold){
                            matched_nodes.push_back(current_node_index);
                        }
                        
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
            
            void add_pathogen(string pathogen_id, string pathogen_sequence,vector<Cds*>cds,int threat =4){
                graph.insert_node(pathogen_id,pathogen_sequence,cds,threat); 
            }
            
            void add_connection(int source_idx,int destination_idx,int similarityScore){
                graph.add_edge(source_idx,destination_idx,similarityScore);
                graph.add_edge(destination_idx,source_idx,similarityScore);
            }  

            vector<int> get_kmer_based_matches(Node* pathogen_node, float kmer_match_threshold,int kmer){
                return graph.get_kmer_based_matches(pathogen_node,kmer_match_threshold,kmer);
            }


            void insert_node_on_dna_similarity_connections(Node* pathogen_node,float kmer_threshold,int kmer_size,float alignment_threshold){
                pathogen_node->node_idx = graph.insert_node(pathogen_node);
                //[1.] Find Nodes in Graph that have high KMer Similarity with new node by traversing all 
                // KMer Based selection is faster (Helpful in rough selection from compleate/large network or dataset).
                vector<int> matched_nodes_indices = get_kmer_based_matches(pathogen_node,kmer_threshold,kmer_size);
                //[2.] Find Local_Aligned Subsequence based Simaliraty (Smith-Waterman Algrothem for DNA Strings)
                //From selected nodes CDS, filter them with high similarity >Threshold using above mentioned Algo.
                GenomoKit::SequenceMatcher sequence_matcher = GenomoKit::SequenceMatcher();
                for(int current_node_index:matched_nodes_indices){
                    if(pathogen_node->node_idx == current_node_index) continue;
                    float similarity_score = -1;

                    //Nodes may have multiple CDS, Hence mean of similarity is computed
                    for(GenomoKit::Cds* cds : graph.node_list[current_node_index]->cds_list){
                        float score = sequence_matcher.local_alignment_similarity(cds->cds_sequence,
                                                              pathogen_node->info.pathogen_sequence);
                                                                                 
                        //Save matched cds for result report                                                           

                        if (similarity_score == -1) similarity_score = score;
                        else similarity_score += score;
                    }
                    similarity_score = similarity_score/ graph.node_list[current_node_index]->cds_list.size();

                //[3.] Add edges between filtered node and new node.
                    if(similarity_score>= alignment_threshold){
                        add_connection(pathogen_node->node_idx,current_node_index,int(similarity_score*10000));
                    }
                }
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
                    file<< "["<<node->info.threatLevel<<","<<"data2"<< "]:\n";
                    file<<node->info.pathogen_sequence<<":\n";
                    
                    bool first = true;
                    for (Cds* cds: node->cds_list) {
                        if(!first)file<<"\n";
                        else first = false;

                        file <<cds->cds_id<<"|"<<cds->cds_sequence;
                    }
                    file << ":\n";

                    first = true;
                    for(Edge* edge : node->edge_list){
                        if(!first)file<<",";
                        else first = false;
                        file <<"["<<edge->destination<<"-"<<edge->weight<<"]";
                    }
                    file<<";\n\n";
                }
        
                file.close();
            }  

            PathogenNetwork(string path_name) {
                fstream file;
                file.open(path_name + ".gpn", ios::in);
                if (!file.is_open()) {
                    cerr << "Failed to open file for loading." << endl;
                    return;
                }
        
                string pathogen_record;
                //Extracting record one by one 
                while (getline(file, pathogen_record,';')) {

                    //Remove front & rear unwanted strings
                    int start = pathogen_record.find_first_not_of('\n');
                    int end = pathogen_record.find_last_not_of('\n');
                    pathogen_record= (start == std::string::npos) ? "" : pathogen_record.substr(start, end - start + 1);
                    if (pathogen_record.empty()) continue;

                    //Retrive diffrent segments of record
                    stringstream record_stream(pathogen_record);
                    string segment;
                    vector<string> segments;
                    while (std::getline(record_stream,segment,':')){
                        segments.push_back(segment);
                    }
                    string data_segment = segments[0];
                    string sequence_segment = segments[1];
                    string cds_segment = segments[2];
                    string edge_segment = segments[3];

                    //Retrive data from data segment
                    int seprator1 = data_segment.find('|');
                    int seprator2 = data_segment.find('|', seprator1 + 1);
                    int current_idx = stoi(data_segment.substr(0, seprator1));
                    string id = data_segment.substr(seprator1 + 1, seprator2 - seprator1 - 1);
                   
                    int seprator3 = data_segment.find('[', seprator2 + 1);
                    int seprator4 = data_segment.find(',', seprator3 + 1);
                    int threat = stoi(data_segment.substr(seprator3+1, seprator4));
                    //Retrive Data from sequence segment
                    string sequence = sequence_segment.substr(1,sequence_segment.size()-1);

                    //Retrive data from cds segment
                    cds_segment = cds_segment.substr(1,sequence_segment.size()-1);
                    stringstream cds_stream = stringstream(cds_segment);
                    string cds_data;
                    vector<GenomoKit::Cds*> pathogen_cds;
                    while(cds_stream>>cds_data){
                        cds_data = cds_data.substr(0,sequence_segment.size()-2);
                        int seperator = cds_data.find('|');
                        string cds_id = cds_data.substr(0, seperator);
                        string cds_sequence = cds_data.substr(seperator + 1, cds_data.size());
                        GenomoKit::Cds *newCds = new  GenomoKit::Cds(cds_id,cds_sequence);
                        pathogen_cds.push_back(newCds);
                    }
                    
                    // Add Node to Graph
                    graph.insert_node(id,sequence,pathogen_cds,threat);

                    //Retrive edges from edge segment
                    edge_segment = edge_segment.substr(1,sequence_segment.size()-1);
                    stringstream edge_stream(edge_segment);
                    string edge;
                    while (std::getline(edge_stream,edge,',')){

                        int seprator = edge.find('-');
                        int destination = stoi(edge.substr(1, seprator-1));
                        int weight = stoi(edge.substr(seprator+1, edge.size() -seprator - 2));
                        
                        //Add Edges in graph
                        graph.add_edge(current_idx,destination,weight);
                    }
                }                 
            }
    };
}
//============================================================================================
//                    End of pathogen_network module from GenomoKit
//============================================================================================
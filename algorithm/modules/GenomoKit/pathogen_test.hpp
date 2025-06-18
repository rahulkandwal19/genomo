//============================================================================================
//      GenomoKit/pathogen_test : provides methods to test dna strings for pathogen
//============================================================================================
#include "pathogen_network.hpp"

namespace GenomoKit{
    
    class ZoonosisTestResult{
        public:
            vector<pair<string,float>> matched_cds;
            vector<vector<string>> paths;
    };
    
    class ZoonosisTest{
        private:
            GenomoKit::PathogenNetwork pathogen_network;
            GenomoKit::Node *test_node; 
            ZoonosisTestResult result;
            
            int insert_node_on_dna_similarity_connections(Node* test_node,float kmer_threshold,int kmer_size,float alignment_threshold){
                test_node->node_idx = pathogen_network.graph.insert_node(test_node);
                //[1.] Find Nodes in Graph that have high KMer Similarity with new node by traversing all 
                // KMer Based selection is faster (Helpful in rough selection from compleate/large network or dataset).
                vector<int> matched_nodes_indices = pathogen_network.get_kmer_based_matches(test_node,kmer_threshold,kmer_size);

                //[2.] Find Local_Aligned Subsequence based Simaliraty (Smith-Waterman Algrothem for DNA Strings)
                //From selected nodes CDS, filter them with high similarity >Threshold using above mentioned Algo.
                vector<int> final_selected_nodes;
                GenomoKit::SequenceMatcher sequence_matcher = GenomoKit::SequenceMatcher();
                for(int current_node_index:matched_nodes_indices){
                    if(test_node->node_idx == current_node_index) continue;
                    float similarity_score = -1;

                    //Nodes may have multiple CDS, Hence mean of similarity is computed
                    for(GenomoKit::Cds* cds : pathogen_network.graph.node_list[current_node_index]->cds_list){
                        float score = sequence_matcher.local_alignment_similarity(cds->cds_sequence,
                                                                                  test_node->info.pathogen_sequence);
                                                           
                        //Save matched cds for result report                                                           
                        if(score>= alignment_threshold){       
                            result.matched_cds.push_back({cds->cds_id,score});
                        }   


                        if (similarity_score == -1) similarity_score = score;
                        else similarity_score += score;
                    }
                    similarity_score = similarity_score/ pathogen_network.graph.node_list[current_node_index]->cds_list.size();
                    cout<<"Align-Avg : "<<pathogen_network.graph.node_list[current_node_index]->info.pathogen_id<<" "<<similarity_score<<endl;  

                    //[3.] Add edges between filtered node and new node.
                    if(similarity_score>= alignment_threshold){
                        pathogen_network.add_connection(test_node->node_idx,current_node_index,similarity_score*100);
                    }
                }

                return test_node->node_idx;
            }

            vector<pair<int,int>> dijkstra(int src){
                priority_queue<vector<int>, vector<vector<int>>, greater<vector<int>>> pq;

                int V = pathogen_network.graph.node_list.size();
                vector<int> dist(V, INT_MAX);

                pq.push({0, src});
                dist[src] = 0;

                while (!pq.empty()){
                    int u = pq.top()[1];
                    pq.pop();

                    for (auto x : pathogen_network.graph.node_list[u]->edge_list){
                        int v = x->destination;
                        int weight = 10000 -  x->weight;

                        if (dist[v] > dist[u] + weight){
                            
                            dist[v] = dist[u] + weight;
                            pq.push({dist[v], v});
                        }
                    }
                }

                vector<pair<int,int>> result;
                for(int i=0;i<V;i++){
                    result.push_back({dist[i],i});
                }

                sort(result.begin(),result.end());

                return result;
        }

        vector<string> dijkstraShortestPath(int src, int dest) {
            int n = pathogen_network.graph.node_list.size();
            vector<int> dist(n, INT_MAX);
            vector<int> parent(n, -1);
            priority_queue<vector<int>, vector<vector<int>>, greater<vector<int>>> pq;
            dist[src] = 0;
            pq.push({src,0});

            while (!pq.empty()) {
                int current = pq.top()[0];
                pq.pop();

                if (current == dest) break;

                for (auto i : pathogen_network.graph.node_list[current]->edge_list) {
                    int weight = 10000 - i->weight;
                    int v = i->destination;

                    int newDist = dist[current] + weight;
                    if (newDist < dist[v]) {
                        dist[v] = newDist;
                        parent[v] = current;
                        pq.push({v, newDist});
                    }
                }
            }

            vector<string> path;
            for (int at = dest; at != -1; at = parent[at]) {
                path.push_back(pathogen_network.graph.node_list[at]->info.pathogen_id);
            }
            reverse(path.begin(), path.end());

            string srcID = pathogen_network.graph.node_list[src]->info.pathogen_id;
            return (path.front() == srcID) ? path : vector<string>{};
        }

        void getPossiblePaths(int src, int k, int threat){
            vector<pair<int,int>> shotestNodes = dijkstra(src);
            vector<int> destination;
            for(int i = 1; i< shotestNodes.size();i++){
                if(destination.size()==k) break;
                int nodeIdx = shotestNodes[i].second;
                if(pathogen_network.graph.node_list[nodeIdx]->info.threatLevel >= threat){
                    destination.push_back(nodeIdx);
                }
            }

            vector<vector<string>> paths;
            for(int i: destination){
                paths.push_back(dijkstraShortestPath(src,i));
            }
            
            this->result.paths =  paths;
        }


        public:

            ZoonosisTest(GenomoKit::PathogenNetwork pathogen_network,string test_sequence){
                int node_index = pathogen_network.graph.node_list.size();
                this->test_node = new GenomoKit::Node(node_index,"TEST001",test_sequence);
                this->pathogen_network = pathogen_network;
            }

            ZoonosisTestResult is_pathogen_zoonotic(float kmer_threshold,int kmer_size,float alignment_threshold, int threat = 4, int maxPaths=1){
                //[1.] Insert new node to correct place in graph network on basis of dna similarity score
                //While saving matching cds in ZoonosisTestResult object
                insert_node_on_dna_similarity_connections(test_node,kmer_threshold,kmer_size,alignment_threshold);

                //[2.] Find Shortest paths to zoonotic nodes
                getPossiblePaths(test_node->node_idx,maxPaths,threat);

                return result;
            }
    };
}

//============================================================================================
//                    End of pathogen_test module from GenomoKit
//============================================================================================
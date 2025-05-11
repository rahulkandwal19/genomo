//============================================================================================
//      GenomoKit/pathogen_test : provides methods to test dna strings for pathogen
//============================================================================================
#include "pathogen_network.hpp"

namespace GenomoKit{
    
    class ZoonosisTestResult{
        public:
            vector<pair<string,float>> matched_cds;
    };
    
    class ZoonosisTest{
        private:
            GenomoKit::PathogenNetwork pathogen_network;
            GenomoKit::Node *test_node; 
            ZoonosisTestResult result;

            int insert_node_on_dna_similarity_connections(Node* test_node,float kmer_threshold,int kmer_size,float alignment_threshold){
                //[1.] Find Nodes in Graph that have high KMer Similarity with new node by traversing all 
                // KMer Based selection is faster (Helpful in rough selection from compleate/large network or dataset).
                vector<int> matched_nodes_indices = pathogen_network.graph.get_kmer_based_matches(test_node,kmer_threshold,kmer_size);

                //[2.] Find Local_Aligned Subsequence based Simaliraty (Smith-Waterman Algrothem for DNA Strings)
                //From selected nodes CDS, filter them with high similarity >Threshold using above mentioned Algo.
                vector<int> final_selected_nodes;
                GenomoKit::SequenceMatcher sequence_matcher = GenomoKit::SequenceMatcher();
                for(int current_node_index:matched_nodes_indices){
                    
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


                //[3.] Add edges between filtered node and new node.
                    if(similarity_score>= alignment_threshold){
                        pathogen_network.add_connection(test_node->node_idx,current_node_index,similarity_score);
                    }
                }

                return test_node->node_idx;
            }

        public:

            ZoonosisTest(GenomoKit::PathogenNetwork pathogen_network,string test_sequence){
                int node_index = pathogen_network.graph.node_list.size();
                this->test_node = new GenomoKit::Node(node_index,"TEST001",test_sequence);
                this->pathogen_network = pathogen_network;
            }

            ZoonosisTestResult is_pathogen_zoonotic(float kmer_threshold,int kmer_size,float alignment_threshold){
                //[1.] Insert new node to correct place in graph network on basis of dna similarity score
                //While saving matching cds in ZoonosisTestResult object
                insert_node_on_dna_similarity_connections(test_node,kmer_threshold,kmer_size,alignment_threshold);

                //[2.] Find Shortest paths to zoonotic nodes

                //[3.] Estimate Zoonosis 

                return result;
            }
    };
}

//============================================================================================
//                    End of pathogen_test module from GenomoKit
//============================================================================================
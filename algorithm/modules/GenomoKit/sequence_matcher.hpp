//============================================================================================
//      GenomoKit/sequrnce_matcher : provides diffrent "DNA" string similarity methods 
//============================================================================================

#include<iostream>
#include<string>
#include<vector>
#include<algorithm>
#include<set>

using namespace std;

namespace GenomoKit{
    class SequenceMatcher{
        private:
            set<string> find_kmers(string sequence,int k){
                set<string> kmers = *new set<string>();
                for(int i=0;i<sequence.size()-k;i++){
                    kmers.insert(sequence.substr(i,k));
                }
                return kmers;
            }
            
        public:
            double find_kmer_similarity(string sequence1, string sequence2,int k=3){
                set<string> kmers_sequence1 = find_kmers(sequence1,k);
                set<string> kmers_sequence2 = find_kmers(sequence2,k);

                set<string> kmers_intersection;
                set_intersection(
                                kmers_sequence1.begin(),kmers_sequence1.end(),
                                kmers_sequence2.begin(),kmers_sequence2.end(),
                                inserter(kmers_intersection,kmers_intersection.begin())
                                );
                
                int set1_count = kmers_sequence1.size();
                int set2_count = kmers_sequence2.size();
                int intersection_count = kmers_intersection.size();
                int union_count = set1_count+set2_count-intersection_count;

                double similarity = (double)intersection_count/(double)union_count;
                return similarity;
            }    
    };
}
//============================================================================================
//                    End of sequrnce_matcher module from GenomoKit
//============================================================================================
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

            vector<string> smith_waterman_local_align(string sequence1, string sequence2){
                int gap_penalty = -3;
                int match_score = 1;
                int mismatch_penalty = -1;

                int n = sequence1.length() + 1;
                int m = sequence2.length() + 1;
                
                vector<vector<int>> scoreMatrix(n, vector<int>(m, 0));
                vector<vector<int>> traceback(n, vector<int>(m, 0));

                int max_i = 0, max_j = 0, max_score = 0;

                for (int i = 1; i < n; i++) {
                    for (int j = 1; j < m; j++) {
                        int match = sequence1[i - 1] == sequence2[j - 1] ? match_score : mismatch_penalty;
                        int score_diag = scoreMatrix[i - 1][j - 1] + match;
                        int score_up = scoreMatrix[i - 1][j] + gap_penalty;
                        int score_left = scoreMatrix[i][j - 1] + gap_penalty;
                        
                        scoreMatrix[i][j] = max({0, score_diag, score_up, score_left});

                        if (scoreMatrix[i][j] > max_score) {
                            max_score = scoreMatrix[i][j];
                            max_i = i;
                            max_j = j;
                        }
                    }
                }

                string aligned_seq1 = "", aligned_seq2 = "";
                int i = max_i, j = max_j;

                while (i > 0 && j > 0 && scoreMatrix[i][j] > 0) {
                    if (scoreMatrix[i][j] == scoreMatrix[i - 1][j - 1] + (sequence1[i - 1] == sequence2[j - 1] ? match_score : mismatch_penalty)) {
                        aligned_seq1 = sequence1[i - 1] + aligned_seq1;
                        aligned_seq2 = sequence2[j - 1] + aligned_seq2;
                        i--; j--;
                    } else if (scoreMatrix[i][j] == scoreMatrix[i - 1][j] + gap_penalty) {
                        aligned_seq1 = sequence1[i - 1] + aligned_seq1;
                        aligned_seq2 = "-" + aligned_seq2;
                        i--;
                    } else {
                        aligned_seq1 = "-" + aligned_seq1;
                        aligned_seq2 = sequence2[j - 1] + aligned_seq2;
                        j--;
                    }
                }

                vector<string> alligned_sequences = {aligned_seq1, aligned_seq2};
                return alligned_sequences;
            }
            
        public:

            // K-Mer based Sequence Similarity is fast method to roughly estimate similarity
            // Utilizing K-Mers of string and Jaccard index. K is variable parameter
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

            //It is a squence alignment based similarity score that utilizes 
            //smith-waterman algorithm to align strings and compute similarity
            //on local subsequence matching using no of matches in aligned sequence.
            double local_alignment_similarity(string sequence1, string sequence2){
                vector<string> aligned_sequences = smith_waterman_local_align(sequence1,sequence2);

                string sequence1_aligned = aligned_sequences[0];
                string sequence2_aligned = aligned_sequences[1];
    
                int match_count = 0;
                int length = sequence1_aligned.length();

                for (int i = 0; i < length; i++) {
                    if (sequence1_aligned[i] == sequence2_aligned[i] && sequence1_aligned[i] != '-') {
                        match_count++;
                    }
                }
                return (double(match_count) / length);
            }
    };
}
//============================================================================================
//                    End of sequrnce_matcher module from GenomoKit
//============================================================================================
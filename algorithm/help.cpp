#include <iostream>
#include<vector>
#include "modules\GenomoKit\pathogen_test.hpp"

using namespace std;



int main(){
    // Loading and saving .gpn (genome pathogen network)
    // no need to write .gpn in load / save (added internaly by GenomoKit)
    
    GenomoKit::PathogenNetwork pn_load = GenomoKit::PathogenNetwork("test");
    pn_load.save("testSave");
 
    // NOTE: Following is demonstration of manual insertion usage

    GenomoKit::SequenceMatcher matcher = GenomoKit::SequenceMatcher();
    GenomoKit::PathogenNetwork pn = GenomoKit::PathogenNetwork();

    cout<<"Node 1 ------------>"<<endl;
    GenomoKit::Cds *p1c1 = new GenomoKit::Cds("C-p1-1","ATGGCGTACGTTAGCAGTACCGTAGCTAGGTAGCGTACGTTACG");
    GenomoKit::Cds *p1c2 = new GenomoKit::Cds("C-p1-2","GACTGATTGACGTA");
    vector<GenomoKit::Cds*> p1_cds;
    p1_cds.push_back(p1c1);
    p1_cds.push_back(p1c2);
    GenomoKit::Node p1 = GenomoKit::Node(0,"P1","GTTAGCACGTAGTTAGCTATGGCGTACGTTAGCAGTACCGTAGCTAGGTAGCGTACGTTACGACTGATTGACGTAACTAGGTCGATGCAGTAGC",p1_cds);
    pn.insert_node_on_dna_similarity_connections(&p1,0.5,3,0.90);

    cout<<"Node 2 ------------>"<<endl;
    GenomoKit::Cds *p2c1 = new GenomoKit::Cds("C-p2-1","ATGCGTACGTTAGCAGTACCGTAGCTAGGTAGCGTACGTTACG");
    GenomoKit::Cds *p2c2 = new GenomoKit::Cds("C-p2-2","GACTGATTGACGTG");
    vector<GenomoKit::Cds*> p2_cds;
    p2_cds.push_back(p2c1);
    p2_cds.push_back(p2c2);
    GenomoKit::Node p2 = GenomoKit::Node(1,"P2","AGCGTACGTTAGCAGTACCGTAGCTAGGTAGCGTACGTTACGACTGATTGACGTGAATTAGCTCGTAGGCTAGTGAGCTGACTCGTAA",p2_cds);
    pn.insert_node_on_dna_similarity_connections(&p2,0.5,3,0.90);

    cout<<"Node 3 ------------>"<<endl;
    GenomoKit::Cds *p3c1 = new GenomoKit::Cds("C-p3-1","ATGGTACTAGCTAGCTAGTACCGGATCGGATCGTACGCTAGTG");
    GenomoKit::Cds *p3c2 = new GenomoKit::Cds("C-p3-2","GACGCTAGTGCGT");
    vector<GenomoKit::Cds*> p3_cds;
    p3_cds.push_back(p3c1);
    p3_cds.push_back(p3c2);
    GenomoKit::Node p3 = GenomoKit::Node(2,"P3","TAGTCGGTACTAGCTAGCTAGTACCGGATCGGATCGTACGCTAGTGCTGAATCGTACGTGACGCTAGTGCGTCAGCTAGTTGA",p3_cds);
    pn.insert_node_on_dna_similarity_connections(&p3,0.5,3,0.90);

    
    cout<<"Node 4 ------------>"<<endl;
    GenomoKit::Cds *p4c1 = new GenomoKit::Cds("C-p4-1","ATGAGCTAGTACGTCGTTGACGCTTAGGATCGCTAGTGATCGT");
    GenomoKit::Cds *p4c2 = new GenomoKit::Cds("C-p4-2","TGCGATAGCTGAC");
    vector<GenomoKit::Cds*> p4_cds;
    p4_cds.push_back(p4c1);
    p4_cds.push_back(p4c2);
    GenomoKit::Node p4 = GenomoKit::Node(3,"P4","CTAGTACGAGCTAGTACGTCGTTGACGCTTAGGATCGCTAGTGATCGTTAACTAGTGCGATAGCTGACGCTAGCGTAGTTGGA",p4_cds);
    pn.insert_node_on_dna_similarity_connections(&p4,0.5,3,0.90);


    cout<<"Node 5 ------------>"<<endl;
    GenomoKit::Cds *p5c1 = new GenomoKit::Cds("C-p5-1","ATGGTACTAGCTAGCTAGTCCGGATCGGATCGTACGCTAGCGC");
    GenomoKit::Cds *p5c2 = new GenomoKit::Cds("C-p5-2","GAGTAGCTGATCG");
    vector<GenomoKit::Cds*> p5_cds;
    p5_cds.push_back(p5c1);
    p5_cds.push_back(p5c2);
    GenomoKit::Node p5 = GenomoKit::Node(4,"P5","ACGTAGCTGACTGGTACTAGCTAGCTAGTCCGGATCGGATCGTACGCTAGCGCTGAGTAGCTGATCGTAGGAATCGTAGCTGA",p5_cds);
    pn.insert_node_on_dna_similarity_connections(&p5,0.5,3,0.90);


    pn.save("testSaveManual");
    
    

    // Demonstration of working and handling of GenomoKit
    GenomoKit::ZoonosisTest zt = GenomoKit::ZoonosisTest(pn_load,"AGTACGCGTAGCTAGTACCGGATCGTACGTTAGCGGATCGGCTAGCTAGTGACTCGTAGGACGTTGACGCTAGTGCGTAGTTGGA");
    GenomoKit::ZoonosisTestResult result = zt.is_pathogen_zoonotic(0.5,3,0.9,1,5);
    
    cout<<"Paths --->"<<endl;
    for(auto path : result.paths){
        for( auto i: path){
            cout<<i<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
    cout<<endl;
    cout<<endl;

    cout<<"Matches --->"<<endl;
    for(auto cds: result.matched_cds){
        cout<<cds.first<< " - " <<cds.second<<endl;
    }

    return 0;
}

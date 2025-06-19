#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include"../GenomoKit/pathogen_network.hpp"
using namespace std;

int main() {
    
    //---------------------------------------------------------------------------
    string inputCSV = "zoonotic.csv";
    string outputGPN = "non_zoonotic";

    float kmer_threshold = 0.50;               // B/W [0,1]
    int kmer_size = 3;                         // >= 1
    float alignment_threshold = 0.70;          // B/W [0,1]
    //---------------------------------------------------------------------------
    

    GenomoKit::PathogenNetwork pn = *new GenomoKit::PathogenNetwork("testData"); 

    ifstream file(inputCSV); 
    if (!file.is_open()) {
        cerr << "Error opening file!" << endl;
        return 1;
    }

    string line;
    int nodeIDX = pn.graph.node_list.size();
    while (getline(file, line)) {
        stringstream ss(line);
        vector<string> fields;
        string field;

        while (getline(ss, field, ',')) {
            fields.push_back(field);
        }

        string UniqueId = fields[0];
        string Pathogen = fields[1];
        string Complete_Genome = fields[2];
        string Partial_Genome_1 = fields[3];
        string Partial_Genome_2 = fields[4];
        string Partial_Genome_3 = fields[5];
        string Animal_Infected = fields[6];
        int Threat_Level = stoi(fields[7]);
        string FoundWhere = fields[8];

        GenomoKit::Cds *c1 = new GenomoKit::Cds("C"+UniqueId+"-1",Partial_Genome_1);
        GenomoKit::Cds *c2 = new GenomoKit::Cds("C"+UniqueId+"-2",Partial_Genome_2);
        GenomoKit::Cds *c3 = new GenomoKit::Cds("C"+UniqueId+"-3",Partial_Genome_3);
        vector<GenomoKit::Cds*> cds;
        cds.push_back(c1);
        cds.push_back(c2);
        cds.push_back(c3);
        GenomoKit::Node node= GenomoKit::Node(nodeIDX,UniqueId,Complete_Genome,cds,Threat_Level);
        pn.insert_node_on_dna_similarity_connections(&node,kmer_threshold,kmer_size,alignment_threshold);
        cout<<"file : " <<nodeIDX+1<<" processed"<<endl;
        nodeIDX++;
    }

    pn.save(outputGPN);
    file.close();
    return 0;
}
#include "crow_all.h"
#include "algorithm/modules/GenomoKit/pathogen_test.hpp"
#include <fstream>
#include <sstream>

using namespace std;
using namespace GenomoKit;

// Safe parsers to avoid stoi/stof crashes
int safe_stoi(const crow::json::rvalue& val, int fallback = 0) {
    try {
        return std::stoi(std::string(val.s()));
    } catch (...) {
        return fallback;
    }
}

float safe_stof(const crow::json::rvalue& val, float fallback = 0.0f) {
    try {
        return std::stof(std::string(val.s()));
    } catch (...) {
        return fallback;
    }
}

int main() {
    crow::SimpleApp app;

    // Serve index.html at root "/"
   // Add this after your existing routes
    CROW_ROUTE(app, "/<path>")([](std::string path) {
    // Build the full file path
    std::string file_path = "webapp/frontend/" + path;
    
    // Check for specific file types
    if (path.find(".html") != std::string::npos ||
        path.find(".js") != std::string::npos ||
        path.find(".css") != std::string::npos) {
        
        std::ifstream file(file_path);
        if (!file.is_open()) {
            return crow::response(404, "File not found");
        }

        std::stringstream buffer;
        buffer << file.rdbuf();
        return crow::response(buffer.str());
    }
    return crow::response(404, "Invalid file type");
});

    // Serve static files like landing.js
    CROW_ROUTE(app, "/pages/landing_page/<string>")([](std::string filename) {
        std::ifstream file("webapp/frontend/pages/landing_page/" + filename);
        if (!file.is_open()) return crow::response(404, "File not found");

        std::stringstream buffer;
        buffer << file.rdbuf();
        return crow::response(buffer.str());
    });

    // Main API for analysis
   
CROW_ROUTE(app, "/calculate-similarity").methods(crow::HTTPMethod::Post)([](const crow::request& req) {
    auto body = crow::json::load(req.body);
    if (!body) return crow::response(400, "Invalid JSON");

    try {
        string dna1 = body["dna1"].s();
        string dna2 = body["dna2"].s();
        int kmerSize = safe_stoi(body["kmerSize"], 3);
        float threshold = safe_stof(body["threshold"], 0.7);
        string algorithm = body["algorithm"].s();

        // Create SequenceMatcher instance
        GenomoKit::SequenceMatcher matcher;
        
        // Calculate scores
        float kmerScore = matcher.find_kmer_similarity(dna1, dna2, kmerSize);
        float alignmentScore = matcher.local_alignment_similarity(dna1, dna2);
        
        // Determine overall similarity based on algorithm choice
        float similarity = (algorithm == "kmer") ? kmerScore : alignmentScore;
        
        // Build response
        crow::json::wvalue res;
        res["similarity"] = similarity;
        res["kmerScore"] = kmerScore;
        res["alignmentScore"] = alignmentScore;
        res["match"] = (similarity >= threshold);
        
        return crow::response{res};
        
    } catch (const std::exception& e) {
        return crow::response(500, string("Error: ") + e.what());
    }
});
    CROW_ROUTE(app, "/analyze").methods(crow::HTTPMethod::Post)([](const crow::request& req) {
        auto body = crow::json::load(req.body);
        if (!body) return crow::response(400, "Invalid JSON");

        try {
            // Defensive input extraction
            string dna = body["dna"].s();
            float kmer = safe_stof(body["kmer"], 0.5);
            float cds = safe_stof(body["cds"], 0.7);
            int stack = safe_stoi(body["stack"], 3);
            int threatlvl = safe_stoi(body["threatlvl"], 2);
            int maxpath = safe_stoi(body["maxpath"], 3);

                CROW_LOG_INFO << "Received parameters:";
                CROW_LOG_INFO << "dna: " << dna;
                CROW_LOG_INFO << "kmer: " << kmer;
                CROW_LOG_INFO << "cds: " << cds;
                CROW_LOG_INFO << "stack: " << stack;
                CROW_LOG_INFO << "threatlvl: " << threatlvl;
                CROW_LOG_INFO << "maxpath: " << maxpath;
        
                CROW_LOG_INFO << "Loading pathogen network from data/opppppp3";
                PathogenNetwork network("data/opppppp3");
                CROW_LOG_INFO << "Loaded pathogen network with " 
              << network.graph.node_list.size() << " nodes";
            ZoonosisTest test(network, dna);
            ZoonosisTestResult result = test.is_pathogen_zoonotic(kmer, stack, cds, threatlvl, maxpath);
            
            // JSON response setup
            crow::json::wvalue res_json;
            std::stringstream summary;

            // CDS matches
            if (result.matched_cds.empty()) {
                summary << "No CDS match found.\n";
            } else {
                summary << "Matched CDS:\n";
                for (const auto& cds_match : result.matched_cds) {
                    summary << "• " << cds_match.first << " → Score: " << cds_match.second << "\n";
                }
            }

            // Paths
            if (!result.paths.empty()) {
                summary << "\nTop Paths to High Threat Pathogens:\n";
                int count = 1;
                for (const auto& path : result.paths) {
                    summary << count++ << ". ";
                    for (size_t i = 0; i < path.size(); ++i) {
                        summary << path[i];
                        if (i != path.size() - 1) summary << " → ";
                    }
                    summary << "\n";
                }
            } else {
                summary << "No high threat paths found.\n";
            }

            // Add parsed paths to JSON
            crow::json::wvalue::list pathList;
            for (const auto& path : result.paths) {
                crow::json::wvalue::list p;
                for (const auto& node : path) p.push_back(node);
                pathList.push_back(std::move(p));
            }

            res_json["summary"] = summary.str();
            res_json["paths"] = std::move(pathList);
            return crow::response{res_json};

        } catch (const std::exception& e) {
            CROW_LOG_ERROR << "Exception: " << e.what();
            return crow::response(500, std::string("Server Error: ") + e.what());
        }
    });

    app.port(18080).multithreaded().run();
}

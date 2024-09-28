#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include "TMVA/DecisionTree.h"
#include "TMVA/MethodBDT.h"
#include "TMVA/Reader.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <fstream>

// Function to convert GBRForest to XML
void GBRForestToXML(const GBRForest* gbrForest, const std::string& outputFileName) {
    std::ofstream xmlFile(outputFileName);
    
    if (!xmlFile.is_open()) {
        std::cerr << "Error: Could not open file " << outputFileName << " for writing." << std::endl;
        return;
    }

    // Begin the XML structure
    xmlFile << "<GBRForest>" << std::endl;

    // Loop over decision trees in the forest
    const std::vector<TMVA::DecisionTree*>& trees = gbrForest->Trees();
    for (size_t i = 0; i < trees.size(); ++i) {
        TMVA::DecisionTree* tree = trees[i];
        xmlFile << "  <DecisionTree id=\"" << i << "\">" << std::endl;

        // Traverse the decision tree and serialize its nodes
        tree->WriteXML(xmlFile, "    ");  // TMVA has a built-in XML writer for decision trees

        xmlFile << "  </DecisionTree>" << std::endl;
    }

    // End the XML structure
    xmlFile << "</GBRForest>" << std::endl;

    xmlFile.close();
    std::cout << "GBRForest successfully exported to " << outputFileName << std::endl;
}

int main(int argc, char** argv) {
    // Load the GBRForest from a file (assuming it is stored in a ROOT file)
    TFile* file = TFile::Open("../../../VertexCompositeAnalysis/VertexCompositeAnalyzer/data/GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt46MassPeak_NoPtYPtErrNHitDLAngle2D_v3.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open GBRForest file." << std::endl;
        return 1;
    }

    GBRForest* gbrForest = dynamic_cast<GBRForest*>(file->Get("GBRForest"));
    if (!gbrForest) {
        std::cerr << "Error: GBRForest not found in file." << std::endl;
        return 1;
    }

    // Export to XML
    GBRForestToXML(gbrForest, "gbrForest.xml");

    // Clean up
    file->Close();
    delete file;

    return 0;
}

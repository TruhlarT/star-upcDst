#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TObject.h>
#include <iostream>
#include <string>
#include <vector>

void mergeROOTFiles() {
    // Path to the system directory where the input ROOT files are stored
    std::string systemDir = "";

    // List of filenames (stored in the nHitsStudy system directory)
    std::vector<std::string> filenames;
    filenames.push_back("nDedxLoose.root");
    filenames.push_back("nFitLoose.root");
    filenames.push_back("nDedxTight.root");
    filenames.push_back("nFitTight.root");
    filenames.push_back("nLoose.root");
    filenames.push_back("nTight.root");

    // Create the output ROOT file named "nHitsStudy.root"
    TFile* outputFile = new TFile("nHitsStudy.root", "RECREATE");

    // Create the main "Embedding" directory in the output ROOT file
    TDirectory* embeddingDir = outputFile->mkdir("Embedding");
    embeddingDir->cd(); // Move into the "Embedding" directory in the ROOT file

    // Loop over input files
    for (size_t i = 0; i < filenames.size(); ++i) {
        // Build the full path to the input file in the system directory
        std::string filePath = systemDir + filenames[i];

        // Open the input file
        TFile* inputFile = TFile::Open(filePath.c_str());
        if (!inputFile || inputFile->IsZombie()) {
            std::cerr << "Error opening file: " << filePath << std::endl;
            continue;
        }

        // Create a subfolder inside the "Embedding" directory with the name of the input file (without ".root")
        std::string subfolderName = filenames[i].substr(0, filenames[i].find(".root")); // Remove ".root"
        TDirectory* subfolderDir = embeddingDir->mkdir(subfolderName.c_str());
        subfolderDir->cd(); // Change to the subfolder inside the ROOT file

        // Get the "Embedding" directory from the input file
        TDirectory* embeddingInputDir = (TDirectory*)inputFile->Get("Embedding");
        if (!embeddingInputDir) {
            std::cerr << "Error: 'Embedding' folder not found in file: " << filePath << std::endl;
            inputFile->Close();
            continue;
        }

        // Loop over all objects in the "Embedding" directory of the input file
        TIter nextKey(embeddingInputDir->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)nextKey())) {
            // Read the object from the input file
            TObject* obj = key->ReadObj();

            // Check if the object is a directory (TDirectory) or another object
            if (obj->InheritsFrom(TDirectory::Class())) {
                // Manually copy all contents of the subdirectory
                TDirectory* subdir = (TDirectory*)obj;
                TDirectory* newSubDir = subfolderDir->mkdir(subdir->GetName()); // Create a matching directory in the output
                newSubDir->cd(); // Change to the new subdirectory

                // Loop over objects in this subdirectory
                TIter subDirKey(subdir->GetListOfKeys());
                TKey* subkey;
                while ((subkey = (TKey*)subDirKey())) {
                    // Read each object and write it to the new directory
                    TObject* subObj = subkey->ReadObj();
                    // Check if the object is a directory (TDirectory) or another object
                    if (subObj->InheritsFrom(TDirectory::Class())) {
                        // Manually copy all contents of the subdirectory
                        TDirectory* subsubdir = (TDirectory*)subObj;
                        TDirectory* newsubSubDir = newSubDir->mkdir(subsubdir->GetName()); // Create a matching directory in the output
                        newsubSubDir->cd(); // Change to the new subdirectory

                        // Loop over objects in this subdirectory
                        TIter subsubDirKey(subsubdir->GetListOfKeys());
                        TKey* subsubkey;
                        while ((subsubkey = (TKey*)subsubDirKey())) {
                            // Read each object and write it to the new directory
                            TObject* subsubObj = subsubkey->ReadObj();
                            newsubSubDir->cd(); // Ensure we're writing in the correct subdirectory
                            cout<<"Writting object "<<subObj->GetName()<<endl;
                            subsubObj->Write(); // Write the object into the output file
                            delete subsubObj; // Clean up
                        }
                    } else {
                        // Otherwise, it's a regular object, so we write it directly
                        newSubDir->cd(); // Make sure we're in the correct output subfolder
                        cout<<"Writting object "<<obj->GetName()<<endl;
                        obj->Write(); // Write the object to the current directory
                    }

                    newSubDir->cd(); // Ensure we're writing in the correct subdirectory
                    cout<<"Writting object "<<subObj->GetName()<<endl;
                    subObj->Write(); // Write the object into the output file
                    delete subObj; // Clean up
                }
            } else {
                // Otherwise, it's a regular object, so we write it directly
                subfolderDir->cd(); // Make sure we're in the correct output subfolder
                cout<<"Writting object "<<obj->GetName()<<endl;
                obj->Write(); // Write the object to the current directory
            }

            delete obj; // Clean up the object after writing it
        }

        // Close the input file
        inputFile->Close();
    }

    // Close the output file
    outputFile->Close();
    std::cout << "Merging completed!" << std::endl;
}
    
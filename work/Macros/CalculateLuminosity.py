# macro to study luminosity in run 17 for CPT2noBBCL
# assume a root file per run

#run with: python CalculateLuminosity.py
# output: all in [pb-1]
# Total number of events in root file: 587164368.0
# Total integrated luminosity: 124.292292917
# Total integrated luminosity from run by run: 120.43179982
# Total corrected integrated luminosity from run by run: 29.4323986326
# Total minimum corrected integrated luminosity from run by run: 28.2280806344
# Total maximum corrected integrated luminosity from run by run: 30.6367166308


import os
import ROOT
import numpy as np

lumi_path = "/gpfs01/star/pwg_tasks/upc02/CEP_input_files/lists/lumForCPT2noBBCL.list"  
file_path = "/gpfs01/star/pwg_tasks/upc02/CEP_input_files/lists/lumPerRun.list" 
#From the table: https://www.star.bnl.gov/protected/common/common2017/trigger2017/lumipp500GeV/lum_pertriggerid_pp2017_500GeV.html
#Get total integrated lumi for RP_CPT2noBBCL = 142.006 Lum [pb-1]
total_lumi_from_Jaime = 142.006
total_events_from_Jaime = 670.845e+6

path_to_rootfiles = "/gpfs01/star/pwg/truhlar/Run17_P20ic/luminosity/" 

def read_file_and_store_values(file_path):
    data_map = {}

    with open(file_path, 'r') as file:
        for line in file:
            values = line.split()
            run = int(float(values[0]))
            lumi = float(values[4]) # [pb-1]
            events = int(float(values[11]))

            time0 = float(values[1])
            time1 = float(values[2])
            prescale = float(values[5])
            livetime = float(values[6])
            timeDiff = time1-time0
            inst_lumi = 0
            if livetime*timeDiff != 0:
                inst_lumi = prescale*lumi*1000000/(livetime*timeDiff); # convert from [pb-1] pico to [mub] mikro

            data_map[run] = (lumi, events, inst_lumi)

    return data_map


lumi_map = read_file_and_store_values(lumi_path)
#print(lumi_map)

total_integrated_luminosity = 0.0 # lumi from Jaime * events in root / events from Jaime
total_int_lum_perrun = 0.0 # run by run method
total_cor_int_lum_perrun = 0.0 # corrected by veto eff
total_cor_int_lum_perrun_min = 0.0 # corrected by veto eff
total_cor_int_lum_perrun_max = 0.0 # corrected by veto eff


# List all files in the root directory
root_files = os.listdir(path_to_rootfiles)
for file in root_files:
    #print(file)
    if file.endswith(".root"):
        root_file = ROOT.TFile.Open(path_to_rootfiles+file)

        # Access the histogram named "AnalysisFlow"
        histogram = root_file.Get("AnalysisFlow")

        # Get the content of the second bin
        events_in_root = histogram.GetBinContent(2)  # Indexing starts from 1, so the second bin [CPT trigger] is at index 2

        # Navigate to the TofQA folder
        # List all histograms within the TofQA folder
        tofqa_folder = root_file.Get("TofQA")
        histogram_names = []
        for key in tofqa_folder.GetListOfKeys():
            obj = key.ReadObj()
            if isinstance(obj, ROOT.TH1):
                _, tmp = obj.GetName().split("_", 1)
                histogram_names.append(tmp)       
        #histogram_names = [int(name) for name in histogram_names]
        run = int(histogram_names[0])        
        #print("Events in run '{}': {}".format(run, events_in_root))
        # Check if the run number exists in the map
        #print run
        if run in lumi_map:
            lumi, events, inst_lumi = lumi_map[run]
            #print("In run {} there are {}/{} events with {} lumi".format(run,events_in_root,events,lumi))
            total_integrated_luminosity += events_in_root
            lum_perrun = lumi * (events_in_root/events) # [pb-1]
            with open(file_path, 'a') as file:  # 'a' mode appends to the file, use 'w' to overwrite existing content
                file.write("{} {}\n".format(run, lum_perrun))
            total_int_lum_perrun += lum_perrun # [pb-1]
            veto_correction = 1.42 * np.exp(-0.01488 * inst_lumi);
            total_cor_int_lum_perrun += lumi * (events_in_root/events )*veto_correction # [pb-1]
            total_cor_int_lum_perrun_min += lumi * (events_in_root/events )*(veto_correction-0.01) # [pb-1]
            total_cor_int_lum_perrun_max += lumi * (events_in_root/events )*(veto_correction+0.01) # [pb-1]

print("Total number of events in root file: {}".format(total_integrated_luminosity))
total_integrated_luminosity = total_lumi_from_Jaime * (total_integrated_luminosity / total_events_from_Jaime ) # [pb-1]

print("Total integrated luminosity: {} [pb-1]".format(total_integrated_luminosity))
print("Total integrated luminosity from run by run: {} [pb-1]".format(total_int_lum_perrun))
print("Total corrected integrated luminosity from run by run: {} [pb-1]".format(total_cor_int_lum_perrun))
print("Total minimum corrected integrated luminosity from run by run: {} [pb-1]".format(total_cor_int_lum_perrun_min))
print("Total maximum corrected integrated luminosity from run by run: {} [pb-1]".format(total_cor_int_lum_perrun_max))


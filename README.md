# CEP Analysis Repository

This repository contains code for performing Central Exclusive Production (CEP) analysis. The code is based on [*star-upcDst*](https://github.com/adamjaro/star-upcDst), a new framework mainly developed by [Jarda Adam](https://github.com/adamjaro) and Tomas Truhlar to simplify analysis that are related to forward and UPC physics. 

The code consists of two parts:
- [*star-upcDst*](https://github.com/adamjaro/star-upcDst): Use MuDst and PicoDst files of one's interest as input, run so-called upcDstMaker to produce a upcDst file with information that are needed for one's analysis only. 
- work directory: Use the result of part 1, the upcDst, to build your personal analysis tree or histograms

## [*star-upcDst*](https://github.com/adamjaro/star-upcDst):

- Make a clean area on RACF and checkout the repository:

<pre><code> git clone https://github.com/TruhlarT/star-upcDst.git </pre></code>

- Go to the main directory star-upcDst:

<pre><code> cd star-upcDst </pre></code>

- Setup the StRoot and build ( already include compiling ) by doing:

<pre><code> ./build-maker.sh </pre></code>

- To perform analysis on upcDst files, one needs to compile the package and setup the link to the library. This can be simply done in a few steps:

<pre><code> mkdir build </pre></code>
<pre><code> cd build </code></pre>
<pre><code> cmake ../ </code></pre>
<pre><code> make </code></pre>

- For the upcDst production, check the original [*star-upcDst*](https://github.com/adamjaro/star-upcDst) repository. But be aware that some changes were made. For example an option to merge old upcDst with PicoDst (global tracks) was added. Also some scripts were modified.  

## work directory:
- This directory contains the code to fully reproduce CEP analysis. 
- The analysis is set by two files "include/AnaConfig.h" and "include/RunDef.h". The first one contains all cuts, paths etc. The second one contains list of booleans determining what part of analysis should be run.
- Firstly, one has to go to the work dir and compile the code:

<pre><code> cd ../work </code></pre>
<pre><code> cmake .</code></pre>
<pre><code> make </code></pre>

### Run the code locally:
- To run the code locally on a single file: 
<pre><code> AnalysisManager /star/data01/pwg_tasks/upc03/pp17/ExtendedUpcDst/18100001.root </code></pre>
- To run the code locally on a list of files: 
<pre><code> AnalysisManager /gpfs01/star/pwg_tasks/upc02/CEP_input_files/lists/extendedUpc.list </code></pre>
- Or to run on the i-th file from list of files: 
<pre><code> AnalysisManager /gpfs01/star/pwg_tasks/upc02/CEP_input_files/lists/extendedUpc.list 111 </code></pre>
- The output will be stored in the work dir as "AnalysisOutput.root"

### Submitting jobs on condor:
- Before the jobs can be submitted, one has to set the output directory. The directory has to be set in "SubmitPlugin.py":
<pre><code> top = "/gpfs01/star/pwg/truhlar/Run17_P20ic" </pre></code>
- And in "PrintStat.py" and "RunMerge.py":
<pre><code> defaultdir = "/gpfs01/star/pwg/truhlar/Run17_P20ic" </pre></code>

- Also one has to create sched directory for the scheduler output
<pre><code> mdkir sched </pre></code>

- To submit the jobs on condor, run:
<pre><code> SubmitPlugin.py qaTest /gpfs01/star/pwg_tasks/upc02/CEP_input_files/lists/extendedUpc.list </pre></code>
- If the list is not given, then the default one from SubmitPlugin.py is used.
- The "qaTest" serves as a tag for the given run of the analysis. Use uniq tag for each analysis.

- One can check the job status using "PrintStat.py":
<pre><code> PrintStat.py qaTest </pre></code>
- If tag is not given, then the last submitted jobs are checked.
- If there are any failed jobs, one can resubmit them by:
<pre><code> PrintStat.py r qaTest </pre></code>


- Each jobs produce an individual output. One can merge the output files by:
<pre><code> RunMerge.py qaTest </pre></code>
- Again, if tag is not given, then the last submitted jobs are merged.

### Create final plots:
- Plot manager was developed to process the output of the main analysis describe above and produce the plots to the thesis and paper.
- To create the plots for the analysis with tag = "qaTest", run:
<pre><code> runPlotManager.py qaTest </pre></code>
- Again, if tag is not given, then the last submitted analysis is processed.
- If one wants to process a local analysis, run:
<pre><code> runPlotManager.py l </pre></code>
- The root file called "FinalPlots.root" is produced in work dir containing all final plots, except the dEdx plot. 
- If the runPlotManager is run over small statics (e.g. a single run), many errors may occur. Ignore them. They originate from failed fits.
- To produce the dEdx plot, one has to run with root4star:
<pre><code> root4star -l </pre></code>
<pre><code> .L bichselG10.C </pre></code>
<pre><code> .x bichselG10.C("qaTest") </pre></code>
- But firstly, the path to the file should be modified:
<pre><code> TString inputFileLocation = "/gpfs01/star/pwg/truhlar/Run17_P20ic/" + input + "/merged/StRP_production_0000.root"; </pre></code>

## Contact:
Tomas Truhlar: <Tomas.Truhlar@fjfi.cvut.cz>









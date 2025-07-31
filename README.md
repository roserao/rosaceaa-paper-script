# rosaceaa-paper-script

Scripts and public datasets used to perform data analysis and generate plots for the paper "Rosace-AA: Enhancing Interpretation of Deep Mutational Scanning Data with Amino Acid Substitution and Position-Specific Insights".


## Folder: DMS_rosace2_realdata

A Snakemake pipeline for running Rosace-AA on public datasets and summarize the results. There are four sets of computational analysis within this folder.

### Subfolder: results

Running three models of Rosace-AA on the cleaned datasets of the data folder ("data/*.rda"), downloaded and processed from the ProteinGym Database, and then compute the variance decomposition results using the contained R script.

#### References
__ProteinGym__:  Pascal Notin, Aaron Kollasch, Daniel Ritter, and et al. Proteingym: Large-scale benchmarks for protein fitness prediction and design. Advances in Neural Information Processing Systems, volume 36, pages 64331–64379. Curran Associates, Inc., 2023.

__CAR11__: Iana Meitlis, Eric J. Allenspach, Bradly M. Bauman, and et al. Multiplexed functional assessment of genetic variants in card11. The American Journal of Human Genetics, 107(6):1029–1043, Dec 2020.

__CAS9__: Jeffrey M. Spencer and Xiaoliu Zhang. Deep mutational scanning of s. pyogenes cas9 reveals important functional domains. Scientific Reports, 7(1):16836, Dec 2017.

__CD19__: Justin R. Klesmith, Lan Wu, Roy R. Lobb, and et al. Fine epitope mapping of the cd19 extracellular domain promotes design. Biochemistry, 58(48):4869–4881, Dec 2019.

__RNC__: Ryan Weeks and Marc Ostermeier. Fitness and functional
landscapes of the e. coli rnase iii gene rnc. Molecular
Biology and Evolution, 40(3):msad047, 02 2023.

__SPG1__: C. Anders Olson, Nicholas C. Wu, and Ren Sun. A comprehensive biophysical description of pairwise epistasis throughout an entire protein domain. Current Biology, 24(22):2643–2651, 2014.


### Subfolder: results_domainase

Running three models of Rosace-AA on the Human Domainome dataset under the "data/domainase-flat" folder. R scripts contain the plotting functions for Figure 2.

#### References

__Domainome__: A Beltran, X Jiang, Y Shen, and et al. Site saturation mutagenesis of 500 human protein domains reveals the contribution of protein destabilization to genetic disease. bioRxiv, 2024.

### Subfolder: results_oct

Running three models of Rosace-AA on the OCT1 and its random or guided downsampling dataset. R scripts contain the downsampling procedures and plotting functions for Figure 3.

#### References

__OCT1__: SW Yee, CB Macdonald, D Mitrovic, and et al. The full spectrum of slc22 oct1 mutations illuminates the bridge between drug transporter biophysics and pharmacogenomics. Molecular Cell, 84(10):1932–1947.e10, May 2024.

### Subfolder: results_met

Running three models of Rosace-AA on the MET kinase domain DMS data under 12 inhibitor and 2 genetic backgrounds. R scripts contain the plotting functions for Figure 4.

#### References

__MET__: GO Estevam, EM Linossi, J Rao, and et al. Mapping kinase domain resistance mechanisms for the met receptor tyrosine kinase via deep mutational scanning. bioRxiv, 2024.

## Folder: DMS_rosace2_simulation

A Snakemake pipeline for running Rosette-AA simulation to validate that the model and software run correctly. The simulation details are specified in the supplementary note of the paper.

We did the simulation with four datasets: 
- Data 1: all variants are classified as neutral (e.g. MET kinase domain DMS without IL3)
- Data 2: only neutral and positive groups are identified (e.g. OCT1 drug cytotoxicity screen)
- Data 3: only neutral and negative groups are identified (e.g. MET kinase domain DMS with IL3)
- Data 4: all three groups can be identified (unpublished data)

For data 2-4, we simulate the experiments with 3 replicates and 3 selection rounds (4 time points), and then apply the 5 scenarios (all random, position-only, 3 types of position and substitution). Then, we compute the fdr and sensitivity and compare it with the Rosace model.

#### References

__OCT1__: SW Yee, CB Macdonald, D Mitrovic, and et al. The full spectrum of slc22 oct1 mutations illuminates the bridge between drug transporter biophysics and pharmacogenomics. Molecular Cell, 84(10):1932–1947.e10, May 2024.

__MET__: Estevam GO, Linossi EM, Macdonald CB, Espinoza CA, Michaud JM, Coyote-Maestas W, et al. Conserved regulatory motifs in the juxtamembrane domain and kinase N-lobe revealed through deep mutational scanning of the MET receptor tyrosine kinase domain. eLife. 2023.




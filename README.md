# DISSERTATION REPOSITORY #

#
# Title: "From perception to action: Characterization of motor function in autism through EEG analysis"
# Institution: Universidade de Aveiro (in collaboration with the Institute of Systems and Robotics – University of Coimbra)
# Author: Filipa Andrade da Câmara
# Year: 2025


This repository has the main materials developed and used throughout the dissertation. 
It includes the original ISR article that served as the data source, as well as MATLAB scripts used 
for the analytical and statistical pipelines described in the study

Each file corresponds to a specific part of the methodology detailed in Chapter 5 and Section 5.3 
of the dissertation.


# 1. ISR_article.pdf  
   - Original publication from the Institute of Systems and Robotics (ISR-UC).  
   - Serves as the empirical basis for the analyses in Chapter 4: “Dataset and Experimental Paradigm (ISR)”.

# 2. Advanced_visualization_pipeline.m  
   - MATLAB script implementing the Advanced Visualization Pipeline (Chapter 5, Subsection 5.2.5).  
   - Generates group-level time–frequency maps and scalp topoplots from `.dattimef` files 
     precomputed in EEGLAB.  
   - Integrates both temporal-spectral and spatial visualization of ERSPs.

# 3. Shapiro_Wilk_test.m  
   - MATLAB script for the normality tests described in Section 5.3: “Statistical Analyses”.  
   - Applies the Shapiro–Wilk test (α = 0.05) to mean power values (8–10 Hz) from central channels (C3, Cz, C4).  
   - Outputs QQ-plots, histograms, and W/p-values supporting the decision between parametric and non-parametric tests.

# 4. tftopo.m  
   - Modified EEGLAB function used by the visualization pipeline.  
   - Separates the channels used for time–frequency computation (selchans) and those used for scalp mapping (mapchans).  
   - Enables combined visualization of central-ROI ERSPs with full-head topographies.

# 5. Appendix
   - Additional figures and plots respectively mentioned throughout the text
   - Examples of '.dattimef' structures
   
2025 Filipa Andrade da Câmara – for academic use only.

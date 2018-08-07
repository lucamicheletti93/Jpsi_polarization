void Compilemyclass(){
  gROOT -> ProcessLineSync(".x /home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/Binning/Binning.cxx+") ;
  gROOT -> ProcessLineSync(".L comparison_weighted_nonweighted_distrib.C++") ;
}

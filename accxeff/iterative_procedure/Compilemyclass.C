void Compilemyclass(){
  gROOT -> ProcessLineSync(".x /home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/Binning/Binning.cxx+") ;
  gROOT -> ProcessLineSync(".L accxeff_iterator.C++") ;
}

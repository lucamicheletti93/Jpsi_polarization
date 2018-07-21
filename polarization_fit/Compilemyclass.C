void Compilemyclass(){
  gROOT -> ProcessLineSync(".x /home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/Binning/Binning.cxx+") ;
  gROOT -> ProcessLineSync(".L helicity_polarization_fit_2D.C++") ;
}

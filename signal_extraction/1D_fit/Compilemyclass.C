void Compilemyclass(){
  gROOT -> ProcessLineSync(".x /home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/Binning/Binning.cxx+") ;
  //gROOT -> ProcessLineSync(".L binned_mass_fit.C++") ;
}

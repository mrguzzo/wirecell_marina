struct PFevalInfo{
  bool flag_NCDelta;
  bool flag_showerMomentum;
  
  Int_t run;                                  // run
  Int_t subrun;                               // subrun
  Int_t event;                                // event 
  Int_t neutrino_type;
  Float_t reco_nuvtxX;                        // reco neutrino vertex X
  Float_t reco_nuvtxY;                        // reco neutrino vertex Y
  Float_t reco_nuvtxZ;                        // reco neutrino vertex Z
  Float_t reco_showervtxX;                    // reco shower vertex X
  Float_t reco_showervtxY;                    // reco shower vertex Y
  Float_t reco_showervtxZ;                    // reco shower vertex Z
  Float_t reco_showerKE;                      // reco shower kinetic energy [MeV]
  Float_t reco_muonvtxX;                      // reco muon vertex X
  Float_t reco_muonvtxY;                      // reco muon vertex Y
  Float_t reco_muonvtxZ;                      // reco muon vertex Z
  Float_t reco_muonMomentum[4];               // reco muon momenta 0=X, 1=Y, 2=Z, 3=time
  Float_t reco_showerMomentum[4];             // reco shower momenta 0=X, 1=Y, 2=Z, 3=time
  
  
  Float_t nuvtx_diff;                         // truth-reco neutrino vertex distance
  Float_t showervtx_diff;                     // truth-reco shower vertex distance
  Float_t muonvtx_diff;                       // truth-reco muon vertex distance
  Float_t truth_showerMomentum[4];            // truth shower momenta 0=X, 1=Y, 2=Z, 4=time
  Float_t truth_corr_nuvtxX;                  // truth sce neutrino vertex X
  Float_t truth_corr_nuvtxY;                  // truth sce neutrino vertex Y
  Float_t truth_corr_nuvtxZ;                  // truth sce neutrino vertex Z
  Float_t truth_corr_showervtxX;              // truth sce shower vertex X
  Float_t truth_corr_showervtxY;              // truth sce shower vertex Y
  Float_t truth_corr_showervtxZ;              // truth sce shower vertex Z
  Float_t truth_showerKE;                     // truth shower kinetic energy
  Float_t truth_corr_muonvtxX;                // truth sce muon vertex X
  Float_t truth_corr_muonvtxY;                // truth sce muon vertex Y
  Float_t truth_corr_muonvtxZ;                // truth sce muon vertex Z
  Float_t truth_muonvtxX;                     // truth muon vertex X
  Float_t truth_muonvtxY;                     // truth muon vertex Y
  Float_t truth_muonvtxZ;                     // truth muon vertex Z
  Float_t truth_muonendX;                     // truth muon end X
  Float_t truth_muonendY;                     // truth muon end Y
  Float_t truth_muonendZ;                     // truth muon end Z
  Float_t truth_muonMomentum[4];              // truth muon momenta 0=X, 1=Y, 2=Z, 3=time
  Float_t truth_nuEnergy;                     // truth neutrino energy [MeV]
  Float_t truth_energyInside;                 
  Float_t truth_electronInside;               
  Int_t truth_nuPdg;                          // truth neutrino pdg
  Bool_t truth_isCC;                          // truth = CC, false = NC
  Float_t truth_vtxX;                         // truth neutrino vertex X
  Float_t truth_vtxY;                         // truth neutrino vertex Y
  Float_t truth_vtxZ;                         // truth neutrino vertex Z
  Float_t truth_nuTime;                       // truth neutrino interaction time [units=??]
  Int_t truth_nuIntType;

  //
  Int_t truth_NprimPio;
  Float_t truth_pio_energy_1;
  Float_t truth_pio_energy_2;
  Float_t truth_pio_angle;
  Int_t truth_NCDelta;
  Float_t reco_protonMomentum[4];
  
}; // PFevalInfo

void loadPFevalInfo(TChain* tree0, PFevalInfo& tagger_info, int flag=1) {
  tagger_info.flag_NCDelta = false;
  tagger_info.flag_showerMomentum = false;
  
  tree0->SetBranchAddress("run", &tagger_info.run);
  tree0->SetBranchAddress("subrun", &tagger_info.subrun);
  tree0->SetBranchAddress("event", &tagger_info.event);
  tree0->SetBranchAddress("neutrino_type", &tagger_info.neutrino_type);
  tree0->SetBranchAddress("reco_nuvtxX", &tagger_info.reco_nuvtxX);
  tree0->SetBranchAddress("reco_nuvtxY", &tagger_info.reco_nuvtxY);
  tree0->SetBranchAddress("reco_nuvtxZ", &tagger_info.reco_nuvtxZ);
  tree0->SetBranchAddress("reco_showervtxX", &tagger_info.reco_showervtxX);
  tree0->SetBranchAddress("reco_showervtxY", &tagger_info.reco_showervtxY);
  tree0->SetBranchAddress("reco_showervtxZ", &tagger_info.reco_showervtxZ);
  tree0->SetBranchAddress("reco_showerKE", &tagger_info.reco_showerKE);
  tree0->SetBranchAddress("reco_muonvtxX", &tagger_info.reco_muonvtxX);
  tree0->SetBranchAddress("reco_muonvtxY", &tagger_info.reco_muonvtxY);
  tree0->SetBranchAddress("reco_muonvtxZ", &tagger_info.reco_muonvtxZ);
  tree0->SetBranchAddress("reco_muonMomentum", &tagger_info.reco_muonMomentum[0]);

  if (tree0->GetBranch("reco_showerMomentum")){
    tagger_info.flag_showerMomentum = true;
    tree0->SetBranchAddress("reco_showerMomentum",&tagger_info.reco_showerMomentum[0]);
    if (flag==1){
      tree0->SetBranchAddress("truth_showerMomentum",&tagger_info.truth_showerMomentum[0]);
    }
  }
  
  if (flag==1){
    tree0->SetBranchAddress("nuvtx_diff", &tagger_info.nuvtx_diff);
    tree0->SetBranchAddress("showervtx_diff", &tagger_info.showervtx_diff);
    tree0->SetBranchAddress("muonvtx_diff", &tagger_info.muonvtx_diff);
    tree0->SetBranchAddress("truth_corr_nuvtxX", &tagger_info.truth_corr_nuvtxX);
    tree0->SetBranchAddress("truth_corr_nuvtxY", &tagger_info.truth_corr_nuvtxY);
    tree0->SetBranchAddress("truth_corr_nuvtxZ", &tagger_info.truth_corr_nuvtxZ);
    tree0->SetBranchAddress("truth_corr_showervtxX", &tagger_info.truth_corr_showervtxX);
    tree0->SetBranchAddress("truth_corr_showervtxY", &tagger_info.truth_corr_showervtxY);
    tree0->SetBranchAddress("truth_corr_showervtxZ", &tagger_info.truth_corr_showervtxZ);
    tree0->SetBranchAddress("truth_showerKE", &tagger_info.truth_showerKE);
    tree0->SetBranchAddress("truth_corr_muonvtxX", &tagger_info.truth_corr_muonvtxX);
    tree0->SetBranchAddress("truth_corr_muonvtxY", &tagger_info.truth_corr_muonvtxY);
    tree0->SetBranchAddress("truth_corr_muonvtxZ", &tagger_info.truth_corr_muonvtxZ);
    tree0->SetBranchAddress("truth_muonvtxX", &tagger_info.truth_muonvtxX);
    tree0->SetBranchAddress("truth_muonvtxY", &tagger_info.truth_muonvtxY);
    tree0->SetBranchAddress("truth_muonvtxZ", &tagger_info.truth_muonvtxZ);
    tree0->SetBranchAddress("truth_muonendX", &tagger_info.truth_muonendX);
    tree0->SetBranchAddress("truth_muonendY", &tagger_info.truth_muonendY);
    tree0->SetBranchAddress("truth_muonendZ", &tagger_info.truth_muonendZ);
    tree0->SetBranchAddress("truth_muonMomentum", &tagger_info.truth_muonMomentum[0]);
    tree0->SetBranchAddress("truth_nuEnergy", &tagger_info.truth_nuEnergy);
    tree0->SetBranchAddress("truth_energyInside", &tagger_info.truth_energyInside);
    tree0->SetBranchAddress("truth_electronInside", &tagger_info.truth_electronInside);
    tree0->SetBranchAddress("truth_nuPdg", &tagger_info.truth_nuPdg);
    tree0->SetBranchAddress("truth_isCC", &tagger_info.truth_isCC);
    tree0->SetBranchAddress("truth_vtxX", &tagger_info.truth_vtxX);
    tree0->SetBranchAddress("truth_vtxY", &tagger_info.truth_vtxY);
    tree0->SetBranchAddress("truth_vtxZ", &tagger_info.truth_vtxZ);
    tree0->SetBranchAddress("truth_nuTime", &tagger_info.truth_nuTime);
    tree0->SetBranchAddress("truth_nuIntType", &tagger_info.truth_nuIntType);

    if (tree0->GetBranch("truth_NCDelta")){
      tagger_info.flag_NCDelta = true;
      tree0->SetBranchAddress("truth_NCDelta",&tagger_info.truth_NCDelta);
      tree0->SetBranchAddress("truth_NprimPio",&tagger_info.truth_NprimPio);
      tree0->SetBranchAddress("truth_pio_energy_1",&tagger_info.truth_pio_energy_1);
      tree0->SetBranchAddress("truth_pio_energy_2",&tagger_info.truth_pio_energy_2);
      tree0->SetBranchAddress("truth_pio_angle",&tagger_info.truth_pio_angle);
      tree0->SetBranchAddress("reco_protonMomentum",&tagger_info.reco_protonMomentum[0]);
    }
    
  }
} // loadPFevalInfo


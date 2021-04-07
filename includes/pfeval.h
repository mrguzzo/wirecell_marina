struct PFevalInfo{
  bool flag_NCDelta;
  bool flag_showerMomentum;
  
  Int_t run;
  Int_t subrun;
  Int_t event;
  Int_t neutrino_type;
  Float_t reco_nuvtxX;
  Float_t reco_nuvtxY;
  Float_t reco_nuvtxZ;
  Float_t reco_showervtxX;
  Float_t reco_showervtxY;
  Float_t reco_showervtxZ;
  Float_t reco_showerKE;
  Float_t reco_muonvtxX;
  Float_t reco_muonvtxY;
  Float_t reco_muonvtxZ;
  Float_t reco_muonMomentum[4];
  Float_t reco_showerMomentum[4];
  
  
  Float_t nuvtx_diff;
  Float_t showervtx_diff;
  Float_t muonvtx_diff;
  Float_t truth_showerMomentum[4];             // 0=X, 1=Y, 2=Z, 4=time
  Float_t truth_corr_nuvtxX;
  Float_t truth_corr_nuvtxY;
  Float_t truth_corr_nuvtxZ;
  Float_t truth_corr_showervtxX;
  Float_t truth_corr_showervtxY;
  Float_t truth_corr_showervtxZ;
  Float_t truth_showerKE;
  Float_t truth_corr_muonvtxX;
  Float_t truth_corr_muonvtxY;
  Float_t truth_corr_muonvtxZ;
  Float_t truth_muonvtxX;
  Float_t truth_muonvtxY;
  Float_t truth_muonvtxZ;
  Float_t truth_muonendX;
  Float_t truth_muonendY;
  Float_t truth_muonendZ;
  Float_t truth_muonMomentum[4];
  Float_t truth_nuEnergy;
  Float_t truth_energyInside;
  Float_t truth_electronInside;
  Int_t truth_nuPdg;
  Bool_t truth_isCC;
  Float_t truth_vtxX;
  Float_t truth_vtxY;
  Float_t truth_vtxZ;
  Float_t truth_nuTime;
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


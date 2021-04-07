struct PFDumpInfo{
    
    Int_t           truth_Ntrack;
    Int_t           truth_id[9243];   //[truth_Ntrack]
    Int_t           truth_pdg[9243];   //[truth_Ntrack]
    Int_t           truth_process[9243];   //[truth_Ntrack]
    Int_t           truth_mother[9243];   //[truth_Ntrack]
    Float_t         truth_startXYZT[9243][4];   //[truth_Ntrack]
    Float_t         truth_endXYZT[9243][4];   //[truth_Ntrack]
    Float_t         truth_startMomentum[9243][4];   //[truth_Ntrack]
    Float_t         truth_endMomentum[9243][4];   //[truth_Ntrack]
    vector<vector<int> > *truth_daughters;

    Int_t           reco_Ntrack;
    Int_t           reco_id[147];   //[reco_Ntrack]
    Int_t           reco_pdg[147];   //[reco_Ntrack]
    Int_t           reco_process[147];   //[reco_Ntrack]
    Int_t           reco_mother[147];   //[reco_Ntrack]
    Float_t         reco_startXYZT[147][4];   //[reco_Ntrack]
    Float_t         reco_endXYZT[147][4];   //[reco_Ntrack]
    Float_t         reco_startMomentum[147][4];   //[reco_Ntrack]
    Float_t         reco_endMomentum[147][4];   //[reco_Ntrack]
    //vector<vector<int> > *reco_daughters = nullptr;
};

void loadPFDumpInfo(TChain* tree0, PFDumpInfo& pfdump_info, int flag=1) {
    tree0->SetBranchAddress("reco_Ntrack", &pfdump_info.reco_Ntrack);
    tree0->SetBranchAddress("reco_id", &pfdump_info.reco_id[0]);
    tree0->SetBranchAddress("reco_pdg", &pfdump_info.reco_pdg[0]);
    tree0->SetBranchAddress("reco_process", &pfdump_info.reco_process[0]);
    tree0->SetBranchAddress("reco_mother", &pfdump_info.reco_mother[0]);
    tree0->SetBranchAddress("reco_startXYZT", &pfdump_info.reco_startXYZT[0][0]);
    tree0->SetBranchAddress("reco_endXYZT", &pfdump_info.reco_endXYZT[0][0]);
    tree0->SetBranchAddress("reco_startMomentum", &pfdump_info.reco_startMomentum[0][0]);
    tree0->SetBranchAddress("reco_endMomentum", &pfdump_info.reco_endMomentum[0][0]);
    //tree0->SetBranchAddress("reco_daughters", &pfdump_info.reco_daughters);

    if( flag==1 ){
        tree0->SetBranchAddress("truth_Ntrack", &pfdump_info.truth_Ntrack);
        tree0->SetBranchAddress("truth_id", &pfdump_info.truth_id[0]);
        tree0->SetBranchAddress("truth_pdg", &pfdump_info.truth_pdg[0]);
        tree0->SetBranchAddress("truth_process", &pfdump_info.truth_process[0]);
        tree0->SetBranchAddress("truth_mother", &pfdump_info.truth_mother[0]);
        tree0->SetBranchAddress("truth_startXYZT", &pfdump_info.truth_startXYZT[0][0]);
        tree0->SetBranchAddress("truth_endXYZT", &pfdump_info.truth_endXYZT[0][0]);
        tree0->SetBranchAddress("truth_startMomentum", &pfdump_info.truth_startMomentum[0][0]);
        tree0->SetBranchAddress("truth_endMomentum", &pfdump_info.truth_endMomentum[0][0]);
        //tree0->SetBranchAddress("truth_daughters", &pfdump_info.truth_daughters[0]);
    }

}
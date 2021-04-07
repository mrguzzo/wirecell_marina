#include "/home/marinaguzzo/Desktop/wirecell_marina/includes/FUNCTIONS.h"

void draw_1d_histogram(TCanvas* c1, TH1F* hist, string xlabel, string ylabel, int line_color, string draw_mode){

    c1->cd();

    hist->GetXaxis()->SetTitle(xlabel.c_str());
    hist->GetYaxis()->SetTitle(ylabel.c_str());

    hist->SetLineColor(line_color);

}

void calculate_distance(){
    
    // --------------------------- //
    // ----- declare classes ----- //
    // --------------------------- //

    readingCodes read;
    wireCellFunctions wcfunc;
    inputInfo info;

    // -------------------------- //
    // ----- selection cuts ----- //
    // -------------------------- //

    bool apply_is_nueCC      = false;
    bool apply_fc            = true;

    int run_number           = 2;
    bool write_run_number    = true;

    // ---------------------- //
    // ----- open files ----- //
    // ---------------------- //

    string filename = "data/nue_overlay_run2.root";
    TChain* T_pot = new TChain("wcpselection/T_pot");             T_pot->Add(filename.c_str());        POTInfo pot_info;            loadPOTInfo(T_pot, pot_info);
    TChain* T_kine = new TChain("wcpselection/T_KINEvars");       T_kine->Add(filename.c_str());       KineInfo kine_info;          loadKineInfo(T_kine, kine_info);
    TChain* T_eval = new TChain("wcpselection/T_eval");           T_eval->Add(filename.c_str());       EvalInfo eval_info;          loadEvalInfo(T_eval, eval_info);
    TChain* T_pfeval = new TChain("wcpselection/T_PFeval");       T_pfeval->Add(filename.c_str());     PFevalInfo pfeval_info;      loadPFevalInfo(T_pfeval, pfeval_info);
    TChain* T_tagger = new TChain("wcpselection/T_BDTvars");      T_tagger->Add(filename.c_str());     TaggerInfo tagger_info;      loadTaggerInfo(T_tagger, tagger_info);
    TChain* T_pfdump = new TChain("wcpselection/T_PFDump");       T_pfdump->Add(filename.c_str());     PFDumpInfo pfdump_info;      loadPFDumpInfo(T_pfdump, pfdump_info);

    // --------------------- //
    // ----- variables ----- //
    // --------------------- //

    int pdg_target = 12;

    int neutron=0, muon=0, kaon=0, pion=0, proton=0, photon=0, electron=0;

    enum {k_neutron, k_muon, k_kaon, k_pion, k_proton, k_photon, k_electron, k_MAX_particleID};
    std::vector<TH1F*> h_nue_particles;         h_nue_particles.resize(k_MAX_particleID); 
    for(unsigned int n=0; n<k_MAX_particleID; n++){
        h_nue_particles.at(n) = new TH1F(Form("particle_multiplicity_nue_%u",n),"",10,0,10);
    }

    // ----------------------------- //
    // ----- loop over entries ----- //
    // ----------------------------- //

    for(unsigned int n=0; n<T_pfdump->GetEntries(); n++){

        T_kine->GetEntry(n);
        T_eval->GetEntry(n);
        T_pfeval->GetEntry(n);
        T_tagger->GetEntry(n);
        T_pfdump->GetEntry(n);

        if( eval_info.truth_nuPdg != pdg_target ) continue;        // skip if it is not nue

        TVector3 truth_nuvtx(pfeval_info.truth_corr_nuvtxX, pfeval_info.truth_corr_nuvtxY, pfeval_info.truth_corr_nuvtxZ);
        TVector3 reco_nuvtx(pfeval_info.reco_nuvtxX, pfeval_info.reco_nuvtxY, pfeval_info.reco_nuvtxZ);

        // --------------------------------------- particle multiplicity

        neutron=0; muon=0; kaon=0; pion=0; proton=0; photon=0; electron=0;

        for(unsigned int i=0; i<pfdump_info.truth_Ntrack; i++){ // loop over daughter particles and count the particle multiplicity

            if( pfdump_info.truth_mother[i]!=0 ) continue; // skip if it's not coming from the neutrino interaction

            if(pfdump_info.truth_pdg[i]==2112) neutron++;                                           // neutron
            if((pfdump_info.truth_pdg[i]==13 || pfdump_info.truth_pdg[i]==-13)) muon++;             // muon
            if((pfdump_info.truth_pdg[i]==321 || pfdump_info.truth_pdg[i]==-321)) kaon++;           // kaon
            if((pfdump_info.truth_pdg[i]==211 || pfdump_info.truth_pdg[i]==-211)) pion++;           // pion
            if((pfdump_info.truth_pdg[i]==2212 || pfdump_info.truth_pdg[i]==-2212)) proton++;       // proton
            if(pfdump_info.truth_pdg[i]==22) photon++;                                              // photon
            if((pfdump_info.truth_pdg[i]==11 || pfdump_info.truth_pdg[i]==-11)) electron++;         // electron
        }

        h_nue_particles.at(k_neutron)->Fill(neutron);
        h_nue_particles.at(k_muon)->Fill(muon);
        h_nue_particles.at(k_kaon)->Fill(kaon);
        h_nue_particles.at(k_pion)->Fill(pion);
        h_nue_particles.at(k_proton)->Fill(proton);
        h_nue_particles.at(k_photon)->Fill(photon);
        h_nue_particles.at(k_electron)->Fill(electron);

        // --------------------------------------- calculate the distance between the neutrino vertex and the electron for electron==2

        if( electron>=2 ){ // we want to study cases with 2 electrons

            for(unsigned int m=0; m<pfdump_info.truth_Ntrack; m++){

                if( pfdump_info.truth_mother[m]!=0 ) continue; // skip if it's not coming from the neutrino interaction

                if((pfdump_info.truth_pdg[m]==11 || pfdump_info.truth_pdg[m]==-11)){

                    TVector3 truth_pos(pfdump_info.truth_startXYZT[m][0], pfdump_info.truth_startXYZT[m][1], pfdump_info.truth_startXYZT[m][2]);
                    TVector3 reco_pos(pfdump_info.reco_startXYZT[m][0], pfdump_info.reco_startXYZT[m][1], pfdump_info.reco_startXYZT[m][2]);

                    std::cout << "event=" << n << "   electron_multiplicity=" << electron << "  truth_dist=" << (truth_pos-truth_nuvtx).Mag() << "   reco_dist=" << (reco_pos-reco_nuvtx).Mag() << std::endl;
                    if( pfeval_info.truth_isCC ) std::cout << "the event number " << n << " is CC.\n";
                
                }
            }
        }

    } // loop over entries

    TCanvas *c1 = new TCanvas();

    // draw_1d_histogram(TCanvas* c1, TH1F* hist, string xlabel, string ylabel, int line_color, string draw_mode)

    draw_1d_histogram(c1, h_nue_particles.at(k_neutron), "Number of Neutrons", )


}
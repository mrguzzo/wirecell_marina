#include "/home/marinaguzzo/Desktop/wirecell_marina/includes/FUNCTIONS.h"

void draw_1d_histogram(TCanvas* c1, TH1F* hist, string xlabel, string ylabel, string filename, bool apply_area_normalisation){

    c1->cd();

    hist->GetXaxis()->SetTitle(xlabel.c_str());
    hist->GetYaxis()->SetTitle(ylabel.c_str());

    if( apply_area_normalisation ){
        hist->Scale(1/hist->GetEntries());
        hist->GetYaxis()->SetRangeUser(0.,1.);
    }

    hist->SetMarkerStyle(8);
    hist->SetMarkerColor(kPink);
    hist->Draw("P");

    c1->SaveAs(filename.c_str());

}

void compare_2_1dhistograms(TCanvas* c1, TH1F* hist1, TH1F* hist2, string xlabel, string ylabel, string filename, bool apply_area_normalisation, string legend1, string legend2){

    c1->cd();

    hist1->GetXaxis()->SetTitle(xlabel.c_str());     // x label
    hist1->GetYaxis()->SetTitle(ylabel.c_str());     // y label
    hist1->SetMarkerStyle(8);                        // hist1 visual style
    hist1->SetMarkerColor(kPink);
    hist1->Draw("P");

    hist2->SetMarkerStyle(8);
    hist2->SetMarkerColor(kBlue);
    hist2->Draw("P same");

    TLegend *leg = new TLegend(0.7020057,0.7707006,0.8696275,0.8747346,NULL,"brNDC");
    leg->AddEntry(hist1,Form("%s (%.0f)",legend1.c_str(), hist1->GetEntries()), "lep");
    leg->AddEntry(hist2,Form("%s (%.0f)",legend2.c_str(), hist2->GetEntries()), "lep");
    leg->SetBorderSize(0);
    leg->Draw("same");

    if( apply_area_normalisation ){
        hist1->Scale(1/hist1->GetEntries());         // apply normalisation
        hist2->Scale(1/hist2->GetEntries());

        //hist1->GetYaxis()->SetRangeUser(0.,1.);      // fix Y range to be [0,1]

        TLatex text_area_normalise;
        text_area_normalise.SetTextColor(kPink);
        text_area_normalise.DrawLatexNDC(0.7048711,0.9214437,"Area Normalized");
    }

    double max = hist1->GetMaximum();
    if( hist2->GetMaximum() > max ) max = hist2->GetMaximum();
    hist1->GetYaxis()->SetRangeUser(0.,max*1.5);

    c1->SaveAs(filename.c_str());

}

void calculating_distance(){
    
    // --------------------------- //
    // ----- declare classes ----- //
    // --------------------------- //

    readingCodes read;
    wireCellFunctions wcfunc;
    inputInfo info;

    read.makeDirectory("plots");

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

    enum neutrino_type {k_nue, k_nuebar, k_MAX_neutrinos};
    int nu_type = 0;

    int neutron=0, muon=0, kaon=0, pion=0, proton=0, photon=0, electron=0;

    enum {k_neutron, k_muon, k_kaon, k_pion, k_proton, k_photon, k_electron, k_MAX_particleID};

    // ---------------------- //
    // ----- histograms ----- //
    // ---------------------- //

    std::vector<TH1F*> h_shower_costheta;            h_shower_costheta.resize(k_MAX_neutrinos);
    std::vector<TH1F*> h_shower_energy;              h_shower_energy.resize(k_MAX_neutrinos);
    std::vector<TH1F*> h_neutrino_energy;            h_neutrino_energy.resize(k_MAX_neutrinos);

    for(unsigned int n=0; n<k_MAX_neutrinos; n++){
        h_shower_costheta.at(n) = new TH1F (Form("shower_costheta_%u",n),"",20,-1,1);
        h_shower_energy.at(n) = new TH1F (Form("shower_energy_%u", n),"",25,0,2500);
        h_neutrino_energy.at(n) = new TH1F (Form("neutrino_energy_%u", n),"",25,0,2500);
    }

    // ---------------------

    std::vector<std::vector<TH1F*>> h_particle_multiplicity; h_particle_multiplicity.resize(k_MAX_neutrinos);

    for(unsigned int n=0; n<k_MAX_neutrinos; n++){
        h_particle_multiplicity.at(n).resize(k_MAX_particleID);
        for(unsigned int i=0; i<k_MAX_particleID; i++){
            h_particle_multiplicity.at(n).at(i) = new TH1F(Form("particle_multiplicity_nue_%u_%u",n,i),"",10,0,10);
        }
    }

    TH1F* h_nue_costheta = new TH1F("","",20,-1,1);

    // ----------------------------- //
    // ----- loop over entries ----- //
    // ----------------------------- //

    for(unsigned int n=0; n<T_pfdump->GetEntries(); n++){

        T_kine->GetEntry(n);
        T_eval->GetEntry(n);
        T_pfeval->GetEntry(n);
        T_tagger->GetEntry(n);
        T_pfdump->GetEntry(n);

        if( eval_info.truth_nuPdg!=12 && eval_info.truth_nuPdg!=-12 ) continue;        // skip if it is not nue/nuebar

        if( eval_info.truth_nuPdg==12 ) nu_type = k_nue;                               // classify the neutrino
        else if( eval_info.truth_nuPdg==-12 ) nu_type = k_nuebar;

        // --------------------------------------- get some initial neutrino information

        TVector3 truth_nuvtx(pfeval_info.truth_corr_nuvtxX, pfeval_info.truth_corr_nuvtxY, pfeval_info.truth_corr_nuvtxZ);
        TVector3 reco_nuvtx(pfeval_info.reco_nuvtxX, pfeval_info.reco_nuvtxY, pfeval_info.reco_nuvtxZ);

        TVector3 truth_shower_momentum(pfeval_info.truth_showerMomentum[0], pfeval_info.truth_showerMomentum[1], pfeval_info.truth_showerMomentum[2]);
        TVector3 truth_shower_momentum_dir = truth_shower_momentum.Unit();

        h_shower_costheta.at(nu_type)->Fill(-read.costheta_wrt_beam(truth_shower_momentum_dir)); // (-) because v_targ_uboone_dir points to the target

        h_shower_energy.at(nu_type)->Fill(pfeval_info.truth_showerKE);

        h_neutrino_energy.at(nu_type)->Fill(pfeval_info.truth_nuEnergy);

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

        h_particle_multiplicity.at(nu_type).at(k_neutron)->Fill(neutron);
        h_particle_multiplicity.at(nu_type).at(k_muon)->Fill(muon);
        h_particle_multiplicity.at(nu_type).at(k_kaon)->Fill(kaon);
        h_particle_multiplicity.at(nu_type).at(k_pion)->Fill(pion);
        h_particle_multiplicity.at(nu_type).at(k_proton)->Fill(proton);
        h_particle_multiplicity.at(nu_type).at(k_photon)->Fill(photon);
        h_particle_multiplicity.at(nu_type).at(k_electron)->Fill(electron);

        // --------------------------------------- calculate the distance between the neutrino vertex and the electron for electron==2, print their energy as well

        /*if( electron>=2 ){ // we want to study cases with 2 electrons

            for(unsigned int m=0; m<pfdump_info.truth_Ntrack; m++){

                if( pfdump_info.truth_mother[m]!=0 ) continue; // skip if it's not coming from the neutrino interaction

                if( (pfdump_info.truth_pdg[m]!=11 && pfdump_info.truth_pdg[m]!=-11) ) continue;  // skip if it's not an electron

                // --- particle distance to the neutrino vertex

                TVector3 truth_pos(pfdump_info.truth_startXYZT[m][0], pfdump_info.truth_startXYZT[m][1], pfdump_info.truth_startXYZT[m][2]);
                TVector3 reco_pos(pfdump_info.reco_startXYZT[m][0], pfdump_info.reco_startXYZT[m][1], pfdump_info.reco_startXYZT[m][2]);

                double truth_distance = (truth_pos-truth_nuvtx).Mag();
                double reco_distance = (reco_pos-reco_nuvtx).Mag();

                std::cout << "evt=" << n << " run=" << pfeval_info.run << " subrun=" << pfeval_info.subrun << std::endl;
                std::cout << "electron_multiplicity=" << electron << std::endl;
                std::cout << "truth_dist=" << (truth_pos-truth_nuvtx).Mag() << "   reco_dist=" << (reco_pos-reco_nuvtx).Mag() << std::endl;
                std::cout << "truth_x=" << truth_pos.X() << "  truth_Y=" << truth_pos.Y() << "  truth_Z=" << truth_pos.Z() << std::endl;

                // --- particle direction

                TVector3 v_particle_startMomentum(pfdump_info.truth_startMomentum[m][0], pfdump_info.truth_startMomentum[m][1], pfdump_info.truth_startMomentum[m][2]);
                TVector3 v_particle_startMomentum_dir = v_particle_startMomentum.Unit();

                double cos_theta = v_particle_startMomentum_dir.Dot(info.v_targ_uboone_dir);

                std::cout << "cos_theta=" << cos_theta << std::endl;

                std::cout << std::endl;

            }

            std::cout << "-----------------------------------" << std::endl;*/
        //}

    } // loop over entries

    // --------------------------- //
    // ----- draw histograms ----- //
    // --------------------------- //

    TCanvas *c1 = new TCanvas();
    gStyle->SetOptStat(0);

    compare_2_1dhistograms(c1, h_particle_multiplicity.at(k_nue).at(k_neutron), h_particle_multiplicity.at(k_nuebar).at(k_neutron), "Number of Neutrons (truth)", "Entries / bin", "plots/multiplicity_neutron.pdf", true, "#nu_{e}", "#bar{#nu}_{e}");
    compare_2_1dhistograms(c1, h_particle_multiplicity.at(k_nue).at(k_muon), h_particle_multiplicity.at(k_nuebar).at(k_muon), "Number of Muons (truth)", "Entries / bin", "plots/multiplicity_muon.pdf", true, "#nu_{e}", "#bar{#nu}_{e}");
    compare_2_1dhistograms(c1, h_particle_multiplicity.at(k_nue).at(k_kaon), h_particle_multiplicity.at(k_nuebar).at(k_kaon), "Number of Kaons (truth)", "Entries / bin", "plots/multiplicity_kaon.pdf", true, "#nu_{e}", "#bar{#nu}_{e}");
    compare_2_1dhistograms(c1, h_particle_multiplicity.at(k_nue).at(k_pion), h_particle_multiplicity.at(k_nuebar).at(k_pion), "Number of Pions (truth)", "Entries / bin", "plots/multiplicity_pion.pdf", true, "#nu_{e}", "#bar{#nu}_{e}");
    compare_2_1dhistograms(c1, h_particle_multiplicity.at(k_nue).at(k_proton), h_particle_multiplicity.at(k_nuebar).at(k_proton), "Number of Protons (truth)", "Entries / bin", "plots/multiplicity_proton.pdf", true, "#nu_{e}", "#bar{#nu}_{e}");
    compare_2_1dhistograms(c1, h_particle_multiplicity.at(k_nue).at(k_photon), h_particle_multiplicity.at(k_nuebar).at(k_photon), "Number of Photons (truth)", "Entries / bin", "plots/multiplicity_photon.pdf", true, "#nu_{e}", "#bar{#nu}_{e}");
    compare_2_1dhistograms(c1, h_particle_multiplicity.at(k_nue).at(k_electron), h_particle_multiplicity.at(k_nuebar).at(k_electron), "Number of Electrons (truth)", "Entries / bin", "plots/multiplicity_electron.pdf", true, "#nu_{e}", "#bar{#nu}_{e}");
    
    compare_2_1dhistograms(c1, h_shower_costheta.at(k_nue), h_shower_costheta.at(k_nuebar), "Shower cos #theta (truth)", "Entries / bin", "plots/shower_costheta.pdf", true, "#nu_{e}", "#bar{#nu}_{e}");
    compare_2_1dhistograms(c1, h_shower_energy.at(k_nue), h_shower_energy.at(k_nuebar), "Shower Energy [MeV] (truth)", "Entries / bin", "plots/shower_energy.pdf", true, "#nu_{e}", "#bar{#nu}_{e}");
    compare_2_1dhistograms(c1, h_neutrino_energy.at(k_nue), h_neutrino_energy.at(k_nuebar), "Neutrino Energy [MeV] (truth)", "Entries / bin", "plots/neutrino_energy.pdf", true, "#nu_{e}", "#bar{#nu}_{e}");

    //draw_1d_histogram(c1, h_nue_costheta, "Truth #nu_{e} cos#theta", "Entries / bin", "plots/nue_costheta.pdf", true);


}
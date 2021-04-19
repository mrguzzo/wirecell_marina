#include "/home/marinaguzzo/Desktop/wirecell_marina/includes/MYLIBRARY.h"


void marina_wirecell(){

    makeHists();

    // ======================================================================================= //
    //                                      OPENING FILES                                      //
    // ======================================================================================= //

    // MC
    TChain* T_pot = new TChain("wcpselection/T_pot"); T_pot->Add("rootfiles/checkout_prodgenie_numi_overlay_run1.root"); POTInfo pot_info; loadPOTInfo(T_pot, pot_info);
    TChain* T_kine = new TChain("wcpselection/T_KINEvars"); T_kine->Add("rootfiles/checkout_prodgenie_numi_overlay_run1.root"); KineInfo kine_info; loadKineInfo(T_kine, kine_info);
    TChain* T_eval = new TChain("wcpselection/T_eval"); T_eval->Add("rootfiles/checkout_prodgenie_numi_overlay_run1.root"); EvalInfo eval_info; loadEvalInfo(T_eval, eval_info);
    TChain* T_pfeval = new TChain("wcpselection/T_PFeval"); T_pfeval->Add("rootfiles/checkout_prodgenie_numi_overlay_run1.root"); PFevalInfo pfeval_info; loadPFevalInfo(T_pfeval, pfeval_info);
    TChain* T_tagger = new TChain("wcpselection/T_BDTvars"); T_tagger->Add("rootfiles/checkout_prodgenie_numi_overlay_run1.root"); TaggerInfo tagger_info; loadTaggerInfo(T_tagger, tagger_info);

    // intrinsic nue 
    TChain* T_pot_nue = new TChain("wcpselection/T_pot"); T_pot_nue->Add("rootfiles/checkout_prodgenie_numi_intrinsic_nue_overlay_run1.root"); POTInfo pot_info_nue; loadPOTInfo(T_pot_nue, pot_info_nue);
    TChain* T_kine_nue = new TChain("wcpselection/T_KINEvars"); T_kine_nue->Add("rootfiles/checkout_prodgenie_numi_intrinsic_nue_overlay_run1.root"); KineInfo kine_info_nue; loadKineInfo(T_kine_nue, kine_info_nue);
    TChain* T_eval_nue = new TChain("wcpselection/T_eval"); T_eval_nue->Add("rootfiles/checkout_prodgenie_numi_intrinsic_nue_overlay_run1.root"); EvalInfo eval_info_nue; loadEvalInfo(T_eval_nue, eval_info_nue);
    TChain* T_pfeval_nue = new TChain("wcpselection/T_PFeval"); T_pfeval_nue->Add("rootfiles/checkout_prodgenie_numi_intrinsic_nue_overlay_run1.root"); PFevalInfo pfeval_info_nue; loadPFevalInfo(T_pfeval_nue, pfeval_info_nue);
    TChain* T_tagger_nue = new TChain("wcpselection/T_BDTvars"); T_tagger_nue->Add("rootfiles/checkout_prodgenie_numi_intrinsic_nue_overlay_run1.root"); TaggerInfo tagger_info_nue; loadTaggerInfo(T_tagger_nue, tagger_info_nue);

    // EXT
    TChain* T_kine_ext = new TChain("wcpselection/T_KINEvars"); T_kine_ext->Add(filename_ext.c_str()); KineInfo kine_info_ext; loadKineInfo(T_kine_ext, kine_info_ext);
    TChain* T_eval_ext = new TChain("wcpselection/T_eval"); T_eval_ext->Add(filename_ext.c_str()); EvalInfo eval_info_ext; loadEvalInfo(T_eval_ext, eval_info_ext, 0);
    TChain* T_tagger_ext = new TChain("wcpselection/T_BDTvars"); T_tagger_ext->Add(filename_ext.c_str()); TaggerInfo tagger_info_ext; loadTaggerInfo(T_tagger_ext, tagger_info_ext);
    TChain* T_pfeval_ext = new TChain("wcpselection/T_PFeval"); T_pfeval_ext->Add(filename_ext.c_str()); PFevalInfo pfeval_info_ext; loadPFevalInfo(T_pfeval_ext, pfeval_info_ext, 0);

    // DATA
    TChain* T_kine_data = new TChain("wcpselection/T_KINEvars"); T_kine_data->Add(filename_data.c_str()); KineInfo kine_info_data; loadKineInfo(T_kine_data, kine_info_data);
    TChain* T_eval_data = new TChain("wcpselection/T_eval"); T_eval_data->Add(filename_data.c_str()); EvalInfo eval_info_data; loadEvalInfo(T_eval_data, eval_info_data, 0);
    TChain* T_tagger_data = new TChain("wcpselection/T_BDTvars"); T_tagger_data->Add(filename_data.c_str()); TaggerInfo tagger_info_data; loadTaggerInfo(T_tagger_data, tagger_info_data);
    TChain* T_pfeval_data = new TChain("wcpselection/T_PFeval");  T_pfeval_data->Add(filename_data.c_str()); PFevalInfo pfeval_info_data;  loadPFevalInfo(T_pfeval_data, pfeval_info_data, 0);

    // ===================================================================================== //
    //                                      GETTING POT                                      //
    // ===================================================================================== //

    // get overlay POT
    for (int n = 0; n < T_pot->GetEntries(); ++n) {
        T_pot->GetEntry(n);
        MCPOT += (float)pot_info.pot_tor875;
    }

    // get intrinsic nue POT
    for (int n = 0; n < T_pot_nue->GetEntries(); ++n) {
        T_pot_nue->GetEntry(n);
        NuePOT += (float)pot_info_nue.pot_tor875;
    }

    // ========================================================================================= //
    //                                      GETTING WEIGHTS                                      //
    // ========================================================================================= //

    wMC = dataPOT/MCPOT; // POT weight   
    wNue = dataPOT/NuePOT; // nue POT weight 
    wEXT = 0.98*2393691.0/2993141.315; // EXT weight
    if (MCPOT == 0) wMC = 1;
    if (NuePOT == 0) wMC = 1;

    std::cout   << "-----------------------------" << "\n"
                << "data POT:     " << dataPOT << "\n"
                << "MC POT:       " << MCPOT << "\n"
                << "Nue POT:      " << NuePOT << "\n"
                << "-----------------------------" << "\n"
                << "MC weight:    " << wMC << "\n"
                << "Nue weight:   " << wNue << "\n"
                << "EXT weight:   " << wEXT << "\n"
                << "data weight:  " << wData << "\n"
                << "-----------------------------" << std::endl;

    // =================================================================================== //
    //                                      MC SAMPLE                                      //
    // =================================================================================== //

    bool isCosmic, isNueCC, isNueBarCC, isNumuCC, isNC, hasPi0, hasChargedPi, hasProton, isOutFV; // topology tagging information
    int nNueCC = 0, nNueBarCC = 0, nNumuCC = 0, nNumuCC_pi0 = 0, nNC = 0, nNC_pi0 = 0, nOutFV = 0, nCosmic = 0, nDirt = 0, nEXT = 0, nData = 0; // counting events
    int nevents = T_eval->GetEntries(); // event loop

    for (int n = 0; n < nevents; ++n) {

        if (n%10000 == 0) std::cout << "[MC] processing event " << n << "/" << nevents << std::endl;

        T_eval->GetEntry(n);
        T_pfeval->GetEntry(n);
        T_kine->GetEntry(n);
        T_tagger->GetEntry(n);

        // ---------------------------- //
        //      STACKED HISTOGRAMS      //
        // ---------------------------- //

        isOutFV = false;
        isCosmic = false;
        isNueCC = false; isNueBarCC = false; isNumuCC = false; isNC = false;
        hasPi0 = false; hasChargedPi = false; hasProton = false;

        // --- classify true information

        if (eval_info.truth_energyInside != 0 && (eval_info.match_completeness_energy/eval_info.truth_energyInside) < 0.1) isCosmic = true;
        if (eval_info.match_completeness_energy/eval_info.truth_energyInside>=0.1 && eval_info.truth_vtxInside==0) isOutFV = true;
        if (pfeval_info.truth_NprimPio > 0) hasPi0 = true;

        if (eval_info.truth_isCC) {
            if (eval_info.truth_nuPdg == 12) isNueCC = true;
            else if (eval_info.truth_nuPdg == -12) isNueBarCC = true;
            else if (abs(eval_info.truth_nuPdg) == 14) isNumuCC = true;
        }
        else { isNC = true; }

        // --- topology

        int top = -1;

        if (isCosmic) { top = cosmic; nCosmic++; }
        else if (isOutFV) { top = outFV; nOutFV++; }
        else if (isNueCC) { top = -1; nNueCC++; } // use intrinsic nue 
        else if (isNueBarCC) { top = -1; nNueBarCC++; }  // use intrinsic nue
        else if (isNumuCC) {
            if (hasPi0) { top = ccnumupi0; nNumuCC_pi0++; }
            else { top = ccnumu; nNumuCC++; }
        }
        else if (isNC) {
            if (hasPi0) { top = ncpi0; nNC_pi0++; }
            else { top = nc;  nNC++; }
        }

        // --- weight

        float genie_weight = (float)(eval_info.weight_cv*eval_info.weight_spline);
        genie_weight=1;

        // --- fill histograms

        if (top > -1) fillStackedHists(eval_info, tagger_info, kine_info, top, wMC*genie_weight);

        // ------------------------- //
        //      ENERGY SPECTRUM      //
        // ------------------------- //

        int nu_type;
        if( pfeval_info.truth_nuPdg == 14 ) nu_type = k_numu;
        else if( pfeval_info.truth_nuPdg == -14 ) nu_type = k_numubar;
        else if( pfeval_info.truth_nuPdg == 12 ) nu_type = k_nue;
        else if( pfeval_info.truth_nuPdg == -12 ) nu_type = k_nuebar;

        if(nu_type==k_numu || nu_type==k_numubar){

            if( !isnan(eval_info.weight_cv) && !isinf(eval_info.weight_cv) ){

                if(eval_info.match_isFC == 1) h_energy_true[k_fc][nu_type]->Fill(pfeval_info.truth_nuEnergy,eval_info.weight_cv);
                else if(eval_info.match_isFC == 0) h_energy_true[k_pc][nu_type]->Fill(pfeval_info.truth_nuEnergy,eval_info.weight_cv);

            }
        }

    } // overlay loop

    // ==================================================================================== //
    //                                      NUE SAMPLE                                      //
    // ==================================================================================== //

    // event loop
    nevents = T_eval_nue->GetEntries();

    for (int n = 0; n < nevents; ++n) {

        if (n%10000 == 0) std::cout << "[Nue] processing event " << n << "/" << nevents << std::endl;

        T_eval_nue->GetEntry(n);
        T_pfeval_nue->GetEntry(n);
        T_kine_nue->GetEntry(n);
        T_tagger_nue->GetEntry(n);

        // ---------------------------- //
        //      STACKED HISTOGRAMS      //
        // ---------------------------- //

        isOutFV = false;
        isCosmic = false;
        isNueCC = false; isNueBarCC = false;

        // --- classify true information

        if (eval_info.truth_energyInside != 0 && (eval_info.match_completeness_energy/eval_info.truth_energyInside) < 0.1) isCosmic = true;
        if (eval_info_nue.match_completeness_energy/eval_info_nue.truth_energyInside>=0.1 && eval_info_nue.truth_vtxInside==0) isOutFV = true;

        if (eval_info_nue.truth_isCC) {
            if (eval_info_nue.truth_nuPdg == 12) isNueCC = true;
            else if (eval_info_nue.truth_nuPdg == -12) isNueBarCC = true;
        }

        // --- topology 

        int top = -1;

        if (isCosmic) { top = cosmic; }
        else if (isOutFV) { top = outFV; }
        else if (isNueCC) { top = ccnue; }
        else if (isNueBarCC) { top = ccnuebar; }

        // --- weight

        float genie_weight = (float)(eval_info_nue.weight_cv*eval_info_nue.weight_spline);
        genie_weight=1;

        // --- fill histograms

        if (top > -1) fillStackedHists(eval_info_nue, tagger_info_nue, kine_info_nue, top, wNue*genie_weight);

        // ------------------------- //
        //      ENERGY SPECTRUM      //
        // ------------------------- //

        int nu_type;
        if (isNueCC) nu_type = k_nue;
        else if (isNueBarCC) nu_type = k_nuebar;

        if (nu_type==k_nue || nu_type==k_nuebar){

            if( !isnan(eval_info_nue.weight_cv) && !isinf(eval_info_nue.weight_cv) ){

                if(eval_info_nue.match_isFC == 1) h_energy_true[k_fc][nu_type]->Fill(pfeval_info_nue.truth_nuEnergy,eval_info_nue.weight_cv);
                else if(eval_info_nue.match_isFC == 0) h_energy_true[k_pc][nu_type]->Fill(pfeval_info_nue.truth_nuEnergy,eval_info.weight_cv);

            }
        }

    } // intrinsic nue loop 

    // ==================================================================================== //
    //                                      EXT SAMPLE                                      //
    // ==================================================================================== //

    nevents = T_eval_ext->GetEntries();

    for (int n = 0; n < nevents; ++n) {

        if (n%10000 == 0) std::cout << "[EXT] processing event " << n << "/" << nevents << std::endl;

        T_eval_ext->GetEntry(n);
        T_kine_ext->GetEntry(n);
        T_tagger_ext->GetEntry(n);

        // ---------------------------- //
        //      STACKED HISTOGRAMS      //
        // ---------------------------- //

        nEXT++;
        fillStackedHists(eval_info_ext, tagger_info_ext, kine_info_ext, ext, wEXT);

    } // EXT loop 

    // ===================================================================================== //
    //                                      DATA SAMPLE                                      //
    // ===================================================================================== //

    nevents = T_eval_data->GetEntries();

    for (int n = 0; n < nevents; ++n) {

        if (n%10000 == 0) std::cout << "[Data] processing event " << n << "/" << nevents << std::endl;

        T_eval_data->GetEntry(n);
        T_kine_data->GetEntry(n);
        T_tagger_data->GetEntry(n);

        // ---------------------------- //
        //      STACKED HISTOGRAMS      //
        // ---------------------------- //

        nData++;
        fillStackedHists(eval_info_data, tagger_info_data, kine_info_data, data, wData);

    } // DATA loop  
    
    drawAllHists();

}
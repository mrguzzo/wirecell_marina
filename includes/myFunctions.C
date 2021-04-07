// Marina Reggiani-Guzzo, April 2021

#include "FUNCTIONS.h"

// -------------------------------------------------------------------------------------- //
// ----- Input information that can be used as global variables throughout the code ----- //
// -------------------------------------------------------------------------------------- //

class inputInfo{

    public:

        const Float_t fidvol = 3;                    // Distance between AV and FV in cm
        const Float_t TPC_xi = -1.55 + fidvol;
        const Float_t TPC_xf = 254.8 - fidvol;
        const Float_t TPC_yi = -115.53 + fidvol;
        const Float_t TPC_yf = 117.47 - fidvol;
        const Float_t TPC_zi = 0.1 + fidvol;
        const Float_t TPC_zf = 1036.9 - fidvol;

        TVector3 v_targ_uboone;
        TVector3 v_targ_uboone_dir;

        inputInfo(); // initialiser
        
};

inputInfo::inputInfo(){

    v_targ_uboone = {-31387.58422, -3316.402543, -60100.2414};
    v_targ_uboone_dir = v_targ_uboone.Unit();

}

// -------------------------------------------------- //
// ----- Functions that I wrote for my analysis ----- //
// -------------------------------------------------- //

class readingCodes{

    inputInfo info;

    public:

        double costheta_wrt_beam(TVector3 vec){
            TVector3 vec_dir = vec.Unit();
            return vec.Dot(info.v_targ_uboone_dir);
        }

        // --------------------------------------------------

        void makeDirectory(string dirname){
            gSystem->Exec(Form("if [ ! -d \"%s\" ]; then echo \"\n%s folder does not exist... creating\"; mkdir -p %s; fi",dirname.c_str(),dirname.c_str(),dirname.c_str()));
        }

        // --------------------------------------------------

        void drawPlots(TCanvas* c1, TH1F* h_nue, TH1F* h_nuebar, string xlabel, string ylabel, string filename,
                       bool apply_fc, bool write_run_number, int run_number, bool area_normalise){

            c1->cd();

            // --- axis

            h_nue->GetYaxis()->SetTitle(Form("%s", ylabel.c_str()));
            h_nue->GetYaxis()->SetLabelSize(0.05);
            h_nue->GetYaxis()->SetTitleSize(0.05);
            h_nue->GetYaxis()->SetTitleOffset(0.9);
            h_nue->GetXaxis()->SetTitle(Form("%s", xlabel.c_str()));
            h_nue->GetXaxis()->SetLabelSize(0.05);
            h_nue->GetXaxis()->SetTitleSize(0.05);
            h_nue->GetXaxis()->SetTitleOffset(1.08);

            if( area_normalise ){
                h_nue->Scale(1/h_nue->GetEntries());
                h_nuebar->Scale(1/h_nuebar->GetEntries());
                h_nue->GetYaxis()->SetRangeUser(0,1.1);
            }

            // --- legend

            TLegend *leg_efficiency = new TLegend(0.6002865,0.7103594,0.8653295,0.8646934,NULL,"brNDC");
            leg_efficiency->AddEntry(h_nue,Form("#nu_{e} (%.0f)",h_nue->GetEntries()), "lep");
            leg_efficiency->AddEntry(h_nuebar,Form("#bar{#nu}_{e} (%.0f)",h_nuebar->GetEntries()), "lep");
            leg_efficiency->SetBorderSize(0);

            h_nue->SetMarkerColor(kPink);
            h_nuebar->SetMarkerColor(kBlue);
            h_nue->SetMarkerStyle(8);
            h_nuebar->SetMarkerStyle(8);
            h_nue->Draw("P");
            h_nuebar->Draw("P same");
            leg_efficiency->Draw("same");

            // --- draw watermark

            c1->cd();

            TLatex text_efficiency;
            if( apply_fc ) text_efficiency.DrawLatexNDC(0.8123209,0.9217759, "FC");
            else text_efficiency.DrawLatexNDC(0.8123209,0.9217759, "PC");

            TLatex text_run_number;
            if( write_run_number ) text_run_number.DrawLatexNDC(0.1174785,0.9260042, Form("Run %u", run_number));
            text_run_number.SetLineColor(kPink);

            TLatex text_area_normalise;
            text_area_normalise.SetTextColor(kPink);
            if( area_normalise ) text_area_normalise.DrawLatexNDC(0.6260745,0.653277,"Area Normalized");

            // --- save pdf

            c1->SaveAs(filename.c_str());

            //c1->Clear();

        }

        // --------------------------------------------------

        void draw2dhist(TCanvas* c1, TH2F* hist, string xlabel, string ylabel, string filename){

            hist->GetXaxis()->SetTitle(xlabel.c_str());
            hist->GetYaxis()->SetTitle(ylabel.c_str());

            c1->cd();
            hist->Draw("colz");
            c1->SaveAs(filename.c_str());
        }

        // --------------------------------------------------

        void distPoint2Point(Float_t x1, Float_t y1, Float_t z1, Float_t x2, Float_t y2, Float_t z2){
            Float_t dist = TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
            std::cout << "x1=" << 
            std::cout << "   dist=" << dist << std::endl;
        }

        // --------------------------------------------------

        void giveMeDaughterParticleEnergy(int entry, PFDumpInfo &pfdump){
            
        }

};

// -------------------------------------------------------------- //
// ----- Functions imported directly from the WireCell team ----- //
// -------------------------------------------------------------- //

class wireCellFunctions{

    inputInfo info;

    public:

        bool inFidVolTrue(EvalInfo &e) {
            return ((e.truth_vtxX >= info.TPC_xi && e.truth_vtxX <= info.TPC_xf) &&
                    (e.truth_vtxY >= info.TPC_yi && e.truth_vtxY <= info.TPC_yf) &&
                    (e.truth_vtxZ >= info.TPC_zi && e.truth_vtxZ <= info.TPC_zf));

        } // inFidVolReco

        // --------------------------------------------------

        bool inFidVolRecoWCP(PFevalInfo &e) {
            return ((e.reco_nuvtxX >= info.TPC_xi && e.reco_nuvtxX <= info.TPC_xf) &&
                    (e.reco_nuvtxY >= info.TPC_yi && e.reco_nuvtxY <= info.TPC_yf) &&
                    (e.reco_nuvtxZ >= info.TPC_zi && e.reco_nuvtxZ <= info.TPC_zf));
        } // inFidVolReco

        // --------------------------------------------------

        Float_t dist(Float_t x1, Float_t y1, Float_t z1, Float_t x2, Float_t y2, Float_t z2) {
            return  sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
        } // dist

        // --------------------------------------------------

        bool is_nueCC(TaggerInfo &tagger_info){
            bool flag = false;
            // default 7.0
            if (tagger_info.numu_cc_flag >=0 && tagger_info.nue_score > 7.0)
                flag = true;
            
            return flag;
        } // is_nueCC

        // --------------------------------------------------

        bool is_numuCC(TaggerInfo &tagger_info){
            bool flag = false;

            if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9)
                flag = true;
            
            return flag;
        } // is_numuCC

        // --------------------------------------------------

        bool is_cc_pi0(KineInfo &kine){
            bool flag = false;

            if (  (kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 ) 
                   && kine.kine_pio_energy_1 > 40 
                   && kine.kine_pio_energy_2 > 25 
                   && kine.kine_pio_dis_1 < 110 
                   && kine.kine_pio_dis_2 < 120 
                   && kine.kine_pio_angle > 0 
                   && kine.kine_pio_angle < 174  
                   && kine.kine_pio_mass > 22 
                   && kine.kine_pio_mass < 300  )

                flag = true;
            
            return flag;
        } // is_cc_pi0

};
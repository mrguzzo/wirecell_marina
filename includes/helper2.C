#include "MYLIBRARY.h"

bool is_nueCC(TaggerInfo& tagger_info);
bool is_cc_pi0(KineInfo& kine);

// ================================================================= //
//                    CHECK THE INFORMATION BELOW                    //
// ================================================================= //

string filename_data = "rootfiles/checkout_data_numi_run1_morestat.root";
string filename_ext = "rootfiles/checkout_data_extnumi_run1.root";
string filename_mc = "rootfiles/checkout_prodgenie_numi_overlay_run1_OFFSETFIXED2.root";
string filename_nue = "rootfiles/checkout_prodgenie_numi_intrinsic_nue_overlay_run1_OFFSETFIXED2.root";

const float dataPOT = 2.064e+20;
const float EXTPOT = 1.1603984e+20;
float MCPOT; // calculated in the code
float NuePOT; // calculated in the code

const float wData = 1;
float wEXT; // calculated in the code
float wNue; // calculated in the code
float wMC; // calculated in the code

int run_number = 1;
const float plotscalefactor = 1.5;         // space above the stacked histogram plotscalefactor*hist->GetMaximum()   

const float fidvol = 3;
const float TPC_xi = -1.55 + fidvol;
const float TPC_xf = 254.8 - fidvol;
const float TPC_yi = -115.53 + fidvol;
const float TPC_yf = 117.47 - fidvol;
const float TPC_zi = 0.1 + fidvol;
const float TPC_zf = 1036.9 - fidvol;

// ================================================ //
//                    HISTOGRAMS                    //
// ================================================ //

enum {ccnue=0, ccnuebar=1, nc=3, ncpi0=2, ccnumu=5, ccnumupi0=4, outFV=6, cosmic=7, dirt=8, ext=9, data=10, k_MAX_topologies};
const char* topologies[k_MAX_topologies] = {"nue", "nuebar",  "NC_pi0", "NC", "numuCC_pi0", "numuCC", "outFV", "cosmic", "dirt", "ext", "data"};
const char* legendLabel[k_MAX_topologies] = {"#nu_{e}", "#bar{#nu}_{e} CC",  "#nu NC #pi_{0}", "#nu NC", "#nu_{#mu} CC #pi_{0}", "#nu_{#mu} CC", "out of FV", "cosmic", "dirt", "EXT", "data"};
const int stackOrder[k_MAX_topologies] = {ccnue, ccnuebar, nc, ncpi0, ccnumu, ccnumupi0, outFV, cosmic, dirt , ext, data};
const int colors[k_MAX_topologies] = {2, 46, 3, 8, 4, 9, 42, 6, 28, 15, 1};
const int lines[k_MAX_topologies] = {2, 46, 3, 8, 4, 9, 42, 6, 28, 15, 1};
const int fills[k_MAX_topologies] = {1001, 1001, 1001, 1001, 1001, 1001, 1001, 3004, 3004, 3004, 1001};

enum {k_pc, k_fc, k_MAX_wirecell};
const char* classification_title[k_MAX_wirecell] = {"PC", "FC"};
enum {k_nue=0, k_nuebar, k_numu, k_numubar, k_MAX_numodes};

TH1F* h_energy_true[k_MAX_wirecell][k_MAX_numodes];

enum TH1F_stacked_vec{
    k_Ereco,
    k_MAX_TH1F_stacked
};
TH1F* h_TH1F_stacked[k_MAX_wirecell][k_MAX_TH1F_stacked][k_MAX_topologies];

/* --------------------------------------- */
void makeHists(){

  // --- stacked histograms
  for(unsigned int i=0; i<k_MAX_wirecell; i++){
    for(unsigned int n=0; n<k_MAX_topologies; n++){
      h_TH1F_stacked[i][k_Ereco][n] = new TH1F(Form("h_Ereco_%u_%u",i,n),"",25,0,2500);
    }
  }

  // --- energy spectrum
  for(unsigned int i=0; i<k_MAX_wirecell; i++){
    for(unsigned int n=0; n<k_MAX_numodes; n++){
      h_energy_true[i][n] = new TH1F(Form("h_energy_%u_%u",i,n),"",25,0,4000);
    }
  }

}
/* --------------------------------------- */
void fillStackedHists(EvalInfo &eval_info, TaggerInfo &tagger_info, KineInfo &kine_info, int top, float weight) {

    if (tagger_info.numu_cc_flag < 0) return; 

    // --- nueCC and FC
    if (eval_info.match_isFC == 1 && is_nueCC(tagger_info)) {
      h_TH1F_stacked[k_fc][k_Ereco][top]->Fill(kine_info.kine_reco_Enu, weight);
    }

    // --- nue and PC
    if (eval_info.match_isFC == 0 && is_nueCC(tagger_info)) {
      h_TH1F_stacked[k_pc][k_Ereco][top]->Fill(kine_info.kine_reco_Enu, weight);
    }

}
/* --------------------------------------- */
void drawStackedHists(TCanvas *c1, int classification, int type, const char* plottitle, const char* savename, int x0, int xf, const char* legendtitle){

    c1->Clear();
    c1->cd();

    // --- legend
    TLegend *legend = new TLegend(0.6,0.5,0.87,0.87,NULL,"brNDC");
    legend->SetBorderSize(0);
    legend->SetNColumns(2);
    for( unsigned int i=0; i<k_MAX_topologies; i++ ){
        if(i==data) legend->AddEntry(h_TH1F_stacked[classification][type][i],Form("%s (%.2f)",legendLabel[i], h_TH1F_stacked[classification][type][i]->Integral(x0,xf)),"le");
        else legend->AddEntry(h_TH1F_stacked[classification][type][i],Form("%s (%.2f)",legendLabel[i], h_TH1F_stacked[classification][type][i]->Integral(x0,xf)),"f");
    }

    // --- create stacked histogram for all topologies but data
    THStack *hs = new THStack("temp_hs",plottitle);
    for( unsigned int i=0; i<k_MAX_topologies-1; i++) {
        h_TH1F_stacked[classification][type][i]->SetFillColor(colors[i]);
        h_TH1F_stacked[classification][type][i]->SetLineColor(lines[i]);
        h_TH1F_stacked[classification][type][i]->SetFillStyle(fills[i]);
        hs->Add(h_TH1F_stacked[classification][type][i]);
    }

    // --- check if there are entries for data so we can calculate the ratio
    if (h_TH1F_stacked[classification][type][data]->GetMaximum()>1){

        // --- make sure everything fits in the canvas
        hs->SetMaximum(1.1*h_TH1F_stacked[classification][type][data]->GetMaximum()); 

        // --- creates temporary plot as the sum of all stacked histograms
        TList *stackHists = hs->GetHists(); 
        TH1* h_pred_error = (TH1*)stackHists->At(0)->Clone();
        h_pred_error->Reset();
        for( int i=0; i<stackHists->GetSize(); i++ ) h_pred_error->Add((TH1*)stackHists->At(i));

        // --- calculates the ratio
        TRatioPlot* rp = new TRatioPlot(h_TH1F_stacked[classification][type][data],h_pred_error);
        rp->SetH1DrawOpt("AXIG");
        rp->SetH2DrawOpt("AXIG");

        std::vector<double> lines = {1.}; // lines to be drawn on the bottom pad
        rp->SetGridlines(lines);

        // --- bottom pad
        rp->Draw();
        rp->GetLowerRefGraph()->SetMaximum(2);
        rp->GetLowerRefGraph()->SetMinimum(0);
        rp->GetLowerRefGraph()->GetYaxis()->SetLabelSize(0.02);
        rp->GetLowerRefGraph()->GetYaxis()->SetTitleSize(0.03);
        rp->GetLowerRefYaxis()->SetTitle("data/(MC + EXT)");

        // --- top pad
        rp->GetUpperPad()->cd();
        hs->Draw("A,HIST");   // stacked histograms

        h_pred_error->SetFillColorAlpha(12, 0.15);
        h_pred_error->SetLineWidth(0);
        h_pred_error->Draw("SAME,A,E2"); // draw uncertainty for the stacked histograms

        h_TH1F_stacked[classification][type][data]->Draw("SAME,A,E"); // data

        legend->Draw("SAME");
    }

    c1->SaveAs(Form("%s", savename));   

}
/* --------------------------------------- */
void compare_2_1dhistograms(TCanvas *c1, TH1F* hist1, TH1F* hist2, string xlabel, string ylabel, string filename, bool apply_area_normalisation, string legend1, string legend2){

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

        TLatex text_area_normalise;
        text_area_normalise.SetTextColor(kPink);
        text_area_normalise.DrawLatexNDC(0.7048711,0.9214437,"Area Normalized");
    }

    double max = hist1->GetMaximum();
    if( hist2->GetMaximum() > max ) max = hist2->GetMaximum();
    hist1->GetYaxis()->SetRangeUser(0.,max*1.1);

    c1->SaveAs(filename.c_str());

}
/* --------------------------------------- */
void drawAllHists(){

  TCanvas* c1 = new TCanvas("c1", "c1",602,189,743,604);

  // --- stacked histograms
  for(unsigned int i=0; i<k_MAX_wirecell; i++){
    drawStackedHists(c1, i, k_Ereco, "; E_{reco} [MeV];Events / bin", Form("plots/Ereco_%s.pdf",classification_title[i]), 0, 25, Form("nueCC, %s", classification_title[i]));
  }

  TCanvas *c2 = new TCanvas();

  // --- energy spectrum
  compare_2_1dhistograms(c2, h_energy_true[k_fc][k_nue],h_energy_true[k_fc][k_nuebar],"True Neutrino Energy [MeV]","Entries / Bin","plots/nuEnergySpectrum_nue_nuebar_FC.pdf",false,"#nu_{e}","#bar{#nu}_{e}");
  compare_2_1dhistograms(c2, h_energy_true[k_pc][k_nue],h_energy_true[k_pc][k_nuebar],"True Neutrino Energy [MeV]","Entries / Bin","plots/nuEnergySpectrum_nue_nuebar_PC.pdf",false,"#nu_{e}","#bar{#nu}_{e}");
  compare_2_1dhistograms(c2, h_energy_true[k_fc][k_numu],h_energy_true[k_fc][k_numubar],"True Neutrino Energy [MeV]","Entries / Bin","plots/nuEnergySpectrum_numu_numubar_FC.pdf",false,"#nu_{#mu}","#bar{#nu}_{#mu}");
  compare_2_1dhistograms(c2, h_energy_true[k_pc][k_numu],h_energy_true[k_pc][k_numubar],"True Neutrino Energy [MeV]","Entries / Bin","plots/nuEnergySpectrum_numu_numubar_PC.pdf",false,"#nu_{#mu}","#bar{#nu}_{#mu}");
  
}
/* --------------------------------------- */

// =============================================== //
//                    FUNCTIONS                    //
// =============================================== //

bool inFidVolTrue(EvalInfo &e) {
  return ((e.truth_vtxX >= TPC_xi && e.truth_vtxX <= TPC_xf) &&
          (e.truth_vtxY >= TPC_yi && e.truth_vtxY <= TPC_yf) &&
          (e.truth_vtxZ >= TPC_zi && e.truth_vtxZ <= TPC_zf));
} // inFidVolReco
/* --------------------------------------- */
bool inFidVolRecoWCP(PFevalInfo &e) {
  return ((e.reco_nuvtxX >= TPC_xi && e.reco_nuvtxX <= TPC_xf) &&
          (e.reco_nuvtxY >= TPC_yi && e.reco_nuvtxY <= TPC_yf) &&
          (e.reco_nuvtxZ >= TPC_zi && e.reco_nuvtxZ <= TPC_zf));
} // inFidVolReco
/* --------------------------------------- */
float dist(float x1, float y1, float z1, float x2, float y2, float z2) {
  return  sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
} // dist
/* --------------------------------------- */
bool is_nueCC(TaggerInfo& tagger_info){
  bool flag = false;
  // default 7.0
  if (tagger_info.numu_cc_flag >=0 && tagger_info.nue_score > 7.0)
    flag = true;
  
  return flag;
} // is_nueCC
/* --------------------------------------- */
bool is_numuCC(TaggerInfo& tagger_info){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9)
    flag = true;
  
  return flag;
} // is_numuCC
/* --------------------------------------- */
bool is_cc_pi0(KineInfo& kine){
  bool flag = false;

  if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 ) && kine.kine_pio_energy_1 > 40 && kine.kine_pio_energy_2 > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174  && kine.kine_pio_mass > 22 && kine.kine_pio_mass < 300)
    flag = true;
  
  return flag;
} // is_cc_pi0
/* --------------------------------------- */
/* --------------------------------------- */
/* --------------------------------------- */
/* --------------------------------------- */
/* --------------------------------------- */
/* --------------------------------------- */
/* --------------------------------------- */
/* --------------------------------------- */

#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TStyle.h>

#include <vector>
#include <cmath>
#include <iostream>

#include "EXP_dNdy.h"

/*------------------------------------------
  工具函数：生成 TGraphErrors
------------------------------------------*/
TGraphErrors* makeGraphError(
    const std::vector<dNdyPoint>& data,
    const char* name
) {
    int n = data.size();
    auto* gr = new TGraphErrors(n);
    gr->SetName(name);
    gr->SetTitle("");

    for (int i = 0; i < n; ++i) {
        double y  = data[i].y;
        double dy = 0.5 * (data[i].y_high - data[i].y_low);

        double val = data[i].value;
        double err = std::sqrt(
            data[i].stat_err * data[i].stat_err +
            data[i].sys_err  * data[i].sys_err
        );

        gr->SetPoint(i, y, val);
        gr->SetPointError(i, dy, err);
    }
    return gr;
}

/*------------------------------------------
  工具函数：统一美化 Graph
------------------------------------------*/
void styleGraph(TGraphErrors* gr, int color, int marker) {
    gr->SetMarkerStyle(marker);
    gr->SetMarkerSize(2.0);
    gr->SetMarkerColor(color);
    gr->SetLineColor(color);
    gr->SetLineWidth(2);
}

/*------------------------------------------
  主函数
------------------------------------------*/
void drawExpData_scaled_dNdDy() {

    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");

    TFile* fout = new TFile("expData.root", "RECREATE");

    /* ---------- 生成 Graph ---------- */
    auto* gr_NA49_17p3 = makeGraphError(
        NA49_PbPb_17p3GeV_005_netB_05Npart_Dy,
        "NA49_17p3GeV"
    );

    auto* gr_BRAHMS_62 = makeGraphError(
        BRAHMS_AuAu_62p4GeV_005_netB_05Npart_Dy,
        "BRAHMS_62p4GeV"
    );

    auto* gr_BRAHMS_200 = makeGraphError(
        BRAHMS_AuAu_200GeV_010_netB_05Npart_Dy,
        "BRAHMS_200GeV"
    );

    auto* gr_STAR_62 = makeGraphError(
        STAR_AuAu_62p4GeV_005_netP_05Npart_Dy,
        "STAR_62p4GeV"
    );

    auto* gr_STAR_200 = makeGraphError(
        STAR_AuAu_200GeV_005_netP_05Npart_Dy,
        "STAR_200GeV"
    );

    /* ---------- 设置样式（颜色 + 形状） ---------- */
    styleGraph(gr_NA49_17p3,  kBlack,    20); // ●
    styleGraph(gr_BRAHMS_62,  kRed+1,    21); // ■
    styleGraph(gr_BRAHMS_200, kBlue+1,   22); // ▲
    styleGraph(gr_STAR_62,    kGreen+2,  33); // ◆
    styleGraph(gr_STAR_200,   kMagenta+2,29); // ★

    /* ---------- Canvas ---------- */
    auto* c1 = new TCanvas("c1", "Scaled net-baryon dN/dy", 1200, 800);
    c1->SetTicks();
    c1->SetLeftMargin(0.13);
    c1->SetBottomMargin(0.12);

    /* ---------- Frame（控制坐标轴） ---------- */
    TH1F* frame = c1->DrawFrame(
        -7.5, 0.0,   // x_min, y_min
         0.5, 0.5    // x_max, y_max（按你的数据可再调）
    );

    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("#Deltay=y-y_{beam}");
    frame->GetYaxis()->SetTitle("(N_{part}/2)dN^{#it{B-#bar{B}}}/d(y-y_{beam})");

    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetLabelSize(0.040);

    /* ---------- 画图 ---------- */
    gr_NA49_17p3 ->Draw("P SAME");
    gr_BRAHMS_62 ->Draw("P SAME");
    gr_BRAHMS_200->Draw("P SAME");
    gr_STAR_62   ->Draw("P SAME");
    gr_STAR_200  ->Draw("P SAME");

    /* ---------- Legend ---------- */
    auto* leg = new TLegend(0.15, 0.60, 0.58, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);

    leg->AddEntry(gr_NA49_17p3,  "NA49 Pb+Pb 17.3 GeV", "p");
    leg->AddEntry(gr_BRAHMS_62,  "BRAHMS Au+Au 62.4 GeV", "p");
    leg->AddEntry(gr_BRAHMS_200, "BRAHMS Au+Au 200 GeV", "p");
    leg->AddEntry(gr_STAR_62,    "STAR Au+Au 62.4 GeV", "p");
    leg->AddEntry(gr_STAR_200,   "STAR Au+Au 200 GeV", "p");

    leg->Draw();

    /* ---------- 写文件 ---------- */
    gr_NA49_17p3 ->Write();
    gr_BRAHMS_62 ->Write();
    gr_BRAHMS_200->Write();
    gr_STAR_62   ->Write();
    gr_STAR_200  ->Write();
    //frame->Write();
    c1->Write();

    fout->Write();
    fout->Close();
}

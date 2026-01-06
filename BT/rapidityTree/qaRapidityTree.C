void qaRapidityTree()
{
    // ===============================
    // 0. I/O
    // ===============================
    TFile* fin = TFile::Open("./rapidityTree/testtree_Au197pnHFB14Au197pnHFB14_62p4_alpha3.00_cent0_5.root");
    if (!fin || fin->IsZombie()) {
        std::cout << "❌ Cannot open file!" << std::endl;
        return;
    }

    TTree *t   = (TTree*)fin->Get("t");

    int    evt, Ncoll;
    double b;
    bool   isProton, isProjectile;
    double y_init, y_final, dy, alpha;

    t->SetBranchAddress("evt",          &evt);
    t->SetBranchAddress("Ncoll",        &Ncoll);
    t->SetBranchAddress("b",            &b);
    t->SetBranchAddress("isProton",     &isProton);
    t->SetBranchAddress("isProjectile", &isProjectile);
    t->SetBranchAddress("y_init",       &y_init);
    t->SetBranchAddress("y_final",      &y_final);
    t->SetBranchAddress("dy",            &dy);
    t->SetBranchAddress("alpha",        &alpha);

    double y_beam = 4.2;   // <<< 修改为你的能量

        TFile *fout = new TFile("qaRapidityTree.root","RECREATE");

    // ===============================
    // Level-0: Tree consistency QA
    // ===============================
    TH2F *h_yinit_yfinal =
        new TH2F("h_yinit_yfinal",
                 ";y_{init};y_{final}",
                 200,-y_beam-1,y_beam+1,
                 200,-y_beam-1,y_beam+1);

    TH2F *h_dy_check =
        new TH2F("h_dy_check",
                 ";y_{final}-y_{init};dy",
                 200,-10,10,
                 200,-10,10);

    TH1F *h_alpha =
        new TH1F("h_alpha",";#alpha",100,0,10);

    // ===============================
    // Level-1: Event geometry QA
    // ===============================
    TH1F *h_b =
        new TH1F("h_b",";b [fm]",100,0,20);

    TH1F *h_Ncoll =
        new TH1F("h_Ncoll",";N_{coll}",100,0,1400);

    TH2F *h_Ncoll_b =
        new TH2F("h_Ncoll_b",
                 ";b;N_{coll}",
                 100,0,20,
                 50,0,50);

    // ===============================
    // Level-2: Particle classification QA
    // ===============================
    TH1F *h_dNdy =
        new TH1F("h_dNdy",";y;counts",
                 500,-y_beam-12,y_beam+12);

    TH1F *h_dNdyA =
        new TH1F("h_dNdyA",";y;counts",
                 500,-y_beam-12,y_beam+12);

    TH1F *h_dNdyB =
        new TH1F("h_dNdyB",";y;counts",
                 500,-y_beam-12,y_beam+12);

    TH1F *h_dNdyAP =
        new TH1F("h_dNdyAP",";y;counts",
                 500,-y_beam-12,y_beam+12);

    TH1F *h_dNdyAN =
        new TH1F("h_dNdyAN",";y;counts",
                 500,-y_beam-12,y_beam+12);

    TH2F *h_y_isProjectile =
        new TH2F("h_y_isProjectile",
                 ";y;isProjectile",
                 200,-y_beam,y_beam,
                 2,-0.5,1.5);

    // ===============================
    // Level-3: Baryon stopping QA
    // ===============================
    TH1F *h_dy =
        new TH1F("h_dy",";#Delta y",200,-10,10);

    TH1F *h_dy_proj =
        new TH1F("h_dy_proj",";#Delta y (projectile)",200,-10,10);

    TProfile *p_dy_Ncoll =
        new TProfile("p_dy_Ncoll",
                     ";N_{coll};<#Delta y>",
                     30,0,30);

    // ===============================
    // Level-4: Global constraint QA
    // ===============================
    TH1F *h_baryon_per_evt =
        new TH1F("h_baryon_per_evt",
                 ";N_{baryon}/event",
                 400,0,400);

    // ===============================
    // Event loop
    // ===============================
    Long64_t nentries = t->GetEntries();

    int current_evt = -1;
    int baryon_count_evt = 0;

    for (Long64_t i = 0; i < nentries; ++i) {
        t->GetEntry(i);

        // ---- event change ----
        if (evt != current_evt) {
            if (current_evt >= 0) {
                h_baryon_per_evt->Fill(baryon_count_evt);
            }
            current_evt = evt;
            baryon_count_evt = 0;

            h_b->Fill(b);
            h_Ncoll->Fill(Ncoll);
            h_Ncoll_b->Fill(b, Ncoll);
        }

        baryon_count_evt++;

        // ---- Level-0 ----
        h_yinit_yfinal->Fill(y_init, y_final);
        h_dy_check->Fill(y_final - y_init, dy);
        h_alpha->Fill(alpha);

        // ---- Level-2 ----
        h_dNdy->Fill(y_final);
        h_y_isProjectile->Fill(y_final, isProjectile);

        if (isProjectile) {
            h_dNdyA->Fill(y_final);

            if (isProton) {
                h_dNdyAP->Fill(y_final);
            } else {
                h_dNdyAN->Fill(y_final);
            }
        } else {
            h_dNdyB->Fill(y_final);
        }

        // ---- Level-3 ----
        h_dy->Fill(dy);

        if (isProjectile) {
            h_dy_proj->Fill(dy);
            p_dy_Ncoll->Fill(Ncoll, dy);
        }
    }

    // last event
    h_baryon_per_evt->Fill(baryon_count_evt);

    // ===============================
    // Output
    // ===============================

    fout->Write();
    fout->Close();

    cout << "QA finished. Output: qaRapidityTree.root" << endl;
}

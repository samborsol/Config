//int colorArr[] = {kRed+1, kOrange+7, kSpring+4, kGreen+3, kAzure+1, kBlue+2, kViolet+5, kViolet-4, kMagenta, kMagenta+2};
int colorArr[] = {kRed+1, kGreen+3, kBlue+2, kMagenta, kMagenta+2};
//int markerArr[] = {kOpenCircle, kOpenSquare, kOpenStar, kOpenTriangleUp, kOpenDiamond, kOpenCross};
int markerArr[] = {kOpenCircle, kFullCircle, kOpenSquare, kFullSquare, kOpenDiamond, kFullDiamond, kOpenCross, kFullCross};
int ncolor = sizeof(colorArr)/sizeof(int);
int nmarker = sizeof(markerArr)/sizeof(int);

void SetHistStyleDefault(TH1 *h, int i, int j) {
//  if (j == kOpenStar || j == kFullStar || j ==kOpenDiamond || j ==kFullDiamond) h->SetMarkerSize(2.000);
//  else if (j == kOpenCross || j == kFullCross) h->SetMarkerSize(1.600);
  if (j>4 && j<6) h->SetMarkerSize(1.800);
  else if (j>=6) h->SetMarkerSize(1.500);
  else h->SetMarkerSize(1.200);

  if (ncolor>i) {
    h->SetMarkerColor(colorArr[i]);
    h->SetLineColor(colorArr[i]);
  } else {
    h->SetMarkerColor(colorArr[i%ncolor]);
    h->SetLineColor(colorArr[i%ncolor]);
  }
  if (nmarker>j) {
    h->SetMarkerStyle(markerArr[j]);
  } else {
    h->SetMarkerStyle(markerArr[j%nmarker]);
  }

  h->GetXaxis()->SetTitleSize(0.048);
  h->GetYaxis()->SetTitleSize(0.048);
  h->GetXaxis()->SetLabelSize(0.048);
  h->GetYaxis()->SetLabelSize(0.048);
}

void SetHistStyleDefault(TGraph *h, int i, int j) {
  if (j == kOpenStar || j == kFullStar || j ==kOpenDiamond || j ==kFullDiamond) h->SetMarkerSize(2.000);
  else if (j == kOpenCross || j == kFullCross) h->SetMarkerSize(1.600);
  else h->SetMarkerSize(1.200);

  if (ncolor>i) {
    h->SetMarkerColor(colorArr[i]);
    h->SetLineColor(colorArr[i]);
  } else {
    h->SetMarkerColor(colorArr[i%ncolor]);
    h->SetLineColor(colorArr[i%ncolor]);
  }
  if (nmarker>j) {
    h->SetMarkerStyle(markerArr[j]);
  } else {
    h->SetMarkerStyle(markerArr[j%nmarker]);
  }

  h->GetXaxis()->SetTitleSize(0.048);
  h->GetYaxis()->SetTitleSize(0.048);
  h->GetXaxis()->SetLabelSize(0.048);
  h->GetYaxis()->SetLabelSize(0.048);
}

void SetHistStyle(TGraph *h, int i, int j, double rmin, double rmax){
  h->SetMinimum(rmin);
  h->SetMaximum(rmax);
  SetHistStyleDefault(h, i, j);
}

void SetHistStyle(TH1 *h, int i, int j, double rmin, double rmax){
  h->GetYaxis()->SetRangeUser(rmin,rmax);
  SetHistStyleDefault(h, i, j);
}

void SetLegendStyle(TLegend* l) {
  l->SetFillColor(0);
  l->SetFillStyle(4000);
  l->SetBorderSize(0);
  l->SetMargin(0.15);
}

void SetStatBox(TPaveStats *p, double x1, double y1, double x2, double y2, int color) {
  p->SetX1NDC(x1);
  p->SetX2NDC(x2);
  p->SetY1NDC(y1);
  p->SetY2NDC(y2);
  p->SetTextColor(color);
  p->SetTextSize(0.035);
  p->SetTextFont(42);
  p->SetBorderSize(0);
}
void SetStatBox(TPaveText *p, double x1, double y1, double x2, double y2, int color) {
  p->SetX1NDC(x1);
  p->SetX2NDC(x2);
  p->SetY1NDC(y1);
  p->SetY2NDC(y2);
  p->SetTextColor(color);
  p->SetTextSize(0.035);
  p->SetTextFont(42);
  p->SetBorderSize(0);
}



Int_t npeaks = 10;

Double_t fpeaks(Double_t *x, Double_t *par) {
  Double_t result = par[0] + par[1]*x[0];
  
  for (Int_t p=0;p<npeaks;p++) {
    Double_t norm  = par[3*p+2];
    Double_t mean  = par[3*p+3];
    Double_t sigma = par[3*p+4];
    result += norm*TMath::Gaus(x[0],mean,sigma);
  }
  return result;
}

void calib_kromek_data(Int_t bfda_chan){
  //Bool_t bDebug=kFALSE;
  Bool_t bDebug = kTRUE;
  if(bDebug)cout<<"________________DEBUG MODE_________________"<<endl;
  //Run matlab in a batch mode to get detector data
  //Int_t bfda_chan=1;
  //gSystem->Setenv("BFDA_CHAN",Form("%d",bfda_chan));
  //TString tset_chan=Form("matlab -nodesktop -nosplash < bfda_dump.m > /dev/null",bfda_chan);
  //matlab < bfda_dump.m > /dev/null
  //gSystem->Exec(tset_chan.Data());

  Float_t lsize=0.2;
  Int_t np=6;
  Float_t slimit=75;
  Float_t slimite=30;
  Float_t cfactor=0.5;
  npeaks = TMath::Abs(np);
  Double_t par[3000];Double_t pare[3000];
  Int_t abin[4096];
  Int_t acount[4096];
  Float_t hvolt;
  Float_t temp;
  TString filename1=Form("det%d.dat",bfda_chan);
  TString DirPath="./";
  ifstream fin1(Form("%s/%s",DirPath.Data(),filename1.Data()));
  if(!fin1) {
    cout << "cannot open file: "<<filename1<<endl;
    return 1;
  }
  fin1>>hvolt>>temp;
  Int_t numf1=0;
  do{
    fin1>>abin[numf1]>>acount[numf1];
    numf1++;
  }while(fin1);
  numf1--;
  gStyle->SetOptStat(0);
  
  TCanvas *c2 = new TCanvas("Kromek Fit","bfda",10,10,1000,900);
  c2->Divide(1,3);
  //First Plot
  c2->cd(1);
  gPad->SetLeftMargin(0.1);
  gPad->SetBottomMargin(0.18);
  TH1F *h;
  TString confname=Form("Kromek_Detector_%d",bfda_chan);
  h=new TH1F(confname.Data(),confname.Data(),numf1,0,numf1);
  for(Int_t i=0;i<numf1;i++){
    Int_t binx=abin[i];
    Float_t content=acount[i];
    h->SetBinContent(binx,content);
  }
  h->Rebin(4);
  h->GetXaxis()->SetRange(18,800);
  h->Draw();

  TH1F *h2 = (TH1F*)h->Clone("h2");
  //Use TSpectrum to find the peak candidates
  TSpectrum *s = new TSpectrum(2*npeaks);
  Int_t nfound = s->Search(h,1.0,"",0.02);//sigma, threshold% to the highest peak
  //  if(bDebug)printf("Found %d candidate peaks to fit\n",nfound);
  //Estimate background using TSpectrum::Background
  TH1 *hb = s->Background(h,10,"same");
  if (np <0) return;
  
  TF1 *fline = new TF1("fline","pol1",100,2000);
  h->Fit("fline","qn");

  par[0] = fline->GetParameter(0);
  par[1] = fline->GetParameter(1);
  npeaks = 0;
  Double_t *xpeaks = s->GetPositionX();
  for(Int_t p=0;p<nfound;p++) {
    Float_t xp = xpeaks[p];
    Int_t bin = h->GetXaxis()->FindBin(xp);
    Float_t yp = h->GetBinContent(bin);
    if (yp-TMath::Sqrt(yp) < fline->Eval(xp)) continue;
    par[3*npeaks+2] = yp;
    par[3*npeaks+3] = xp;
    par[3*npeaks+4] = 30;
    npeaks++;
  }
  for(Int_t i=0;i<npeaks;i++){
    if(bDebug)cout<<i+1<<" "<<par[3*i+3]<<endl;
  }
 
  if(bDebug)printf("Found %d useful peaks to fit\n",npeaks);
  if(bDebug)printf("Now fitting: Be patient\n");

  
  TF1 *fitc = new TF1("fitc",fpeaks,100,3000,2+3*npeaks);
  fitc->SetLineColor(2);
  //we may have more than the default 25 parameters
  TVirtualFitter::Fitter(h2,3*npeaks);
  fitc->SetParameters(par);
  fitc->SetNpx(2000);
  h2->Fit("fitc","QR");
  h2->GetXaxis()->CenterTitle(1);h2->GetYaxis()->CenterTitle(1);
  h2->GetYaxis()->SetTitle("Counts");
  h2->GetXaxis()->SetTitle("ADC Bin");
  h2->GetYaxis()->SetNdivisions(205);
  h2->GetXaxis()->SetNdivisions(310);
  h2->GetXaxis()->SetTitleSize(lsize/3);
  h2->GetYaxis()->SetTitleSize(lsize/3);
  h2->GetXaxis()->SetLabelSize(lsize/3);
  h2->GetYaxis()->SetLabelSize(lsize/3);
  h2->GetYaxis()->SetTitleOffset(0.7);
 
  Float_t chanval[4];
  Float_t sigval[4];
  Float_t enval[4]={122.06065,661.657,1173.228,1332.492};
  Int_t parind[4];
  //Int_t parind[5]={6,3,18,12,15};
  //find the right indeces

  Int_t index[100];
  Float_t cval[100];Float_t csigm[100];
  for(Int_t i=0;i<npeaks;i++){
    Int_t indx=3*i+3;
    cval[i]= fitc->GetParameter(indx);
    csigm[i]= fitc->GetParameter(indx+1);
  }
  TMath::Sort(npeaks,cval,index,0);
  if(bDebug)printf("Sorted %d peaks after fit\n",npeaks);
  for(Int_t i=0;i<npeaks;i++){Int_t ind=index[i];
    if(bDebug)cout<<i+1<<" "<<cval[ind]<<endl;
  }
  
  for(Int_t j=0;j<4;j++){
    for(Int_t i=0;i<npeaks;i++){
      Int_t ind=index[i];
      if(cval[ind]*cfactor<enval[j]+slimit&&cval[ind]*cfactor>enval[j]-slimit&&csigm[ind]<70.){
	//cout<<"found: en:"<<enval[j]<<" chan:"<<cval[ind]<<endl;
	parind[j]=3*ind+3;continue;
      }
    }
  }
  for(Int_t i=0;i<4;i++){fitc->GetParameter(parind[i]);
    chanval[i]=fitc->GetParameter(parind[i]);
    sigval[i]=fitc->GetParameter(parind[i]+1);
    if(bDebug)cout<<"mean:"<<chanval[i]<<" sig:"<<sigval[i]<<endl;
  }

  TPaveText *doc = new TPaveText(0.7,0.5,0.9,0.9,"NDC");
  doc->SetBorderSize(1);
  doc->SetTextAlign(12);
  doc->SetTextFont(72);
  doc->SetFillColor(10);
  for(Int_t i=0;i<4;i++){
    TString ttext=Form("p%d=%7.1f   s%d=%5.1f",i+1,chanval[i],i+1,sigval[i]);
    doc->AddText(ttext.Data());
  }
  h2->GetListOfFunctions()->Add(doc);
  h2->DrawCopy();
  c2->Update();

 //Second Plot
  c2->cd(2);
  gPad->SetLeftMargin(0.1);
  gPad->SetBottomMargin(0.18);

  TGraph *mgr=new TGraph(4,chanval,enval);
  TH1F *hfr = gPad->DrawFrame(100,0,3000,1500.);
  hfr->GetXaxis()->CenterTitle(1);hfr->GetYaxis()->CenterTitle(1);
  hfr->GetYaxis()->SetTitle("Energy (keV)");
  hfr->GetXaxis()->SetTitle("ADC Bin");
  hfr->GetXaxis()->SetTitleSize(lsize/3);
  hfr->GetYaxis()->SetTitleSize(lsize/3);
  hfr->GetXaxis()->SetLabelSize(lsize/3);
  hfr->GetYaxis()->SetLabelSize(lsize/3);
  hfr->GetYaxis()->SetTitleOffset(0.7);
  hfr->GetYaxis()->SetNdivisions(205);
  hfr->GetXaxis()->SetNdivisions(310);
  mgr->SetMarkerColor(2); 
  mgr->SetMarkerStyle(20);
  mgr->SetMarkerSize(0.9);
  mgr->Draw("p");
  Double_t par4[3];
  TF1 *f4 = new TF1("f4","[0]+[1]*x+[2]*x*x",100,3000);
  f4->SetParameters(10,0.5,1e-5);
  f4->SetParNames("a","b","c");
  f4->SetLineColor(1);
  f4->SetLineWidth(1);
  mgr->Fit("f4","QR");
  f4->GetParameters(par4);
  f4->Draw("Same");
  TPaveText *doc1 = new TPaveText(0.75,0.2,0.9,0.55,"NDC");
  doc1->SetBorderSize(1);
  doc1->SetTextAlign(12);
  doc1->SetTextFont(72);
  doc1->SetFillColor(10);
  doc1->AddText(Form("a=%5.2f",par4[0]));
  doc1->AddText(Form("b=%5.4f",par4[1]));
  doc1->AddText(Form("c=%5.4e",par4[2]));
  hfr->GetListOfFunctions()->Add(doc1);
  c2->Update();
  
  //Save parameters to the text file
  TString outgfname=Form("det%d.par",bfda_chan);
  FILE *ogFile = fopen(outgfname.Data(),"w+"); 
  fprintf(ogFile,"%9.6f %9.6f %8.4e\n",par4[0],par4[1],par4[2]);
  fclose(ogFile);

  //Third Plot
  c2->cd(3);
  gPad->SetLeftMargin(0.1);
  gPad->SetBottomMargin(0.18);
  TH1F *h3;
  Int_t numf2=numf1; //number of bins in energy spectrum
  TString confname1=Form("BFDAen_Detector_%d",bfda_chan);
  h3=new TH1F(confname1.Data(),confname1.Data(),numf2,0,numf2);//We want energy range up to 2 MeV
  for(Int_t i=0;i<numf2;i++){
    Float_t en=f4->Eval(abin[i]);
    Int_t binx=h3->GetXaxis()->FindBin(en);
    Float_t content=acount[i];
    h3->SetBinContent(binx,content);
  }
  h3->Rebin(4);
  h3->GetXaxis()->SetRange(20,380);
  h3->Draw();
  
  //Fit the spectrum with known parameters
  TH1F *h4 = (TH1F*)h3->Clone("h4");
  //Use TSpectrum to find the peak candidates
  TSpectrum *s1 = new TSpectrum(2*npeaks);
  nfound = s1->Search(h3,1.0,"",0.02);//sigma, threshold% to the highest peak
  if(bDebug)printf("Found %d candidate peaks to fit\n",nfound);
  //Estimate background using TSpectrum::Background
  //TH1 *hb_1 = s1->Background(h3,20,"same");
  if (np <0) return;
  TF1 *fline1 = new TF1("fline1","pol1",0,1500);
  h3->Fit("fline1","qRn0");

  pare[0] = fline1->GetParameter(0);
  pare[1] = fline1->GetParameter(1);
  npeaks = 0;
  xpeaks = s1->GetPositionX();
  for(Int_t p=0;p<nfound;p++) {
    Float_t xp = xpeaks[p];
    Int_t bin = h4->GetXaxis()->FindBin(xp);
    Float_t yp = h4->GetBinContent(bin);
    if (yp-TMath::Sqrt(yp) < fline1->Eval(xp)) continue;
    pare[3*npeaks+2] = yp;
    pare[3*npeaks+3] = xp;
    pare[3*npeaks+4] = 20;
    npeaks++;
  }
  for(Int_t i=0;i<npeaks;i++){
    //    if(bDebug)cout<<i+1<<" "<<pare[3*i+3]<<endl;
  }
  
  if(bDebug)printf("Found %d useful peaks to fit\n",npeaks);
  if(bDebug)printf("Now fitting: Be patient\n");
  TF1 *fite = new TF1("fite",fpeaks,10,1500,2+3*npeaks);
  fite->SetLineColor(2);
  TVirtualFitter::Fitter(h4,3*npeaks);
  fite->SetParameters(pare);
  fite->SetNpx(500);
  h4->Fit("fite","QR");
  h4->GetXaxis()->CenterTitle(1);h4->GetYaxis()->CenterTitle(1);
  h4->GetYaxis()->SetTitle("Counts");
  h4->GetXaxis()->SetTitle("Energy (keV)");
  h4->GetYaxis()->SetNdivisions(205);
  h4->GetXaxis()->SetNdivisions(310);
  h4->GetXaxis()->SetTitleSize(lsize/3);
  h4->GetYaxis()->SetTitleSize(lsize/3);
  h4->GetXaxis()->SetLabelSize(lsize/3);
  h4->GetYaxis()->SetLabelSize(lsize/3);
  h4->GetYaxis()->SetTitleOffset(0.7);


  Float_t fenval[4];
  Float_t fsigval[4];
  Float_t sfenval[4];
  Float_t sfsigval[4];
  Int_t parinde[4];
  //Int_t parinde[4]={6,3,18,9,12};
  Float_t cvale[100];Float_t csigme[100];
  for(Int_t i=0;i<npeaks;i++){
    Int_t indx=3*i+3;
    cvale[i]= fite->GetParameter(indx);
    csigme[i]= fite->GetParameter(indx+1);
  }
  TMath::Sort(npeaks,cvale,index,0);
  if(bDebug)printf("Sorted %d peaks after fit\n",npeaks);
  for(Int_t i=0;i<npeaks;i++){Int_t ind=index[i];
    if(bDebug)cout<<i+1<<" "<<cvale[ind]<<endl;
  }
  
  for(Int_t j=0;j<4;j++){
    for(Int_t i=0;i<npeaks;i++){
      Int_t ind=index[i];
      if(cvale[ind]<enval[j]+slimite&&cvale[ind]>enval[j]-slimite&&csigme[ind]<70.){
	//cout<<"found: en:"<<enval[j]<<" chan:"<<cval[ind]<<endl;
	parinde[j]=3*ind+3;continue;
      }
    }
  }

  for(Int_t i=0;i<4;i++){
    fenval[i]=fite->GetParameter(parinde[i]);sfenval[i]=fite->GetParError(parinde[i]);
    fsigval[i]=fite->GetParameter(parinde[i]+1); sfsigval[i]=fite->GetParError(parinde[i]+1);
    cout<<"mean: found:"<<fenval[i]<<" expected:"<<f4->Eval(chanval[i])<<" sig:"<<fabs(fsigval[i])<<endl;
  }
  
  TPaveText *doc2 = new TPaveText(0.65,0.4,0.9,0.9,"NDC");
  doc2->SetBorderSize(1);
  doc2->SetTextAlign(12);
  doc2->SetTextFont(72);
  doc2->SetFillColor(10);
  doc2->AddText(Form("HV=%5.1f V     T=%4.2f^{o}C",hvolt,temp));
  for(Int_t i=0;i<4;i++){
    TString ttext=Form("p%d=%7.1f(%2.1f)   s%d=%5.1f(%2.1f)",i+1,fenval[i],sfenval[i],i+1,fabs(fsigval[i]),fabs(sfsigval[i]));
    doc2->AddText(ttext.Data());
  }
  h4->GetListOfFunctions()->Add(doc2);
  h4->DrawCopy();
  c2->Update();


  TString pngoutf=Form("%s.png",confname.Data());
  c2->Print(pngoutf.Data());
  
}

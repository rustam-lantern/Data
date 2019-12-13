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
Bool_t bDebug=kTRUE;
Bool_t bPlot=kTRUE;
Bool_t bPlot1=kTRUE;
//bPlot=kFALSE;
//bPlot1=kFALSE;
void calib_data_fusion(TString filename1,Int_t nnpar){
  if(bDebug)cout<<"________________DEBUG MODE_________________"<<endl;

  
  Int_t np=nnpar;
  Int_t ipeaks=nnpar;
  Float_t enval[10];
  //Float_t slimit=85;
  Float_t slimit=80;
  Float_t slimite=70;

  Float_t slimit_high;
  Float_t cfactor=0.5;
  //Float_t cfactor=3.5;
  //Float_t bin_edge=900;
  Float_t bin_edge=3000;
  npeaks = TMath::Abs(np);
  Double_t par[3000];
  Double_t parc[3000];
  Double_t pare[3000];
  Double_t parec[3000];
  Float_t chanval[10];
  Float_t sigval[10];
  Double_t lsize=0.2; Float_t small = 0.02;

  if(nnpar==4){
    enval[0]=122.136;
    enval[1]=661.66;
    enval[2]=1173.228;
    enval[3]=1332.492;
    slimit_high=70;
  }
  if(nnpar==6){
    enval[0]=122.136;
    enval[1]=511.0;
    enval[2]=661.66;
    enval[3]=1274.5;
    enval[4]=1461.0;
    enval[5]=2614.0;
    // enval[5]=1785.5;
    slimit_high=70;
  }
  if(nnpar==7){
    enval[0]=122.136;
    enval[1]=511.0;
    enval[2]=661.66;
    enval[3]=1274.5;
    enval[4]=1461.0;
    enval[5]=1785.5;
    enval[6]=2614.0;
    slimit_high=200;
  }
  Int_t parind[10];
  //Int_t ipeaks=7;
  //find indeces
  Int_t index[100];
  Float_t cval[100];Float_t csigm[100];

  Int_t abin[4096];
  Int_t acount[4096];
  TString ts;
  TString fname;
  Int_t ncounts;
  Float_t rtime;
  Float_t ltime;
  Double_t a2;
  Double_t a1;
  Double_t a0;

  //TString filename1=Form("irss_det%d.dat",irss_chan);
  TString DirPath="./";
  ifstream fin1(Form("%s/%s",DirPath.Data(),filename1.Data()));
  if(!fin1) {
    cout << "cannot open file: "<<filename1<<endl;
    return 1;
  }
  fin1>>ts>>fname>>ncounts>>rtime>>ltime>>a2>>a1>>a0;
  Int_t numf1=0;
  do{
    fin1>>abin[numf1]>>acount[numf1];
    numf1++;
  }while(fin1);
  numf1--;
  gStyle->SetOptStat(0);
  TString dname=filename1;

  TCanvas *c2 = new TCanvas("Data fusion Fit","Data_fusion",10,10,1000,900);
  c2->Divide(1,3);
  //First Plot
  c2->cd(1);
  gPad->SetLeftMargin(0.1);
  gPad->SetBottomMargin(0.18);
  TH1F *h;
  TString confname=filename1;
  h=new TH1F(confname.Data(),confname.Data(),numf1,0,numf1);
  for(Int_t i=0;i<numf1;i++){
    Int_t binx=abin[i]+1;
    Float_t content=acount[i];
    h->SetBinContent(binx,content);
  }
  //h->Rebin(4);
  h->GetXaxis()->SetRange(25,bin_edge);
  //  h->GetXaxis()->SetRange(150,1000);
  gPad->SetLogy();
  h->Draw();

  //TH1F *h2 = (TH1F*)h->Clone("h2");
  //Use TSpectrum to find the peak candidates
  // TSpectrum *s = new TSpectrum(2*npeaks);
  //Int_t nfound = s->Search(h,10,"",0.003);//sigma, threshold% to the highest peak
  //if(bDebug)printf("Found %d candidate peaks to fit\n",nfound);
  //Estimate background using TSpectrum::Background
  //TH1 *hb = s->Background(h,2,"same");
  //if (hb) c2->Update();
  //if (np <0) return;
  
  TH1F *h2 = (TH1F*)h->Clone("h2");
  //Use TSpectrum to find the peak candidates
   TSpectrum *s = new TSpectrum(2*npeaks);
  Int_t nfound = s->Search(h,10,"",0.003);//sigma, threshold% to the highest peak
  if(bDebug)printf("Found %d candidate peaks to fit\n",nfound);
  //Estimate background using TSpectrum::Background
  TH1 *hb = s->Background(h,15,"same");
  if (hb) c2->Update();
  if (np <0) return;

 
  TF1 *fline = new TF1("fline","pol1",25,bin_edge);
  h->Fit("fline","qn");
  
  par[0] = fline->GetParameter(0);
  par[1] = fline->GetParameter(1);
  npeaks = 0;
  Double_t *xpeaks = s->GetPositionX();
  for(Int_t p=0;p<nfound;p++) {
    Double_t xp = xpeaks[p];
    Int_t bin = h->GetXaxis()->FindBin(xp);
    Float_t yp = h->GetBinContent(bin);
    //if (yp-TMath::Sqrt(yp) < fline->Eval(xp)) continue;
    par[3*npeaks+2] = yp;
    par[3*npeaks+3] = xp;
    //par[3*npeaks+4] = 2;
    //if(xp>300){
    //  par[3*npeaks+4] = 10;
      //par[3*npeaks+4] = 8;
    //}
    if(xp>600&&xp<800){
      //cout<<"check:"<<xp<<" "<<yp<<" estimate: "<<yp-TMath::Sqrt(yp)<<" back: "<<fline->Eval(xp)<<endl;
      //par[3*npeaks+2] = 50;
      //par[3*npeaks+2] = fline->Eval(xp);
      //par[3*npeaks+3] = 689;
    }

   if(xp<50){
      par[3*npeaks+4] = 1.5;
    }
    else par[3*npeaks+4] = 6;

    if(xp>300){
      par[3*npeaks+4] = 8.0;
      // par[3*npeaks+4] = 8.0;
    }

    if(xp>300&&xp<360){
      par[3*npeaks+4] = 5.0;
      //par[3*npeaks+2] = fline->Eval(xp); 
 
      // par[3*npeaks+4] = 8.0;
    }
    if(xp>600){
      par[3*npeaks+4] = 12.0;
      // par[3*npeaks+4] = 8.0;
    }

    npeaks++;
  }
  for(Int_t i=0;i<npeaks;i++){
    parc[3*i+2]=par[3*i+2];
    parc[3*i+3]=par[3*i+3];
    parc[3*i+4]=0;//set to zero intentionally
    if(bDebug)cout<<i+1<<" starting pars:"<<par[3*i+3]<<endl;
  }
 
  if(bDebug)printf("Found %d useful peaks to fit\n",npeaks);
  if(bDebug)printf("Now fitting: Be patient\n");
  
  if(bPlot1){
  
  //for(Int_t i=0;i<npeaks;i++){
    //if(par[3*i+3]>600){
      //cout<<i<<" "<<par[3*i+2]<<" "<<par[3*i+3]<<" "<<par[3*i+4]<<endl;
  // }
  // }
  Double_t lsize=0.2; Float_t small = 0.02;
  
  TF1 *fitc = new TF1("fitc",fpeaks,25,bin_edge,2+3*npeaks);
  fitc->SetLineColor(2);
  //we may have more than the default 25 parameters
  TVirtualFitter::Fitter(h2,8+3*npeaks);
  fitc->SetParameters(par);
  fitc->SetNpx(500);
  
  
  for(Int_t i=0;i<npeaks;i++){
    //cout<<i<<" "<<par[3*i+2]<<" "<<par[3*i+3]<<" "<<par[3*i+4]<<endl;
    if(par[3*i+3]<50){
      //fitc->SetParLimits(3*i+4,1.0,4.0);
    }
    //else fitc->SetParLimits(3*i+4,3.0,30.0);
    // elseif(par[3*i+3]>300&&par[3*i+3]<400){
    //  fitc->SetParLimits(3*i+4,2,10);
    // }
    //else fitc->SetParLimits(3*i+4,2.0,20.0);
  }

  

  h2->Fit("fitc","QR");
  //h2->Fit("fitc","R");

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
  
  
 
  for(Int_t i=0;i<npeaks;i++){
    Int_t indx=3*i+3;
    cval[i]= fitc->GetParameter(indx);
    csigm[i]= fitc->GetParameter(indx+1);
  }
  TMath::Sort(npeaks,cval,index,0);
  if(bDebug)printf("Sorted %d peaks after fit\n",npeaks);

  for(Int_t i=0;i<npeaks;i++){Int_t ind=index[i];
    if(bDebug)cout<<"fitted: "<<i+1<<" "<<cval[ind]<<" "<<csigm[ind]<<endl;
  }


  for(Int_t j=0;j<ipeaks;j++){
    for(Int_t i=0;i<npeaks;i++){
      Int_t ind=index[i];
      Float_t slimit_check=slimit;
      if(cval[ind]>2000./cfactor){
	slimit_check=slimit_high;
      }
      //cout<<j<<" "<<j<<" "<<enval[j]<<" "<<cval[ind]<<" "<<cval[ind]*cfactor<<endl;
      if(cval[ind]*cfactor<enval[j]+slimit_check&&cval[ind]*cfactor>enval[j]-slimit_check&&csigm[ind]<70.){
	Float_t cfactor_calc=enval[j]/cval[ind];
	cout<<"found: en:"<<enval[j]<<" chan:"<<cval[ind]<<" ratio: "<<cfactor<<endl;
	parind[j]=3*ind+3;continue;
       }
    }
  }
  

  for(Int_t i=0;i<ipeaks;i++){
    fitc->GetParameter(parind[i]);
    chanval[i]=fitc->GetParameter(parind[i]);
    sigval[i]=fitc->GetParameter(parind[i]+1);
    if(i==0&&chanval[i]<0&&sigval[i]<1.0){
      cout<<"bad peak:"<<chanval[0]<<" "<<parc[3]<<endl;
      chanval[0]=parc[3];
      sigval[0]=parc[4];
    }
    if(i==0&&sigval[i]>10.0){
      cout<<"bad peak:"<<chanval[0]<<" "<<parc[3]<<endl;
      chanval[0]=parc[3];
      sigval[0]=parc[4];
    }
   
    if(bDebug)cout<<"mean:"<<chanval[i]<<" sig:"<<sigval[i]<<endl;
  }
  
  //chanval[0]=38.5;
  
  //sigval[0]=0;
    chanval[0]=par[3];

  TPaveText *doc = new TPaveText(0.7,0.5,0.9,0.9,"NDC");
  doc->SetBorderSize(1);
  doc->SetTextAlign(12);
  doc->SetTextFont(72);
  doc->SetFillColor(10);
  for(Int_t i=0;i<ipeaks;i++){
    TString ttext=Form("p%d=%7.1f   s%d=%5.1f",i+1,chanval[i],i+1,sigval[i]);
    doc->AddText(ttext.Data());
  }
  h2->GetListOfFunctions()->Add(doc);
  h2->DrawCopy();
  c2->Update();
}

  if(bPlot){
  //_____________________________________________________________________________________//
 //Second Plot
  c2->cd(2);
  gPad->SetLeftMargin(0.1);
  gPad->SetBottomMargin(0.18);

  TGraph *mgr=new TGraph(ipeaks,chanval,enval);
  mgr->Print();

  
  TH1F *hfr = gPad->DrawFrame(5,0,3000,2700.);
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
  TF1 *f4 = new TF1("f4","[0]+[1]*x+[2]*x*x",30,3000);
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
  //TString outgfname=Form("irss_det%d.par",irss_chan);
  //TString dname1=filename1;
  //if(dname1.Contains("-fe80"))dname1.Remove(dname1.Index("-fe80"));
  TString outgfname=Form("%s.par",dname.Data());
  FILE *ogFile = fopen(outgfname.Data(),"w+"); 
  //fprintf(ogFile,"%9.6f %9.6f %8.4e\n",par4[0],par4[1],par4[2]);
  fprintf(ogFile,"%8s %8.2f %8.6f %10.8f \n",dname.Data(),par4[0],par4[1],par4[2]);
  fprintf(ogFile,"        <Cal0>%10.8f</Cal0>\n",par4[0]);
  fprintf(ogFile,"        <Cal1>%10.8f</Cal1>\n",par4[1]);
  fprintf(ogFile,"        <Cal2>%10.8f</Cal2>\n",par4[2]);


  fclose(ogFile);
 
  //_________________________________________________________________________________________________//
 //Third Plot
  c2->cd(3);
  gPad->SetLeftMargin(0.1);
  gPad->SetBottomMargin(0.18);
  gPad->SetLogy();
  gPad->SetGrid();
  TH1F *h3;
  Int_t numf2=numf1; //number of bins in energy spectrum
  TString confname1=Form("%s energy",dname.Data());
  h3=new TH1F(confname1.Data(),confname1.Data(),200,0,2500);
  for(Int_t i=0;i<numf2;i++){
    Float_t en=f4->Eval(abin[i]+1);
    Int_t binx=h3->GetXaxis()->FindBin(en);
    Float_t content=acount[i];
    h3->SetBinContent(binx,content);
  }
  h3->Draw();
  
  
  //Fit the spectrum with known parameters
  TH1F *h4 = (TH1F*)h3->Clone("h4");
  //Use TSpectrum to find the peak candidates
  TSpectrum *s1 = new TSpectrum(2*npeaks);
  Int_t nfound = s1->Search(h3,2.0,"",0.002);//sigma, threshold% to the highest peak
  if(bDebug)printf("Found %d candidate peaks to fit\n",nfound);
  //Estimate background using TSpectrum::Background
  TH1 *hb = s1->Background(h3,2,"same");
  //if (hb) c2->Update();
  if (ipeaks <0) return;
  TF1 *fline1 = new TF1("fline1","pol1",50,2500);
  h3->Fit("fline1","qRn0");
  //  h4->Fit("fline1","qn");

  pare[0] = fline1->GetParameter(0);
  pare[1] = fline1->GetParameter(1);
  npeaks = 0;
  Double_t *speaks = s1->GetPositionX();
  for(Int_t p=0;p<nfound;p++) {
    Double_t xp = speaks[p];
    Int_t bin = h4->GetXaxis()->FindBin(xp);
    Float_t yp = h4->GetBinContent(bin);
    //if (yp-TMath::Sqrt(yp) < fline1->Eval(xp)) continue;
    pare[3*npeaks+2] = yp;
    pare[3*npeaks+3] = xp;
    pare[3*npeaks+4] = 10;
    if(xp>1200){
      pare[3*npeaks+4] = 30;
    }
     if(xp<200){
       //pare[3*npeaks+2] = 115;
       pare[3*npeaks+4] = 8;
    }
    npeaks++;
  }
  for(Int_t i=0;i<npeaks;i++){
    parec[3*i+2] =  pare[3*i+2];
    parec[3*i+3] =  pare[3*i+3];
    parec[3*i+4] =  0;

    if(pare[3*i+3]>2000){
      //cout<<i<<" "<<pare[3*i+2]<<" "<<pare[3*i+3]<<" "<<pare[3*i+4]<<endl;
    }
    if(bDebug)cout<<i+1<<" starting: "<<pare[3*i+3]<<" "<<pare[3*i+4]<<endl;
  }
  
  
  if(bDebug)printf("Found %d useful peaks to fit\n",npeaks);
  if(bDebug)printf("Now fitting: Be patient\n");
  TF1 *fite = new TF1("fite",fpeaks,50,3000,2+3*npeaks);
  fite->SetLineColor(2);
  TVirtualFitter::Fitter(h4,3*npeaks);
  fite->SetParameters(pare);
  fite->SetNpx(1000);
  
  //for(Int_t i=2;i<4;i++){
  //  fite->SetParLimits(3*(npeaks-i)+1,18.0,40.0);
  // }
  // for(Int_t i=0;i<2;i++){
    // fite->SetParLimits(3*(npeaks-i)+1,18.0,50.0);
  // }

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
  //h4->GetXaxis()->SetRange(5,2000);

  Float_t fenval[10];
  Float_t fsigval[10];
  Float_t sfenval[10];
  Float_t sfsigval[10];
  Int_t parinde[10];
  Float_t cvale[100];Float_t csigme[100];
  for(Int_t i=0;i<npeaks;i++){
    Int_t indx=3*i+3;
    cvale[i]= fite->GetParameter(indx);
    csigme[i]= fite->GetParameter(indx+1);
    
  }
  TMath::Sort(npeaks,cvale,index,0);
  if(bDebug)printf("Sorted %d peaks after fit\n",npeaks);
  for(Int_t i=0;i<npeaks;i++){
    Int_t ind=index[i];
    if(bDebug)cout<<i+1<<" "<<cvale[ind]<<" "<< csigme[ind]<<endl;
  }
  
  for(Int_t j=0;j<ipeaks;j++){
    for(Int_t i=0;i<npeaks;i++){
      Int_t ind=index[i];
      if(cvale[ind]<enval[j]+slimite&&cvale[ind]>enval[j]-slimite&&csigme[ind]<70.){
	cout<<"found2: en:"<<enval[j]<<" chan:"<<cvale[ind]<<endl;
	parinde[j]=3*ind+3;continue;
      }
    }
  }

  for(Int_t i=0;i<ipeaks;i++){
    fenval[i]=fite->GetParameter(parinde[i]);//sfenval[i]=fite->GetParError(parinde[i]);
    fsigval[i]=fite->GetParameter(parinde[i]+1); //sfsigval[i]=fite->GetParError(parinde[i]+1);
    if(i==0&&fsigval[i]<1.0){
      Int_t indx=3*i+3;
      fenval[i]=parec[indx];
      fsigval[i]=parec[indx+1];
    }
    cout<<"mean2: found2:"<<fenval[i]<<" expected:"<<f4->Eval(chanval[i])<<" sig:"<<fabs(fsigval[i])<<endl;
  }
  

  TPaveText *doc2 = new TPaveText(0.12,0.2,0.37,0.57,"NDC");
  doc2->SetBorderSize(1);
  doc2->SetTextAlign(12);
  doc2->SetTextFont(72);
  doc2->SetFillColor(10);
  
  for(Int_t i=0;i<ipeaks;i++){
    // TString ttext=Form("p%d=%7.1f(%2.1f)   s%d=%5.1f(%2.1f)",i+1,fenval[i],sfenval[i],i+1,fabs(fsigval[i]),fabs(sfsigval[i]));
    TString ttext=Form("p%d=%7.1f   s%d=%5.1f",i+1,fenval[i],i+1,fabs(fsigval[i]));
      if(fabs(fsigval[i])<100&&fabs(fsigval[i])>1){
      doc2->AddText(ttext.Data());
    }
  }
  h4->GetListOfFunctions()->Add(doc2);
  h4->DrawCopy();
  c2->Update();


  //TString psoutf=Form("%s/Det_%d_%s_%s.ps",OutDir.Data(),cap_chan,HV.Data(),idate.Data());

  //TString gifoutf=Form("%s_%s.gif",dname.Data(),GetTimeStamp().Data());
  TString gifoutf=Form("%s.gif",dname.Data());
  c2->Print(gifoutf.Data());


  }  

}


TString GetTimeStamp()
{
  TDatime dates; 
  TString stDate=dates.AsSQLString();
  TString sTime=stDate;
  TString sDate=stDate;
  Int_t spos=sDate.Index(":");
  sDate.Remove(spos-3);sTime.Remove(0,spos-2);
  sDate.ReplaceAll("-","_");stDate=sDate;
  spos=stDate.Index("_");
  stDate.Remove(spos);sDate.Remove(0,spos+1);
  sDate=Form("%s_%s",sDate.Data(),stDate.Data());
  sTime.ReplaceAll(":","_");
  TString timestr=Form("%s_%s",sDate.Data(),sTime.Data());
  return timestr.Data();
}


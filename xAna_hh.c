// example code to run Bulk Graviton->ZZ->ZlepZhad selections on electron-channel

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <TH1D.h>
#include <TFile.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TGraph.h>

using namespace std;
void xAna_hh(std::string inputFile){

  //get TTree from file ...
  TreeReader data(inputFile.data());

  Long64_t nTotal=0;
  Long64_t nPass[20]={0};

  TH1F* h_SD=new TH1F("h_SD","",100,0,200);
  TH1F* h_nVtx=new TH1F("h_nVtx","",100,0,100);
  Double_t nVtxk[71]; //nVtx of nPass[3]
  Double_t nVtxt[71]; //total nVtx after classified and changed into array form
  Double_t nVtxu[71]; //unit of nVtx
  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){

    if (jEntry % 10000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());

    data.GetEntry(jEntry);
    nTotal++;

    //0. has a good vertex
    Int_t nVtx        = data.GetInt("nVtx");
    if(nVtx<1)continue;
    nPass[0]++;
    h_nVtx->Fill(nVtx);
      //classify nVtx, turns into array form
      for (int j=1;j<=70;j++)
        {
            if ((j-0.5)<=nVtx && nVtx<=(j+0.5)) nVtxt[j]++;
        }
      
    int nFATJet         = data.GetInt("FATnJet");
    const int nJets=nFATJet;
    TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Float_t*  fatjetSDmass = data.GetPtrFloat("FATjetSDmass");
    Int_t*   nSubSoftDropJet = data.GetPtrInt("FATnSubSDJet");
    vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV", nFATJet);
    vector<float>   *subjetSDPx  =  data.GetPtrVectorFloat("FATsubjetSDPx", nFATJet);
    vector<float>   *subjetSDPy  =  data.GetPtrVectorFloat("FATsubjetSDPy", nFATJet);
    vector<float>   *subjetSDPz  =  data.GetPtrVectorFloat("FATsubjetSDPz", nFATJet);
    vector<float>   *subjetSDE   =  data.GetPtrVectorFloat("FATsubjetSDE", nFATJet);
    vector<bool>    &passFatJetLooseID = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    

    int nSubBTag[2]={0}; // check only the leading two fat jets 
    int nGoodFatJet=0;
    for(int ij=0; ij<nJets; ij++)
      {
    	
	TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
	if(thisJet->Pt()<30)continue;
	if(fabs(thisJet->Eta())>2.5)continue;
	if(fatjetSDmass[ij]<95 || fatjetSDmass[ij]>145)continue;
	if(!passFatJetLooseID[ij])continue;


      	for(int is=0; is < nSubSoftDropJet[ij]; is++){
	  if(subjetSDCSV[ij][is] < 0.605)continue;
	  if(nGoodFatJet<2)
	  nSubBTag[nGoodFatJet]++;
      	}

	nGoodFatJet++;

      }
  
    // if each fat jet has at least one subjet btag
    if(nSubBTag[0]>0 && nSubBTag[1]>0)nPass[1]++;

    // if one of the fat jets has at least two subjet btags
    if((nSubBTag[0]>1 && nSubBTag[1]>0) || 
       (nSubBTag[0]>0 && nSubBTag[1]>1))nPass[2]++;

    // if both fat jets have at least two subjet btags
    if(nSubBTag[0]>1 && nSubBTag[1]>1)
        {nPass[3]++;
            // the amount of nVtx
            for(int k=1; k<=70; k++)
            {
                if ((k-0.5)<=nVtx && nVtx<=(k+0.5)) nVtxk[k]++;
    }
  } // end of loop over entries
 
    Double_t efficiency[71];
    for (int l=1;l<=70;l++)
    {
        //nVtxu is the x-axis of the nVtx (unit)
        nVtxu[0]=0;
        nVtxu[l]=l;
        //counting efficiency, denominator cannot be negative
        efficiency[0]=0;
        if (nVtxt[l]!=0) efficiency[l]=(nVtxk[l]/nVtxt[l])*1.0;
    }
  
    //
    TGraph *efu = new TGraph(71,nVtxu,efficiency);
    efu->Draw("AC*");
  }
  std::cout << "nTotal    = " << nTotal << std::endl;
  for(int i=0;i<20;i++)
    if(nPass[i]>0)
      std::cout << "nPass[" << i << "]= " << nPass[i] << std::endl;
    
    
  
}

//TH1F::Divide//count uncertenty
//bynomion uncertenty
//poisson distribution average/error
//HW:p.19(rapidity in 2 coordinate)(eta=-1/2 tan(theta/2)), p.36, p.37
//natural unit
//Statistics(poisson guassion binomial distribution)
//error


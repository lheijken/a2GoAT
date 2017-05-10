#include "PPi0Dal.h"

PPi0Dal::PPi0Dal()
{ 
    CreateHistograms();
}

PPi0Dal::~PPi0Dal()
{
}

Bool_t	PPi0Dal::Init()
{
	cout << "Initialising physics analysis..." << endl;
	cout << "--------------------------------------------------" << endl << endl;

	if(!InitBackgroundCuts()) return kFALSE;
	if(!InitTargetMass()) return kFALSE;
	if(!InitTaggerChannelCuts()) return kFALSE;
	if(!InitTaggerScalers()) return kFALSE;
	cout << "--------------------------------------------------" << endl;
	return kTRUE;
}

Bool_t	PPi0Dal::Start()
{
    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();

    TraverseValidEvents();

    return kTRUE;
}

void	PPi0Dal::ProcessEvent()
{
    //- My own PID.. maybe
    //-- First CB tracks
    Track_neuCB.clear(); Track_neuCB.resize(0);
    Track_chaCB.clear(); Track_chaCB.resize(0);
    for(Int_t i=0; i<GetTracks()->GetNTracks(); i++){
        if(GetTracks()->HasCB(i)){
            FillTime(*GetTracks(),i,time_CB);
            //FillTimeCut(*GetTracks(),i,time_cut_CB);
            Double_t Ecry=GetTracks()->GetClusterEnergy(i);
            Double_t Epid=GetTracks()->GetVetoEnergy(i);
            Double_t Emwpc0=GetTracks()->GetMWPC0Energy(i);
            Double_t Emwpc1=GetTracks()->GetMWPC1Energy(i);
            if(Ecry > 0.0 && !(Epid > 0.0 || Emwpc0 > 0.0 || Emwpc1 > 0.0)){
                Track_neuCB.push_back(i);
                FillTime(*GetTracks(),i,time_CBneu);
            }
            if(Ecry > 0.0 && Epid > 0.0){
                EdE_CBPID->Fill(Ecry,Epid);
                EdEdcTh_CBPID->Fill(Ecry,Epid*sin(GetTracks()->GetTheta(i)));
                Track_chaCB.push_back(i);
            }
        }
    }
    //-- Then TAPS tracks
    Track_neuTA.clear(); Track_neuTA.resize(0);
    Track_chaTA.clear(); Track_chaTA.resize(0);
    for(Int_t i=0; i<GetTracks()->GetNTracks(); i++){
        if(GetTracks()->HasTAPS(i)){
            Double_t Ecry=GetTracks()->GetClusterEnergy(i);
            Double_t Evet=GetTracks()->GetVetoEnergy(i);
            if(Ecry > 0.0 && !(Evet > 0.0)){
                Track_neuTA.push_back(i);
            }
            if(Ecry > 0.0 && Evet > 0.0){
                EdE_TAVet->Fill(Ecry,Evet);
                Double_t thet=GetTracks()->GetTheta(i);
                EdEdcTh_TAVet->Fill(Ecry*sin(thet),Evet*sin(thet));
                Track_chaTA.push_back(i);
            }
        }
    }

    //- Just checking stuff
    Int_t ind1, ind2, ind3, ind4;
    Double_t PIDeCut = 1.05;
    if(Track_neuCB.size()==2){
        ind1 = Track_neuCB[0]; ind2 = Track_neuCB[1];
        IM2g->Fill((GetTracks()->GetVector(ind1)+GetTracks()->GetVector(ind2)).M());
    }
    if(Track_neuCB.size()==1 && Track_chaCB.size()==2){
        ind1 = Track_neuCB[0]; ind2 = Track_chaCB[0];
        ind3 = Track_chaCB[1];
        IMeeg[0]->Fill((GetTracks()->GetVector(ind1)+GetTracks()->GetVector(ind2)+GetTracks()->GetVector(ind3)).M());
        if(GetTracks()->GetVetoEnergy(ind2) < PIDeCut && GetTracks()->GetVetoEnergy(ind3) < PIDeCut){
            IMeeg[1]->Fill((GetTracks()->GetVector(ind1)+GetTracks()->GetVector(ind2)+GetTracks()->GetVector(ind3)).M());
        }
    }
    if(Track_neuCB.size()==1 && Track_chaCB.size()==3){
        ind1 = Track_neuCB[0]; ind2 = Track_chaCB[0];
        ind3 = Track_chaCB[1]; ind4 = Track_chaCB[2];
        Int_t inde[3]; Int_t Nre=0;
        for(Int_t i=0; i<3; i++){
            if(GetTracks()->GetVetoEnergy(Track_chaCB[i]) < PIDeCut){
                inde[Nre]=Track_chaCB[i];
                Nre++;
            }
        }
        if(Nre == 2){
            IMeeg[2]->Fill((GetTracks()->GetVector(ind1)+GetTracks()->GetVector(inde[0])+GetTracks()->GetVector(inde[1])).M());
        }
    }

}

void	PPi0Dal::ProcessScalerRead()
{
	// Fill Tagger Scalers
	FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal);
}

Bool_t	PPi0Dal::Write()
{
    // Write all GH1's and TObjects defined in this class
    return GTreeManager::Write();
}

void    PPi0Dal::CreateHistograms()
{
    time_CB 	= new GH1("time_CB", 	"time_CB", 	1401, -700, 700);
    time_CBneu 	= new GH1("time_CBneu", 	"time_CBneu", 	1401, -700, 700);
    time_cut_CB 	= new GH1("time_cut_CB", 	"time_cut_CB", 	1401, -700, 700);

    IM2g = new TH1F("IM2g","IM2g",1000,-0.5,999.5);
    IMeeg[0] = new TH1F("IMeeg1N2C","IMeeg1N2C",1000,-0.5,999.5);
    IMeeg[1] = new TH1F("IMeeg1N2Ccut","IMeeg1N2Ccut",1000,-0.5,999.5);
    IMeeg[2] = new TH1F("IMeeg1N3Ccut","IMeeg1N3Ccut",1000,-0.5,999.5);

    EdE_CBPID = new TH2F("EdE_CBPID","EdE_CBPID;E_{CB};dE_{PID}",200,0.,800.,200,0.,10.);
    EdEdcTh_CBPID = new TH2F("EdEdcTh_CBPID","EdEdcTh_CBPID;E_{CB};dE_{PID}*sin(#theta)",200,0.,800.,200,0.,10.);
    EdE_TAVet = new TH2F("EdE_TAVet","EdE_TAVet;E_{TAPS};dE_{Veto}",200,0.,800.,200,0.,10.);
    EdEdcTh_TAVet = new TH2F("EdEdcTh_TAVe","EdEdcTh_TAVe;E_{TAPS}*sin(#theta);dE_{Veto}*sin(#theta)",200,0.,800.,200,0.,10.);

    CombChTra[0] = new TH1F("CombChTra_1CBneu","CombChTra_1CBneu",10,0.,10.);
    CombChTra[0]->GetXaxis()->SetBinLabel(2,"3CB"); CombChTra[0]->GetXaxis()->SetBinLabel(3,"2CB(2e)1TA");
    CombChTra[0]->GetXaxis()->SetBinLabel(4,"2CB(1e1p)1TA"); CombChTra[0]->GetXaxis()->SetBinLabel(5,"1CB(e)2TA");
    CombChTra[0]->GetXaxis()->SetBinLabel(6,"1CB(p)2TA");

    TaggerAccScal = new TH1D("TaggerAccScal","TaggerAccScal",352,0,352);

}

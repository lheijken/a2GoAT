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
    Track_neuCB.clear();
    Track_neuCB.resize(0);
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
                EdE_CB->Fill(Ecry,Epid);
            }
        }
    }

    if(Track_neuCB.size()>1){
        for(Int_t i=0; i<(Track_neuCB.size()-1); i++){
            for(Int_t j=1; j<(Track_neuCB.size()); j++){
                IM2g->Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M());
            }
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

    EdE_CB = new TH2F("EdE_CB","EdE_CB;E_{CB};E_{PID}",100,0.,1000.,100,0.,10.);

    TaggerAccScal = new TH1D("TaggerAccScal","TaggerAccScal",352,0,352);

}

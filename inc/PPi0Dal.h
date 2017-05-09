#ifndef __PPi0Dal_h__
#define __PPi0Dal_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "PPhysics.h"

class	PPi0Dal  : public PPhysics
{
private:
    GH1*	time_CB;
    GH1*    time_CBneu;
    GH1*	time_cut_CB;

    TH1*    IM2g;

    TH2*    EdE_CBPID;
    TH2*    EdEdcTh_CBPID;

    TH1*	TaggerAccScal;

    //- Chosen reconstructed tracks
    std::vector<int> Track_neuCB;
    std::vector<int> Track_chaCB;

protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
    virtual Bool_t  Write();
    virtual void    CreateHistograms();
			
public:
    PPi0Dal();
    virtual ~PPi0Dal();
    virtual Bool_t  Init();

};
#endif


//
// export LD_LIBRARY_PATH=<bmf-library-installation-place>:${LD_LIBRARY_PATH}
//
// root -l
// root [] gSystem->Load("libbmf.so"); gSystem->AddIncludePath("-I../include");
// root [] .x mfield2eve.C+
//

#include <math.h>

#include <TEveManager.h>
#include <TEveArrow.h>
#include <TEveLine.h>

#include <TVector2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TPad.h>
#include <TLine.h>

#include <BeastMagneticField.h>

// 63126 = 510 x 126; [cm], [T];
#define _MAP_ZMIN_ (-500.0)
#define _MAP_ZDIM_     501

#define _MAP_RMIN_     0.0
#define _MAP_RDIM_     126

#define _COORD_CFF_    1.0
#define _FIELD_CFF_    1.0
#define _STEP_         2.0

// NB: prefer to use {z,r} 2D coordinate sequence in this code;
static double BR[_MAP_ZDIM_][_MAP_RDIM_], BZ[_MAP_ZDIM_][_MAP_RDIM_];

// Z -> in [cm] here, sorry;
#define _KICK_ZMIN_  160.0
#define _KICK_ZDIM_    100
#define _KICK_ZSTEP_   1.0

// Theta -> in [degree]; from 0 to 40 in 1 degree steps;
#define _KICK_TMIN_      0
#define _KICK_TDIM_     40
#define _KICK_TSTEP_   1.0

// Use this momentum for normalization; does not really matter;
#define _P0_          10.0

// ----------------------------------------------------------------------------------------------

void mfield2eve( void )
{
  double kicksigs[_KICK_TDIM_];

  gStyle->SetOptStat(0);

  TEveManager::Create();

  auto bmf = new BeastMagneticField("../data/mfield.4col.dat");

  {
    // Assume IP at (0,0); may check how beam size along Z axis spoils the picture?;
    double z0 = 0.0;

    // Up to eta~1 with a step of 1 degree;
    for(unsigned is=0; is<_KICK_TDIM_; is++) {
      double theta = is*_KICK_TSTEP_, slope = tan(theta*TMath::Pi()/180.);

      // Between 150cm and 250cm (?) with a step of 1cm;
      double accu = 0.0, avg = 0.0, avg_sqr = 0.0;
      for(unsigned iz=0; iz<_KICK_ZDIM_; iz++) {
	double z = _KICK_ZMIN_ + iz*_KICK_ZSTEP_, r = (z-z0)*slope, br, bz;

	int ret = bmf->GetFieldValue(r, z, br, bz); assert(!ret);
	TVector2 n(1.0, slope), b(bz, br);
	double norm = (b.Norm(n)).Mod();
	// Assume ~420mrad kick for 1 GeV/c and 1.3T*m B*dl (HERMES); may want to cross-check!;
	double kick = (420.0/_P0_) * sqrt(1.0+slope*slope)*_KICK_ZSTEP_ * b.Mod()*tan(b.DeltaPhi(n)) / (1.3*100.0);
	accu += kick;
	
	avg     += accu;
	avg_sqr += accu*accu;
      } //for iz

      avg /= _KICK_ZDIM_; avg_sqr /= _KICK_ZDIM_;

      double dsp = sqrt(avg_sqr - avg*avg);
      kicksigs[is] = dsp;
      printf("   %3d -> <kick> %7.2f [mrad] -> sigma %7.2f [mrad]\n", is, avg, dsp); 
    } //for is
  }

  {
    double qstep = 1.0, alen = qstep*0.8;

    for(unsigned ir=0; ir<25; ir++) {
      unsigned counter = 0;
      double rr = ir*5.0, zz = 0.0;
      
      for( ; ; ) {
	double B[2];

	double bz, br;
	// NB: prefer to regularize GetField() on Z-axis;
	if (!bmf->GetFieldValue(ir ? rr : 0.0, zz, br, bz)) break;
	
	double norm = sqrt(br*br+bz*bz);
	
	TEveArrow *ea = new TEveArrow(0.0, alen*br/norm, alen*bz/norm, 0.0, rr, zz);
	int color = (zz >= _KICK_ZMIN_ && zz <= _KICK_ZMIN_ + _KICK_ZDIM_*_KICK_ZSTEP_) ? kRed : kOrange;
	ea->SetMainColor(color);
	ea->SetPickable(kFALSE);
	ea->SetTubeR(1.50);
	ea->SetConeR(1.00);
	ea->SetConeL(1.00);
	gEve->AddElement(ea);

	zz += (bz/norm)*qstep; rr += (br/norm)*qstep;
	
	if (zz > 300.0) break;
	if (counter++ > 1000) break;
      } //for inf
    } //for ir

    for(unsigned iq=0; iq<5; iq++) {
      double th = (iq+1)*5*TMath::Pi()/180.;

      auto pt = new TEveLine(2);
      pt->SetMainColor(kWhite);
      pt->SetPickable(kFALSE);
      pt->SetPoint(0,   0., 0.,   0.);
      pt->SetPoint(1,   0., 300.*tan(th), 300.);
      gEve->AddElement(pt);
    } //for iq
  }

  gEve->FullRedraw3D(kTRUE);

  {
    TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 800, 500);
    TH1D *dspH = new TH1D("", "", _KICK_TDIM_, 0.0, _KICK_TDIM_*_KICK_TSTEP_);

    double p0 = 30.0;

    // Up to eta~1 with a step of 1 degree;
    for(unsigned is=0; is<_KICK_TDIM_; is++) {
      double theta = is*_KICK_TSTEP_;
      int ibin = (int)floor((theta-_KICK_TMIN_)/_KICK_TSTEP_);

      if (ibin >= 0 && ibin < _KICK_TDIM_)
	dspH->SetBinContent(is+1, (_P0_/p0)*kicksigs[ibin]);
    } //for is

    dspH->SetMinimum( 0.0);
    dspH->SetMaximum(10.0);
    dspH->GetXaxis()->SetTitle("Scattered hadron polar angle, [degree]");
    dspH->GetXaxis()->SetTitleFont(52);
    dspH->GetXaxis()->SetTitleSize(0.05);
    dspH->GetXaxis()->SetLabelFont(52);
    dspH->GetXaxis()->SetLabelSize(0.04);
    dspH->GetXaxis()->SetTitleOffset(0.90);
    dspH->GetYaxis()->SetTitle("Direction variation, [mrad]");
    dspH->GetYaxis()->SetTitleFont(52);
    dspH->GetYaxis()->SetTitleSize(0.05);
    dspH->GetYaxis()->SetLabelFont(52);
    dspH->GetYaxis()->SetLabelSize(0.04);
    dspH->GetYaxis()->SetTitleOffset(0.80);
    gPad->SetGrid();
    dspH->SetLineWidth(2);
    dspH->SetFillColor(kGreen+1);
    dspH->Draw();

    TLine *line = new TLine(25., 0., 25., 10.);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->SetLineStyle(5);
    line->Draw();
  }
} // mfield2eve()

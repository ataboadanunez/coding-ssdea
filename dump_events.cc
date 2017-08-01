#include <iostream>
#include <fstream>
#include "IoSd.h"
#include "IoSdData.h"
#include "Ec.h"

#include "TFile.h"
#include "TH1F.h"

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
// Reads a CDAS sd_* file and dumps some information from UUB
//
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

using namespace std;


/**
 * \brief Export Single Histogram into ASCII file
 */
inline
Bool_t
SingleExportAscii(TH1F* hist, TString filename, TString folder = "./", TString separator = "  ")
{
  Int_t i, j;
  Double_t xcenter, xwidth;
  Bool_t success = kFALSE;
  filename = folder + filename;
  ofstream file_out(filename);
  file_out << "# Output " << hist->ClassName() << ": " << hist->GetName() << " (" << hist->GetTitle() << ")\n";
  if (hist->GetDimension() == 1) {
    file_out << "# BinCenter" << separator << "Content" << separator << "BinHalfWidth" << separator << "Error\n";
    for (i = 1; i <= hist->GetNbinsX(); i++)
      file_out << hist->GetBinCenter(i) << separator << hist->GetBinContent(i) << separator << hist->GetBinWidth(i) / 2 << separator << hist->GetBinError(i) << endl;
    if (i > 1)
      success = kTRUE;
  } else if (hist->GetDimension() == 2) {
    file_out << "# xBinCenter" << separator << "yBinCenter" << separator << "Content" << separator << "xBinHalfWidth" << separator << "yBinHalfWidth" << separator << "Error" << endl;
    for (i = 1; i <= hist->GetNbinsX(); i++) {
      xcenter = hist->GetXaxis()->GetBinCenter(i);
      xwidth = hist->GetXaxis()->GetBinWidth(i) / 2;
      for (j = 1; j <= hist->GetNbinsY(); j++)
        file_out << xcenter << separator << hist->GetYaxis()->GetBinCenter(j) << separator << hist->GetBinContent(i, j) << separator << xwidth << separator << hist->GetYaxis()->GetBinWidth(j) / 2 << separator << hist->GetBinError(i, j) << endl;
      if (j > 1)
        file_out << endl; // produce a blank line after each set of Ybins for a certain Xbin. Gnuplot likes this.
    }
    if (i > 1)
      success = kTRUE;
  }
  file_out.close();
  if (success == kTRUE)
    cout << "*** TRemiHistExport: Histogram " << hist->GetName() << " written to " << filename << endl;
  return success;
}


int main(int argc, char* argv[])
{
  if (argc == 1) {
    cout << "Will dump UUB traces into files of the structure eventid_statid.txt" << endl;
    cout << "Usage: " << argv[0] << " <files>" << endl;
    return 1;
  }

  // options of what will be dumped
  const bool dump_traces = true;
  const bool dump_histos = true;

  IoSd input(argc - 1, argv + 1); // open all files

  int IsUUB[2000];
  int nbuub = 0;

  for (int i = 0; i < 2000; i++)
    IsUUB[i] = 0;

  int error = 0;
  int nevt = 0;
  const int ntot = input.LastEvent() - input.FirstEvent();

  // get first event id in current files
  unsigned int firstid;
  unsigned int firstgps;
  for (EventPos pos = input.FirstEvent(); pos < input.LastEvent(); pos = input.NextEvent()) { // loop on all events
    TEcEvent event(pos);
    firstid = event.Id;
    firstgps = event.trigger().second();
    break;
  }

  ostringstream fnoutuub;
  ostringstream fnoutold;

  fnoutuub  << "./eventinfo_uub_" << firstgps << "_" << firstid << ".txt";
  ofstream evout_uub(fnoutuub.str());

  fnoutold  << "./eventinfo_old_" << firstgps << "_" << firstid << ".txt";
  ofstream evout_old(fnoutold.str());

  for (EventPos pos = input.FirstEvent(); pos < input.LastEvent(); pos = input.NextEvent()) { // loop on all events

    TEcEvent event(pos);

    unsigned int gpssecond = event.trigger().second();
    unsigned int eventid = event.Id;

    for (unsigned int i = 0 ; i < event.fCalibStations.size(); i++) { // loop on all fCalibstations

      unsigned int statid = event.fCalibStations[i].Id;

      // station is a UUB
      if (event.fCalibStations[i].IsUUB) {

        cout << "# Event " << eventid << " Station " << statid << " is a UUB." << endl;
        cout << "# Data error " << event.fCalibStations[i].Error << " (256 means has FADC traces)" << endl;

        if (event.fCalibStations[i].Error == 256) { // error codes for UUB are usual ones + 256: 0+256 for no error

          evout_uub << gpssecond << " " << eventid << " " << statid << " ";

          for (unsigned int l = 0; l < kIoSd::NPMT; l++) {

            // evout_uub << event.fCalibStations[i].Calib->VemPeak[l] << " ";
            // evout_uub << event.fCalibStations[i].Calib->VemCharge[l] << " ";

            evout_uub << event.fCalibStations[i].fPmt[l].fCalibratedState << " ";
            evout_uub << event.fCalibStations[i].fPmt[l].fHighGainSat << " ";
            evout_uub << event.fCalibStations[i].fPmt[l].fLowGainSat << " ";
            evout_uub << event.fCalibStations[i].fPmt[l].fVemPeak << " " << event.fCalibStations[i].fPmt[l].fVemCharge << " ";
            evout_uub << event.fCalibStations[i].fPmt[l].fPeakInVEM << " " << event.fCalibStations[i].fPmt[l].fSigInVEM << " ";

          }

          // evout_uub << event.fCalibStations[i].SigOverPeak() << " " << event.fCalibStations[i].fStartBin;
          evout_uub << "\n";

          if (dump_traces) {

            ostringstream fnout;
            fnout  << "./eventdumps/" << gpssecond << "_" << eventid << "_" << statid << ".txt";

            ofstream fout(fnout.str());

            cout << "Writing out to file " << fnout.str() << endl;

            // write out traces
            for (int k = 0; k < event.fCalibStations[i].UFadc->NSample; k++) {
              for (int j = 0; j < 5; j++)
                fout << event.fCalibStations[i].UFadc->GetValue(j, 0, k) << "  " << event.fCalibStations[i].UFadc->GetValue(j, 1, k) << "  ";
              fout << "\n";
            }

            fout.close();

          }

          if (dump_histos) {

            ostringstream calibout;

            // write out calibration histograms for each PMT
            for (unsigned int j = 0; j < 4; j++) {

              if (event.fCalibStations[i].HCharge(j) == NULL)
                break;

              calibout.str("");
              calibout.clear();

              calibout << gpssecond << "_" << eventid << "_" << statid << "_" << j << ".calib";
              SingleExportAscii(event.fCalibStations[i].HCharge(j),
                                calibout.str(), "./calibrations/");
            }

          }

        }

        if (statid < 2000 && !IsUUB[statid]) {
          IsUUB[statid] = 1;
          nbuub++;
        }

      }

      if ((statid == 1764) || (statid == 1739) || (statid == 1733)) {

        evout_old << gpssecond << " " << eventid << " " << statid << " ";

        for (unsigned int l = 0; l < kIoSd::NPMT; l++) {

          // evout_old << event.fCalibStations[i].Calib->VemPeak[l] << " ";
          // evout_old << event.fCalibStations[i].Calib->VemCharge[l] << " ";

          evout_old << event.fCalibStations[i].fPmt[l].fCalibratedState << " ";
          evout_old << event.fCalibStations[i].fPmt[l].fHighGainSat << " ";
          evout_old << event.fCalibStations[i].fPmt[l].fLowGainSat << " ";
          evout_old << event.fCalibStations[i].fPmt[l].fVemPeak << " " << event.fCalibStations[i].fPmt[l].fVemCharge << " ";
          evout_old << event.fCalibStations[i].fPmt[l].fPeakInVEM << " " << event.fCalibStations[i].fPmt[l].fSigInVEM << " ";

        }

        // evout_old << event.fCalibStations[i].SigOverPeak() << " " << event.fCalibStations[i].fStartBin;
        evout_old << "\n";

        // if (dump_traces) {

        //   ostringstream fnout;
        //   fnout  << "./eventdumps/" << gpssecond << "_" << eventid << "_" << statid << ".txt";

        //   ofstream fout(fnout.str());

        //   cout << "Writing out to file " << fnout.str() << endl;

        //   // write out traces
        //   for (int k = 0; k < event.fCalibStations[i].UFadc->NSample; k++) {
        //     for (int j = 0; j < 5; j++)
        //       fout << event.fCalibStations[i].UFadc->GetValue(j, 0, k) << "  " << event.fCalibStations[i].UFadc->GetValue(j, 1, k) << "  ";
        //     fout << "\n";
        //   }

        //   fout.close();

        // }

        if (dump_histos) {

          ostringstream calibout;

          // write out calibration histograms for each PMT
          for (unsigned int j = 0; j < 3; j++) {

            if (event.fCalibStations[i].HCharge(j) == NULL)
              break;

            calibout.str("");
            calibout.clear();

            calibout << gpssecond << "_" << eventid << "_" << statid << "_" << j << ".calib";
            SingleExportAscii(event.fCalibStations[i].HCharge(j),
                              calibout.str(), "./calibrations_old/");
          }

        }

      }
    }

    nevt++;
    if (nevt % 100 == 0) cerr << "# read " << nevt << " events out of " << ntot << endl;

  }

  evout_old.close();
  evout_uub.close();

  input.Close();
  cerr << "# Found " << nbuub << " stations with UUB";
  if (nbuub) cerr << ":";
  cerr << endl;
  for (int i = 0; i < 2000; i++)
    if (IsUUB[i]) cerr << "#   " << i << endl;
  return 0;
};

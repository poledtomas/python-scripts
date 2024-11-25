#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
double with_gap_calculation(double QxA, double QyA, double QxB, double QyB, int MA, int MB) {
    double numerator = QxA * QxB + QyA * QyB;
    double denominator = static_cast<double>(MA) * MB;
    return numerator / denominator;
}
void vn_SP_pT_pid(const char* direct, int NoF, int NoE, double order, int pid, string centrality)
{
  cout << endl << endl << "Calculating v_" << (int)order << "{2} (pT)..." << endl;
  cout << "Processing events from directory: " << direct << endl;

  // Cuts
  const double etaCut = 0.8;
  const double ptMinCut = 0.9;
  const double ptMaxCut = 4.5;
  const int ptBins = 12;
  const int eventStep = 1; // number of events in one super-event (to increase statistics)
  const int NS = 10; // number of subsamples to estimate the error
  const double dpt = (ptMaxCut-ptMinCut)/ptBins;

  // Maximum number of particles - length of arrays
  //const int NP = 60000;

  double px, py, pz, m, E;
  double ele;
  int id;
  int npart, nevents = 0;
  int npar, maxpar;

  // arrays for calculation v_n
  double vn2[ptBins] = {0.0}, vn4[ptBins] = {0.0};
  double dn2[ptBins] = {0.0}, dn2_nom[ptBins] = {0.0}, dn2_denom[ptBins] = {0.0};
  double dn4[ptBins] = {0.0}, dn4_nom[ptBins] = {0.0}, dn4_denom[ptBins] = {0.0};
  double cn2 = 0.0, cn2_nom = 0.0, cn2_denom = 0.0;
  double cn4 = 0.0, cn4_nom = 0.0, cn4_denom = 0.0;
  double vn2_err[ptBins] = {0.0}, vn4_err[ptBins] = {0.0};
  std::pair<double, double> eta_range_A = {2.8, 5.1};
std::pair<double, double> eta_range_C = {-0.8, 0.8};
std::pair<double, double> eta_range_B = {-3.7, -1.7};

  // Loop over files
  for (int ifls = 1; ifls < NoF + 1; ifls++) {
    FILE *infile;
    char line[500];
    char delims[] = " ,\n\t";
    char *strtokresult = NULL;
    char *pars[12];

    string filename;
    string buffer;
    filename.append(direct);
    filename.append(to_string(ifls));
    filename.append("_fin.oscar");

    infile = fopen(filename.c_str(), "r");
    if (infile == NULL)
    {
      cout << "Warning: Missing file #" << ifls << endl;
    }
    else
    {
      fgets(line,500,infile);
      fgets(line,500,infile);
      fgets(line,500,infile);
      // Loop over events in file
      for (int iev = 0; iev < NoE; iev += eventStep)
      {
        int MA =0, MB =0;
        int RFP = 0, POI[ptBins] = {0}; // RFP = M (all particles), POI = m (particles in given ptBin)
        double cumulant2 = 0.0, diff_cumulant2[ptBins] = {0.0};
        double cumulant4 = 0.0, diff_cumulant4[ptBins] = {0.0};
        double QxA = 0.0, QyA = 0.0;
        double QxB = 0.0, QyB = 0.0;
        double ux[ptBins] = {0.0}, uy[ptBins] = {0.0};
       
        complex<double> Q, Q2, q[ptBins], q2[ptBins];

        for (int k = 0; k < eventStep; k++) 
        {
          fgets(line,500,infile);
          strtokresult = strtok(line, delims);
          npar = 0;
          maxpar = 5;
          while( (strtokresult != NULL) && (npar < maxpar) )
          {
            pars[npar]= strtokresult;
            strtokresult = strtok( NULL, delims );
            npar += 1;
          }
          int npart = atoi(pars[4]);
          if (npart > 0)
          {
            nevents++;
            // Loop over particles
            for (int i = 0; i < npart; i++)
            {
              fgets(line,500,infile);
              strtokresult = strtok(line, delims);
              npar = 0;
              maxpar = 12;
              while( (strtokresult != NULL) && (npar < maxpar) )
              {
                pars[npar]= strtokresult;
                strtokresult = strtok( NULL, delims );
                npar += 1;
              }
              m  = atof(pars[4]);
              E  = atof(pars[5]);
              px = atof(pars[6]);
              py = atof(pars[7]);
              pz = atof(pars[8]);
              id = atoi(pars[9]);

              float pt = sqrt(px*px+py*py);
              float pabs = sqrt(px*px+py*py+pz*pz);
              float eta = 0.5*log((pabs+pz)/(pabs-pz));
              float phi = atan2(py, px);
              float rap = 0.5*log((E+pz)/(E-pz));

              if (pt>ptMinCut && pt<ptMaxCut){
                int ptBin = (pt-ptMinCut)/dpt;
                if (eta > 2.8 && eta < 5.1 && abs(id) != pid) {
                  QxA += std::cos(order * phi);
                  QyA += std::sin(order * phi);
                  MA += 1;
                } 
              else if (eta > -3.7 && eta < -1.7 && abs(id) != pid) {
                  QxB += std::cos(order * phi);
                  QyB += std::sin(order * phi);
                  MB += 1;
                } 
              else if (abs(id) == pid && eta > -0.8 && eta < 0.8) {
                ux[ptBin] += std::cos(order * phi);
                uy[ptBin] += std::sin(order * phi);
                POI[ptBin] += 1;
              }

              }

              
            }
          }
          fgets(line,500,infile);
        }

        if (MA+MB > 0)
        {
  
          cumulant2 =  with_gap_calculation(QxA, QyA, QxB, QyB, MA, MB);
          
          cn2_nom += MA*MB*cumulant2;
          cn2_denom += MA*MB;
          
          for (int ipt = 0; ipt < ptBins; ipt++)
          {
            if (POI[ipt] > 0)
            {
              RFP=MA +MB;
              Q = complex<double>(QxA+QxB, QyA+QyB);
              q[ipt] = complex<double>(ux[ipt], uy[ipt]);
              diff_cumulant2[ipt] = real(q[ipt]*conj(Q)-complex<double>(POI[ipt],0))/(POI[ipt]*RFP-POI[ipt]);
              
              dn2_nom[ipt] += (POI[ipt]*RFP - POI[ipt]) * diff_cumulant2[ipt];
              dn2_denom[ipt] += POI[ipt]*RFP - POI[ipt];

            }
          }
        }
      }
    }
    if ((ifls)%((int)NoF/20) == 0) cout << (int)100*ifls/NoF << "%" << endl;
  }

  cout << "Total number of events: " << nevents << endl;

  // final calculation of v_n
  cn2 = cn2_nom / cn2_denom;

  // Write results into the text file (append)
  ofstream fout;
  string file="vn_SP_pT_pair_";
  string dat= ".dat";
  string mezera = "_";
  
  string hadron = to_string(pid);
  string filename = file+hadron+mezera+centrality+dat;
  fout.open(filename, ofstream::app);

  fout << direct<<"\t"<<pid<<endl;
  fout << nevents<<endl;

  for (int ipt = 0; ipt < ptBins; ipt++)
  {
    dn2[ipt] = dn2_nom[ipt] / dn2_denom[ipt];
    vn2[ipt] = dn2[ipt] / sqrt(cn2);


    double ptBin = (ipt+0.5)*dpt + ptMinCut;
  

    cout << ptBin << "\t" << vn2[ipt]  << endl;
    fout << ptBin << "\t" << vn2[ipt]  << endl;
  }
  fout << endl;
  fout.close();

  cout << "Results have been written to:  "<<filename << endl;
  
}
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

void count(const char* direct, int NoF, int NoE, int pid, string centrality)
{
  cout << endl << endl << "Calculating v_" << "{2} (pT)..." << endl;
  cout << "Processing events from directory: " << direct << endl;

  // Cuts
  const double etaCut = 0.8;
  const double ptMinCut = 0.9;
  const double ptMaxCut = 4.5;

  // Maximum number of particles - length of arrays
  //const int NP = 60000;

  double px, py, pz, m, E;
  double ele;
  int id;
  int npart, nevents = 0;
  int npar, maxpar;
  int count_particle =0;
  int eventStep=1;

  
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
            

              if (pt>ptMinCut && pt<ptMaxCut && id ==pid){
                count_particle +=1;
              
              }

              
            }
          }
          fgets(line,500,infile);
        }

        
      }
    }
    //if ((ifls)%((int)NoF/20) == 0) cout << (int)100*ifls/NoF << "%" << endl;
  }

  cout << "Total number of events: " << nevents << endl;

  ofstream fout;
  string file="count_";
  string dat= ".dat";
  string mezera = "_";
  
  string hadron = to_string(pid);
  string filename = file+hadron+mezera+centrality+dat;
  fout.open(filename, ofstream::app);

  fout << direct<<"\t"<<pid<<endl;
  fout << nevents<<endl;

 

  cout << nevents << "\t" << count_particle  << endl;
  fout << nevents << "\t" <<count_particle  << endl;
  
  fout << endl;
  fout.close();

  cout << "Results have been written to:  "<<filename << endl;
  
}
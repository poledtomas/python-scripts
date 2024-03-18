import os
import math as mt
import sys
from decimal import Decimal
def spectrum(directory,output_file, NoF, NoE):
    print(f"\n\nCalculating dN/dpt...")
    print(f"Processing events from directory: {directory}")

    # Cuts
    eventStep = 1
    ptMin = 0.2
    ptMax = 2.0
    yCut = 0.1

    nBins = 50
    dpt = (ptMax-ptMin)/nBins
    
    dNpi=[0]*nBins
    dNK=[0]*nBins
    dNp=[0]*nBins
    dND=[0]*nBins
    dNantiD=[0]*nBins

    dNpi_ev=[0]*nBins
    dNK_ev=[0]*nBins
    dNp_ev=[0]*nBins
    dND_ev=[0]*nBins
    dNantiD_ev=[0]*nBins 
    nevents = 0
    
    # Loop over files
    for ifls in range(1, NoF + 1):
        filename = os.path.join(directory, f"{ifls}_fin.oscar")
        
        if not os.path.exists(filename):
            print(f"Warning: Missing file #{ifls}")
            continue
        if os.path.getsize(filename)<200:
            continue

        with open(filename, "r") as infile:
            for _ in range(3):  # Skip the first three lines
                next(infile)

            

            # Loop over events in file
            for iev in range(0, NoE, eventStep):
                for _ in range(eventStep):
                    line = infile.readline()
                    pars = line.strip().split(" ")
                    npart =int(pars[4])
                    
                   
                    if npart > 0:
                        nevents += 1

                        # Loop over particles
                        for _ in range(npart):
                            line = infile.readline()
                            pars = line.strip().split(" ")
                            
                            m  = float(pars[4])
                            E  = (float(pars[5]))
                            px = float(pars[6])
                            py = float(pars[7])
                            pz = float(pars[8])
                            id = float(pars[9])
                            
                            pt = mt.sqrt(px*px+py*py)
                           
                            if (id == 211 or id == 321 or id == 2212 or id == 1000010020 or id == -1000010020):
                                if E != abs(pz):
                                    

                                    rap = 0.5*mt.log((E+pz)/(E-pz))
                                    ptBin = int((pt - ptMin) / dpt)
                                    if (abs(rap) < yCut and pt > ptMin and pt < ptMax):
                                        if (id == 211):
                                            dNpi_ev[ptBin]+= 1
                                        elif (id == 321):
                                            dNK_ev[ptBin]+= 1
                                        elif (id == 2212):
                                            dNp_ev[ptBin]+= 1
                            
                                        elif (id == 1000010020):
                                            print(id)
                                            dND_ev[ptBin]+= 1
                                        elif (id == -1000010020):
                                            dNantiD_ev[ptBin]+= 1 
                        for ibin in range(nBins):
                            dNpi[ibin] += dNpi_ev[ibin]
                            dNK[ibin] += dNK_ev[ibin]
                            dNp[ibin] += dNp_ev[ibin]
                            dND[ibin] += dND_ev[ibin] 
                            dNantiD[ibin] += dNantiD[ibin] 
                    infile.readline()

        if NoF: 
            if (ifls % (int(NoF / 20)) == 0):
                print(f"{int(100 * ifls / NoF)}%")
               

  
    with open("spectrum_"+output_file+".dat", "w") as file:
        file.write("%s \n"%directory)
        for i in range(nBins):
            pt= ptMin + (i+0.5)*dpt

            dNpi[i] /= nevents
            dNK[i] /= nevents
            dNp[i] /= nevents
            dND[i] /= nevents
            dNantiD[i] /= nevents

            dNpi[i] /= (2*mt.pi*dpt*2*yCut*(ptMin + (i+0.5)*dpt))
            dNK[i] /= (2*mt.pi*dpt*2*yCut*(ptMin + (i+0.5)*dpt))
            dNp[i] /= (2*mt.pi*dpt*2*yCut*(ptMin + (i+0.5)*dpt))
            dND[i] /= (2*mt.pi*dpt*2*yCut*(ptMin + (i+0.5)*dpt))
            dNantiD[i] /= (2*mt.pi*dpt*2*yCut*(ptMin + (i+0.5)*dpt))
             
            file.write("%s\t"%round(pt,2))
            file.write("%s\t"%dNpi[i])
            file.write("%s\t"%dNK[i])
            file.write("%s\t"%dNp[i])
            file.write("%s\t"%dND[i])
            file.write("%s\n"%dNantiD[i])
            
            print(str(dNpi[i]) + "\t"+ str(dNp[i])+"\t"+str(dNK[i])+"\t"+str(dND[i])+"\t"+str(dNantiD[i]))
spectrum("/storage/praha1/home/poledto1/hydro/hybrid/sampler.out/lhc2760-10-20-3000-deuteron/","llhc2760-10-20-3000-deuteron", NoF=100, NoE=100)
#spectrum(str(sys.argv[1]),str(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]))
#vn_C_pT("/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc5020-10-20-norm135","lhc5020-10-20-norm135", NoF=10, NoE=10, order=2)

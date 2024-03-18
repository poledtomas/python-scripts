import os
import math as mt
import sys

def dndeta(directory,output_file, NoF, NoE):
    print(f"\n\nCalculating dN/deta...")
    print(f"Processing events from directory: {directory}")

    # Cuts
    eventStep = 1
    etaMin = -5.0
    etaMax =  5.0
    nBins = 50
    deta = (etaMax-etaMin)/nBins
    dN = [0]*nBins
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
                            
                            px = float(pars[6])
                            py = float(pars[7])
                            pz = float(pars[8])
                            ele=float(pars[11])
                            
                            
                            eta = 0.5*mt.log((mt.sqrt(px*px+py*py+pz*pz)+pz)/(mt.sqrt(px*px+py*py+pz*pz)-pz))
            
                            
                            if (ele != 0 and eta < etaMax and eta > etaMin):
                                etaBin = int((eta - etaMin) / deta)
                                dN[etaBin] += 1
                                
                    infile.readline()

        if NoF: 
            if (ifls % (int(NoF / 20)) == 0):
                print(f"{int(100 * ifls / NoF)}%")
    
                
   
    with open("dndeta_"+output_file+".dat", "w") as file:
        print("Results have been written to: %s.dat"%output_file)
        file.write("%s \n"%directory)
        for i in range(nBins):
            eta= etaMin + (i+0.5)*deta

            dN[i] /= (nevents*deta)
            print(round(eta,2),dN[i])
            
            
  
#dndeta("/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/diplomka/lhc5020-20-30-norm135/","lhc5020-20-30-norm135", NoF=10, NoE=100)
dndeta(str(sys.argv[1]),str(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]))
#vn_C_pT("/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc5020-10-20-norm135","lhc5020-10-20-norm135", NoF=10, NoE=10, order=2)

import os
import math as mt
import sys

def count_pid(directory, NoF, NoE, pid):
    print(f"\n\nCalculating dN/deta...")
    print(f"Processing events from directory: {directory}")

    # Cuts
    eventStep = 1
    etaMin = -5.0
    etaMax =  5.0
    nevents = 0
    particle = 0
    antiparticle = 0
    
    # Loop over files
    count =0
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
                            id = float(pars[9])
                            
                            
                            eta = 0.5*mt.log((mt.sqrt(px*px+py*py+pz*pz)+pz)/(mt.sqrt(px*px+py*py+pz*pz)-pz))
            
                            if (ele != 0 and eta < etaMax and eta > etaMin):
                                count = count +1
                            if (ele != 0 and eta < etaMax and eta > etaMin and id ==pid):
                              
                               particle = particle+1 
                            if (ele != 0 and eta < etaMax and eta > etaMin and id == -pid):
                                antiparticle = antiparticle+1    
                    infile.readline()

        if NoF: 
            if (ifls % (int(NoF / 20)) == 0):
                print(f"{int(100 * ifls / NoF)}%")
    
    print(f"Number of particles is: {count}. Number of deuterons is: {particle}. Number of anti-deuterons is: {antiparticle}")
    
                
   
   
            
 
#dndeta("/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/diplomka/lhc5020-20-30-norm135/","lhc5020-20-30-norm135", NoF=10, NoE=100)
count_pid(str(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]))
#count_pid("/storage/praha1/home/poledto1/hydro/hybrid/sampler.out/lhc2760-10-20-3000-deuteron", NoF=100, NoE=100, pid=1000010020)

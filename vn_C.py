import os
import math as mt
import sys


def ave_2particle_correlation(Qxn, Qyn, M):
    Qn_sq = Qxn * Qxn + Qyn * Qyn

    numerator = Qn_sq - M
    
    denominator = M * (M - 1)
    return numerator / denominator

def ave_4particle_correlation(Qxn, Qyn, Qx2n, Qy2n, M):
    Qn_sq = Qxn * Qxn + Qyn * Qyn
    Q2n_sq = Qx2n * Qx2n + Qy2n * Qy2n
    Qn_isq = Qxn * Qxn - Qyn * Qyn

    first = Qn_sq * Qn_sq + Q2n_sq - 2 * (Qx2n * Qn_isq + 2 * Qy2n * Qxn * Qyn)
    second = 2 * (2 * (M - 2) * Qn_sq - M * (M - 3))

    numerator = first - second
    denominator = M * (M - 1) * (M - 2) * (M - 3)

    return numerator / denominator


def vn_C(directory,output_file, NoF, NoE, order):
    print(f"\n\nCalculating v_{int(order)}{{2}} and v_{int(order)}{{4}} (pT)...")
    print(f"Processing events from directory: {directory}")

    # Cuts
    eventStep = 1  # number of events in one super-event (to increase statistics)
    etaCut = 1.6
    ptMinCut = 0.2
    ptMaxCut = 5.0
    
    cn2_nom=0
    cn2_denom=0
    cn4_nom=0
    cn4_denom=0

    
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

            nevents = 0

            # Loop over events in file
            for iev in range(0, NoE, eventStep):
                Qx=Qy=Qx2=Qy2=0
                RFP = 0
                
                cumulant2=0
                cumulant4=0
                

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
                            E  = float(pars[5])
                            px = float(pars[6])
                            py = float(pars[7])
                            pz = float(pars[8])
                            id = float(pars[9])
                            ele=float(pars[11])
                            
                            
                            pt = float(mt.sqrt(px*px+py*py))
                            pabs = float(mt.sqrt((px*px+py*py+pz*pz)))
                            eta = float(0.5*(mt.log((pabs+pz)/(pabs-pz))))
                            phi = float(mt.atan2(py, px))
                            
                            if (ele != 0 and abs(eta) < etaCut and pt > ptMinCut and pt < ptMaxCut):
                                RFP += 1
                                Qx += mt.cos(order * phi)
                                Qy += mt.sin(order * phi)
                                Qx2 += mt.cos(2*order * phi)
                                Qy2 += mt.sin(2*order * phi)
                                
                        infile.readline()

                if RFP>0:
                    
                    cumulant2 = ave_2particle_correlation(Qx,Qy,RFP)
                    cn2_nom += RFP*(RFP-1)*cumulant2
                    cn2_denom += RFP*(RFP-1)
                    

                    cumulant4 = ave_4particle_correlation(Qx,Qy,Qx2,Qy2,RFP)
                    cn4_nom += RFP*(RFP-1)*(RFP-2)*(RFP-3)*cumulant4
                    cn4_denom += RFP*(RFP-1)*(RFP-2)*(RFP-3)

                                   
                   
                            
    cn2 = cn2_nom / cn2_denom
    cn4 = cn4_nom / cn4_denom -2 * cn2 * cn2
    
    vn2=mt.sqrt(cn2)
    vn4=pow(-cn4,0.25)
    
    print(f"v_{int(order)}{{2}}: ",vn2)
    print(f"v_{int(order)}{{4}}: ",vn4)
    with open(output_file+"ele.dat", "w") as file:
        print("Results have been written to: %s.dat"%output_file)
        file.write("%s \n"%directory)
        file.write("%s \t"%vn2)
        file.write("%s \t"%vn4)
        

#vn_C(str(sys.argv[1]),str(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]))
vn_C("/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc5020-00-05-norm135","test", NoF=20, NoE=10, order=2)

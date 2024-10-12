import os
import numpy as np
import math as mt
import sys


def dndeta_new(directory, output_file, NoF, NoE):
    print("\n\nCalculating dN/deta...")
    print(f"Processing events from directory: {directory}")

    # Cuts
    eventStep = 1
    etaMin = -5.0
    etaMax = 5.0
    nBins = 50
    deta = (etaMax - etaMin) / nBins
    dN = np.zeros(nBins)
    nevents = 0

    # Loop over files
    for ifls in range(1, NoF + 1):
        filename = os.path.join(directory, f"{ifls}_fin.oscar")

        if not os.path.exists(filename):
            print(f"Warning: Missing file #{ifls}")
            continue
        if os.path.getsize(filename) < 200:
            continue

        with open(filename, "r") as infile:
            for _ in range(3):
                next(infile)

            for iev in range(0, NoE, eventStep):
                for _ in range(eventStep):
                    line = infile.readline().strip().split(" ")
                    npart = int(line[4])

                    if npart > 0:
                        nevents += 1

                        for _ in range(npart):
                            line = infile.readline().strip().split(" ")
                            px = float(line[6])
                            py = float(line[7])
                            pz = float(line[8])
                            ele = float(line[11])
                            p = mt.sqrt(px**2 + py**2 + pz**2)
                            eta = (
                                0.5 * mt.log((p + pz) / (p - pz)) if p - pz != 0 else 0
                            )
                            if ele != 0 and etaMin < eta < etaMax:
                                etaBin = int((eta - etaMin) / deta)
                                dN[etaBin] += 1
                    infile.readline()
        if NoF and (ifls % (int(NoF / 20)) == 0):
            print(f"{int(100 * ifls / NoF)}%")

    dN /= nevents * deta

    output_path = f"dndeta_{output_file}.dat"
    with open(output_path, "w") as file:
        print(f"Results have been written to: {output_path}")
        file.write(f"{directory}\n")

        for i in range(nBins):
            eta_center = etaMin + (i + 0.5) * deta
            file.write(f"{round(eta_center, 2)}\t{dN[i]}\n")
            print(round(eta_center, 2), dN[i])


# dndeta("/path/to/directory", "output_file_name", NoF=10, NoE=100)
dndeta_new(str(sys.argv[1]), str(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))

import math as mt
import sys

import numpy as np


def dndy(directory, output_file, NoF, NoE):
    print("\n\nCalculating dndy...")
    print(f"Processing events from directory: {directory}")

    # Cuts
    yMin = -5.0
    yMax = 5.0
    nBins = 50
    dy = (yMax - yMin) / nBins
    dN = np.zeros(nBins)
    nevents = 0
    # Loop over files

    for ifls in range(1, NoF + 1):
        filename = f"{directory}{ifls}_fin.oscar"

        try:
            with open(filename, "r") as infile:
                # Skip first 3 lines
                next(infile)
                next(infile)
                next(infile)

                # Loop over events in the file
                for iev in range(NoE):
                    line = infile.readline().strip()
                    pars = line.split()

                    if len(pars) < 5:
                        continue

                    npart = int(pars[4])

                    if npart > 0:
                        nevents += 1

                        # Loop over particles
                        for i in range(npart):
                            line = infile.readline().strip()
                            pars = line.split()

                            if len(pars) < 12:
                                continue
                            E = float(pars[5])
                            px = float(pars[6])
                            py = float(pars[7])
                            pz = float(pars[8])

                            ele = int(pars[11])

                            # Calculate rapidity double rap = 0.5*log((E[i]+pz[i])/(E[i]-pz[i]));
                            p = mt.sqrt(px**2 + py**2 + pz**2)
                            y = 0.5 * mt.log((E + pz) / (E - pz)) if p - pz != 0 else 0
                            if ele != 0 and yMin < y < yMax:
                                yBin = int((y - yMin) / dy)
                                dN[yBin] += 1

                    infile.readline()

        except FileNotFoundError:
            print(f"Warning: Missing file #{ifls}")

    dN /= nevents * dy

    output_path = f"dndy_{output_file}.dat"
    with open(output_path, "w") as file:
        print(f"Results have been written to: {output_path}")
        file.write(f"{directory}\n")

        for i in range(nBins):
            y_center = yMin + (i + 0.5) * dy
            file.write(f"{round(y_center, 2)}\t{dN[i]}\n")
            print(round(y_center, 2), dN[i])


dndy(str(sys.argv[1]), str(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))

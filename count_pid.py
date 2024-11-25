import math as mt
import sys


def count_pid(directory, NoF, NoE, pid):
    print("\n\nCalculating sum of particles...")
    print(f"Processing events from directory: {directory}")

    # Cuts
    etaMin = -5.0
    etaMax = 5.0
    nevents = 0
    # Loop over files
    count = 0
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

                            px = float(pars[6])
                            py = float(pars[7])
                            pz = float(pars[8])
                            id = int(pars[9])
                            ele = int(pars[11])

                            # Calculate eta
                            p = mt.sqrt(px**2 + py**2 + pz**2)
                            eta = (
                                0.5 * mt.log((p + pz) / (p - pz)) if p - pz != 0 else 0
                            )
                            if ele != 0 and etaMin < eta < etaMax:
                                if abs(id) == pid:
                                    count += 1

                    infile.readline()

        except FileNotFoundError:
            print(f"Warning: Missing file #{ifls}")

    print(f"Number of particles == {pid} is {count}")


count_pid(
    str(sys.argv[1]),
    int(sys.argv[2]),
    int(sys.argv[3]),
    int(sys.argv[4]),
)

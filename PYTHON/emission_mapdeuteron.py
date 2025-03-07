import math
import sys


def emission_mapdeuteron(direct, output_file, NoF, NoE, pid, phi_min, phi_max):
    print("\n\n Emmission map...")
    print(f"Processing events from directory: {direct}")

    ptMinCut = 1.9
    ptMaxCut = 2.1
    nevents = 0

    eventStep = 1

    with open("emission_mapdeuteron" + output_file + str(pid) + ".dat", "a") as file:
        # Loop over files
        for ifls in range(1, NoF + 1):
            filename = f"{direct}/{ifls}_fin.oscar"

            try:
                with open(filename, "r") as infile:
                    # Skip header lines
                    for _ in range(3):
                        infile.readline()

                    # Loop over events in file
                    for iev in range(1, NoE + 1):
                        for _ in range(eventStep):
                            line = infile.readline().strip().split()

                            if len(line) <= 5:
                                npart = int(line[4])

                            else:
                                continue

                            if npart > 0:
                                nevents += 1

                                # Loop over particles
                                for _ in range(npart):
                                    line = infile.readline().strip().split()
                                    x = float(line[1])
                                    y = float(line[2])
                                    z = float(line[3])
                                    px = float(line[6])
                                    py = float(line[7])
                                    id = int(line[9])

                                    pt = math.sqrt(px**2 + py**2)
                                    phi = math.atan2(y, x)

                                    if (ptMinCut < pt < ptMaxCut) and (id == pid):
                                        phi = phi % (2 * math.pi)

                                        if phi_min < phi < phi_max:
                                            file.write(
                                                str(x)
                                                + "\t"
                                                + str(y)
                                                + "\t"
                                                + str(z)
                                                + "\n"
                                            )

                            infile.readline()

            except FileNotFoundError:
                print(f"Warning: Missing file #{ifls}")

        print(f"Total number of events: {nevents}")
    file.close()


emission_mapdeuteron(
    str(sys.argv[1]),
    str(sys.argv[2]),
    int(sys.argv[3]),
    int(sys.argv[4]),
    int(sys.argv[5]),
    float(sys.argv[6]),
    float(sys.argv[7]),
)

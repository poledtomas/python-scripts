import math
import sys
import os


def map(direct, output_file, NoF, NoE, pid):
    print(f"\n\n Emission map of  {int(pid)} (pT)...")
    print(f"Processing events from directory: {direct}")

    ptMinCut = 1.9
    ptMaxCut = 2.1
    nevents = 0
    eventStep = 1
    t = []
    x = []
    y = []
    z = []

    # Loop over files
    for ifls in range(1, NoF + 1):
        filename = f"{direct}{ifls}_fin.oscar"
        try:
            with open(filename, "r") as infile:
                # Skip header lines
                for _ in range(3):
                    infile.readline()

                # Loop over events in file
                for iev in range(1, NoE + 1):
                    for _ in range(eventStep):
                        line = infile.readline().strip().split()

                        if 0 < len(line) and len(line) <= 5:
                            npart = int(line[4])
                        else:
                            continue

                        if npart > 0:
                            nevents += 1

                            # Loop over particles
                            for _ in range(npart):
                                line = infile.readline().strip().split()
                                if len(line) > 5:
                                    time = float(line[13])
                                    xi = float(line[1])
                                    yi = float(line[2])
                                    zi = float(line[3])
                                    px = float(line[6])
                                    py = float(line[7])
                                    id = int(line[9])

                                    pt = math.sqrt(px**2 + py**2)

                                    if (ptMinCut < pt < ptMaxCut) and (id == pid):
                                        # Adjust phi to start from -pi/8
                                        t.append(time)
                                        x.append(xi)
                                        y.append(yi)
                                        z.append(zi)
                                else:
                                    continue
                        infile.readline()
            infile.close()
        except FileNotFoundError:
            print(f"Warning: Missing file #{ifls}")

    print(f"Total number of events: {nevents}")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_filename = os.path.join(script_dir, f"time_map_all_{pid}_{output_file}.dat")
    print(output_filename)

    with open(output_filename, "a") as fout:
        fout.write(f"{direct}\t{pid}\n")
        fout.write(f"{nevents}\n")
        for i in range(len(x)):
            fout.write(f"{t[i]}\t{x[i]}\t{y[i]}\t{z[i]}\n")
        fout.close()


map(
    str(sys.argv[1]),
    str(sys.argv[2]),
    int(sys.argv[3]),
    int(sys.argv[4]),
    int(sys.argv[5]),
)

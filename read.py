def neco(direct, NoF, NoE, order):
    print(f"\n\nCalculating v_{int(order)}{{2}} (pT)...")
    print(f"Processing events from directory: {direct}")

    nevents = 0
    eventStep = 1
    # Loop over files
    for ifls in range(1, NoF + 1):
        filename = f"{direct}/{ifls}_fin.oscar"

        try:
            with open(filename, "r") as infile:
                # Skip header lines
                for _ in range(3):
                    infile.readline()

                # Loop over events in file
                for iev in range(0, NoE, eventStep):
                    for _ in range(eventStep):
                        line = infile.readline().strip().split()
                        npart = int(line[4])

                        if npart > 0:
                            nevents += 1

                            # Loop over particles
                            for _ in range(npart):
                                line = infile.readline().strip().split()

                                print(line)

                        infile.readline()

        except FileNotFoundError:
            print(f"Warning: Missing file #{ifls}")

        if ifls % (NoF // 20) == 0:
            print(f"{100 * ifls // NoF}%")


if __name__ == "__main__":
    neco(
        "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-30-40-deuteron/",
        3000,
        500,
        2,
    )

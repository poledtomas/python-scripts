import math


def neco(direct, NoF, NoE, order):
    print(f"\n\nCalculating v_{int(order)}{{2}} (pT)...")
    print(f"Processing events from directory: {direct}")

    # Cuts
    etaCut = 0.8
    ptMinCut = 0.2
    ptMaxCut = 3.0
    ptBins = 14
    eventStep = 1  # Number of events in one super-event (to increase statistics)
    dpt = (ptMaxCut - ptMinCut) / ptBins

    # Arrays for calculation of v_n
    vn2 = [0.0] * ptBins
    dn2 = [0.0] * ptBins
    dn2_nom = [0.0] * ptBins
    dn2_denom = [0.0] * ptBins
    cn2 = 0.0
    cn2_nom = 0.0
    cn2_denom = 0.0

    nevents = 0

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
                    RFP = 0
                    POI = [0] * ptBins
                    cumulant2 = 0.0
                    diff_cumulant2 = [0.0] * ptBins
                    Qx, Qy, Qx2, Qy2 = 0.0, 0.0, 0.0, 0.0
                    qx = [0.0] * ptBins
                    qy = [0.0] * ptBins
                    qx2 = [0.0] * ptBins
                    qy2 = [0.0] * ptBins

                    for _ in range(eventStep):
                        line = infile.readline().strip().split()

                        if len(line) < 5:
                            continue
                        else:
                            npart = int(line[4])
                        if npart > 0:
                            nevents += 1

                            # Loop over particles
                            for _ in range(npart):
                                line = infile.readline().strip().split()
                                px = float(line[6])
                                py = float(line[7])
                                pz = float(line[8])
                                id = int(line[9])

                                pt = math.sqrt(px**2 + py**2)
                                pabs = math.sqrt(px**2 + py**2 + pz**2)
                                eta = 0.5 * math.log((pabs + pz) / (pabs - pz))
                                phi = math.atan2(py, px)

                                if (
                                    abs(eta) < etaCut
                                    and ptMinCut < pt < ptMaxCut
                                    and abs(id) == 211
                                ):
                                    RFP += 1
                                    ptBin = int((pt - ptMinCut) / dpt)
                                    POI[ptBin] += 1
                                    Qx += math.cos(order * phi)
                                    Qy += math.sin(order * phi)
                                    Qx2 += math.cos(2 * order * phi)
                                    Qy2 += math.sin(2 * order * phi)
                                    qx[ptBin] += math.cos(order * phi)
                                    qy[ptBin] += math.sin(order * phi)
                                    qx2[ptBin] += math.cos(2 * order * phi)
                                    qy2[ptBin] += math.sin(2 * order * phi)

                        infile.readline()

                    if RFP > 1:
                        Q = complex(Qx, Qy)
                        # Calculation of cumulants
                        cumulant2 = (abs(Q) ** 2 - RFP) / (RFP * (RFP - 1))

                        cn2_nom += RFP * (RFP - 1) * cumulant2
                        cn2_denom += RFP * (RFP - 1)

                        for ipt in range(ptBins):
                            if POI[ipt] > 0:
                                q = complex(qx[ipt], qy[ipt])
                                diff_cumulant2[ipt] = (
                                    q * Q.conjugate() - complex(POI[ipt], 0)
                                ).real / (POI[ipt] * RFP - POI[ipt])

                                dn2_nom[ipt] += (
                                    POI[ipt] * RFP - POI[ipt]
                                ) * diff_cumulant2[ipt]
                                dn2_denom[ipt] += POI[ipt] * RFP - POI[ipt]

        except FileNotFoundError:
            print(f"Warning: Missing file #{ifls}")

    print(f"Total number of events: {nevents}")

    # Final calculation of v_n
    cn2 = cn2_nom / cn2_denom
    print(cn2)

    for ipt in range(ptBins):
        dn2[ipt] = dn2_nom[ipt] / dn2_denom[ipt]
        vn2[ipt] = dn2[ipt] / math.sqrt(cn2)
        ptBin = (ipt + 0.5) * dpt + ptMinCut

        print(round(ptBin, 2), vn2[ipt])


if __name__ == "__main__":
    neco(
        "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-30-40-deuteron/",
        3,
        500,
        2,
    )

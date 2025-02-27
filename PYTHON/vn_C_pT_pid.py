import math as mt


def def_ave_2partucle_correlation(pnx, pny, Qx, Qy, mq, mp, M):
    diff_cumulant = float((pnx * Qx + pny * Qy - mq) / (mp * M - mq))

    return diff_cumulant


def def_ave_4particle_correlation(
    pxn, pyn, qx2n, qy2n, Qxn, Qyn, Qx2n, Qy2n, M, mp, mq
):
    qxn = pxn
    qyn = pyn

    Qn_sq = Qxn * Qxn + Qyn * Qyn

    one = (
        pxn * Qxn * Qxn * Qxn
        + pxn * Qxn * Qyn * Qyn
        + pyn * Qxn * Qxn * Qyn
        + pyn * Qyn * Qyn * Qyn
    )
    two = qx2n * Qxn * Qxn - qx2n * Qyn * Qyn + 2 * qy2n * Qxn * Qyn
    three = pxn * Qxn * Qx2n - pyn * Qyn * Qx2n + pxn * Qyn * Qy2n + pyn * Qxn * Qy2n
    four = 2 * M * (pxn * Qxn + pyn * Qyn)
    five = 2 * mq * Qn_sq
    six = 7 * (qxn * Qxn + qyn * Qyn)
    seven = Qxn * qxn + Qyn * qyn
    eight = qx2n * Qx2n + qy2n * Qy2n
    nine = 2 * (pxn * Qxn + pyn * Qyn)
    ten = 2 * mq * M
    eleven = 6 * mq

    numerator = (
        one - two - three - four - five + six - seven + eight + nine + ten - eleven
    )
    denominator = (mp * M - 3 * mq) * (M - 1) * (M - 2)

    return numerator / denominator


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


def four_particle_correlation(Qxn, Qyn, Qx2n, Qy2n, M):
    Qn = complex(Qxn, Qyn)
    Qn2 = complex(Qx2n, Qy2n)

    term1 = (
        abs(Qn * Qn * Qn.conjugate() * Qn.conjugate())
        + abs(Qn2 * Qn2.conjugate())
        - 2 * (Qn2 * Qn.conjugate() * Qn.conjugate()).real
    )
    term2 = 2 * (M - 2) * abs(Qn * Qn.conjugate()) - M * (M - 3)
    numerator = term1 - 2 * term2
    denominator = M * (M - 1) * (M - 2) * (M - 3)
    return numerator / denominator


def vn_C_pT_pid(directory, output_file, NoF, NoE, order, pid):
    print(f"\n\nCalculating v_{int(order)}{{2}} and v_{int(order)}{{4}} (pT)...")
    print(f"Processing events from directory: {directory}")

    # Cuts
    etaCut = 1.6
    ptMinCut = 0.2
    ptMaxCut = 3.0
    ptBins = 14
    dpt = (ptMaxCut - ptMinCut) / ptBins
    POI = [0] * ptBins
    dn2 = [0] * ptBins
    vn2 = [0] * ptBins

    dn2_nom = [0] * ptBins
    dn2_denom = [0] * ptBins

    cn2_nom = 0
    cn2_denom = 0

    for ifls in range(1, NoF + 1):
        filename = f"{directory}{ifls}_fin.oscar"

        try:
            with open(filename, "r") as infile:
                # Skip first 3 lines
                next(infile)
                next(infile)
                next(infile)
                nevents = 0

                # Loop over events in the file
                for iev in range(NoE):
                    line = infile.readline().strip()
                    pars = line.split()

                    Qx = Qy = 0
                    RFP = 0
                    POI = [0] * ptBins
                    qx = [0] * ptBins
                    qy = [0] * ptBins

                    cumulant2 = 0

                    cummulant = [0] * ptBins

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
                            id = float(pars[9])
                            ele = float(pars[11])

                            pt = float(mt.sqrt(px * px + py * py))
                            pabs = float(mt.sqrt((px * px + py * py + pz * pz)))
                            eta = float(0.5 * (mt.log((pabs + pz) / (pabs - pz))))
                            phi = float(mt.atan2(py, px))

                            if (
                                ele != 0
                                and abs(eta) < etaCut
                                and pt > ptMinCut
                                and pt < ptMaxCut
                                and id == pid
                            ):
                                RFP += 1
                                ptBin = int((pt - ptMinCut) / dpt)
                                POI[ptBin] += 1
                                Qx += mt.cos(order * phi)
                                Qy += mt.sin(order * phi)
                                qx[ptBin] += mt.cos(order * phi)
                                qy[ptBin] += mt.sin(order * phi)

                if RFP > 0:
                    cumulant2 = ave_2particle_correlation(Qx, Qy, RFP)
                    cn2_nom += RFP * (RFP - 1) * cumulant2
                    cn2_denom += RFP * (RFP - 1)

                    for ipt in range(ptBins):
                        if POI[ipt] > 0:
                            cummulant[ipt] = def_ave_2partucle_correlation(
                                qx[ipt], qy[ipt], Qx, Qy, POI[ipt], POI[ipt], RFP
                            )

                            dn2_nom[ipt] += (POI[ipt] * RFP - POI[ipt]) * cummulant[ipt]
                            dn2_denom[ipt] += POI[ipt] * RFP - POI[ipt]

        except FileNotFoundError:
            print(f"Warning: Missing file #{ifls}")

    cn2 = cn2_nom / cn2_denom

    print(f"v_{int(order)}{{2}}: ", mt.sqrt(cn2))

    with open(output_file + "_" + str(pid) + "_ele.dat", "w") as file:
        print("Results have been written to: %s.dat" % output_file)
        file.write("%s \n" % directory)
        for ipt in range(ptBins):
            dn2[ipt] = dn2_nom[ipt] / dn2_denom[ipt]
            vn2[ipt] = dn2[ipt] / mt.sqrt(cn2)

            ptBin = (ipt + 0.5) * dpt + ptMinCut
            file.write("%s\t" % round(ptBin, 2))
            file.write("%s\t" % vn2[ipt])
            file.write("%s\t" % cn2)

            print(round(ptBin, 2), vn2[ipt])


# Main execution block for command-line arguments
if __name__ == "__main__":
    vn_C_pT_pid(
        "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-20-30-deuteron/",
        "lhc2760-20-30-deuteron",
        100,
        100,
        2,
        2212,
    )
# vn_C_pT_pid(
#    str(sys.argv[1]),
#    str(sys.argv[2]),
#    int(sys.argv[3]),
#    int(sys.argv[4]),
#    int(sys.argv[5]),
#    int(sys.argv[6]),
# )

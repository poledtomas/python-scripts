import numpy as np
import math as mt
import os
import sys


def def_ave_2partucle_correlation(pnx, pny, Qx, Qy, mq, mp, M):
    diff_cumulant = (pnx * Qx + pny * Qy - mq) / (mp * M - mq)
    return diff_cumulant


def def_ave_4particle_correlation(
    pxn, pyn, qx2n, qy2n, Qxn, Qyn, Qx2n, Qy2n, M, mp, mq
):
    Qn_sq = Qxn**2 + Qyn**2

    one = pxn * (Qxn**3 + Qxn * Qyn**2) + pyn * (Qxn**2 * Qyn + Qyn**3)
    two = qx2n * (Qxn**2 - Qyn**2) + 2 * qy2n * Qxn * Qyn
    three = pxn * (Qxn * Qx2n - pyn * Qyn * Qx2n) + (pxn * Qyn + pyn * Qxn) * Qy2n
    four = 2 * M * (pxn * Qxn + pyn * Qyn)
    five = 2 * mq * Qn_sq
    six = 7 * (pxn * Qxn + pyn * Qyn)
    seven = Qxn * pxn + Qyn * pyn
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
    Qn_sq = Qxn**2 + Qyn**2
    numerator = Qn_sq - M
    denominator = M * (M - 1)
    return numerator / denominator


def ave_4particle_correlation(Qxn, Qyn, Qx2n, Qy2n, M):
    Qn_sq = Qxn**2 + Qyn**2
    Q2n_sq = Qx2n**2 + Qy2n**2
    Qn_isq = Qxn**2 - Qyn**2

    first = Qn_sq**2 + Q2n_sq - 2 * (Qx2n * Qn_isq + 2 * Qy2n * Qxn * Qyn)
    second = 2 * (2 * (M - 2) * Qn_sq - M * (M - 3))

    numerator = first - second
    denominator = M * (M - 1) * (M - 2) * (M - 3)
    return numerator / denominator


def vn_C_pT_pid_new(directory, output_file, NoF, NoE, order, pid):
    print(f"\n\nCalculating v_{int(order)}{{2}} and v_{int(order)}{{4}} (pT)...")
    print(f"Processing events from directory: {directory}")

    eventStep = 1
    etaCut = 1.6
    ptMinCut = 0.2
    ptMaxCut = 3.0
    ptBins = 14
    dpt = (ptMaxCut - ptMinCut) / ptBins

    POI = np.zeros(ptBins)
    dn2_nom = np.zeros(ptBins)
    dn2_denom = np.zeros(ptBins)
    dn4_nom = np.zeros(ptBins)
    dn4_denom = np.zeros(ptBins)

    cn2_nom = 0
    cn2_denom = 0
    cn4_nom = 0
    cn4_denom = 0

    # Loop over files
    for ifls in range(1, NoF + 1):
        filename = os.path.join(directory, f"{ifls}_fin.oscar")
        if not os.path.exists(filename) or os.path.getsize(filename) < 200:
            print(f"Warning: Missing or small file #{ifls}")
            continue

        with open(filename, "r") as infile:
            for _ in range(3):  # Skip first three lines
                next(infile)

            nevents = 0

            # Loop over events
            for iev in range(0, NoE, eventStep):
                Qx = Qy = Qx2 = Qy2 = 0
                RFP = 0
                POI = np.zeros(ptBins)
                qx = np.zeros(ptBins)
                qy = np.zeros(ptBins)
                qx2 = np.zeros(ptBins)
                qy2 = np.zeros(ptBins)

                for _ in range(eventStep):
                    line = infile.readline().strip().split(" ")
                    npart = int(line[4])

                    if npart > 0:
                        nevents += 1

                        # Loop over particles
                        for _ in range(npart):
                            line = infile.readline().strip().split(" ")
                            m, E, px, py, pz, id, ele = map(
                                float,
                                [
                                    line[4],
                                    line[5],
                                    line[6],
                                    line[7],
                                    line[8],
                                    line[9],
                                    line[11],
                                ],
                            )

                            pt = mt.sqrt(px**2 + py**2)
                            pabs = mt.sqrt(px**2 + py**2 + pz**2)
                            eta = 0.5 * mt.log((pabs + pz) / (pabs - pz))
                            phi = mt.atan2(py, px)

                            if (
                                ele != 0
                                and abs(eta) < etaCut
                                and ptMinCut < pt < ptMaxCut
                                and id == pid
                            ):
                                RFP += 1
                                ptBin = int((pt - ptMinCut) / dpt)
                                POI[ptBin] += 1
                                Qx += np.cos(order * phi)
                                Qy += np.sin(order * phi)
                                Qx2 += np.cos(2 * order * phi)
                                Qy2 += np.sin(2 * order * phi)
                                qx[ptBin] += np.cos(order * phi)
                                qy[ptBin] += np.sin(order * phi)
                                qx2[ptBin] += np.cos(2 * order * phi)
                                qy2[ptBin] += np.sin(2 * order * phi)

                        infile.readline()

                if RFP > 0:
                    cumulant2 = ave_2particle_correlation(Qx, Qy, RFP)
                    cn2_nom += RFP * (RFP - 1) * cumulant2
                    cn2_denom += RFP * (RFP - 1)

                    cumulant4 = ave_4particle_correlation(Qx, Qy, Qx2, Qy2, RFP)
                    cn4_nom += RFP * (RFP - 1) * (RFP - 2) * (RFP - 3) * cumulant4
                    cn4_denom += RFP * (RFP - 1) * (RFP - 2) * (RFP - 3)

                    for ipt in range(ptBins):
                        if POI[ipt] > 0:
                            commulant = def_ave_2partucle_correlation(
                                qx[ipt], qy[ipt], Qx, Qy, POI[ipt], POI[ipt], RFP
                            )
                            cumulantfour = def_ave_4particle_correlation(
                                qx[ipt],
                                qy[ipt],
                                qx2[ipt],
                                qy2[ipt],
                                Qx,
                                Qy,
                                Qx2,
                                Qy2,
                                RFP,
                                POI[ipt],
                                POI[ipt],
                            )

                            dn2_nom[ipt] += (POI[ipt] * RFP - POI[ipt]) * commulant
                            dn2_denom[ipt] += POI[ipt] * RFP - POI[ipt]

                            dn4_nom[ipt] += (
                                (POI[ipt] * RFP - 3 * POI[ipt])
                                * (RFP - 1)
                                * (RFP - 2)
                                * cumulantfour
                            )
                            dn4_denom[ipt] += (
                                (POI[ipt] * RFP - 3 * POI[ipt]) * (RFP - 1) * (RFP - 2)
                            )

        if NoF:
            print(f"{int(100 * ifls / NoF)}%")

    cn2 = cn2_nom / cn2_denom
    cn4 = cn4_nom / cn4_denom - 2 * cn2 * cn2
    vn_4 = (-cn4) ** 0.25

    print(f"v_{int(order)}{{2}}: ", mt.sqrt(cn2))
    print(f"v_{int(order)}{{4}}: ", vn_4)

    with open(f"{output_file}_{pid}_ele.dat", "w") as file:
        print(f"Results written to: {output_file}_{pid}_ele.dat")
        file.write(f"{directory}\n")
        for ipt in range(ptBins):
            dn2 = dn2_nom[ipt] / dn2_denom[ipt]
            vn2 = dn2 / mt.sqrt(cn2)

            dn4 = dn4_nom[ipt] / dn4_denom[ipt] - 2 * dn2 * cn2
            vn4 = -dn4 / (-cn4) ** 0.75

            ptBin = (ipt + 0.5) * dpt + ptMinCut
            file.write(f"{round(ptBin, 2)}\t{vn2}\t{vn4}\t{cn2}\t{cn4}\n")
            print(round(ptBin, 2), vn2, vn4)


vn_C_pT_pid_new(
    str(sys.argv[1]),
    str(sys.argv[2]),
    int(sys.argv[3]),
    int(sys.argv[4]),
    int(sys.argv[5]),
    int(sys.argv[6]),
)

import math
import random
import sys


def vn_C_pT(direct, output, NoF, NoE, order):
    print(f"\n\nCalculating v_{int(order)}{{2}} (pT)...")
    print(f"Processing events from directory: {direct}")

    # Cuts
    etaCut = 0.8
    ptMinCut = 0.2
    ptMaxCut = 3.0
    ptBins = 14
    eventStep = 1  # number of events in one super-event (to increase statistics)
    NS = 10  # number of subsamples to estimate the error
    dpt = (ptMaxCut - ptMinCut) / ptBins

    px = py = pz = m = E = 0.0
    id_val = 0
    nevents = 0

    # arrays for calculation v_n
    vn2 = [0.0] * ptBins
    vn4 = [0.0] * ptBins
    dn2 = [0.0] * ptBins
    dn2_nom = [0.0] * ptBins
    dn2_denom = [0.0] * ptBins
    dn4 = [0.0] * ptBins
    dn4_nom = [0.0] * ptBins
    dn4_denom = [0.0] * ptBins
    cn2 = cn2_nom = cn2_denom = 0.0
    cn4 = cn4_nom = cn4_denom = 0.0
    vn2_err = [0.0] * ptBins
    vn4_err = [0.0] * ptBins

    # arrays for calculation error of v_n
    cn2_sub = [0.0] * NS
    cn2_nom_sub = [0.0] * NS
    cn2_denom_sub = [0.0] * NS
    cn4_sub = [0.0] * NS
    cn4_nom_sub = [0.0] * NS
    cn4_denom_sub = [0.0] * NS

    dn2_sub = [[0.0] * ptBins for _ in range(NS)]
    dn2_nom_sub = [[0.0] * ptBins for _ in range(NS)]
    dn2_denom_sub = [[0.0] * ptBins for _ in range(NS)]
    dn4_sub = [[0.0] * ptBins for _ in range(NS)]
    dn4_nom_sub = [[0.0] * ptBins for _ in range(NS)]
    dn4_denom_sub = [[0.0] * ptBins for _ in range(NS)]
    vn2_sub = [[0.0] * ptBins for _ in range(NS)]
    vn4_sub = [[0.0] * ptBins for _ in range(NS)]

    def safe_next_line(file_obj) -> str:
        line = file_obj.readline()
        if not line:
            return ""
        return line

    # Loop over files
    for ifls in range(1, NoF + 1):
        filename = f"{direct}{ifls}_fin.oscar"
        try:
            with open(filename, "r") as infile:
                # Skip headers
                safe_next_line(infile)
                safe_next_line(infile)
                safe_next_line(infile)

                # Loop over events in file
                for _ in range(0, NoE, eventStep):
                    RFP = 0
                    POI = (
                        [0] * ptBins
                    )  # RFP = M (all particles), POI = m (particles in given ptBin)
                    cumulant2 = 0.0
                    diff_cumulant2 = [0.0] * ptBins
                    cumulant4 = 0.0
                    diff_cumulant4 = [0.0] * ptBins
                    Qx = Qy = 0.0
                    Qx2 = Qy2 = 0.0
                    qx = [0.0] * ptBins
                    qy = [0.0] * ptBins
                    qx2 = [0.0] * ptBins
                    qy2 = [0.0] * ptBins

                    for _k in range(eventStep):
                        line = safe_next_line(infile)
                        if not line:
                            break
                        parts = line.split()
                        if len(parts) < 5:
                            continue
                        npart = int(parts[4])
                        if npart > 0:
                            nevents += 1
                            # Loop over particles
                            for _i in range(npart):
                                line = safe_next_line(infile)
                                if not line:
                                    break
                                pars = line.split()
                                if len(pars) < 12:
                                    continue
                                m = float(pars[4])
                                E = float(pars[5])
                                px = float(pars[6])
                                py = float(pars[7])
                                pz = float(pars[8])
                                id_val = int(pars[9])

                                pt = math.sqrt(px * px + py * py)
                                pabs = math.sqrt(px * px + py * py + pz * pz)
                                eta = 0.5 * math.log((pabs + pz) / (pabs - pz))
                                phi = math.atan2(py, px)
                                rap = 0.5 * math.log((E + pz) / (E - pz))

                                if (
                                    abs(eta) < etaCut
                                    and pt > ptMinCut
                                    and pt < ptMaxCut
                                ):
                                    RFP += 1
                                    ptBin = int((pt - ptMinCut) / dpt)
                                    if 0 <= ptBin < ptBins:
                                        POI[ptBin] += 1
                                        Qx += math.cos(order * phi)
                                        Qy += math.sin(order * phi)
                                        Qx2 += math.cos(2 * order * phi)
                                        Qy2 += math.sin(2 * order * phi)
                                        qx[ptBin] += math.cos(order * phi)
                                        qy[ptBin] += math.sin(order * phi)
                                        qx2[ptBin] += math.cos(2 * order * phi)
                                        qy2[ptBin] += math.sin(2 * order * phi)

                        # skip end-of-event line
                        safe_next_line(infile)

                    if RFP > 0:
                        Q = complex(Qx, Qy)
                        Q2 = complex(Qx2, Qy2)
                        # calculation of cumulants
                        cumulant2 = (abs(Q) * abs(Q) - RFP) / (RFP * (RFP - 1))
                        cumulant4 = (
                            (abs(Q) ** 4)
                            + (abs(Q2) ** 2)
                            - 2 * (Q2 * Q.conjugate() * Q.conjugate()).real
                            - 4 * (RFP - 2) * abs(Q) * abs(Q)
                            + 2 * RFP * (RFP - 3)
                        ) / (RFP * (RFP - 1) * (RFP - 2) * (RFP - 3))

                        cn2_nom += RFP * (RFP - 1) * cumulant2
                        cn2_denom += RFP * (RFP - 1)
                        cn4_nom += (RFP * (RFP - 1) * (RFP - 2) * (RFP - 3)) * cumulant4
                        cn4_denom += RFP * (RFP - 1) * (RFP - 2) * (RFP - 3)

                        # distribute events to NS samples for error calculation
                        r = random.randrange(NS)
                        cn2_nom_sub[r] += RFP * (RFP - 1) * cumulant2
                        cn2_denom_sub[r] += RFP * (RFP - 1)
                        cn4_nom_sub[r] += (
                            RFP * (RFP - 1) * (RFP - 2) * (RFP - 3)
                        ) * cumulant4
                        cn4_denom_sub[r] += RFP * (RFP - 1) * (RFP - 2) * (RFP - 3)

                        for ipt in range(ptBins):
                            if POI[ipt] > 0:
                                q = complex(qx[ipt], qy[ipt])
                                q2 = complex(qx2[ipt], qy2[ipt])
                                diff_cumulant2[ipt] = (
                                    q * Q.conjugate() - complex(POI[ipt], 0)
                                ).real / (POI[ipt] * RFP - POI[ipt])

                                term1 = q * Q * Q.conjugate() * Q.conjugate()
                                term2 = q2 * Q.conjugate() * Q.conjugate()
                                term3 = q * Q * Q2.conjugate()
                                term4 = complex((9 - 2 * RFP), 0) * q * Q.conjugate()
                                term5 = 2 * POI[ipt] * abs(Q) * abs(Q)
                                term6 = Q * q.conjugate()
                                term7 = q2 * Q2.conjugate()
                                diff_cumulant4[ipt] = (
                                    term1
                                    - term2
                                    - term3
                                    + term4
                                    - term5
                                    - term6
                                    + term7
                                    + complex(2 * POI[ipt] * RFP - 6 * POI[ipt], 0)
                                ).real / (
                                    (POI[ipt] * RFP - 3 * POI[ipt])
                                    * (RFP - 1)
                                    * (RFP - 2)
                                )

                                dn2_nom[ipt] += (
                                    POI[ipt] * RFP - POI[ipt]
                                ) * diff_cumulant2[ipt]
                                dn2_denom[ipt] += POI[ipt] * RFP - POI[ipt]

                                dn4_nom[ipt] += (
                                    (POI[ipt] * RFP - 3 * POI[ipt])
                                    * (RFP - 1)
                                    * (RFP - 2)
                                    * diff_cumulant4[ipt]
                                )
                                dn4_denom[ipt] += (
                                    (POI[ipt] * RFP - 3 * POI[ipt])
                                    * (RFP - 1)
                                    * (RFP - 2)
                                )

                                dn2_nom_sub[r][ipt] += (
                                    POI[ipt] * RFP - POI[ipt]
                                ) * diff_cumulant2[ipt]
                                dn2_denom_sub[r][ipt] += POI[ipt] * RFP - POI[ipt]

                                dn4_nom_sub[r][ipt] += (
                                    (POI[ipt] * RFP - 3 * POI[ipt])
                                    * (RFP - 1)
                                    * (RFP - 2)
                                    * diff_cumulant4[ipt]
                                )
                                dn4_denom_sub[r][ipt] += (
                                    (POI[ipt] * RFP - 3 * POI[ipt])
                                    * (RFP - 1)
                                    * (RFP - 2)
                                )

            if NoF >= 20 and (ifls) % (NoF // 20) == 0:
                print(f"{int(100 * ifls / NoF)}%")
        except FileNotFoundError:
            print(f"Warning: Missing file #{ifls}")

    print(f"Total number of events: {nevents}")

    # final calculation of v_n
    cn2 = cn2_nom / cn2_denom if cn2_denom != 0 else 0.0
    cn4 = cn4_nom / cn4_denom - 2 * cn2 * cn2 if cn4_denom != 0 else 0.0

    for j in range(NS):
        cn2_sub[j] = cn2_nom_sub[j] / cn2_denom_sub[j] if cn2_denom_sub[j] != 0 else 0.0
        cn4_sub[j] = (
            cn4_nom_sub[j] / cn4_denom_sub[j] - 2 * cn2_sub[j] * cn2_sub[j]
            if cn4_denom_sub[j] != 0
            else 0.0
        )

    for ipt in range(ptBins):
        dn2[ipt] = dn2_nom[ipt] / dn2_denom[ipt] if dn2_denom[ipt] != 0 else 0.0
        vn2[ipt] = dn2[ipt] / math.sqrt(cn2) if cn2 > 0 else 0.0

        dn4[ipt] = (
            dn4_nom[ipt] / dn4_denom[ipt] - 2 * dn2[ipt] * cn2
            if dn4_denom[ipt] != 0
            else 0.0
        )
        vn4[ipt] = -dn4[ipt] / ((-cn4) ** 0.75) if cn4 < 0 else 0.0

        vn2_mean = vn2_sd = 0.0
        vn4_mean = vn4_sd = 0.0
        for j in range(NS):
            dn2_sub[j][ipt] = (
                dn2_nom_sub[j][ipt] / dn2_denom_sub[j][ipt]
                if dn2_denom_sub[j][ipt] != 0
                else 0.0
            )
            vn2_sub[j][ipt] = (
                dn2_sub[j][ipt] / math.sqrt(cn2_sub[j]) if cn2_sub[j] > 0 else 0.0
            )

            dn4_sub[j][ipt] = (
                dn4_nom_sub[j][ipt] / dn4_denom_sub[j][ipt]
                - 2 * dn2_sub[j][ipt] * cn2_sub[j]
                if dn4_denom_sub[j][ipt] != 0
                else 0.0
            )
            vn4_sub[j][ipt] = (
                -dn4_sub[j][ipt] / ((-cn4_sub[j]) ** 0.75) if cn4_sub[j] < 0 else 0.0
            )

            vn2_mean += vn2_sub[j][ipt]
            vn2_sd += vn2_sub[j][ipt] * vn2_sub[j][ipt]
            vn4_mean += vn4_sub[j][ipt]
            vn4_sd += vn4_sub[j][ipt] * vn4_sub[j][ipt]

        vn2_mean /= NS
        vn4_mean /= NS

        vn2_err[ipt] = math.sqrt(
            max(vn2_sd / NS - vn2_mean * vn2_mean, 0.0)
        ) / math.sqrt(NS)
        vn4_err[ipt] = math.sqrt(
            max(vn4_sd / NS - vn4_mean * vn4_mean, 0.0)
        ) / math.sqrt(NS)

    # Write results into the text file (append)
    with open(f"v{int(order)}_{output}_cumulants_pT.dat", "a") as fout:
        fout.write(f"{direct}\t{int(order)}\n")
        for ipt in range(ptBins):
            ptBin = (ipt + 0.5) * dpt + ptMinCut
            print(f"{ptBin}\t{vn2[ipt]}\t{vn2_err[ipt]}\t{vn4[ipt]}\t{vn4_err[ipt]}")
            fout.write(
                f"{ptBin}\t{vn2[ipt]}\t{vn2_err[ipt]}\t{vn4[ipt]}\t{vn4_err[ipt]}\n"
            )
        fout.write("\n")

    print(f"Results have been written to 'v{int(order)}_{output}_cumulants_pT.dat'")


vn_C_pT(
    str(sys.argv[1]),
    str(sys.argv[2]),
    int(sys.argv[3]),
    int(sys.argv[4]),
    float(sys.argv[5]),
)

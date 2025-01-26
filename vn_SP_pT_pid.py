import math
import sys
import os


def with_gap_calculation(QxA, QyA, QxB, QyB, MA, MB):
    numerator = QxA * QxB + QyA * QyB
    denominator = MA * MB
    return numerator / denominator


def scalar_product(ux, uy, Qx, Qy, MA, MB):
    numerator = ux * Qx + uy * Qy
    denominator = MA + MB
    return numerator / denominator


def vn_SP_pT_pid(direct, output_file, NoF, NoE, order, pid):
    print(f"\n\nCalculating v_{int(order)}{{2}} (pT)...")
    print(f"Processing events from directory: {direct}")

    ptMinCut = 0.0
    ptMaxCut = 5.0
    ptBins = 15
    dpt = (ptMaxCut - ptMinCut) / ptBins
    nevents = 0

    eventStep = 1
    vn2 = [0.0] * ptBins
    dn2_nom = [0.0] * ptBins
    dn2_denom = [0.0] * ptBins

    cn2 = 0.0

    cn2_nomAB = 0
    cn2_denomAB = 0
    eta_range_A = (2.8, 5.1)
    eta_range_C = (-0.8, 0.8)
    eta_range_B = (-3.7, -1.7)

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
                    MA = MB = 0
                    POI = [0] * ptBins
                    diff_cumulant2 = [0.0] * ptBins

                    QxA, QxB, QyA, QyB = (
                        0,
                        0,
                        0,
                        0,
                    )

                    ux = [0] * ptBins
                    uy = [0] * ptBins

                    for _ in range(eventStep):
                        line = infile.readline().strip().split()
                        # print(filename,line)

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
                                    # nevents +=1
                                    # print(line)
                                    px = float(line[6])
                                    py = float(line[7])
                                    pz = float(line[8])
                                    id = int(line[9])

                                    pt = math.sqrt(px**2 + py**2)
                                    pabs = math.sqrt(px**2 + py**2 + pz**2)
                                    eta = 0.5 * math.log((pabs + pz) / (pabs - pz))
                                    phi = math.atan2(py, px)

                                    if ptMinCut < pt < ptMaxCut:
                                        ptBin = int((pt - ptMinCut) / dpt)

                                        if 0 <= ptBin < ptBins:
                                            if (
                                                eta_range_A[0] < eta < eta_range_A[1]
                                                and abs(id) != pid
                                            ):
                                                QxA += math.cos(order * phi)
                                                QyA += math.sin(order * phi)
                                                MA += 1

                                            elif (
                                                eta_range_B[0] < eta < eta_range_B[1]
                                                and abs(id) != pid
                                            ):
                                                QxB += math.cos(order * phi)
                                                QyB += math.sin(order * phi)
                                                MB += 1

                                            elif (
                                                abs(id) == pid
                                                and eta_range_C[0]
                                                < eta
                                                < eta_range_C[1]
                                            ):
                                                ux[ptBin] += math.cos(order * phi)
                                                uy[ptBin] += math.sin(order * phi)
                                                POI[ptBin] += 1
                                else:
                                    continue

                        infile.readline()

                    if MA > 0 and MB > 0:
                        gap = with_gap_calculation(QxA, QyA, QxB, QyB, MA, MB)
                        cn2_nomAB += MA * MB * gap
                        cn2_denomAB += MA * MB
                        for ipt in range(ptBins):
                            if POI[ipt] > 0:
                                RFP = MA + MB
                                Q = complex(QxA + QxB, QyA + QyB)
                                q = complex(ux[ipt], uy[ipt])

                                diff_cumulant2[ipt] = (
                                    q * Q.conjugate() - complex(POI[ipt], 0)
                                ).real / (POI[ipt] * RFP - POI[ipt])

                                dn2_nom[ipt] += (
                                    POI[ipt] * RFP - POI[ipt]
                                ) * diff_cumulant2[ipt]
                                dn2_denom[ipt] += POI[ipt] * RFP - POI[ipt]

            infile.close()
        except FileNotFoundError:
            print(f"Warning: Missing file #{ifls}")

    print(f"Total number of events: {nevents}")

    cn2 = cn2_nomAB / cn2_denomAB
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_filename = os.path.join(script_dir, f"vn_SP_pT_pair_{pid}_{output_file}.dat")
    print(output_filename)

    # output_filename = f"vn_SP_pT_pair_{pid}_{output_file}.dat"
    with open(output_filename, "a") as fout:
        fout.write(f"{direct}\t{int(order)}\t{pid}\n")
        fout.write(f"{nevents}\n")
        for i in range(ptBins):
            ptBin = (i + 0.5) * dpt + ptMinCut
            vn2[ipt] = (dn2_nom[i] / dn2_denom[i]) / math.sqrt((cn2))
            print(round(ptBin, 2), vn2[ipt])
            fout.write(str(round(ptBin, 2)) + "\t" + str(vn2[ipt]) + "\n")
        fout.close()


vn_SP_pT_pid(
    str(sys.argv[1]),
    str(sys.argv[2]),
    int(sys.argv[3]),
    int(sys.argv[4]),
    int(sys.argv[5]),
    int(sys.argv[6]),
)

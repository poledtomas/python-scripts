import math
import sys
import os


def with_gap_calculation(QxA, QyA, QxB, QyB, MA, MB):
    numerator = QxA * QxB + QyA * QyB
    denominator = MA * MB
    return numerator / denominator


def vn_SP_pT_pid_err(direct, output_file, NoF, NoE, order, pid):
    print(f"\n\nCalculating v_{int(order)}{{2}} (pT)...")
    print(f"Processing events from directory: {direct}")

    ptMinCut = 0.0
    ptMaxCut = 5.0
    ptBins = 15
    dpt = (ptMaxCut - ptMinCut) / ptBins
    nevents = 0

    vn2 = [0.0] * ptBins
    vn2_error = [0.0] * ptBins
    dn2_nom = [0.0] * ptBins
    dn2_denom = [0.0] * ptBins
    dn2_nom_sq = [0.0] * ptBins  # Squared sum for error calculation

    cn2_nomAB = 0
    cn2_denomAB = 0

    eta_range_A = (2.8, 5.1)
    eta_range_C = (-0.8, 0.8)
    eta_range_B = (-3.7, -1.7)

    for ifls in range(1, NoF + 1):
        filename = f"{direct}{ifls}_fin.oscar"
        try:
            with open(filename, "r") as infile:
                for _ in range(3):
                    infile.readline()

                for iev in range(1, NoE + 1):
                    MA = MB = 0
                    POI = [0] * ptBins
                    diff_cumulant2 = [0.0] * ptBins

                    QxA = QyA = QxB = QyB = 0
                    ux = [0] * ptBins
                    uy = [0] * ptBins

                    line = infile.readline().strip().split()
                    if len(line) > 0 and len(line) >= 5:
                        npart = int(line[4])
                    else:
                        continue

                    if npart > 0:
                        nevents += 1
                        for _ in range(npart):
                            line = infile.readline().strip().split()
                            if len(line) > 5:
                                px = float(line[6].strip().replace("\x00", ""))
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
                                            and eta_range_C[0] < eta < eta_range_C[1]
                                        ):
                                            ux[ptBin] += math.cos(order * phi)
                                            uy[ptBin] += math.sin(order * phi)
                                            POI[ptBin] += 1
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
                                dn2_nom_sq[ipt] += (
                                    (POI[ipt] * RFP - POI[ipt]) * diff_cumulant2[ipt]
                                ) ** 2
        except FileNotFoundError:
            print(f"Warning: Missing file #{ifls}")

    cn2 = cn2_nomAB / cn2_denomAB
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_filename = os.path.join(script_dir, f"vn_SP_pT_pair_{pid}_{output_file}.dat")

    with open(output_filename, "a") as fout:
        fout.write(f"{direct}\t{int(order)}\t{pid}\n")
        fout.write(f"{nevents}\n")
        for i in range(ptBins):
            ptBin = (i + 0.5) * dpt + ptMinCut
            vn2[i] = (dn2_nom[i] / dn2_denom[i]) / math.sqrt(cn2)
            vn2_error[i] = (
                math.sqrt(
                    (1.0 / (4 * dn2_denom[i]))
                    * (dn2_nom_sq[i] / dn2_denom[i] - (dn2_nom[i] / dn2_denom[i]) ** 2)
                )
                / math.sqrt(nevents)
                if nevents > 0
                else 0
            )
            fout.write(f"{round(ptBin, 2)}\t{vn2[i]}\t{vn2_error[i]}\n")
            print(round(ptBin, 2), vn2[i], vn2_error[i])


def main():
    vn_SP_pT_pid_err(
        str(sys.argv[1]),
        str(sys.argv[2]),
        int(sys.argv[3]),
        int(sys.argv[4]),
        int(sys.argv[5]),
        int(sys.argv[6]),
    )


if __name__ == "__main__":
    main()

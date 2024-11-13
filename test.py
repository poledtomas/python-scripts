import math
import os


def def_ave_2partucle_correlation(pnx, pny, Qx, Qy, mq, mp, M):
    diff_cumulant = float((pnx * Qx + pny * Qy - mq) / (mp * M - mq))

    return diff_cumulant


def ave_2particle_correlation(uxn, uyn, QxB, QyB, M):
    numerator = uxn * QxB + uyn * QyB - 1
    denominator = M - 1

    return numerator / denominator


def with_gap_calculation(QxA, QyA, QxB, QyB, MA, MB):
    numerator = QxA * QxB + QyA * QyB
    denominator = MA * MB
    return numerator / denominator


def vn_SP_pT_pid(directory, output_file, NoF, NoE, order, pid):
    """
    Calculates vn using the Scalar Product (SP) Method for transverse momentum bins.
    """
    print(
        f"Processing events from directory: {directory} using Scalar Product (SP) method"
    )

    etaCut = 0.8
    etaAmax = 5.1
    etaAmin = 2.8
    etaBmax = -1.7
    etaBmin = -3.7
    ptMinCut = 0.2
    ptMaxCut = 3.0
    ptBins = 14
    dpt = (ptMaxCut - ptMinCut) / ptBins
    nevents = 0

    dn2_nom = [0] * ptBins
    dn2_denom = [0] * ptBins
    denom = nom = 0
    # Loop over files
    for ifls in range(1, NoF + 1):
        filename = os.path.join(directory, f"{ifls}_fin.oscar")
        if not os.path.exists(filename):
            print(f"Warning: Missing file #{ifls}")
            continue

        with open(filename, "r") as infile:
            # Skip headers
            for _ in range(3):
                infile.readline()

            # Process each event
            for iev in range(NoE):
                line = infile.readline().strip()
                pars = line.split()
                Qx = Qy = QxA = QyA = QxB = QyB = 0
                M = MA = MB = 0
                POI = ux = uy = [0] * ptBins

                commulant = [0] * ptBins

                if len(pars) < 5:
                    continue

                npart = int(pars[4])

                if npart > 0:
                    nevents += 1

                    for i in range(npart):
                        line = infile.readline().strip()
                        pars = line.split()

                        if len(pars) < 12:
                            continue

                        px = float(pars[6])
                        py = float(pars[7])
                        pz = float(pars[8])
                        id = int(pars[9])

                        pabs = math.sqrt(px * px + py * py + pz * pz)
                        pt = math.sqrt(px * px + py * py)
                        eta = 0.5 * math.log((pabs + pz) / (pabs - pz))
                        phi = math.atan2(py, px)

                        if abs(eta) < etaCut and pt > ptMinCut and pt < ptMaxCut:
                            if abs(id) == pid:
                                ptBin = int((pt - ptMinCut) / dpt)
                                POI[ptBin] += 1
                                ux[ptBin] += math.cos(order * phi)
                                uy[ptBin] += math.sin(order * phi)
                            else:
                                Qx += math.cos(order * phi)
                                Qy += math.sin(order * phi)
                                M += 1
                        if (
                            etaAmin < eta
                            and eta < etaAmax
                            and pt > ptMinCut
                            and pt < ptMaxCut
                        ):
                            QxA += math.cos(order * phi)
                            QyA += math.sin(order * phi)
                            MA += 1

                        elif (
                            etaBmin < eta
                            and eta < etaBmax
                            and pt > ptMinCut
                            and pt < ptMaxCut
                        ):
                            QxB += math.cos(order * phi)
                            QyB += math.sin(order * phi)
                            MB += 1

                if M > 0:
                    de = with_gap_calculation(QxA, QyA, QxB, QyB, MA, MB)
                    denom += (MA + MB) * de
                    nom += MA + MB
                    for ipt in range(ptBins):
                        if POI[ipt] > 0:
                            commulant[ipt] = def_ave_2partucle_correlation(
                                ux[ipt], uy[ipt], Qx, Qy, POI[ipt], POI[ipt], M
                            )
                            dn2_nom[ipt] += (POI[ipt] * M - POI[ipt]) * commulant[ipt]
                            dn2_denom[ipt] += POI[ipt] * M - POI[ipt]

    vn = [0] * ptBins
    for ipt in range(ptBins):
        vn[ipt] = (dn2_nom[ipt] / dn2_denom[ipt]) / math.sqrt(denom / nom)
        print(vn[ipt])


# Main execution block for command-line arguments
if __name__ == "__main__":
    vn_SP_pT_pid(
        "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-20-30-deuteron/",
        "lhc2760-20-30-deuteron",
        10,
        100,
        2,
        2212,
    )
    # if len(sys.argv) != 7:
    #    print(
    #        "Usage: python script.py <directory> <output_file> <NoF> <NoE> <order> <pid>"
    #    )
    # else:
    #    vn_SP_pT_pid(
    #        sys.argv[1],
    #        sys.argv[2],
    #        int(sys.argv[3]),
    #        int(sys.argv[4]),
    #        int(sys.argv[5]),
    #       int(sys.argv[6]),
    #    )

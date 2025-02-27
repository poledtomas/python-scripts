import os
import numpy as np
import math
from scipy.special import iv
import sys


def vn_EP_pT_pid(directory, output_file, NoF, NoE, order, pid):
    print(f"Processing events from directory: {directory}")

    # Cuts and parameters
    etaCut = 1.0
    ptMinCut = 0.2
    ptMaxCut = 3.0
    nBins = 14
    dpt = (ptMaxCut - ptMinCut) / nBins

    vn_obs = np.zeros(nBins)
    sd1 = np.zeros(nBins)
    vnerr = np.zeros(nBins)
    nevents = np.zeros(nBins, dtype=int)

    Rn = 0.0

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

            # Loop over events
            for iev in range(NoE):
                line = infile.readline().strip()
                npart = int(float(line.split()[4]))

                if npart > 0:
                    Qx, Qy = 0.0, 0.0
                    E_ev, px_ev, py_ev, pz_ev, id_ev, ele_ev = [], [], [], [], [], []

                    # First loop over particles to calculate the flow vector
                    for _ in range(npart):
                        line = infile.readline().strip().split()
                        E, px, py, pz = (
                            float(line[5]),
                            float(line[6]),
                            float(line[7]),
                            float(line[8]),
                        )
                        id, ele = int(line[9]), int(line[11])

                        pabs = math.sqrt(px**2 + py**2 + pz**2)
                        pt = math.sqrt(px**2 + py**2)
                        eta = 0.5 * math.log((pabs + pz) / (pabs - pz))
                        phi = math.atan2(py, px)

                        if (
                            abs(eta) < etaCut
                            and ptMinCut < pt < ptMaxCut
                            and abs(id) == pid
                        ):
                            Qx += pt * math.cos(order * phi)
                            Qy += pt * math.sin(order * phi)

                        id_ev.append(id)
                        E_ev.append(E)
                        px_ev.append(px)
                        py_ev.append(py)
                        pz_ev.append(pz)
                        ele_ev.append(ele)

                    # Second loop over particles to calculate vn
                    _vn_obs = np.zeros(nBins)
                    QxA, QxB, QyA, QyB = 0.0, 0.0, 0.0, 0.0
                    _nvn = np.zeros(nBins, dtype=int)

                    for i in range(npart):
                        pt = math.sqrt(px_ev[i] ** 2 + py_ev[i] ** 2)
                        pabs = math.sqrt(px_ev[i] ** 2 + py_ev[i] ** 2 + pz_ev[i] ** 2)
                        eta = 0.5 * math.log((pabs + pz_ev[i]) / (pabs - pz_ev[i]))
                        phi = math.atan2(py_ev[i], px_ev[i])

                        if (
                            ptMinCut < pt < ptMaxCut
                            and abs(eta) < etaCut
                            and ele_ev[i] != 0
                            and abs(id_ev[i]) == pid
                        ):
                            _Qx = Qx - pt * math.cos(order * phi)
                            _Qy = Qy - pt * math.sin(order * phi)
                            psin = math.atan2(_Qy, _Qx) / order
                            ptBin = int((pt - ptMinCut) // dpt)
                            _vn_obs[ptBin] += math.cos(order * phi) * math.cos(
                                order * psin
                            ) + math.sin(order * phi) * math.sin(order * psin)
                            _nvn[ptBin] += 1

                        if abs(eta) < etaCut and ptMinCut < pt < ptMaxCut:
                            if i % 2 == 0:
                                QxA += pt * math.cos(order * phi)
                                QyA += pt * math.sin(order * phi)
                            else:
                                QxB += pt * math.cos(order * phi)
                                QyB += pt * math.sin(order * phi)

                    psi2A = math.atan2(QyA, QxA) / order
                    psi2B = math.atan2(QyB, QxB) / order
                    Rn += math.cos(order * (psi2A - psi2B))

                    for i in range(nBins):
                        if _nvn[i] > 0:
                            nevents[i] += 1
                            _vn_obs[i] /= _nvn[i]
                        vn_obs[i] += _vn_obs[i]
                        sd1[i] += _vn_obs[i] ** 2

                # Skip the last line in event
                infile.readline()

    print(f"Number of events: {nevents[0]}")

    for i in range(nBins):
        vn_obs[i] /= nevents[i]
        vnerr[i] = math.sqrt(sd1[i] / nevents[i] - vn_obs[i] ** 2)

    Rn = math.sqrt(Rn / nevents[0])
    print(f"Rn^sub = {Rn}")

    # Calculate ksi using iterative method
    pf = math.sqrt(math.pi) / (2.0 * math.sqrt(2.0))
    ksiMin, ksiMax = 0.0, 20.0

    while ksiMax - ksiMin > 0.01:
        ksi = 0.5 * (ksiMin + ksiMax)
        R = (
            pf
            * ksi
            * math.exp(-0.25 * ksi**2)
            * (iv(0, 0.25 * ksi**2) + iv(1, 0.25 * ksi**2))
        )
        if R > Rn:
            ksiMax = ksi
        else:
            ksiMin = ksi

    ksi = math.sqrt(2) * 0.5 * (ksiMin + ksiMax)
    print(f"ksi = {ksi}")

    Rn = (
        pf
        * ksi
        * math.exp(-0.25 * ksi**2)
        * (iv(0, 0.25 * ksi**2) + iv(1, 0.25 * ksi**2))
    )

    # Write results into the text file
    output_filename = f"vn_EP_pT_pair_{pid}_{output_file}.dat"
    with open(output_filename, "a") as fout:
        fout.write(f"{directory}\t{int(order)}\t{pid}\n")

        for ipt in range(nBins):
            vn = vn_obs[ipt] / Rn
            vnerr[ipt] = vnerr[ipt] / Rn / math.sqrt(nevents[ipt])
            ptBin = ptMinCut + (ipt + 0.5) * dpt
            print(f"{ptBin:.3f}\t{vn:.3f}\t{vnerr[ipt]:.3f}")
            fout.write(f"{ptBin:.3f}\t{vn:.3f}\t{vnerr[ipt]:.3f}\n")

    print("Results have been written to the output file.")


vn_EP_pT_pid(
    str(sys.argv[1]),
    str(sys.argv[2]),
    int(sys.argv[3]),
    int(sys.argv[4]),
    int(sys.argv[5]),
    int(sys.argv[6]),
)

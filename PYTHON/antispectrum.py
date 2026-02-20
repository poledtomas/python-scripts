import math
import os
import re
import sys
from typing import List, Optional

_TOKEN_SPLIT_RE = re.compile(r"[\s,]+")


def _tokenize(line: str) -> List[str]:
    line = line.strip()
    if not line:
        return []
    return [tok for tok in _TOKEN_SPLIT_RE.split(line) if tok]


def _safe_float(token: str) -> Optional[float]:
    try:
        return float(token)
    except (TypeError, ValueError):
        return None


def _safe_int(token: str) -> Optional[int]:
    try:
        return int(token)
    except (TypeError, ValueError):
        return None


def antispectrum(directory, output_file, NoF, NoE):

    print("\n\nCalculating dN/dp_t...")
    print(f"Processing events from directory: {directory}")

    ptMin = 0.2
    ptMax = 2.0
    yCut = 0.1
    nBins = 36
    pi = math.pi
    dpt = (ptMax - ptMin) / nBins

    dNpi = [0.0] * nBins
    dNK = [0.0] * nBins
    dNp = [0.0] * nBins

    dNpi_err_acc = [0.0] * nBins
    dNK_err_acc = [0.0] * nBins
    dNp_err_acc = [0.0] * nBins

    nevents = 0
    progress_step = max(1, NoF // 20) if NoF > 0 else 1

    for ifls in range(1, NoF + 1):
        filename = f"{directory}{ifls}_fin.oscar"
        print(f"Processing file: {filename}")

        if not os.path.exists(filename):
            print(f"Warning: Missing file #{ifls}")
            if ifls % progress_step == 0:
                print(f"{int(100 * ifls / NoF)}%")
            continue

        try:
            with open(filename, "r", encoding="utf-8", errors="replace") as f:
                # přesně jako 3x fgets v C
                for _ in range(3):
                    if f.readline() == "":
                        break

                for _iev in range(NoE):
                    header = f.readline()
                    if header == "":
                        break

                    parts = _tokenize(header)

                    # očekáváme hlavičku ve tvaru
                    # "# event <ievent> ensemble <iensemble> out <npart>"
                    if len(parts) < 2 or parts[0] != "#" or parts[1] != "event":
                        continue

                    # počet částic je poslední číslo na řádku (např. 1106)
                    npart_val = _safe_int(parts[-1])
                    if npart_val is None or npart_val <= 0:
                        continue

                    npart = npart_val
                    nevents += 1

                    dNpi_ev = [0.0] * nBins
                    dNK_ev = [0.0] * nBins
                    dNp_ev = [0.0] * nBins

                    for i in range(npart):
                        line = f.readline()
                        if line == "":
                            break

                        pars = _tokenize(line)
                        if len(pars) < 10:
                            continue

                        # indexy odpovídají pars[4..9] v C
                        m = _safe_float(pars[4])
                        E = _safe_float(pars[5])
                        px = _safe_float(pars[6])
                        py = _safe_float(pars[7])
                        pz = _safe_float(pars[8])
                        pid = _safe_int(pars[9])

                        if None in (m, E, px, py, pz, pid):
                            continue

                        pt = math.sqrt(px * px + py * py)

                        denom = E - pz
                        num = E + pz
                        if denom <= 0.0 or num <= 0.0:
                            continue

                        rap = 0.5 * math.log(num / denom)

                        if abs(rap) < yCut and pt > ptMin and pt < ptMax:
                            ptBin = int((pt - ptMin) / dpt)
                            if 0 <= ptBin < nBins:
                                if pid == -211:
                                    dNpi_ev[ptBin] += 1.0
                                if pid == -321:
                                    dNK_ev[ptBin] += 1.0
                                if pid == -2212:
                                    dNp_ev[ptBin] += 1.0

                    _ = f.readline()

                    for ibin in range(nBins):
                        dNpi[ibin] += dNpi_ev[ibin]
                        dNK[ibin] += dNK_ev[ibin]
                        dNp[ibin] += dNp_ev[ibin]
                        dNpi_err_acc[ibin] += dNpi_ev[ibin] * dNpi_ev[ibin]
                        dNK_err_acc[ibin] += dNK_ev[ibin] * dNK_ev[ibin]
                        dNp_err_acc[ibin] += dNp_ev[ibin] * dNp_ev[ibin]

        except OSError as e:
            print(f"Warning: Could not read file #{ifls}: {e}")

        if ifls % progress_step == 0:
            print(f"{int(100 * ifls / NoF)}%")

    print(f"Total number of events: {nevents}")

    if nevents <= 0:
        print("No events found; nothing written.")
        return

    with open(output_file + ".dat", "a", encoding="utf-8") as fout:
        fout.write(f"{directory}\n")
        fout.write("pt\tpi\terr_pi\tK\terr_K\tp\terr_p\n")

        for i in range(nBins):
            mean_pi = dNpi[i] / nevents
            mean_k = dNK[i] / nevents
            mean_p = dNp[i] / nevents

            # sqrt( <x^2> - <x>^2 ) / sqrt(N)
            var_pi = dNpi_err_acc[i] / nevents - mean_pi * mean_pi
            var_k = dNK_err_acc[i] / nevents - mean_k * mean_k
            var_p = dNp_err_acc[i] / nevents - mean_p * mean_p

            var_pi = max(0.0, var_pi)
            var_k = max(0.0, var_k)
            var_p = max(0.0, var_p)

            err_pi = math.sqrt(var_pi) / math.sqrt(nevents)
            err_k = math.sqrt(var_k) / math.sqrt(nevents)
            err_p = math.sqrt(var_p) / math.sqrt(nevents)

            pt_center = ptMin + (i + 0.5) * dpt
            norm = 2.0 * pi * dpt * 2.0 * yCut * pt_center

            mean_pi /= norm
            mean_k /= norm
            mean_p /= norm
            err_pi /= norm
            err_k /= norm
            err_p /= norm

            fout.write(
                f"{pt_center}\t{mean_pi}\t{err_pi}\t{mean_k}\t{err_k}\t{mean_p}\t{err_p}\n"
            )

        fout.write("\n")

    print(f"Results have been written to '{output_file}'")


antispectrum(str(sys.argv[1]), str(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))

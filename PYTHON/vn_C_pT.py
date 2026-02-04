#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
import random
import re
import sys
from pathlib import Path
from typing import List, Optional, Sequence, Set

_DELIM_RE = re.compile(r"[\s,]+")


def _tokens(line: str) -> List[str]:
    line = line.strip()
    if not line:
        return []
    return [t for t in _DELIM_RE.split(line) if t]


def _safe_readline(f) -> Optional[str]:
    line = f.readline()
    if line == "":
        return None
    return line


def _is_probably_directory_prefix(direct: str) -> bool:
    try:
        return Path(direct).is_dir()
    except OSError:
        return False


def _is_probably_single_file(direct: str) -> bool:
    try:
        p = Path(direct)
        return p.is_file() and p.name.endswith("_fin.oscar")
    except OSError:
        return False


def _input_filename(direct: str, ifls: int) -> str:
    if _is_probably_single_file(direct):
        return direct
    if _is_probably_directory_prefix(direct):
        return str(Path(direct) / f"{ifls}_fin.oscar")
    return f"{direct}{ifls}_fin.oscar"


def _autodetect_nof_from_dir(direct: str) -> Optional[int]:
    path = Path(direct)
    if not path.is_dir():
        return None
    indices: Set[int] = set()
    for candidate in path.glob("*_fin.oscar"):
        m = re.match(r"^(\d+)_fin\.oscar$", candidate.name)
        if m:
            try:
                indices.add(int(m.group(1)))
            except ValueError:
                pass
    if not indices:
        return None
    return max(indices)


def vn_C_pT(
    direct: str,
    NoF: int,
    NoE: Optional[int],
    order: float,
    *,
    out_path: str = "vn_cumulants_pT.dat",
    seed: Optional[int] = 0,
) -> None:
    print(f"\n\nCalculating v_{int(order)}{{2}} (pT)...")
    print(f"Processing events from directory: {direct}")

    # Cuts
    etaCut = 0.8
    ptMinCut = 0.2
    ptMaxCut = 3.0
    ptBins = 14
    eventStep = 1  # number of events in one super-event
    NS = 10  # number of subsamples to estimate the error
    dpt = (ptMaxCut - ptMinCut) / ptBins

    # arrays for calculation v_n
    vn2 = [0.0] * ptBins
    vn4 = [0.0] * ptBins

    dn2 = [0.0] * ptBins
    dn2_nom = [0.0] * ptBins
    dn2_denom = [0.0] * ptBins

    dn4 = [0.0] * ptBins
    dn4_nom = [0.0] * ptBins
    dn4_denom = [0.0] * ptBins

    cn2_nom = 0.0
    cn2_denom = 0.0
    cn4_nom = 0.0
    cn4_denom = 0.0

    vn2_err = [0.0] * ptBins
    vn4_err = [0.0] * ptBins

    # arrays for calculation error of v_n (subsamples)
    cn2_nom_sub = [0.0] * NS
    cn2_denom_sub = [0.0] * NS
    cn4_nom_sub = [0.0] * NS
    cn4_denom_sub = [0.0] * NS

    dn2_nom_sub = [[0.0] * ptBins for _ in range(NS)]
    dn2_denom_sub = [[0.0] * ptBins for _ in range(NS)]
    dn4_nom_sub = [[0.0] * ptBins for _ in range(NS)]
    dn4_denom_sub = [[0.0] * ptBins for _ in range(NS)]

    rng = random.Random(seed) if seed is not None else random.Random()

    nevents = 0
    progress_step = max(1, NoF // 20) if NoF > 0 else 1

    for ifls in range(1, NoF + 1):
        filename = _input_filename(direct, ifls)

        try:
            with open(filename, "r", encoding="utf-8", errors="replace") as infile:
                # skip headers
                for _ in range(3):
                    if _safe_readline(infile) is None:
                        break

                if NoE is None:
                    iev_iter = iter(int, 1)  # infinite
                else:
                    iev_iter = range(0, NoE, eventStep)

                for _iev in iev_iter:
                    RFP = 0
                    POI = [0] * ptBins
                    diff_cumulant2 = [0.0] * ptBins
                    diff_cumulant4 = [0.0] * ptBins

                    Qx = 0.0
                    Qy = 0.0
                    Qx2 = 0.0
                    Qy2 = 0.0
                    qx = [0.0] * ptBins
                    qy = [0.0] * ptBins
                    qx2 = [0.0] * ptBins
                    qy2 = [0.0] * ptBins

                    eof = False
                    for _k in range(eventStep):
                        header = _safe_readline(infile)
                        if header is None:
                            eof = True
                            break

                        parts = _tokens(header)
                        if len(parts) < 5:
                            # malformed header; try continue to next
                            continue

                        try:
                            npart = int(parts[4])
                        except ValueError:
                            continue

                        if npart > 0:
                            nevents += 1

                            for _i in range(npart):
                                pline = _safe_readline(infile)
                                if pline is None:
                                    eof = True
                                    break

                                p = _tokens(pline)
                                if len(p) < 9:
                                    continue

                                try:
                                    E = float(p[5])
                                    px = float(p[6])
                                    py = float(p[7])
                                    pz = float(p[8])
                                except ValueError:
                                    continue

                                pt = math.sqrt(px * px + py * py)
                                if not (ptMinCut < pt < ptMaxCut):
                                    continue

                                pabs = math.sqrt(px * px + py * py + pz * pz)
                                denom_eta = pabs - pz
                                numer_eta = pabs + pz
                                if denom_eta <= 0.0 or numer_eta <= 0.0:
                                    continue

                                eta = 0.5 * math.log(numer_eta / denom_eta)
                                if abs(eta) >= etaCut:
                                    continue

                                phi = math.atan2(py, px)

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

                        # C++ reads one extra line after each event
                        _safe_readline(infile)

                        if eof:
                            break

                    if eof and NoE is None:
                        # End of file in EOF-mode
                        break

                    if RFP <= 1:
                        continue

                    Q = complex(Qx, Qy)
                    Q2 = complex(Qx2, Qy2)

                    absQ2 = abs(Q) ** 2
                    cumulant2 = (absQ2 - RFP) / (RFP * (RFP - 1))

                    cn2_nom += RFP * (RFP - 1) * cumulant2
                    cn2_denom += RFP * (RFP - 1)

                    # only well-defined for RFP >= 4
                    cumulant4 = None
                    if RFP >= 4:
                        term = (
                            (abs(Q) ** 4)
                            + (abs(Q2) ** 2)
                            - 2.0 * (Q2 * Q.conjugate() * Q.conjugate()).real
                            - 4.0 * (RFP - 2) * absQ2
                            + 2.0 * RFP * (RFP - 3)
                        )
                        denom = RFP * (RFP - 1) * (RFP - 2) * (RFP - 3)
                        cumulant4 = term / denom
                        cn4_nom += denom * cumulant4
                        cn4_denom += denom

                    # distribute events to NS samples for error calculation
                    r = rng.randrange(NS)
                    cn2_nom_sub[r] += RFP * (RFP - 1) * cumulant2
                    cn2_denom_sub[r] += RFP * (RFP - 1)
                    if cumulant4 is not None:
                        denom4 = RFP * (RFP - 1) * (RFP - 2) * (RFP - 3)
                        cn4_nom_sub[r] += denom4 * cumulant4
                        cn4_denom_sub[r] += denom4

                    # differential cumulants
                    for ipt in range(ptBins):
                        m = POI[ipt]
                        if m <= 0:
                            continue

                        if RFP <= 1:
                            continue

                        q = complex(qx[ipt], qy[ipt])
                        q2 = complex(qx2[ipt], qy2[ipt])

                        denom2 = m * RFP - m
                        if denom2 <= 0:
                            continue

                        diff_cumulant2[ipt] = (
                            (q * Q.conjugate()) - complex(m, 0.0)
                        ).real / denom2

                        dn2_nom[ipt] += denom2 * diff_cumulant2[ipt]
                        dn2_denom[ipt] += denom2
                        dn2_nom_sub[r][ipt] += denom2 * diff_cumulant2[ipt]
                        dn2_denom_sub[r][ipt] += denom2

                        if RFP < 4:
                            continue

                        denom4 = (m * RFP - 3 * m) * (RFP - 1) * (RFP - 2)
                        if denom4 == 0:
                            continue

                        term1 = q * Q * Q.conjugate() * Q.conjugate()
                        term2 = q2 * Q.conjugate() * Q.conjugate()
                        term3 = q * Q * Q2.conjugate()
                        term4 = complex((9 - 2 * RFP), 0.0) * q * Q.conjugate()
                        term5 = complex(2 * m * absQ2, 0.0)
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
                            + complex((2 * m * RFP - 6 * m), 0.0)
                        ).real / denom4

                        dn4_nom[ipt] += denom4 * diff_cumulant4[ipt]
                        dn4_denom[ipt] += denom4
                        dn4_nom_sub[r][ipt] += denom4 * diff_cumulant4[ipt]
                        dn4_denom_sub[r][ipt] += denom4

        except FileNotFoundError:
            print(f"Warning: Missing file #{ifls}")
        except OSError as e:
            print(f"Warning: Could not read file #{ifls}: {e}")

        if ifls % progress_step == 0:
            print(f"{int(100 * ifls / NoF)}%")

    print(f"Total number of events: {nevents}")

    if cn2_denom == 0.0:
        raise RuntimeError("No valid events for cn2 (cn2_denom=0).")

    cn2 = cn2_nom / cn2_denom

    if cn4_denom == 0.0:
        cn4 = float("nan")
    else:
        cn4 = cn4_nom / cn4_denom - 2.0 * cn2 * cn2

    cn2_sub = [float("nan")] * NS
    cn4_sub = [float("nan")] * NS
    for j in range(NS):
        if cn2_denom_sub[j] > 0.0:
            cn2_sub[j] = cn2_nom_sub[j] / cn2_denom_sub[j]
        if cn4_denom_sub[j] > 0.0 and not math.isnan(cn2_sub[j]):
            cn4_sub[j] = (
                cn4_nom_sub[j] / cn4_denom_sub[j] - 2.0 * cn2_sub[j] * cn2_sub[j]
            )

    # final calculation of v_n
    for ipt in range(ptBins):
        if dn2_denom[ipt] > 0.0:
            dn2[ipt] = dn2_nom[ipt] / dn2_denom[ipt]
        else:
            dn2[ipt] = float("nan")

        if dn4_denom[ipt] > 0.0:
            dn4_raw = dn4_nom[ipt] / dn4_denom[ipt]
            dn4[ipt] = dn4_raw - 2.0 * dn2[ipt] * cn2
        else:
            dn4[ipt] = float("nan")

        if cn2 > 0.0 and not math.isnan(dn2[ipt]):
            vn2[ipt] = dn2[ipt] / math.sqrt(cn2)
        else:
            vn2[ipt] = float("nan")

        if (not math.isnan(cn4)) and cn4 < 0.0 and (not math.isnan(dn4[ipt])):
            vn4[ipt] = -dn4[ipt] / ((-cn4) ** 0.75)
        else:
            vn4[ipt] = float("nan")

        # subsample error estimate
        vn2_vals: List[float] = []
        vn4_vals: List[float] = []
        for j in range(NS):
            if cn2_denom_sub[j] <= 0.0 or not (cn2_sub[j] > 0.0):
                continue

            if dn2_denom_sub[j][ipt] > 0.0:
                dn2_j = dn2_nom_sub[j][ipt] / dn2_denom_sub[j][ipt]
                vn2_j = dn2_j / math.sqrt(cn2_sub[j])
                vn2_vals.append(vn2_j)

            if (
                cn4_denom_sub[j] > 0.0
                and (not math.isnan(cn4_sub[j]))
                and cn4_sub[j] < 0.0
            ):
                if dn4_denom_sub[j][ipt] > 0.0:
                    dn4_j_raw = dn4_nom_sub[j][ipt] / dn4_denom_sub[j][ipt]
                    dn4_j = (
                        dn4_j_raw - 2.0 * dn2_j * cn2_sub[j]
                    )  # uses dn2_j from above when available
                    vn4_j = -dn4_j / ((-cn4_sub[j]) ** 0.75)
                    vn4_vals.append(vn4_j)

        if len(vn2_vals) >= 2:
            mean = sum(vn2_vals) / len(vn2_vals)
            sd = sum(v * v for v in vn2_vals)
            vn2_err[ipt] = math.sqrt(
                max(0.0, sd / len(vn2_vals) - mean * mean)
            ) / math.sqrt(len(vn2_vals))
        else:
            vn2_err[ipt] = float("nan")

        if len(vn4_vals) >= 2:
            mean = sum(vn4_vals) / len(vn4_vals)
            sd = sum(v * v for v in vn4_vals)
            vn4_err[ipt] = math.sqrt(
                max(0.0, sd / len(vn4_vals) - mean * mean)
            ) / math.sqrt(len(vn4_vals))
        else:
            vn4_err[ipt] = float("nan")

    # Write results into the text file (append)
    out_file = Path(out_path)
    with out_file.open("a", encoding="utf-8") as fout:
        fout.write(f"{direct}\t{int(order)}\n")
        for ipt in range(ptBins):
            ptBinCenter = (ipt + 0.5) * dpt + ptMinCut
            print(
                f"{ptBinCenter}\t{vn2[ipt]}\t{vn2_err[ipt]}\t{vn4[ipt]}\t{vn4_err[ipt]}"
            )
            fout.write(
                f"{ptBinCenter}\t{vn2[ipt]}\t{vn2_err[ipt]}\t{vn4[ipt]}\t{vn4_err[ipt]}\n"
            )
        fout.write("\n")

    print(f"Results have been written to '{out_path}'")


def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Compute v_n{2}, v_n{4} vs pT from OSCAR files (Python port of vn_C_pT.c)."
    )
    p.add_argument(
        "direct",
        help=(
            "Buď prefix pro soubory (použije se f'{direct}{i}_fin.oscar'), "
            "nebo přímo adresář s '1_fin.oscar', '2_fin.oscar', ..., "
            "nebo přímo cesta k jednomu '*_fin.oscar'."
        ),
    )
    p.add_argument(
        "--nof",
        type=int,
        default=None,
        help="Počet souborů (NoF). Pokud nezadáš a direct je adresář, autodetekuje se z '*_fin.oscar'.",
    )
    p.add_argument(
        "--noe",
        type=int,
        default=None,
        help="Počet eventů na soubor (NoE). Pokud nezadáš, čte se až do EOF.",
    )
    p.add_argument(
        "--order", type=float, required=True, help="Harmonický řád n (např. 2, 3, ...)"
    )
    p.add_argument(
        "--out",
        default="vn_cumulants_pT.dat",
        help="Output file to append to (default: vn_cumulants_pT.dat)",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Seed pro náhodné rozdělení eventů do subsamplů (default: 0). Použij --seed -1 pro nedeterministické.",
    )
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_argparser().parse_args(argv)

    nof = args.nof
    if _is_probably_single_file(args.direct):
        nof = 1
    if nof is None:
        nof = _autodetect_nof_from_dir(args.direct)
        if nof is None:
            print(
                "Error: '--nof' není zadané a nepodařilo se autodetekovat z adresáře.",
                file=sys.stderr,
            )
            return 2

    seed: Optional[int]
    if args.seed == -1:
        seed = None
    else:
        seed = args.seed

    try:
        vn_C_pT(args.direct, nof, args.noe, args.order, out_path=args.out, seed=seed)
    except RuntimeError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

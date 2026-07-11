"""
Izvoz PRAVIH instanci (iz .mat) u jednostavan, JEZIK-NEUTRALAN tekstualni format.

Cilj: da iste podatke moze da procita i Python i buduci C++/MATLAB port, bez
zavisnosti od scipy/.mat. Format je obican tekst (lako se parsira svuda):

    NAME <ime>
    N_CUSTOMERS <n>
    N_STATIONS <s>
    N_VEHICLES <M>
    Q <kapacitet>        R <potrosnja>     DMAX <domet>
    SPEED <v>            TMAX <t>          PSTART <p>
    SERVICE <svc>        REFUEL <ref>
    CAP <eta_1> ... <eta_s>
    NODES                # index type x y   (type: d/c/f)
    0 d <x> <y>
    1 c <x> <y>
    ...
    EOF

Rastojanja se NE upisuju (euklidska su, identicna originalnoj distance_table
do ~1e-15, sto je potvrdjeno) -> svaki jezik ih preracuna iz koordinata.

Koriscenje:
    python export_instances.py --instances-dir /tmp/METS-Algorithm/Instances \\
        --out instances_txt
"""

import os
import glob
import argparse
from mets_loader import load_mat_instance, instance_set_name


def export_instance(inst, path):
    types = (['d'] + ['c'] * inst.n_customers + ['f'] * inst.n_stations)
    with open(path, 'w') as f:
        f.write(f"NAME {inst.name}\n")
        f.write(f"N_CUSTOMERS {inst.n_customers}\n")
        f.write(f"N_STATIONS {inst.n_stations}\n")
        f.write(f"N_VEHICLES {inst.n_vehicles}\n")
        f.write(f"Q {inst.tank_capacity:g}    R {inst.consumption_rate:g}    DMAX {inst.d_max:g}\n")
        f.write(f"SPEED {inst.speed:g}    TMAX {inst.t_max:g}    PSTART {inst.p_start:g}\n")
        f.write(f"SERVICE {inst.service_time_customer:g}    REFUEL {inst.refueling_time:g}\n")
        f.write("CAP " + " ".join(str(int(c)) for c in inst.station_capacity) + "\n")
        f.write("NODES\n")
        for i in range(inst.n_total):
            f.write(f"{i} {types[i]} {inst.coords[i,0]:.6f} {inst.coords[i,1]:.6f}\n")
        f.write("EOF\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--instances-dir', default='benchmark_instances')
    ap.add_argument('--out', default='instances_txt')
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    n = 0
    for p in sorted(glob.glob(os.path.join(args.instances_dir, "*.mat"))):
        base = os.path.basename(p)
        if base.startswith('jd'):
            continue  # Beijing instances - preskoci za sada
        try:
            inst = load_mat_instance(p)
        except Exception as e:
            print(f"  preskocen {base}: {e}")
            continue
        out = os.path.join(args.out, f"{instance_set_name(p)}_{base.split('_')[1].replace('.mat','')}.txt")
        export_instance(inst, out)
        n += 1
    print(f"Izvezeno {n} instanci u {args.out}/")


if __name__ == "__main__":
    main()

# Archived measurement run `2026-10-03-sandbox-cpu`

## Environment

```
date_utc=2026-10-03T18:06:15Z
python=3.12.0
warp=1.17.0
numpy=2.5.3
cpu=13th Gen Intel(R) Core(TM) i9-13900K
logical_cpus=32
kernel=6.18.40.1-microsoft-standard-WSL2
git_sha=e8e35c4d9769efe58c7af74099aeb187bb996c2a
scripts_dirty=no
step_timeout_s=300
repeats=3
script_sha256:
  c3fc1ac92597faca53fd8b297eef8863356b1b41931960a15aa82cc1e45e071b  /home/wahln/aiprojects/ion-mc-worktrees/v3-001/validation/scripts/warp-architecture/compile_probe.py
  496d0546f52922ed12287a3e5033ee6ab31e974b4a0126c00ad665a5acd878cc  /home/wahln/aiprojects/ion-mc-worktrees/v3-001/validation/scripts/warp-architecture/philox_lib.py
  2ae4ac324f356b6afc6de650f2ae7be8150eefb4dcf62bc9c1ef619709646cbc  /home/wahln/aiprojects/ion-mc-worktrees/v3-001/validation/scripts/warp-architecture/precision.py
  32953749aff468c7bac36f5d8cf990284dfff4abbe8a7ff641ec90bb7d02717d  /home/wahln/aiprojects/ion-mc-worktrees/v3-001/validation/scripts/warp-architecture/probe_pyscope.py
  e17c1868b728bbec084220abef619a72fd6b8845d448529b922b2f7505913db1  /home/wahln/aiprojects/ion-mc-worktrees/v3-001/validation/scripts/warp-architecture/pyprec.py
  6e443c81072ebecb8c8149233b7df6038c836b9e0cf16b6d34b5dd206043cf21  /home/wahln/aiprojects/ion-mc-worktrees/v3-001/validation/scripts/warp-architecture/pyprec2.py
  934ed88358e4d5b2369c82359e35119759c939e1d69192dac12b18d847cef87f  /home/wahln/aiprojects/ion-mc-worktrees/v3-001/validation/scripts/warp-architecture/pyscope_cost.py
  9402533200486b331d27439e57f25ff6862bb26cb23a472a395e57890165249c  /home/wahln/aiprojects/ion-mc-worktrees/v3-001/validation/scripts/warp-architecture/recur.py
  590a92867c169ca650d3c33f0a7efdbbd574ee26c27f6f4546373a5fd4da57c3  /home/wahln/aiprojects/ion-mc-worktrees/v3-001/validation/scripts/warp-architecture/run_toy.py
  ead8fc2c0707dc4b2011f975d8e118d565d6cac4a048d2680ac8c9a3b5f6400d  /home/wahln/aiprojects/ion-mc-worktrees/v3-001/validation/scripts/warp-architecture/structt.py
  5faeeae3905d56a2046be11ce790dfcbea080753066aa656b9c950d30ec24390  /home/wahln/aiprojects/ion-mc-worktrees/v3-001/validation/scripts/warp-architecture/summarize.py
  8fb03f84235762599c2f7fea40d7f080c4741afdf3fbdeacd56df87b2254a8d8  /home/wahln/aiprojects/ion-mc-worktrees/v3-001/validation/scripts/warp-architecture/toy_transport.py
  513905222427c2331dec912141a9a9937639120f4d5d7fccdae65b863b721be5  /home/wahln/aiprojects/ion-mc-worktrees/v3-001/validation/scripts/warp-architecture/run_all.sh
  590895e72c942ad7b5dd53ef95232bff5a4233d269be181a197a0254b2306df1  /home/wahln/aiprojects/ion-mc-worktrees/v3-001/validation/scripts/warp-architecture/../rng/philox.py
  496d0546f52922ed12287a3e5033ee6ab31e974b4a0126c00ad665a5acd878cc  /home/wahln/aiprojects/ion-mc-worktrees/v3-001/validation/scripts/warp-architecture/../rng/philox_lib.py
  1256d5527935da9481debe827d1495673df0ca356c84d5bfe390011f1f998cc5  /home/wahln/aiprojects/ion-mc-worktrees/v3-001/validation/scripts/warp-architecture/../rng/rng_overlap.py
  4883d7b76491491819ecdb43fea58eab85d4722394cc58e4a38a56857c54ea7e  /home/wahln/aiprojects/ion-mc-worktrees/v3-001/validation/scripts/warp-architecture/../rng/rng_seed_dupes.py
```

## Steps

| Step | Command | Exit | Result lines (verbatim, warnings removed) |
|---|---|---|---|
| 01-r1-compile-cold-backward-on | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python compile_probe.py 1` | 0 | enable_backward=True: CPU cold compile 1.40s (fresh cache /home/wahln/.local/share/ionmc-experiment/experiment-v3/tmp/claude-1000/ionmc-warp-cold-ldxsorlb, warp 1.17.0) |
| 01-r2-compile-cold-backward-on | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python compile_probe.py 1` | 0 | enable_backward=True: CPU cold compile 1.37s (fresh cache /home/wahln/.local/share/ionmc-experiment/experiment-v3/tmp/claude-1000/ionmc-warp-cold-ukjw2jyj, warp 1.17.0) |
| 01-r3-compile-cold-backward-on | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python compile_probe.py 1` | 0 | enable_backward=True: CPU cold compile 1.39s (fresh cache /home/wahln/.local/share/ionmc-experiment/experiment-v3/tmp/claude-1000/ionmc-warp-cold-w1mbbhra, warp 1.17.0) |
| 02-r1-compile-cold-backward-off | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python compile_probe.py 0` | 0 | enable_backward=False: CPU cold compile 1.16s (fresh cache /home/wahln/.local/share/ionmc-experiment/experiment-v3/tmp/claude-1000/ionmc-warp-cold-q36_2e3_, warp 1.17.0) |
| 02-r2-compile-cold-backward-off | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python compile_probe.py 0` | 0 | enable_backward=False: CPU cold compile 1.14s (fresh cache /home/wahln/.local/share/ionmc-experiment/experiment-v3/tmp/claude-1000/ionmc-warp-cold-_2wyz9x9, warp 1.17.0) |
| 02-r3-compile-cold-backward-off | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python compile_probe.py 0` | 0 | enable_backward=False: CPU cold compile 1.15s (fresh cache /home/wahln/.local/share/ionmc-experiment/experiment-v3/tmp/claude-1000/ionmc-warp-cold-dj2ywk5a, warp 1.17.0) |
| 03-r1-toy-f32-1thread | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python run_toy.py f32 20000 1` | 0 | prec=f32 N=20000 simulated=20000 threads=1 wall=0.26s hist/s=78419 edep/primary=148.15MeV particles=23844 truncated=0 escaped=0 stack_overflow=0 steps/primary=160 init=0.14s load_module=0.00s host_peak_MiB=319 |
| 03-r2-toy-f32-1thread | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python run_toy.py f32 20000 1` | 0 | prec=f32 N=20000 simulated=20000 threads=1 wall=0.26s hist/s=76738 edep/primary=148.15MeV particles=23844 truncated=0 escaped=0 stack_overflow=0 steps/primary=160 init=0.15s load_module=0.00s host_peak_MiB=319 |
| 03-r3-toy-f32-1thread | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python run_toy.py f32 20000 1` | 0 | prec=f32 N=20000 simulated=20000 threads=1 wall=0.26s hist/s=76187 edep/primary=148.15MeV particles=23844 truncated=0 escaped=0 stack_overflow=0 steps/primary=160 init=0.14s load_module=0.00s host_peak_MiB=319 |
| 04-r1-toy-f64-1thread | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python run_toy.py f64 20000 1` | 0 | prec=f64 N=20000 simulated=20000 threads=1 wall=0.34s hist/s=58307 edep/primary=148.15MeV particles=23852 truncated=0 escaped=0 stack_overflow=0 steps/primary=160 init=0.14s load_module=0.01s host_peak_MiB=351 |
| 04-r2-toy-f64-1thread | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python run_toy.py f64 20000 1` | 0 | prec=f64 N=20000 simulated=20000 threads=1 wall=0.35s hist/s=56721 edep/primary=148.15MeV particles=23852 truncated=0 escaped=0 stack_overflow=0 steps/primary=160 init=0.15s load_module=0.00s host_peak_MiB=352 |
| 04-r3-toy-f64-1thread | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python run_toy.py f64 20000 1` | 0 | prec=f64 N=20000 simulated=20000 threads=1 wall=0.35s hist/s=57966 edep/primary=148.15MeV particles=23852 truncated=0 escaped=0 stack_overflow=0 steps/primary=160 init=0.14s load_module=0.00s host_peak_MiB=352 |
| 05-r1-toy-f32-4threads | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python run_toy.py f32 80000 4` | 0 | prec=f32 N=80000 simulated=80000 threads=4 wall=0.27s hist/s=299423 edep/primary=148.16MeV particles=95286 truncated=0 escaped=0 stack_overflow=0 steps/primary=160 init=0.14s load_module=0.00s host_peak_MiB=483 |
| 05-r2-toy-f32-4threads | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python run_toy.py f32 80000 4` | 0 | prec=f32 N=80000 simulated=80000 threads=4 wall=0.27s hist/s=296824 edep/primary=148.16MeV particles=95286 truncated=0 escaped=0 stack_overflow=0 steps/primary=160 init=0.14s load_module=0.00s host_peak_MiB=483 |
| 05-r3-toy-f32-4threads | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python run_toy.py f32 80000 4` | 0 | prec=f32 N=80000 simulated=80000 threads=4 wall=0.27s hist/s=301643 edep/primary=148.16MeV particles=95286 truncated=0 escaped=0 stack_overflow=0 steps/primary=160 init=0.14s load_module=0.00s host_peak_MiB=483 |
| 06-r1-toy-f32-8threads | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python run_toy.py f32 160000 8` | 0 | prec=f32 N=160000 simulated=160000 threads=8 wall=0.27s hist/s=591856 edep/primary=148.18MeV particles=190309 truncated=0 escaped=0 stack_overflow=0 steps/primary=160 init=0.15s load_module=0.00s host_peak_MiB=648 |
| 06-r2-toy-f32-8threads | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python run_toy.py f32 160000 8` | 0 | prec=f32 N=160000 simulated=160000 threads=8 wall=0.30s hist/s=538851 edep/primary=148.18MeV particles=190309 truncated=0 escaped=0 stack_overflow=0 steps/primary=160 init=0.14s load_module=0.00s host_peak_MiB=648 |
| 06-r3-toy-f32-8threads | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python run_toy.py f32 160000 8` | 0 | prec=f32 N=160000 simulated=160000 threads=8 wall=0.28s hist/s=579494 edep/primary=148.18MeV particles=190309 truncated=0 escaped=0 stack_overflow=0 steps/primary=160 init=0.14s load_module=0.00s host_peak_MiB=648 |
| 07-r1-toy-f32-16threads | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python run_toy.py f32 320000 16` | 0 | prec=f32 N=320000 simulated=320000 threads=16 wall=0.39s hist/s=830374 edep/primary=148.17MeV particles=380911 truncated=0 escaped=0 stack_overflow=0 steps/primary=160 init=0.14s load_module=0.00s host_peak_MiB=978 |
| 07-r2-toy-f32-16threads | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python run_toy.py f32 320000 16` | 0 | prec=f32 N=320000 simulated=320000 threads=16 wall=0.38s hist/s=848205 edep/primary=148.17MeV particles=380911 truncated=0 escaped=0 stack_overflow=0 steps/primary=160 init=0.13s load_module=0.00s host_peak_MiB=978 |
| 07-r3-toy-f32-16threads | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python run_toy.py f32 320000 16` | 0 | prec=f32 N=320000 simulated=320000 threads=16 wall=0.38s hist/s=846536 edep/primary=148.17MeV particles=380911 truncated=0 escaped=0 stack_overflow=0 steps/primary=160 init=0.15s load_module=0.00s host_peak_MiB=978 |
| 08-precision | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python precision.py` | 0 | == float32 sequential accumulation into one voxel (worst case: hottest Bragg-peak voxel) ==<br>n=   100000 rel.err f32 sequential=4.65e-06  (MC rel. stat. err of mean ~ 1.6e-03)<br>n=  1000000 rel.err f32 sequential=1.56e-05  (MC rel. stat. err of mean ~ 5.1e-04)<br>n=  3000000 rel.err f32 sequential=1.33e-03  (MC rel. stat. err of mean ~ 2.9e-04)<br>n= 10000000 rel.err f32 sequential=3.16e-03  (MC rel. stat. err of mean ~ 1.6e-04)<br>== same, split over 40 batch accumulators, then f64 sum of batches ==<br>n=10000000 40 batches: rel.err=3.45e-07<br>== position/direction accumulation over 300 steps: float32 vs float64 ==<br>final \|p32-p64\| over 2000 paths: median=3.54e-05 mm max=1.70e-04 mm; float32 ulp at 300 mm=3.1e-05 mm<br>== energy near end of range: E from residual range R (Bragg-Kleeman R=a E^p, water) ==<br>E= 150.0 MeV ds=1e+00 mm: dE64=5.3070e-01 dE32=5.3072e-01 rel.err=3.7e-05<br>E= 150.0 MeV ds=1e-02 mm: dE64=5.2998e-03 dE32=5.3101e-03 rel.err=1.9e-03<br>E= 150.0 MeV ds=1e-04 mm: dE64=5.2997e-05 dE32=7.6294e-05 rel.err=4.4e-01<br>E=  10.0 MeV ds=1e+00 mm: dE64=5.4802e+00 dE32=5.4802e+00 rel.err=2.2e-07<br>E=  10.0 MeV ds=1e-02 mm: dE64=4.2713e-02 dE32=4.2714e-02 rel.err=2.8e-05<br>E=  10.0 MeV ds=1e-04 mm: dE64=4.2643e-04 dE32=4.2725e-04 rel.err=1.9e-03<br>E=   1.0 MeV ds=1e-02 mm: dE64=2.8257e-01 dE32=2.8257e-01 rel.err=1.3e-07<br>E=   1.0 MeV ds=1e-04 mm: dE64=2.5134e-03 dE32=2.5134e-03 rel.err=3.8e-06<br>E=   0.1 MeV ds=1e-04 mm: dE64=1.5753e-02 dE32=1.5753e-02 rel.err=2.5e-07 |
| 09-rng-seed-dupes | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python ../rng/rng_seed_dupes.py` | 0 | seeds 1,2: identical start states between batches (N=1000000): 260<br>seeds 1000,1001: identical start states between batches (N=1000000): 260<br>seeds 42,43: identical start states between batches (N=1000000): 260<br>distinct pcg(i) for i<1e6: 1000000 |
| 10-rng-overlap-1e5x1000 | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python ../rng/rng_overlap.py 100000 1000` | 0 | N=100000 L=1000 draws=1.000e+08 fraction_of_2^32=0.023 repeated_states=1.143e+06 (1.143%) t=1.5s |
| 11-rng-overlap-1e6x1000 | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python ../rng/rng_overlap.py 1000000 1000` | 0 | N=1000000 L=1000 draws=1.000e+09 fraction_of_2^32=0.233 repeated_states=1.073e+08 (10.734%) t=14.8s |
| 12-rng-overlap-1e6x2000 | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python ../rng/rng_overlap.py 1000000 2000` | 0 | N=1000000 L=2000 draws=2.000e+09 fraction_of_2^32=0.466 repeated_states=3.999e+08 (19.997%) t=29.7s |
| 13-r1-philox-kat-and-cost | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python ../rng/philox.py` | 0 | Random123 KAT pass: True<br>python adapter KAT pass: True<br>kernel==python on 20000 random blocks: True; python adapter 5.7 us/block<br>philox: 2.07 ns/uniform (Warp CPU, 1 thread), mean=0.49995<br>warp-pcg: 2.33 ns/uniform (Warp CPU, 1 thread), mean=0.50004<br>all philox checks passed |
| 13-r2-philox-kat-and-cost | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python ../rng/philox.py` | 0 | Random123 KAT pass: True<br>python adapter KAT pass: True<br>kernel==python on 20000 random blocks: True; python adapter 5.4 us/block<br>philox: 1.99 ns/uniform (Warp CPU, 1 thread), mean=0.49995<br>warp-pcg: 2.30 ns/uniform (Warp CPU, 1 thread), mean=0.50004<br>all philox checks passed |
| 13-r3-philox-kat-and-cost | `/home/wahln/aiprojects/ion-mc-worktrees/v3-001/.venv/bin/python ../rng/philox.py` | 0 | Random123 KAT pass: True<br>python adapter KAT pass: True<br>kernel==python on 20000 random blocks: True; python adapter 5.6 us/block<br>philox: 2.10 ns/uniform (Warp CPU, 1 thread), mean=0.49995<br>warp-pcg: 2.30 ns/uniform (Warp CPU, 1 thread), mean=0.50004<br>all philox checks passed |

## Statistics of repeated steps (successful repeats only)

| Step | Metric | Repeats | Median | Min | Max |
|---|---|---|---|---|---|
| 01-compile-cold-backward-on | cold compile s | 3 | 1.39 | 1.37 | 1.4 |
| 02-compile-cold-backward-off | cold compile s | 3 | 1.15 | 1.14 | 1.16 |
| 03-toy-f32-1thread | hist/s | 3 | 76738 | 76187 | 78419 |
| 04-toy-f64-1thread | hist/s | 3 | 57966 | 56721 | 58307 |
| 05-toy-f32-4threads | hist/s | 3 | 299423 | 296824 | 301643 |
| 06-toy-f32-8threads | hist/s | 3 | 579494 | 538851 | 591856 |
| 07-toy-f32-16threads | hist/s | 3 | 846536 | 830374 | 848205 |
| 13-philox-kat-and-cost | philox ns/uniform | 3 | 2.07 | 1.99 | 2.1 |
| 13-philox-kat-and-cost | warp-pcg ns/uniform | 3 | 2.3 | 2.3 | 2.33 |

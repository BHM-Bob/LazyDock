'''
Date: 2026-09-04
Description: abfe-status: scan an ABFE run output tree and report the progress
             of every FEP window leg (ligand/complex x vdw/coul/bonded):
             finished windows / total, running step, elapsed wall time.
             Read-only; never touches the simulation files.
'''
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Union

PathLike = Union[os.PathLike, str, bytes]


def _read_mdp_param(mdp: Path, key: str) -> Optional[str]:
    """Fetch one parameter value from a GROMACS mdp file (best-effort)."""
    try:
        for line in mdp.read_text(errors='ignore').splitlines():
            line = line.strip()
            if line.startswith(key) and '=' in line:
                return line.split('=', 1)[1].strip()
    except OSError:
        pass
    return None


def _step_state(step_dir: Path, step: str) -> str:
    """Return state of one window step: done / running / pending / errored."""
    finished = step_dir / f'{step}.finished'
    if finished.exists():
        return 'done'
    # step has output (gro for everything after min; min itself) -> running
    out = step_dir / f'{step}.gro'
    if out.exists() or (step_dir / f'{step}.tpr').exists():
        return 'running'
    return 'pending'


def _window_state(win_dir: Path) -> tuple:
    """Scan one lambda window dir; return (state, current_step).

    state: done (prod .finished) / active (some step running) /
           queued (all chain steps done except prod, prod not started yet) /
           pending (nothing started yet).
    current_step: name of the step currently being executed / None
    """
    steps = ['00_min', '01_nvt', '02_npt', '03_npt_norest', 'prod']
    # prod finished -> whole window done
    if (win_dir / 'prod' / 'prod.finished').exists():
        return 'done', None

    # chain steps: 00_min completes on gro (no .finished), others on .finished
    for step in steps[:4]:
        if not _step_completed(win_dir, step):
            sd = win_dir / step
            # queued: only grompp artifacts (mdp) exist, mdrun never started
            has_run = (sd / f'{step}.tpr').exists() or (sd / f'{step}.log').exists()
            return ('active' if has_run else 'queued'), step

    # all prep steps done; now prod
    if (win_dir / 'prod' / 'prod.tpr').exists() or (win_dir / 'prod' / 'prod.gro').exists():
        return 'active', 'prod'
    return 'queued', 'prod'


def _step_completed(win_dir: Path, step: str) -> bool:
    """Whether a step has completed (engine completion rules).

    00_min (minimize) has no .finished marker: completes when 00_min.gro
    exists; every other step completes when its .finished exists.
    """
    sd = win_dir / step
    if step == '00_min':
        return (sd / '00_min.gro').exists()
    return (sd / f'{step}.finished').exists()


# files written during the actual mdrun; static grompp artifacts (mdout.mdp,
# *.mdp) are excluded from mtime-span estimates
_RUN_OUTPUT_SUFFIXES = ('.tpr', '.log', '.gro', '.cpt', '.edr', '.trr', '.xvg')


def _read_log_wall(log_path: Path) -> Optional[float]:
    """Parse mdrun log 'Time: <core> <wall> <pct>' line -> WALL seconds.

    mdrun writes this line only when it finishes normally, so a running
    step never has it (falls back to the mtime-span estimate instead).
    """
    try:
        for line in log_path.read_text(errors='ignore').splitlines():
            if line.strip().startswith('Time:'):
                parts = line.split()
                if len(parts) >= 3:
                    return float(parts[2])   # parts[1] is core t, parts[2] wall t
    except (OSError, ValueError):
        pass
    return None


def _step_mtime_span(step_dir: Path) -> Optional[float]:
    """max(mtime) - min(mtime) over mdrun output files only (seconds).

    Only files written during the actual mdrun count (tpr/log/gro/cpt/
    edr/trr/xvg), so files grompp created earlier (*.mdp, mdout.mdp)
    cannot inflate the span. For a running step this equals
    now - mdrun start; for a finished step without 'Time:' in its log
    (00_min steep) it equals tpr -> last output.
    """
    ts = []
    try:
        for f in step_dir.iterdir():
            if f.is_file() and f.name.endswith(_RUN_OUTPUT_SUFFIXES):
                ts.append(f.stat().st_mtime)
    except OSError:
        return None
    if len(ts) < 2:
        return None
    return max(ts) - min(ts)


def _step_walltime(win_dir: Path, step: str) -> Optional[float]:
    """Wall time of one step.

    - finished step whose log has a 'Time:' line: use it (authoritative).
    - finished step without 'Time:' (00_min steep has no Time line):
      use the mtime span of the step directory (tpr first written ->
      final output).
    - running step: mtime span of the step directory (keeps growing while
      mdrun writes output; the newest file approximates "now").
    Returns None when the step has no usable output at all.
    """
    sd = win_dir / step
    if not sd.is_dir():
        return None
    log = sd / f'{step}.log'
    if log.exists():
        w = _read_log_wall(log)
        if w is not None:
            return w
    return _step_mtime_span(sd)


def _window_walltime(win_dir: Path) -> Optional[float]:
    """Total wall time of one window = sum over its steps.

    Each step's wall time comes from its mdrun log 'Time:' when available
    (01_nvt+ after min), otherwise from the step directory mtime span
    (00_min steep steps, and running steps where log has no Time line yet).
    """
    steps = ['00_min', '01_nvt', '02_npt', '03_npt_norest', 'prod']
    total = 0.0
    any_contribution = False
    for step in steps:
        w = _step_walltime(win_dir, step)
        if w is not None:
            total += w
            any_contribution = True
    return total if any_contribution else None


def _fmt_hms(seconds: float) -> str:
    h, rem = divmod(int(seconds), 3600)
    m, s = divmod(rem, 60)
    if h >= 100:
        return f'{h}h{m:02d}m'
    return f'{h}h{m:02d}m{s:02d}s'


def _running_windows_detail(win_dirs: List[Path], sim_root: Path) -> List[dict]:
    """List windows that are currently running, with their active step."""
    detail = []
    for wd in win_dirs:
        st, step = _window_state(wd)
        if st != 'active':
            continue
        detail.append({
            'window': wd.name,
            'step': step or '',
        })
    return detail


def scan_abfe_progress(out_root: PathLike) -> List[dict]:
    """Scan an ABFE output tree and collect status rows.

    Returns a list of dicts with keys:
        ligand, leg, lam_type, state, n_finished, n_total, current_step, elapsed
    """
    out_root = Path(out_root)
    rows = []
    if not out_root.exists():
        return rows
    # locate ligand/replica roots: {out}/<ligand>/<replica>
    roots = []
    try:
        top_dirs = [p for p in out_root.iterdir() if p.is_dir()]
    except OSError:
        return rows
    for lig_dir in sorted(top_dirs):
        if lig_dir.name.startswith('.'):
            continue
        try:
            rep_dirs = [p for p in lig_dir.iterdir() if p.is_dir()]
        except OSError:
            continue
        for rep_dir in sorted(rep_dirs):
            if rep_dir.name.startswith('.'):
                continue
            roots.append((lig_dir.name, rep_dir.name, lig_dir, rep_dir))
    if not roots:
        # maybe out_root is exactly the ligand dir
        for rep_dir in sorted(p for p in out_root.iterdir() if p.is_dir()):
            if rep_dir.name.startswith('.'):
                continue
            roots.append((out_root.name, rep_dir.name, out_root, rep_dir))

    for lig_name, rep_name, lig_dir, rep_dir in roots:
        for sys_type in ('ligand', 'complex'):
            sim_root = rep_dir / sys_type / 'fep' / 'simulation'
            if not sim_root.is_dir():
                continue
            lam_types = ('vdw', 'coul') if sys_type == 'ligand' else ('vdw', 'coul', 'bonded')
            for lam in lam_types:
                win_dirs = sorted(sim_root.glob(f'{lam}.*'),
                                  key=lambda p: int(p.name.split('.')[1]))
                if not win_dirs:
                    continue
                n_fin = 0
                n_running = 0
                n_queued = 0
                win_walltimes = []
                for wd in win_dirs:
                    st, _ = _window_state(wd)
                    if st == 'done':
                        n_fin += 1
                    elif st == 'active':
                        n_running += 1
                    elif st == 'queued':
                        n_queued += 1
                    wt = _window_walltime(wd)
                    if wt is not None:
                        win_walltimes.append(wt)
                # overall leg state from the "most advanced" window
                leg_done = n_fin == len(win_dirs)
                leg_state = ('done' if leg_done
                             else ('running' if n_running else
                                   ('queued' if n_queued else 'pending')))
                # running windows detail (window + active step)
                running_detail = _running_windows_detail(win_dirs, sim_root)
                # elapsed: sum of per-window wall times (each window = sum of
                # its step mdrun log times); avg = mean per-window wall time
                elapsed = _fmt_hms(sum(win_walltimes)) if win_walltimes else None
                avg = (_fmt_hms(sum(win_walltimes) / len(win_walltimes))
                       if win_walltimes else None)
                rows.append({
                    'ligand': lig_name,
                    'leg': f'{sys_type} {lam}',
                    'state': leg_state,
                    'n_finished': n_fin,
                    'n_total': len(win_dirs),
                    'n_running': n_running,
                    'elapsed': elapsed,      # sum of per-window wall times
                    'avg': avg,              # mean per-window wall time
                    'running': running_detail,
                })
    return rows


def render_status(rows: List[dict]) -> str:
    """Render status rows as a rich table."""
    if not rows:
        return '(no ABFE progress data found under this directory)'
    from rich.console import Console
    from rich.table import Table
    console = Console()
    table = Table(title='ABFE run progress')
    table.add_column('ligand')
    table.add_column('leg', style='cyan')
    table.add_column('state')
    table.add_column('finished')
    table.add_column('running')
    table.add_column('elapsed (sum)', justify='right')
    table.add_column('avg', justify='right')
    for r in rows:
        state = r['state']
        style = 'green' if state == 'done' else ('yellow' if state == 'running' else 'white')
        table.add_row(r['ligand'], r['leg'],
                      f'[{style}]{state}[/{style}]',
                      f"{r['n_finished']}/{r['n_total']}",
                      str(r['n_running']), r['elapsed'] or '-', r.get('avg') or '-')
    console.print(table)
    # running windows detail section
    running_rows = [(r['ligand'], r['leg'], w['window'], w['step'])
                    for r in rows for w in r.get('running', [])]
    if running_rows:
        det = Table(title='running windows (window -> active step)')
        det.add_column('ligand')
        det.add_column('leg', style='cyan')
        det.add_column('window')
        det.add_column('step', style='yellow')
        for lig, leg, win, step in running_rows:
            det.add_row(lig, leg, win, step)
        console.print('')
        console.print(det)
    return ''
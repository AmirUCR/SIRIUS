import os, sys, subprocess, pathlib, tarfile, shutil, signal, argparse

CACHE_BASE = pathlib.Path(os.path.expanduser("~/.cache/sirius"))
ENV_VERSION = "py27env-8"  # bump if you change the tarball
CACHE_DIR = CACHE_BASE / ENV_VERSION

def _cache_dirs():
    """Return (current_cache, all_versioned_caches)."""
    cur = CACHE_DIR
    all_dirs = []
    if CACHE_BASE.exists():
        for p in CACHE_BASE.iterdir():
            if p.is_dir() and p.name.startswith("py27env-"):
                all_dirs.append(p)
    return cur, all_dirs

def purge_cache(all_versions: bool = False) -> int:
    cur, all_dirs = _cache_dirs()
    targets = all_dirs if all_versions else [cur]
    for p in targets:
        try:
            shutil.rmtree(p)
        except FileNotFoundError:
            pass
    return 0

def purge_cache_main():
    parser = argparse.ArgumentParser(description="Purge SIRIUS cached Python 2.7 environment")
    parser.add_argument("--all", action="store_true",
                        help="remove all cached py27 environments (all versions)")
    args = parser.parse_args()
    rc = purge_cache(all_versions=args.all)
    print("[SIRIUS] Cache purged." + (" (all versions)" if args.all else ""))
    return rc

def _gc_old_caches():
    """Delete old py27env-* caches that are not the current ENV_VERSION."""
    _, all_dirs = _cache_dirs()
    for p in all_dirs:
        if p.name != ENV_VERSION:
            try:
                shutil.rmtree(p)
            except Exception:
                pass

def ensure_py27_env(pkg_root: pathlib.Path) -> pathlib.Path:
    env_python = CACHE_DIR / "bin" / "python"
    if env_python.exists():
        return env_python

    # Prefer .xz (smaller), fall back to .gz
    candidates = [
        pkg_root / "resources" / "py27_env.tar.xz",
        pkg_root / "resources" / "py27_env.tar.gz",
    ]
    tarball = next((p for p in candidates if p.exists()), None)
    if not tarball:
        sys.exit(f"[SIRIUS] Missing tarball under {pkg_root}/resources")

    if CACHE_DIR.exists():
        shutil.rmtree(CACHE_DIR)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    import tarfile, os, subprocess

    def _is_within(base_dir: str, member_name: str) -> bool:
        # purely lexical (no symlink resolution)
        base = os.path.abspath(base_dir)
        dest = os.path.abspath(os.path.normpath(os.path.join(base_dir, member_name)))
        return dest == base or dest.startswith(base + os.sep)

    with tarfile.open(tarball, "r:*") as tf:
        for m in tf.getmembers():
            # reject absolute paths or '..' traversal lexically
            if not _is_within(str(CACHE_DIR), m.name) or os.path.isabs(m.name):
                sys.exit("[SIRIUS] Unsafe path in py27 archive; aborting.")
        tf.extractall(CACHE_DIR)

    # Fix symlinks/paths created by conda-pack
    cu = CACHE_DIR / "bin" / "conda-unpack"
    if cu.exists():
        try:
            subprocess.check_call([str(cu)], cwd=str(CACHE_DIR))
        except Exception:
            pass

    # Ensure a runnable python
    p = CACHE_DIR / "bin" / "python"
    if not p.exists():
        for q in (CACHE_DIR / "bin").glob("python2*"):
            if q.exists():
                (CACHE_DIR / "bin" / "python").symlink_to(q.name)
                p = CACHE_DIR / "bin" / "python"
                break
    if not p.exists():
        sys.exit("[SIRIUS] Failed to unpack bundled Python 2.7 environment (no bin/python found).")
    p.chmod(p.stat().st_mode | 0o111)
    return p

def _die(msg: str):
    sys.stderr.write(f"[SIRIUS] {msg}\n")
    return 1

def main():
    if len(sys.argv) > 1 and sys.argv[1] in {"--purge-cache", "purge-cache"}:
        # Allow: `sirius --purge-cache` and `sirius purge-cache`
        all_flag = "--all" in sys.argv[2:]
        return purge_cache(all_versions=all_flag)

    _gc_old_caches()  # auto-clean stale caches on normal runs
    
    root             = pathlib.Path(__file__).resolve().parent
    py2              = ensure_py27_env(root)
    binary           = root / "build" / "sirius"
    lib_dir          = root / "lib"
    diversifier      = root / "GeneDiversifier" / "sequenceDiversifier.py"
    default_host_csv = root / "GeneDiversifier" / "data" / "human_codon_usage.csv"

    if not diversifier.exists():
        sys.exit(f"[SIRIUS] Missing GeneDiversifier at {diversifier}")
    if not binary.exists():
        sys.exit(f"[SIRIUS] Missing binary at {binary} (wheel not built or CMake failed)")

    env = os.environ.copy()
    env["GENE_DIVERSIFIER_PY"]     = str(diversifier)
    env["GENE_DIVERSIFIER_PYTHON"] = str(py2)
    env["SIRIUS_DEFAULT_HOST_CSV"] = str(default_host_csv)
    env["LD_LIBRARY_PATH"] = str(lib_dir) + (os.pathsep + env["LD_LIBRARY_PATH"] if "LD_LIBRARY_PATH" in env else "")

    # Start child in its own process group so signals hit the whole tree
    try:
        proc = subprocess.Popen(
            [str(binary)] + sys.argv[1:],
            env=env,
            preexec_fn=os.setsid  # new PGID for child
        )
        try:
            return proc.wait()
        except KeyboardInterrupt:
            # Forward SIGINT to the child's process group, exit cleanly
            try:
                os.killpg(proc.pid, signal.SIGINT)
            except ProcessLookupError:
                pass  # already exited
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                # Escalate if it hangs
                try:
                    os.killpg(proc.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
            return 130
    except KeyboardInterrupt:
        # Ctrl-C before/after spawn — exit quietly
        return 130

if __name__ == "__main__":
    sys.exit(main())
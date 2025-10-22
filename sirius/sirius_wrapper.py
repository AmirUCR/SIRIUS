import os, sys, subprocess, pathlib, tarfile, shutil, signal

def ensure_py27_env(pkg_root: pathlib.Path) -> pathlib.Path:
    """Return path to bundled py27 python; unpack tarball if needed."""
    env_dir = pkg_root / "resources" / "py27_env"
    py2     = env_dir / "bin" / "python"
    if py2.exists():
        return py2

    tarball = pkg_root / "resources" / "py27_env.tar.gz"
    if not tarball.exists():
        sys.exit(f"[SIRIUS] Missing tarball: {tarball}")

    # Try unpacking into site-packages first
    try:
        env_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tarball, "r:gz") as tf:
            tf.extractall(env_dir)
        if py2.exists():
            return py2
    except Exception:
        pass  # fall back to user cache

    # Fallback: unpack to user cache (handles read-only site-packages)
    cache_dir = pathlib.Path(os.path.expanduser("~/.cache/sirius/py27_env"))
    cache_dir.mkdir(parents=True, exist_ok=True)
    # Clean stale cache if needed
    if any(cache_dir.iterdir()):
        shutil.rmtree(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tarball, "r:gz") as tf:
        tf.extractall(cache_dir)
    py2_cache = cache_dir / "bin" / "python"
    if not py2_cache.exists():
        sys.exit("[SIRIUS] Failed to unpack bundled Python 2.7 environment.")
    return py2_cache

def _die(msg: str):
    sys.stderr.write(f"[SIRIUS] {msg}\n")
    return 1

def main():
    root             = pathlib.Path(__file__).resolve().parent
    py2              = ensure_py27_env(root)
    binary           = root / "build" / "sirius"
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
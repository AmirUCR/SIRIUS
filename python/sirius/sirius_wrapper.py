import os, sys, subprocess, pathlib, tarfile, shutil, signal

CACHE_BASE = pathlib.Path(os.path.expanduser("~/.cache/sirius"))
ENV_VERSION = "py27env-1"  # bump if you change the tarball
CACHE_DIR = CACHE_BASE / ENV_VERSION

def ensure_py27_env(pkg_root: pathlib.Path) -> pathlib.Path:
    env_python = CACHE_DIR / "bin" / "python"
    if env_python.exists():
        return env_python

    tarball = pkg_root / "resources" / "py27_env.tar.gz"
    if not tarball.exists():
        sys.exit(f"[SIRIUS] Missing tarball: {tarball}")

    if CACHE_DIR.exists():
        shutil.rmtree(CACHE_DIR)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tarball, "r:gz") as tf:
        tf.extractall(CACHE_DIR)

    env_python = CACHE_DIR / "bin" / "python"
    if not env_python.exists():
        sys.exit("[SIRIUS] Failed to unpack bundled Python 2.7 environment.")
    return env_python


def _die(msg: str):
    sys.stderr.write(f"[SIRIUS] {msg}\n")
    return 1

def main():
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
# setup.py
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
import subprocess, pathlib, sys, shutil, os
from glob import glob

class CMakeBuild(build_ext):
    def run(self):
        root = pathlib.Path(__file__).resolve().parent
        src_dir = root / "cpp"
        build_temp = root / "build_temp"
        build_temp.mkdir(exist_ok=True)

        # Configure + build
        subprocess.check_call(["cmake", str(src_dir)], cwd=build_temp)
        subprocess.check_call(["cmake", "--build", ".", "--config", "Release"], cwd=build_temp)

        # Binary produced by CMake
        binary_src = build_temp / "sirius"
        if not binary_src.exists():
            sys.exit(f"[BUILD][FATAL] Missing binary: {binary_src}")

        # Prepare package staging dir
        build_py = self.get_finalized_command("build_py")
        self.run_command("build_py")
        pkg_root = pathlib.Path(build_py.build_lib) / "sirius"
        (pkg_root / "build").mkdir(parents=True, exist_ok=True)
        shutil.copy2(binary_src, pkg_root / "build" / "sirius")

        # Copy OR-Tools shared libs into sirius/lib/
        ortools_lib_dir = root / "cpp" / "ortools" / "lib"
        pkg_lib_dir = pkg_root / "lib"
        pkg_lib_dir.mkdir(parents=True, exist_ok=True)

        # Copy libortools and any nearby dependent .so’s (abseil, protobuf, etc.)
        matched = False
        for p in glob(str(ortools_lib_dir / "libortools.so*")):
            shutil.copy2(p, pkg_lib_dir / pathlib.Path(p).name)
            matched = True
        if not matched:
            sys.exit(f"[BUILD][FATAL] No libortools.so* found under {ortools_lib_dir}")

        # Optionally include other .so dependencies from the same dir
        for p in glob(str(ortools_lib_dir / "*.so*")):
            name = os.path.basename(p)
            if name.startswith("libortools.so"):
                continue
            shutil.copy2(p, pkg_lib_dir / name)

setup(
    name="sirius",
    version="1.0",
    packages=find_packages(where="python", include=["sirius", "sirius.*"]),
    package_dir={"": "python"},
    include_package_data=True,
    package_data={
        "sirius": [
            "build/sirius",
            "resources/py27_env.tar.gz",
            "lib/*.so*",
        ],
        "sirius.GeneDiversifier": ["*.py", "data/*.csv"],
    },
    ext_modules=[Extension("sirius_dummy", sources=[])],
    cmdclass={"build_ext": CMakeBuild},
)

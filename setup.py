from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
import subprocess, pathlib, shutil, os, re, stat, sys

root = pathlib.Path(__file__).resolve().parent
long_description = (root / "README.md").read_text()

class CMakeBuild(build_ext):
    def build_extension(self, ext):
        src_dir    = root / "cpp"
        build_dir  = pathlib.Path(self.build_temp) / "cmake-build"
        build_dir.mkdir(parents=True, exist_ok=True)

        # === Configure CMake ===
        ortools_args = [
            "-DBUILD_DEPS=ON",
            "-DUSE_COINOR=OFF", "-DUSE_HIGHS=OFF", "-DUSE_SCIP=OFF",
            "-DUSE_BOP=OFF", "-DUSE_GLOP=OFF", "-DUSE_MATH_OPT=OFF",
            "-DUSE_PDLP=OFF", "-DBUILD_TESTING=OFF",
            "-DBUILD_SAMPLES=OFF", "-DBUILD_EXAMPLES=OFF",
            "-DUSE_DOTNET_8=OFF",
        ]
        py = sys.executable
        cmake_args = [
            "-DCMAKE_BUILD_TYPE=Release",
            "-DCMAKE_POSITION_INDEPENDENT_CODE=ON",
            "-DPYBIND11_FINDPYTHON=ON",
            f"-DPython_EXECUTABLE={py}",
            f"-DPython3_EXECUTABLE={py}",
            f"-DPYTHON_EXECUTABLE={py}",
            "-DPython3_FIND_VIRTUALENV=FIRST",
        ]
        for var in ("CMAKE_C_FLAGS","CMAKE_CXX_FLAGS",
                    "CMAKE_EXE_LINKER_FLAGS","CMAKE_MODULE_LINKER_FLAGS","CMAKE_SHARED_LINKER_FLAGS"):
            if os.environ.get(var):
                cmake_args.append(f"-D{var}={os.environ[var]}")

        print("[BUILD] Configuring CMake…")
        subprocess.check_call(["cmake", "-S", str(src_dir), "-B", str(build_dir)] + ortools_args + cmake_args)

        print("[BUILD] Building…")
        jobs = os.environ.get("CMAKE_BUILD_PARALLEL_LEVEL", "4")
        subprocess.check_call(["cmake", "--build", str(build_dir), "--config", "Release", "--parallel", jobs, "--verbose"])

        # === 1) Install the pybind extension into the ABI-tagged path ===
        ext_path = pathlib.Path(self.get_ext_fullpath(ext.name))
        ext_path.parent.mkdir(parents=True, exist_ok=True)
        so_candidates = list(build_dir.glob("_sirius*.so"))
        if not so_candidates:
            raise RuntimeError("Could not find _sirius*.so in build dir")
        shutil.copyfile(so_candidates[0], ext_path)
        os.chmod(ext_path, os.stat(ext_path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

        # === 2) Stage the native executable + deps into build_lib/sirius ===
        pkg_stage = pathlib.Path(self.build_lib) / "sirius"
        bin_out   = pkg_stage / "build"
        lib_out   = pkg_stage / "lib"
        bin_out.mkdir(parents=True, exist_ok=True)
        lib_out.mkdir(parents=True, exist_ok=True)

        exe = build_dir / "sirius"  # thanks to CMAKE_RUNTIME_OUTPUT_DIRECTORY
        if not exe.exists():
            raise RuntimeError(f"did not find built executable at {exe}")
        target_exe = bin_out / "sirius"
        shutil.copy2(exe, target_exe)

        # bundle non-system .so used by the exe
        def is_system(path: str) -> bool:
            return path.startswith(("/lib", "/usr/lib", "/usr/lib64"))
        ldd_out = subprocess.check_output(["ldd", str(exe)], text=True, stderr=subprocess.STDOUT)
        for line in ldd_out.splitlines():
            m = re.search(r"=>\s*(/[^ ]+)", line)
            if m:
                p = m.group(1)
                if not is_system(p):
                    dst = lib_out / pathlib.Path(p).name
                    if not dst.exists():
                        shutil.copy2(p, dst)

        strip = shutil.which("strip")
        if strip:
            for path in [ext_path, target_exe] + list((lib_out).glob("*.so*")):
                try:
                    subprocess.check_call([strip, "--strip-unneeded", str(path)])
                except Exception:
                    pass

        # RPATH so the exe finds ../lib at runtime
        patchelf = shutil.which("patchelf")
        if patchelf:
            subprocess.check_call([patchelf, "--set-rpath", "$ORIGIN/../lib:$ORIGIN", str(target_exe)])
        else:
            print("[BUILD] Warning: patchelf not found; wrapper will fallback to LD_LIBRARY_PATH.")

setup(
    name="sirius-bio",
    version="1.4",
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(where="python", include=["sirius", "sirius.*"]),
    package_dir={"": "python"},
    include_package_data=True,
    package_data={
        "sirius": [
            "resources/py27_env.tar.xz",
            "build/sirius",
            "lib/*",
        ],
        "sirius.GeneDiversifier": ["data/*.csv"],
    },
    cmdclass={"build_ext": CMakeBuild},
    # dummy extension just to trigger build_ext; CMake builds the real .so
    ext_modules=[Extension("sirius._sirius", sources=["cpp/src/pybind_module.cpp"])],
    entry_points={"console_scripts": [
        "sirius = sirius.sirius_wrapper:main",
        "sirius-bio = sirius.sirius_wrapper:main"
        ]},
)

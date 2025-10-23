from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
import subprocess, pathlib, sys, shutil, os
from glob import glob

class CMakeBuild(build_ext):
    def run(self):
        root = pathlib.Path(__file__).resolve().parent
        src_dir = root / "cpp"
        build_temp = root / "build_temp"

        if build_temp.exists():
            shutil.rmtree(build_temp)

        build_temp.mkdir(exist_ok=True)

        # Configure + build
        print("[BUILD] Configuring CMake...")
        subprocess.check_call(["cmake", str(src_dir), 
                               "-DBUILD_DEPS=ON", 
                               "-DUSE_COINOR=OFF", 
                               "-DUSE_HIGHS=OFF", 
                               "-DUSE_SCIP=OFF", 
                               "-DBUILD_SAMPLES=OFF", 
                               "-DBUILD_EXAMPLES=OFF", 
                               "-DUSE_DOTNET_8=OFF"], cwd=build_temp)

        print("[BUILD] Building with CMake...")
        subprocess.check_call(["cmake", "--build", ".", "--config", "Release"], cwd=build_temp)

        # Binary produced by CMake
        binary_src = build_temp / "sirius"
        if not binary_src.exists():
            sys.exit(f"[BUILD][FATAL] Missing binary: {binary_src}")

        # Prepare package staging dir
        build_py = self.get_finalized_command("build_py")
        self.run_command("build_py")
        pkg_root = pathlib.Path(build_py.build_lib) / "sirius"
        
        # Copy executable
        exe_staging_dir = pkg_root / "build"
        exe_staging_dir.mkdir(parents=True, exist_ok=True)
        exe_path = exe_staging_dir / "sirius"
        print(f"[BUILD] Copying executable to {exe_path}")
        shutil.copy2(binary_src, exe_path)

        # Copy libraries
        ortools_lib_dir = root / "cpp" / "ortools" / "lib"
        pkg_lib_dir = pkg_root / "lib"
        pkg_lib_dir.mkdir(parents=True, exist_ok=True)
        print(f"[BUILD] Copying libraries to {pkg_lib_dir}")

        matched = False 
        # grab unversioned, SONAME, and fully-versioned files
        for pat in ("libortools.so", "libortools.so.*"):
            for p in glob(str(ortools_lib_dir / pat)):
                print(f"[BUILD] Copying {p}")
                shutil.copy2(p, pkg_lib_dir / os.path.basename(p))
                matched = True
        if not matched:
            sys.exit(f"[BUILD][FATAL] No libortools.so* found under {ortools_lib_dir}")

        # also copy sibling deps (absl/protobuf, etc.) if present
        for p in glob(str(ortools_lib_dir / "*.so*")):
            name = os.path.basename(p)
            if name.startswith("libortools.so"):
                continue
            print(f"[BUILD] Copying dependency {p}")
            shutil.copy2(p, pkg_lib_dir / name) 

        # Set the RPATH of the executable so auditwheel can find
        # the libraries we just copied. This is only needed on Linux.
        if sys.platform == "linux":
            print(f"[BUILD] Setting RPATH on {exe_path}")
            try:
                # The executable is in <pkg_root>/build/sirius
                # The libraries are in <pkg_root>/lib/
                # The relative path from the exe's dir to the lib dir is ../lib
                # $ORIGIN is a linker token meaning "the directory of the executable"
                rpath = "$ORIGIN/../lib"
                subprocess.check_call([
                    "patchelf",
                    "--set-rpath",
                    rpath,
                    str(exe_path)
                ])
                print(f"[BUILD] RPATH set to {rpath}")
            except (subprocess.CalledProcessError, FileNotFoundError) as e:
                print(f"[BUILD][FATAL] patchelf failed. This is required for manylinux wheels.")
                print(f"[BUILD][FATAL] Please install patchelf (e.g., 'yum install patchelf')")
                print(f"[BUILD][FATAL] Error: {e}")
                sys.exit(1)

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
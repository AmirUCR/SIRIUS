from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
import subprocess, pathlib, sys, shutil

class CMakeBuild(build_ext):
    def run(self):
        root = pathlib.Path(__file__).resolve().parent
        build_temp = root / "build_temp"
        build_temp.mkdir(exist_ok=True)

        subprocess.check_call(["cmake", str(root)], cwd=build_temp)
        subprocess.check_call(["cmake", "--build", ".", "--config", "Release"], cwd=build_temp)

        binary_src = root / "sirius" / "build" / "sirius"
        if not binary_src.exists():
            sys.exit(f"[BUILD][FATAL] Missing binary: {binary_src}")

        build_py = self.get_finalized_command("build_py")
        self.run_command("build_py")  # ensures build/lib exists
        pkg_build_dir = pathlib.Path(build_py.build_lib) / "sirius" / "build"
        pkg_build_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(binary_src), str(pkg_build_dir / "sirius"))

setup(
    name="sirius",
    version="1.0.0",
    packages=find_packages(include=["sirius", "sirius.*"]),
    include_package_data=True,
    package_data={
        "sirius": ["build/sirius", "resources/py27_env.tar.gz"],
        "sirius.GeneDiversifier": ["*.py", "data/*.csv"],
    },
    ext_modules=[Extension("sirius_dummy", sources=[])],
    cmdclass={"build_ext": CMakeBuild},
)

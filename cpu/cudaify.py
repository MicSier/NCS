#!/usr/bin/env python3
"""Generate a g++-compilable copy of NCS (header + examples) for the CPU emulation.

Usage:  python3 cpu/cudaify.py [OUT_DIR]      (default: build-cpu)

Rewrites, in NCS.h, gl.h, examples/**/*.cu and cpu/tests/**/*.cu:
  kernel<<<grid, block[, shmem[, stream]]>>>(args)
      -> cpu_launch(cpu_cfg(grid, block[, ...]), kernel, args)
  extern __shared__ T name[];
      -> T* name = ncs_cpu::dyn_shared<T>();
The directory layout is mirrored so relative includes ("../../NCS.h") still work.
Compile the result with:  g++ -std=c++17 -O3 -fopenmp -x c++ -Icpu <file>.cu
"""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "build-cpu").resolve()

# `name <<<cfg>>> (` -- tolerant of the "<< <" / ">> >" spacing Visual Studio produces.
LAUNCH = re.compile(r"(\b[A-Za-z_]\w*)\s*<<\s*<(.*?)>>\s*>\s*\(\s*(\)?)", re.S)
EXTERN_SHARED = re.compile(r"extern\s+__shared__\s+([\w:]+)\s+(\w+)\s*\[\s*\]\s*;")


def rewrite(text: str) -> str:
    def launch(m):
        kernel, cfg, no_args = m.group(1), " ".join(m.group(2).split()), m.group(3)
        head = f"cpu_launch(cpu_cfg({cfg}), {kernel}"
        return head + ")" if no_args else head + ", "

    text = LAUNCH.sub(launch, text)
    text = EXTERN_SHARED.sub(r"\1* \2 = ncs_cpu::dyn_shared<\1>();", text)
    return text


def main():
    sources = [ROOT / "NCS.h", ROOT / "gl.h",
               *sorted((ROOT / "examples").rglob("*.cu")),
               *sorted((ROOT / "cpu" / "tests").rglob("*.cu"))]
    for src in sources:
        dst = OUT / src.relative_to(ROOT)
        dst.parent.mkdir(parents=True, exist_ok=True)
        text = src.read_text(encoding="utf-8-sig")  # some examples carry a UTF-8 BOM
        dst.write_text(rewrite(text), encoding="utf-8")
    print(f"cudaify: wrote {len(sources)} files to {OUT}")


if __name__ == "__main__":
    main()

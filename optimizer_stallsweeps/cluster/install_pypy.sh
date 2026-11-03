#!/bin/bash
# Install PyPy (a faster Python for pure-Python number crunching like this
# optimiser) into ~/pypy-dbf, from the PyPy project's own download server.
# No root needed.  Pick the current Python-3.11 build for Linux x86_64 at
# https://pypy.org/download.html and pass its name, e.g.
#     ./install_pypy.sh pypy3.11-v7.3.20-linux64
# Then:  DBF_PYPY=1 sbatch slurm/check.sbatch   (does PyPy match the PC?)
set -e
name="${1:-pypy3.11-v7.3.20-linux64}"
url="https://downloads.python.org/pypy/${name}.tar.bz2"
dest="$HOME/pypy-dbf"
tmp=$(mktemp -d)
echo "downloading $url"
curl -fL --retry 3 -o "$tmp/pypy.tar.bz2" "$url"
mkdir -p "$dest"
tar -xjf "$tmp/pypy.tar.bz2" -C "$dest" --strip-components=1
rm -rf "$tmp"
"$dest/bin/pypy3" -c 'import platform; print("installed", platform.python_implementation(), platform.python_version())'

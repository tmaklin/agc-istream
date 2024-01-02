#!/bin/sh
set -euxo pipefail
agc_path=$1
g++ agc_istream.cpp -o agc_istream -I $agc_path/src/core -std=c++17 -lzstd -lagc -lz -L $agc_path

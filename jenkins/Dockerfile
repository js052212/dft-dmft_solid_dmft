# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:unstable-ubuntu-clang

# Install missing system dependencies
RUN apt-get update && apt-get install -y meson ninja-build python3-setuptools python3-skimage python3-tomli

ENV OMP_NUM_THREADS=1

# Install upstream TRIQS apps that solid_dmft depends on
# cthyb (also needed at import time by python/solid_dmft/dmft_tools/solvers/hubbardI_interface.py)
RUN cd / && git clone -b unstable --depth 1 https://github.com/TRIQS/cthyb.git cthyb.src \
    && mkdir -p cthyb.build && cd cthyb.build \
    && cmake ../cthyb.src \
    && make -j4 install

# dft_tools
RUN cd / && git clone -b unstable --depth 1 https://github.com/TRIQS/dft_tools.git dft_tools.src \
    && mkdir -p dft_tools.build && cd dft_tools.build \
    && cmake ../dft_tools.src \
    && make -j4 install

# hubbardI
RUN cd / && git clone -b unstable --depth 1 https://github.com/TRIQS/hubbardI.git hubbardI.src \
    && mkdir -p hubbardI.build && cd hubbardI.build \
    && cmake ../hubbardI.src \
    && make install

# Hartree-Fock solver
RUN cd / && git clone -b unstable --depth 1 https://github.com/TRIQS/hartree_fock.git hartree.src \
    && mkdir -p hartree.build && cd hartree.build \
    && cmake ../hartree.src \
    && make install

# maxent
RUN cd / && git clone -b unstable --depth 1 https://github.com/TRIQS/maxent.git maxent.src \
    && mkdir -p maxent.build && cd maxent.build \
    && cmake ../maxent.src \
    && make install

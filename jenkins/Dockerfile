# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:unstable-ubuntu-clang

# Install missing system dependencies
RUN apt-get update && apt-get install -y meson ninja-build python3-setuptools python3-skimage python3-tomli python3-sphinx-autorun

ENV OMP_NUM_THREADS=1

# The upstream TRIQS apps that solid_dmft imports at runtime (cthyb, dft_tools,
# dftkit, hubbardI, hartree_fock, maxent) are no longer cloned and built here.
# They are consumed from their kept unstable installs on the shared /mnt/triqs
# mount via PYTHONPATH, set up in jenkins/Jenkinsfile.

FROM flow123d/install-gnu:4.0.1 as installer
FROM flow123d/ci-gnu:4.0.0a01_0b2509

# Install system packages and GUI dependencies
RUN sudo apt-get update && apt install -y \
    python3 python3-pip python3-venv \
    curl libglu1-mesa redis git

# Set up Python virtual environment
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

RUN sudo python3 -m pip install --upgrade numpy>=1.23.3
# install mpi4py using mpi in base container
COPY --from=installer /.mpiccpath /.mpiccpath
RUN env MPICC=$(cat /.mpiccpath) python3 -m pip install --no-cache-dir --ignore-installed  mpi4py

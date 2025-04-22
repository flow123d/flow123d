# This Dockerfile generates benchmark meshes
FROM ubuntu:latest

# Install Python and venv
RUN apt-get update && apt-get install -y python3 python3-pip python3-venv

# Create virtual environment and activate it
RUN python3 -m venv /opt/venv

# Upgrade pip in the virtual environment
RUN /opt/venv/bin/pip install --upgrade pip

# Add venv to PATH so later RUN commands use it
ENV PATH="/opt/venv/bin:$PATH"

RUN python3 -m pip install gmsh
COPY gmsh.sh /usr/local/bin/gmsh
RUN chmod u+x /usr/local/bin/gmsh

CMD gmsh

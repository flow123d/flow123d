# Stage 1: Build stage (source_image)
ARG source_image
FROM ${source_image} as source_image

# Stage 2: Production image
ARG base_image
FROM ${base_image}

# Build-time metadata as defined at http://label-schema.org
ARG flow_version
ARG git_hash
ARG build_date
ARG source_location
ARG flow_install_location=/opt/flow123d

LABEL org.label-schema.build-date="$build_date" \
      org.label-schema.name="Flow123d" \
      org.label-schema.description="Transport Processes in Fractured Media" \
      org.label-schema.url="http://flow123d.github.io/" \
      org.label-schema.vcs-ref="$git_hash" \
      org.label-schema.vcs-url="https://github.com/flow123d/flow123d" \
      org.label-schema.vendor="Technical University of Liberec" \
      org.label-schema.version="$flow_version" \
      org.label-schema.schema-version="1.0"

# Copy entire installation from source image to /opt/flow123d
COPY --from=source_image ${source_location} ${flow_install_location} 

# ---------------------------------------------------------
# Install shared libraries (non recursively) into /usr/local/lib
# ---------------------------------------------------------
# Copy only files from ${flow_install_location}/lib (for example, *.so)
RUN cp ${flow_install_location}/lib/*.so /usr/local/lib/ && \
    rm -f ${flow_install_location}/lib/*.so

# ---------------------------------------------------------
# Install Python extension modules (e.g., libpy.so) into Python’s platlib directory
# ---------------------------------------------------------
RUN cp ${flow_install_location}/lib/python/*.so \
       $(python -c "import sysconfig; print(sysconfig.get_path('platlib'))") && \
    rm -f ${flow_install_location}/lib/python/*.so

# ---------------------------------------------------------
# Install pure Python code (e.g., flow123d package) into Python’s purelib directory
# ---------------------------------------------------------
RUN cp ${flow_install_location}/lib/python/src/*.py \
       $(python -c "import sysconfig; print(sysconfig.get_path('purelib'))")

# Copy additional local files (e.g., user-installed pip packages)
COPY --from=source_image /root/.local /root/.local

# Update PATH to include binaries from the local installation and user-local
ENV PATH="/root/.local/bin:${flow_install_location}/bin:${PATH}"

# Optionally, save the image version info for runtime use:
# RUN echo ${flow_version} > /.dockerversion
      

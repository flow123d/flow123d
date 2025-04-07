# This Dockerfile will install flow123d into a flow123d/install image
ARG source_image
ARG base_image
FROM ${source_image} as source_image

ARG base_image
FROM ${base_image}

ARG base_image
ARG source_location
ARG flow_version
ARG flow_install_location=/opt/flow123d
ARG git_hash
ARG build_date

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
       $(python3 -c "import sysconfig; print(sysconfig.get_path('platlib'))") && \
    rm -f ${flow_install_location}/lib/python/*.so

# ---------------------------------------------------------
# Install pure Python code (e.g., flow123d package) into Python’s purelib directory
# ---------------------------------------------------------
RUN cp -r ${flow_install_location}/lib/python/py123d \
       $(python3 -c "import sysconfig; print(sysconfig.get_path('purelib'))")

# Copy additional local files (e.g., user-installed pip packages)
#COPY --from=source_image /root/.local /root/.local


# Update systemwide dynamic linker setting.
# Add additional library directories for the dynamic linker
RUN echo "/usr/local/lib" > /etc/ld.so.conf.d/flow123d.conf && \
    for dir in /usr/local/*; do \
      if [ -d "$dir/lib" ]; then \
         echo "$dir/lib" >> /etc/ld.so.conf.d/flow123d.conf; \
      fi; \
    done && ldconfig
ENV LD_LIBRARY_PATH=/usr/local/lib:/usr/lib:${LD_LIBRARY_PATH}

# Update PATH to include binaries from the local installation and user-local
ENV PATH="/root/.local/bin:${flow_install_location}/bin:${PATH}"

# Set path to flow123d root dir
ENV FLOW123D_DIR="/opt/flow123d"

# ---------------------------------------------------------
# Install Python dependencies from requirements.txt in the build context
# ---------------------------------------------------------
COPY requirements.txt /tmp/requirements.txt
RUN pip3 install --no-cache-dir -r /tmp/requirements.txt && rm /tmp/requirements.txt

#Optionally, save the image version info for runtime use:
#RUN echo ${flow_version} > /.dockerversion
      

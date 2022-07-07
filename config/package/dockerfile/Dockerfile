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

# Build-time metadata as defined at http://label-schema.org
LABEL org.label-schema.build-date="$build_date" \
      org.label-schema.name="Flow123d" \
      org.label-schema.description="Transport Processes in Fractured Media" \
      org.label-schema.url="http://flow123d.github.io/" \
      org.label-schema.vcs-ref="$git_hash" \
      org.label-schema.vcs-url="https://github.com/flow123d/flow123d" \
      org.label-schema.vendor="Technical University of Liberec" \
      org.label-schema.version="$flow_version" \
      org.label-schema.schema-version="1.0"
      
      
COPY --from=source_image /${source_location} ${flow_install_location} 

# update PATH
ENV PATH="${flow_install_location}/bin:${PATH}"

# save the image version
RUN echo ${flow_version} > /.dockerversion

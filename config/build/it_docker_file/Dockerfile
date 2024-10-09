# This Dockerfile will install flow123d into a flow123d/install image
ARG source_image

FROM ${source_image} as source_image

ARG flow_version
ARG flow_install_location=/opt/flow123d
ARG git_hash
ARG build_date

LABEL org.label-schema.build-date="$build_date" \
      org.label-schema.name="Integration flow123d image" \
      org.label-schema.vcs-ref="$git_hash" \
      org.label-schema.version="$flow_version" \
      org.label-schema.schema-version="1.0"
      
      
COPY --from=source_image /opt/flow123d /opt/flow123d

COPY --from=source_image /root/.local /root/.local

ENV PATH="/root/.local/bin:${flow_install_location}/bin:${PATH}"

WORKDIR /opt/flow123d

CMD ["./config/build/project/src/linux/fterm.sh", "rel"]
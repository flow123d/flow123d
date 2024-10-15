ARG source_image
FROM ${source_image} as source_image

ARG flow_version
ARG flow_install_location=/opt/flow123d
ARG git_hash
ARG build_date

LABEL org.label-schema.build-date=${build_date} \
      org.label-schema.name="Integration flow123d image" \
      org.label-schema.vcs-ref=${git_hash} \
      org.label-schema.version=${flow_version} \
      org.label-schema.schema-version="1.0"

FROM source_image as final

COPY --from=source_image /opt/flow123d ${flow_install_location}
COPY --from=source_image /root/.local /root/.local

ENV PATH="/root/.local/bin:${flow_install_location}/bin:${PATH}"

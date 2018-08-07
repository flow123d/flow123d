# This template will alter given image by creating new username
# and by executing a custom setup.sh script

FROM flow123d/@docker_image@


COPY . /tmp/
USER root
RUN /tmp/setup.sh

# change user from root(0:0) to a current user
USER @uid@:@gid@

# when running container specify --login flag
# this will automatically load bash.bashrc and profile files in /etc/
CMD ["/bin/bash", "--login"]

# Script that mimicks the release build on the Jenkins CI. 
# Can hlep to debug some build errors happening during release build.

flow_repo=`pwd`/flow123d_repo
git clone --depth=100 -b master https://github.com/flow123d/flow123d.git ${flow_repo}
${flow_repo}/bin/fterm rel_gnu -- --privileged -di --name contgnurelease --volume ${flow_repo}:/opt/flow123d/flow123d
#docker exec contgnurelease git clone -b master https://github.com/flow123d/flow123d.git /opt/flow123d/flow123d
docker exec contgnurelease cp /opt/flow123d/flow123d/config/config-jenkins-docker-release.cmake /opt/flow123d/flow123d/config.cmake
docker exec contgnurelease make -C /opt/flow123d/flow123d -j4 all

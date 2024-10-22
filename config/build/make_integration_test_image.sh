#!/bin/bash
# Usage:
# ./make_integration_test_image.sh <environment> <release_tag> [push]

# Stop on first error
set -e
set -x

environment=$1              # gnu | intel
release_tag=$2              # version
docker_command=$3           # push

if [ "${docker_command}" == "push" ]
then
    docker_command="docker push"
else
    docker_command="echo docker push"
fi


flow_repo_host="$( cd "$( dirname "${BASH_SOURCE[0]}" )"/../.. && pwd )"
cd ${flow_repo_host}

source_image="stepanmoc/ci-${environment}:${release_tag}"
target_image="stepanmoc/it-${environment}:${release_tag}"

if [[ "$(docker images -q ${target_image} 2> /dev/null)" != "" ]]; then
    docker rmi -f ${target_image}
fi

docker pull ${source_image}

docker build --build-arg source_image=${source_image} \
             --build-arg flow_version=${release_tag} \
             --build-arg git_hash=$(git rev-parse --short HEAD) \
             --build-arg build_date=$(date -u +"%Y-%m-%dT%H:%M:%SZ") \
             --tag ${target_image} \
             ${flow_repo_host}/config/build/it_docker_file

container_id=$(docker run -d ${target_image})

docker cp ${flow_repo_host}/tests ${container_id}:/opt/flow123d/bin/tests
docker cp ${flow_repo_host}/config/build/run_tests_with_it_image.sh ${container_id}:/opt/flow123d/bin

docker commit ${container_id} ${target_image}

docker stop ${container_id}
docker rm ${container_id}

${docker_command} ${target_image}

echo "Docker image '${target_image}' has been built and pushed."

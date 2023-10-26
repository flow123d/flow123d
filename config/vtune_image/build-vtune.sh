#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
dev_ver=`cat ${SCRIPT_DIR}/../build/image_tag`
default_img="flow123d/flow-dev-gnu-rel:${dev_ver}"
src_img=${1:-$default_img}
docker build ${SCRIPT_DIR} --build-arg src_img=${src_img} --tag flow123d/flow-dev-gnu-vtune:${dev_ver}

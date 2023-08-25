#!/bin/bash
# Script will import docker image into system

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

function check_image {
    if [[ "$(docker images -q $1 2> /dev/null)" != "" ]]; then
        read -r -p "Image with name $1 already exists. Do you want to remove this image? (y/n): " response
        echo
        if [[ $response =~ ^[Yy]$ ]]; then
            docker rmi -f $1
        fi
    fi
}

# get image path and import into to machine
echo "Pulling docker image '@IMAGE_TAG@'"
check_image "@IMAGE_TAG@"
docker pull "@IMAGE_TAG@"
docker images --format "{{.Repository}}:{{.Tag}}" | grep "@IMAGE_TAG@"

if [ $? -eq 0 ]; then
    echo "Installation of Flow123d finished successfully."
    echo "Run Flow123d using script fterm.sh or flow123d.sh in bin folder."
    echo "You can start by printing the version of Flow123d:"
    echo "  bin/fterm.sh run --version"
else
    echo "Error during installation"
    exit 1
fi

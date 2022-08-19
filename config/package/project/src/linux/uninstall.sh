#!/bin/bash
# Script will import remove docker image from the system

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


# get image path and import into to machine
echo "Removing docker images '@IMAGE_TAG@' and @IMAGE_TAG@:user"
docker rmi -f "@IMAGE_TAG@"
docker rmi -f "@IMAGE_TAG@:user"

echo "Removing install dir: $CWD"
rm -rfI "$CWD"
#!/bin/bash

set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# This could fail if the image is not there, so I temporarily set +e
set +e
OLDIMAGE=`docker images | awk '{print $1, $2, $3}' | grep '^layer-raman-ir-tool\ latest\ ' | awk '{print $3}'`
set -e

# To build container
docker build -t layer-raman-ir-tool "$DIR/.."

# Check the new image name
NEWIMAGE=`docker images | awk '{print $1, $2, $3}' | grep '^layer-raman-ir-tool\ latest\ ' | awk '{print $3}'`

LAYER_CONTAINER=`docker ps --filter="name=layer-raman-ir-tool-instance" -q`
if [ "$LAYER_CONTAINER" != "" ]
then
    docker kill "$LAYER_CONTAINER"
fi  

# If it built correctly, remove the old image (otherwise it remains there as `<none>`)
# Note: I need to do it after killing the corresponding container
# Also, avoid to remove it if nothing had changed, caching was used, and so
# the new one is the same as the old one
if [ "$OLDIMAGE" != "" -a "$OLDIMAGE" != "$NEWIMAGE" ]
then
    docker rmi "$OLDIMAGE"
fi

# To launch container
docker run -d -p 8091:80 --rm --name=layer-raman-ir-tool-instance layer-raman-ir-tool

# Pass '-n' to avoid opening a new browser window
if [ "$1" != "-n" ]
then
    # Give it a second to let apache start
    sleep 1
    python -c "import webbrowser; webbrowser.open('http://localhost:8091')"
    echo "Browser opened at http://localhost:8091"
    echo "Pass -n to avoid opening it"
else
    echo "You can access the webservice at:"
    echo "http://localhost:8091"
fi

echo ""
echo "You can kill the service with:"
echo "docker kill layer-raman-ir-tool-instance"

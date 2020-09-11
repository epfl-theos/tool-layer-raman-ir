#!/bin/bash

set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# To build container
docker build -t layer-raman-tool "$DIR/.."

LAYER_CONTAINER=`docker ps --filter="name=layer-raman-tool-instance" -q`
if [ "$LAYER_CONTAINER" != "" ]
then
    docker kill "$LAYER_CONTAINER"
fi  

# To launch container
docker run -d -p 8091:80 --rm --name=layer-raman-tool-instance layer-raman-tool

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
echo "docker kill layer-raman-tool-instance"

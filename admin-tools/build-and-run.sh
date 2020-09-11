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

echo "You can access the webservice at:"
echo "http://localhost:8091"
echo ""
echo "You can kill the service with:"
echo "docker kill layer-raman-tool-instance"

python -c "import webbrowser; webbrowser.open('http://localhost:8091')"

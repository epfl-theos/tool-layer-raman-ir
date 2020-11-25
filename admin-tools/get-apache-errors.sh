#!/bin/bash

if [ "$1" == "--reload" ]
then
    while true
    do
        docker exec -it layer-raman-ir-tool-instance tail -f /var/log/apache2/error.log
        echo "Reloading in 1 sec, press CTRL+C to stop..."
        sleep 1
    done
else
    docker exec -it layer-raman-ir-tool-instance tail -f /var/log/apache2/error.log
fi

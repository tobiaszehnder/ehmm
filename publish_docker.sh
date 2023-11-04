#!/bin/bash

SCRIPT_DIR=$(cd $(dirname $0); pwd)

docker build -t bigdatainbiomedicine/inspect-ehmm $SCRIPT_DIR && \
    docker push bigdatainbiomedicine/inspect-ehmm
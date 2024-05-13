FROM ubuntu:20.04

MAINTAINER Jared Lewis <jared.lewis@climate-resource.com>

ARG DEBIAN_FRONTEND=noninteractive
ARG CMAQ_VERSION="5.0.2"
ARG CMAQ_BIN="Linux2_armgfort"
ENV TZ=Etc/UTC

RUN apt-get update && \
    apt-get install -y build-essential gfortran m4 csh git jq wget libmpich-dev libnetcdff-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/cmaq
COPY templates /opt/cmaq/templates/
COPY scripts/build_ioapi.sh /opt/cmaq/scripts/build_ioapi.sh

RUN ls /opt/cmaq && CMAQ_VERSION=${CMAQ_VERSION} BIN=${CMAQ_BIN} bash /opt/cmaq/scripts/build_ioapi.sh

ENTRYPOINT ["/bin/bash"]
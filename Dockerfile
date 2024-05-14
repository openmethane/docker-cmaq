FROM ubuntu:20.04

MAINTAINER Jared Lewis <jared.lewis@climate-resource.com>

ARG DEBIAN_FRONTEND=noninteractive
ARG CMAQ_VERSION="5.0.2"
ARG BIN="Linux2_armgfort"
ENV TZ=Etc/UTC

RUN apt-get update && \
    apt-get install -y build-essential gfortran m4 csh git jq wget libmpich-dev libnetcdff-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/cmaq
COPY templates /opt/cmaq/templates/

COPY scripts/build_00_ioapi.sh /opt/cmaq/scripts/build_00_ioapi.sh
RUN BIN=${BIN} bash /opt/cmaq/scripts/build_00_ioapi.sh

COPY scripts/build_10_cmaq.sh /opt/cmaq/scripts/build_10_cmaq.sh
RUN CMAQ_VERSION=${CMAQ_VERSION} BIN=${BIN} bash /opt/cmaq/scripts/build_10_cmaq.sh

ENTRYPOINT ["/bin/bash"]
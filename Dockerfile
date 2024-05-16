FROM ubuntu:20.04

MAINTAINER Jared Lewis <jared.lewis@climate-resource.com>

ARG DEBIAN_FRONTEND=noninteractive
ENV CMAQ_VERSION="5.0.2"
ENV ARCH="aarch64"
ENV BIN="Linux2_${ARCH}gfort"
ENV TZ=Etc/UTC

RUN apt-get update && \
    apt-get install -y build-essential gfortran m4 csh git jq wget libopenmpi-dev libnetcdff-dev unzip nano && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/cmaq
COPY templates/ioapi /opt/cmaq/templates/ioapi
COPY scripts/build_00_ioapi.sh /opt/cmaq/scripts/build_00_ioapi.sh
RUN bash /opt/cmaq/scripts/build_00_ioapi.sh

COPY templates/CMAQ /opt/cmaq/templates/CMAQ
COPY scripts/build_10_cmaq.sh /opt/cmaq/scripts/build_10_cmaq.sh
RUN bash /opt/cmaq/scripts/build_10_cmaq.sh

COPY scripts/build_11_cmaq_adj.sh /opt/cmaq/scripts/build_11_cmaq_adj.sh
COPY src/cmaq_adj /opt/cmaq/cmaq_adj
RUN bash /opt/cmaq/scripts/build_11_cmaq_adj.sh

ENTRYPOINT ["/bin/bash"]
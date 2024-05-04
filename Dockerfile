FROM ghcr.io/climate-resource/wrf:4.5.1

MAINTAINER Jared Lewis <jared.lewis@climate-resource.com>

ARG DEBIAN_FRONTEND=noninteractive
ARG CMAQ_VERSION="5.0.2"
ENV TZ=Etc/UTC

WORKDIR /opt/cmaq
COPY scripts templates /opt/cmaq/

RUN CMAQ_VERSION=${CMAQ_VERSION} bash /opt/cmaq/scripts/build_cmaq.sh

ENTRYPOINT ["/bin/bash"]
FROM ubuntu:16.04


RUN apt-get update && \
    apt-get install -y g++ make libargtable2-dev libgsl-dev

ADD . /source
RUN cd /source &&\
    make &&\
    install -m 0755 TG_3D /usr/bin &&\
    install -m 0755 docker-entrypoint.sh /usr/bin &&\
    mkdir -p /tg

WORKDIR /tg

ENTRYPOINT ["/usr/bin/docker-entrypoint.sh"]

FROM ubuntu:latest

RUN apt-get update && apt-get install -y gcc g++ make zlib1g-dev

ADD . /src
WORKDIR /src
RUN make all

FROM ubuntu:latest
COPY --from=0 /src/bwbble /usr/bin/bwbble

RUN apt-get update && apt-get install -y libgomp1

WORKDIR /
ENTRYPOINT ["/usr/bin/bwbble"]
CMD []
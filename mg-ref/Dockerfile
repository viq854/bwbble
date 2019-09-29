FROM ubuntu:latest

RUN apt-get update && apt-get install -y gcc g++ make zlib1g-dev

ADD . /src
WORKDIR /src
RUN make all

FROM ubuntu:latest
COPY --from=0 /src/comb /usr/bin/comb
COPY --from=0 /src/data_prep /usr/bin/data_prep
COPY --from=0 /src/sam_pad /usr/bin/sam_pad
ADD ./help.sh /usr/bin/help

RUN apt-get update && apt-get install -y libgomp1

WORKDIR /
CMD ["/usr/bin/help"]
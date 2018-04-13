FROM ubuntu:trusty
LABEL app="omegaMap"
LABEL description="Estimating diversifying selection and functional constraint in the presence of recombination"
LABEL CITATION="https://dx.doi.org/10.1534/genetics.105.044917"
LABEL maintainer="Daniel Wilson"
LABEL build-type="From source"
RUN apt-get -yqq update
RUN apt-get -yqq install make g++
ENV MKDIR /tmp/omegaMap
RUN mkdir $MKDIR
COPY . $MKDIR
WORKDIR $MKDIR
RUN make
RUN mv $MKDIR/omegaMap /usr/bin/
RUN rm $MKDIR/*.o
WORKDIR /home/ubuntu
ENTRYPOINT ["/usr/bin/omegaMap"]

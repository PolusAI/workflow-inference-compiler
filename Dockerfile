FROM ubuntu:22.04

# Install miniconda
RUN apt-get update && apt-get install -y wget
RUN MINICONDA="Miniconda3-latest-Linux-x86_64.sh" && \
    wget --quiet https://repo.continuum.io/miniconda/$MINICONDA && \
    chmod +x $MINICONDA && \
    ./$MINICONDA -b -p /miniconda && \
    rm -f $MINICONDA
ENV PATH /miniconda/bin:$PATH

# Install wic
RUN apt-get install -y git
RUN git clone --recursive https://github.com/PolusAI/workflow-inference-compiler.git
#RUN conda create --name wic
#RUN conda activate wic
# The above command prints
# CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'.
# It still prints that even if we run `conda init bash` first.
# But this is a Docker image; we don't necessarily need to additionally isolate
# wic within a conda environment. Let's just install it globally!
RUN cd workflow-inference-compiler && ./conda_devtools.sh
RUN cd workflow-inference-compiler && pip install -e ".[test]"

ADD Dockerfile .

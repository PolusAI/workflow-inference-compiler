#!/bin/bash -e
# Unlike the default github action runners which spin up a brand new machine every time,
# self-hosted runners are not necessarily isolated. In particular, the docker cache
# needs to be manually updated by explicitly performing docker pull/build commands.
# NOTE: For now use the following explicit list. In the future, consider using the
# cwl-utils library to recursively search for the dockerPull: tags within each workflow.

sudo docker build --no-cache --pull -f Dockerfile_debian -t polusai/wic_debian ..
sudo docker build --no-cache --pull -f Dockerfile_debian_pypy -t polusai/wic_debian_pypy ..
sudo docker build --no-cache --pull -f Dockerfile_ubuntu -t polusai/wic_ubuntu ..
sudo docker build --no-cache --pull -f Dockerfile_ubuntu_pypy -t polusai/wic_ubuntu_pypy ..

sudo docker build --no-cache --pull -f Dockerfile_fedora -t polusai/wic_fedora ..
sudo docker build --no-cache --pull -f Dockerfile_fedora_pypy -t polusai/wic_fedora_pypy ..
sudo docker build --no-cache --pull -f Dockerfile_redhat -t polusai/wic_redhat ..
sudo docker build --no-cache --pull -f Dockerfile_redhat_pypy -t polusai/wic_redhat_pypy ..
sudo docker build --no-cache --pull -f Dockerfile_amazon -t polusai/wic_amazon ..
sudo docker build --no-cache --pull -f Dockerfile_amazon_pypy -t polusai/wic_amazon_pypy ..

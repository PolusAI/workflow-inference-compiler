#!/bin/bash -e
# Unlike the default github action runners which spin up a brand new machine every time,
# self-hosted runners are not necessarily isolated. In particular, the docker cache
# needs to be manually updated by explicitly performing docker pull commands.
# NOTE: For now use the following explicit list. In the future, consider using the
# cwl-utils library to recursively search for the dockerPull: tags within each workflow.
docker pull jakefennick/autodock_vina
docker pull jakefennick/data
docker pull jakefennick/scripts
docker pull jakefennick/biosimspace
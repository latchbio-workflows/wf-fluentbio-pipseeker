# DO NOT CHANGE
from 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:fe0b-main

workdir /tmp/docker-build/work/

shell [ \
    "/usr/bin/env", "bash", \
    "-o", "errexit", \
    "-o", "pipefail", \
    "-o", "nounset", \
    "-o", "verbose", \
    "-o", "errtrace", \
    "-O", "inherit_errexit", \
    "-O", "shift_verbose", \
    "-c" \
]
env TZ='Etc/UTC'
env LANG='en_US.UTF-8'

arg DEBIAN_FRONTEND=noninteractive

run curl -L https://fbs-public.s3.us-east-2.amazonaws.com/public-pipseeker-releases/pipseeker-v3.2.0/pipseeker-v3.2.0-linux.tar.gz -o pipseeker.tar.gz &&\
    tar -xzvf pipseeker.tar.gz &&\
    mv pipseeker*/pipseeker /bin/ &&\
    chmod +x /bin/pipseeker

run apt-get install unzip

# Latch SDK
# DO NOT REMOVE
run pip install latch==2.36.10
run mkdir /opt/latch

# Copy workflow data (use .dockerignore to skip files)
copy . .latch/* /root/


# Latch workflow registration metadata
# DO NOT CHANGE
arg tag
# DO NOT CHANGE
env FLYTE_INTERNAL_IMAGE $tag

workdir /root

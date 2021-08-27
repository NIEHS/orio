# Build and test docker image
# This docker file can be used to mount the source directory and do a make init and make test in a neutral environment

FROM centos/python-38-centos7

USER root
# Create a group and user to run our app
ARG APP_USER=appuser
RUN groupadd -r ${APP_USER} && useradd --no-log-init -r -g ${APP_USER} ${APP_USER}

# Install packages needed to run your application
RUN yum update -y

RUN yum -y install libpcre3 mime-support postgresql-client nano less gcc

RUN pip install requests
RUN pip install clint
RUN pip install pymysql

WORKDIR /orio

CMD tail -f /dev/null

# docker build -t niehs/orio:latest .
# docker run  -v `pwd`:/orio -t niehs/orio:latest
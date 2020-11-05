FROM materialscloud/tools-barebone:1.1.0

LABEL maintainer="Giovanni Pizzi <giovanni.pizzi@epfl.ch>"

COPY ./requirements.txt /home/app/code/requirements.txt
RUN pip3 install -U 'pip>=10' setuptools wheel

# install packages as normal user (app, provided by passenger)
USER app
WORKDIR /home/app/code
# Install pinned versions of packages
RUN pip3 install --user -r requirements.txt

# Go back to root. Also, it should remain as user root for startup
USER root
# Copy various files: configuration, user templates, the actual python code, ...

COPY ./config.yaml /home/app/code/webservice/static/config.yaml
COPY ./user_templates/ /home/app/code/webservice/templates/user_templates/
COPY ./compute/ /home/app/code/webservice/compute/
COPY ./user_static/ /home/app/code/webservice/user_static/

# Set proper permissions on files just copied
RUN chmod -R o+rX /home/app/code/webservice/
RUN chown -R app:app /home/app/code/webservice/

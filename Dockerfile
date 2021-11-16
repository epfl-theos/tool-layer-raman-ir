FROM dokku/tools-barebone

LABEL maintainer="Giovanni Pizzi <giovanni.pizzi@epfl.ch>"

# Copy various files: configuration, user templates, the actual python code, ...

COPY ./config.yaml /home/app/code/webservice/static/config.yaml
COPY ./user_templates/ /home/app/code/webservice/templates/user_templates/
COPY ./compute/ /home/app/code/webservice/compute/
COPY ./user_static/ /home/app/code/webservice/user_static/

# Set proper permissions on files just copied
RUN chmod -R o+rX /home/app/code/webservice/
RUN chown -R app:app /home/app/code/webservice/

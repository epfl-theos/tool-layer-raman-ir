FROM tools-barebone

MAINTAINER Giovanni Pizzi <giovanni.pizzi@epfl.ch>

ENV HOME /root
COPY ./app_requirements.txt $HOME/app_requirements.txt
RUN pip3 install -r $HOME/app_requirements.txt
ENV HOME /home/app

COPY ./config.yaml /home/app/code/webservice/static/config.yaml
COPY ./user_static/ /home/app/code/webservice/user_static/
COPY ./user_templates/ /home/app/code/webservice/templates/user_templates/
COPY ./compute/ /home/app/code/webservice/compute/

# Set proper permissions
RUN chown -R app:app /home/app/code


#### Add custom-tool's specific commands here below

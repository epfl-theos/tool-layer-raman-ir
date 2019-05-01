FROM tools-barebone

MAINTAINER Giovanni Pizzi <giovanni.pizzi@epfl.ch>

COPY ./config.yaml /home/app/code/webservice/static/config.yaml
COPY ./user_templates/ /home/app/code/webservice/templates/user_templates/
COPY ./compute/ /home/app/code/webservice/compute/
COPY ./user_static/ /home/app/code/webservice/user_static/

WORKDIR $HOME
COPY ./app_requirements.txt app_requirements.txt
RUN pip3 install -r $HOME/app_requirements.txt

# Set proper permissions
RUN chown -R app:app $HOME

#### Add custom-tool's specific commands here below
